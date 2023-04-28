/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/clusterization/clusterization_algorithm.hpp"
#include "traccc/clusterization/spacepoint_formation.hpp"
#include "traccc/cuda/clusterization/clusterization_algorithm.hpp"
#include "traccc/cuda/seeding/seeding_algorithm.hpp"
#include "traccc/cuda/seeding/track_params_estimation.hpp"
#include "traccc/cuda/utils/stream.hpp"
#include "traccc/efficiency/seeding_performance_writer.hpp"
#include "traccc/io/read_cells.hpp"
#include "traccc/io/read_digitization_config.hpp"
#include "traccc/io/read_geometry.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/io/../../../src/csv/read_cells.hpp"
#include "traccc/options/common_options.hpp"
#include "traccc/options/full_tracking_input_options.hpp"
#include "traccc/options/handle_argument_errors.hpp"
#include "traccc/performance/collection_comparator.hpp"
#include "traccc/performance/container_comparator.hpp"
#include "traccc/performance/timer.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"

#include "traccc/geometry/geometry.hpp"
#include "traccc/geometry/digitization_config.hpp"
// VecMem include(s).
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/host_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>

//Root Reader
#include <TROOT.h>
#include <TChain.h>
#include "IDreader.h"
#include "DataList.h"

// System include(s).
#include <exception>
#include <iomanip>
#include <iostream>

namespace po = boost::program_options;

traccc::io::cell_reader_output read_cells_from_htt_event(DataList *data, IDreader *idr, const traccc::geometry *geom, const traccc::digitization_config *dcfg, vecmem::memory_resource *mr)
{
    traccc::cell_collection_types::host result_cells;
    traccc::cell_module_collection_types::host result_modules;
    
    std::vector<unsigned int> cellCounts;
    cellCounts.reserve(5000);

    if (mr != nullptr) {
        result_modules = traccc::cell_module_collection_types::host{0, mr};
    } else {
        result_modules = traccc::cell_module_collection_types::host(0);
    }
    result_modules.reserve(5000);

    // Create a cell collection, which holds on to a flat list of all the cells
    // and the position of their respective cell counter & module.
    std::vector<std::pair<traccc::io::csv::cell, unsigned int>> allCells;
    allCells.reserve(50000);

    for(int ihit=0; ihit<data->m_Hits_ ; ihit++){
        if(data->m_Hits_m_hitType[ihit]!=0) continue;
        if(data->m_Hits_m_detType[ihit]!=1) continue;

        uint atlasid = data->m_Hits_m_identifierHash[ihit];
        unsigned long int geometry_id = idr->Give50muPixelDetectorID(atlasid,1);
        int channel0 = data->m_Hits_m_phiIndex[ihit];
        int channel1 = data->m_Hits_m_etaIndex[ihit];
        double timestamp = 0;
        double value = data->m_Hits_m_ToT[ihit];
        //traccc::cell current_cell {channel0, channel1, value, timestamp, geometry_id};
        traccc::io::csv::cell current_cell {geometry_id, ihit, channel0, channel1, timestamp, value};

        auto rit = std::find_if(result_modules.rbegin(), result_modules.rend(),
            [&current_cell](const traccc::cell_module& mod) {
                return mod.module == current_cell.geometry_id;
            });

        if (rit == result_modules.rend()) {
            // Add new cell and new cell counter if a new module is found
            const traccc::cell_module mod = get_module(current_cell, geom, dcfg);
            allCells.push_back({current_cell, result_modules.size()});
            result_modules.push_back(mod);
            cellCounts.push_back(1);
        } else {
            // Add a new cell and update cell counter if repeat module is found
            const unsigned int pos =
                std::distance(result_modules.begin(), rit.base()) - 1;
            allCells.push_back({current_cell, pos});
            ++(cellCounts[pos]);
        }
    }

    // Transform the cellCounts vector into a prefix sum for accessing
    // positions in the result vector.
    std::partial_sum(cellCounts.begin(), cellCounts.end(), cellCounts.begin());

    // The total number cells.
    const unsigned int totalCells = allCells.size();

    // Construct the result collection.
    if (mr != nullptr) {
        result_cells = traccc::cell_collection_types::host{totalCells, mr};
    } else {
        result_cells = traccc::cell_collection_types::host(totalCells);
    }

    // Member "-1" of the prefix sum vector
    unsigned int nCellsZero = 0;
    // Fill the result object with the read csv cells
    for (unsigned int i = 0; i < totalCells; ++i) {
        const traccc::io::csv::cell& c = allCells[i].first;

        // The position of the cell counter this cell belongs to
        const unsigned int& counterPos = allCells[i].second;

        unsigned int& prefix_sum_previous =
            counterPos == 0 ? nCellsZero : cellCounts[counterPos - 1];
        result_cells[prefix_sum_previous++] = traccc::cell{
            c.channel0, c.channel1, c.value, c.timestamp, counterPos};
    }

    if (cellCounts.size() == 0) {
        return {result_cells, result_modules};
    }
    /* This is might look a bit overcomplicated, and could be made simpler by
     * having a copy of the prefix sum vector before incrementing its value when
     * filling the vector. however this seems more efficient, but requires
     * manually setting the 1st & 2nd modules instead of just the 1st.
     */
    const auto comp = [](const traccc::cell& c1, const traccc::cell& c2) {
        return c1.channel1 < c2.channel1;
    };
        
    // Sort the cells belonging to the first module.
    std::sort(result_cells.begin(), result_cells.begin() + nCellsZero, comp);
    // Sort the cells belonging to the second module.
    std::sort(result_cells.begin() + nCellsZero,
              result_cells.begin() + cellCounts[0], comp);

    // Sort cells belonging to all other modules.
    for (unsigned int i = 1; i < cellCounts.size() - 1; ++i) {
        std::sort(result_cells.begin() + cellCounts[i - 1],
                  result_cells.begin() + cellCounts[i], comp);
    }

    return {result_cells, result_modules};
}

int seq_run(const traccc::full_tracking_input_config& i_cfg,
            const traccc::HTTSim_options& common_opts, bool run_cpu) {

    //testing
    int nevents = common_opts.events;

    // Read the surface transforms
    auto surface_transforms = traccc::io::read_geometry(i_cfg.detector_file);

    // Read the digitization configuration file
    auto digi_cfg =
        traccc::io::read_digitization_config(i_cfg.digitization_config_file);

    // Output stats
    uint64_t n_cells = 0;
    uint64_t n_modules = 0;
    // uint64_t n_clusters = 0;
    uint64_t n_measurements = 0;
    uint64_t n_spacepoints = 0;
    uint64_t n_spacepoints_cuda = 0;
    uint64_t n_seeds = 0;
    uint64_t n_seeds_cuda = 0;

    // Memory resources used by the application.
    vecmem::host_memory_resource host_mr;
    vecmem::cuda::host_memory_resource cuda_host_mr;
    vecmem::cuda::device_memory_resource device_mr;
    traccc::memory_resource mr{device_mr, &cuda_host_mr};

    traccc::cuda::stream stream;

    vecmem::cuda::async_copy copy{stream.cudaStream()};
    
    traccc::clusterization_algorithm ca(host_mr);
    traccc::spacepoint_formation sf(host_mr);
    
    traccc::cuda::clusterization_algorithm ca_cuda( mr, copy, stream, common_opts.target_cells_per_partition);
    //traccc::cuda::seeding_algorithm sa_cuda(mr, copy, stream);
    traccc::cuda::track_params_estimation tp_cuda(mr, copy, stream);

    traccc::performance::timing_info elapsedTimes;

    // Read from root file
    TChain *myC = new TChain("HTTEventTree");
    myC->Add((common_opts.input_directory.c_str()));
    DataList *data = new DataList(myC);
    const std::string geometry_file_name = traccc::io::data_directory() + i_cfg.detector_file;
    IDreader *idr = new IDreader(geometry_file_name.c_str());

    //Gen number of entries and number of events from root file
    Long64_t nentries = data->fChain->GetEntriesFast();
    Int_t nfiles = (nevents < nentries ? nevents : nentries);

    for (Long64_t jentry=0; jentry<nfiles;jentry++) {
        traccc::io::cell_reader_output alt_read_out_per_event;
        traccc::clusterization_algorithm::output_type measurements_per_event;
        traccc::spacepoint_formation::output_type spacepoints_per_event;
        traccc::spacepoint_collection_types::buffer spacepoints_cuda_buffer(0, *mr.host);

        //Load the current entry tree from data into ientry
        Long64_t ientry = data->LoadTree(jentry);
        if (ientry < 0) break;
        data->fChain->GetEntry(jentry); 
        
        alt_read_out_per_event = 
            read_cells_from_htt_event(data, idr, &surface_transforms, &digi_cfg, &cuda_host_mr);
        //Create the event data holders (vectors)
        traccc::cell_collection_types::host cells_events_host = alt_read_out_per_event.cells;
        traccc::cell_module_collection_types::host module_events_host = alt_read_out_per_event.modules;

        /*for(int ihit=0; ihit<data->m_Hits_ ; ihit++){
            if(data->m_Hits_m_hitType[ihit]!=0) continue;
            if(data->m_Hits_m_detType[ihit]!=1) continue;

            uint atlasid = data->m_Hits_m_identifierHash[ihit];
            unsigned long int geometry_id = idr->Give50muPixelDetectorID(atlasid,1);
            int channel0 = data->m_Hits_m_phiIndex[ihit];
            int channel1 = data->m_Hits_m_etaIndex[ihit];
            double timestamp = 0;
            double value = data->m_Hits_m_ToT[ihit];

            traccc::cell current_cell {channel0, channel1, value, timestamp, geometry_id};
            traccc::io::csv::cell csv_cell {geometry_id, ihit, channel0, channel1, timestamp, value};
            traccc::cell_module current_module = traccc::io::csv::get_module(csv_cell, &surface_transforms, &digi_cfg);

            cells_events_host.push_back(current_cell);
            module_events_host.push_back(current_module);
        }*/

        //Allocate GPU memory and copy the cell and module data from the host type to buffer type
        traccc::cell_collection_types::buffer cells_events_buffer (cells_events_host.size(), mr.main);
        //std::cout<<"Event "<< jentry <<" has "<<cells_events_host.size()<<" cells."<<std::endl;
        //std::cout<<"Event "<< jentry <<" has "<<module_events_host.size()<<" modules."<<std::endl;
        copy(vecmem::get_data(cells_events_host), cells_events_buffer);
        traccc::cell_module_collection_types::buffer module_events_buffer (module_events_host.size(), mr.main);
        copy(vecmem::get_data(module_events_host), module_events_buffer);

        spacepoints_cuda_buffer = ca_cuda(cells_events_buffer, module_events_buffer).first;
        stream.synchronize();
        
        if (run_cpu) {

            /*-----------------------------
                Clusterization (cpu)
            -----------------------------*/
            {
                traccc::performance::timer t("Clusterization  (cpu)",  elapsedTimes);
                measurements_per_event = ca(cells_events_host, module_events_host);
            // stop measuring clusterization cpu timer
            }
            /*---------------------------------
                Spacepoint formation (cpu)
            ---------------------------------*/
            {
                traccc::performance::timer t("Spacepoint formation  (cpu)", elapsedTimes);
                spacepoints_per_event = sf(measurements_per_event, module_events_host);
            // stop measuring spacepoint formation cpu timer
            }
        }
        // Statistics
        n_cells += cells_events_host.size();
        n_spacepoints += spacepoints_per_event.size();
        n_spacepoints_cuda += spacepoints_cuda_buffer.size();

    }

    std::cout << "==> Statistics ... " << std::endl;
    std::cout << "- read    " << n_spacepoints << " spacepoints from "
              << n_modules << " modules" << std::endl;
    std::cout << "- created        " << n_cells << " cells" << std::endl;
    std::cout << "- created        " << n_measurements << " measurements     "
              << std::endl;
    std::cout << "- created        " << n_spacepoints << " spacepoints     "
              << std::endl;
    std::cout << "- created (cuda) " << n_spacepoints_cuda
              << " spacepoints     " << std::endl;
    std::cout << "==>Elapsed times...\n" << elapsedTimes << std::endl;

    delete data;
    delete idr;
    delete myC;
    return 0;
}

// The main routine
//
int main(int argc, char* argv[]) {
    // Set up the program options
    po::options_description desc("Allowed options");

    // Add options
    desc.add_options()("help,h", "Give some help with the program's options");
    traccc::HTTSim_options common_opts(desc);
    traccc::full_tracking_input_config full_tracking_input_cfg(desc);
    desc.add_options()("run_cpu", po::value<bool>()->default_value(false),
                       "run cpu tracking as well");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    // Check errors
    traccc::handle_argument_errors(vm, desc);

    // Read options
    common_opts.read(vm);
    full_tracking_input_cfg.read(vm);
    auto run_cpu = vm["run_cpu"].as<bool>();

    std::cout << "Running " << argv[0] << " "
              << full_tracking_input_cfg.detector_file << " "
              << common_opts.input_directory << " " << common_opts.events
              << std::endl;

    return seq_run(full_tracking_input_cfg, common_opts, run_cpu);
}
