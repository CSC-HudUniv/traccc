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
#include "traccc/options/common_options.hpp"
#include "traccc/options/full_tracking_input_options.hpp"
#include "traccc/options/handle_argument_errors.hpp"
#include "traccc/performance/collection_comparator.hpp"
#include "traccc/performance/container_comparator.hpp"
#include "traccc/performance/timer.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"

// VecMem include(s).
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/host_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>

//Root Reader
#include <TROOT.h>
#include <TChain.h>
//#include "traccc/clusterization/IDreader.h"
#include "traccc/clusterization/DataList.h"

// System include(s).
#include <exception>
#include <iomanip>
#include <iostream>

namespace po = boost::program_options;

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

    traccc::clusterization_algorithm ca(host_mr);
    traccc::spacepoint_formation sf(host_mr);
    traccc::seeding_algorithm sa(host_mr);
    traccc::track_params_estimation tp(host_mr);

    traccc::cuda::stream stream;

    vecmem::cuda::async_copy copy{stream.cudaStream()};

    traccc::cuda::clusterization_algorithm ca_cuda( mr, copy, stream, common_opts.target_cells_per_partition);
    traccc::cuda::seeding_algorithm sa_cuda(mr, copy, stream);
    traccc::cuda::track_params_estimation tp_cuda(mr, copy, stream);

    // performance writer
    traccc::seeding_performance_writer sd_performance_writer(
        traccc::seeding_performance_writer::config{});
    if (i_cfg.check_performance) {
        sd_performance_writer.add_cache("CPU");
        sd_performance_writer.add_cache("CUDA");
    }

    traccc::performance::timing_info elapsedTimes;

    // Read from root file
    TChain *myC = new TChain("HTTEventTree");
    myC->Add((common_opts.input_directory.c_str()));
    DataList *data = new DataList(myC);
    
    //Gen number of entries and number of events from root file
    Long64_t nentries = data->fChain->GetEntriesFast();
    Int_t nfiles = (nevents < nentries ? nevents : nentries);


    for (Long64_t jentry=0; jentry<nfiles;jentry++) {
        
        
        traccc::clusterization_algorithm::output_type measurements_per_event;
        traccc::spacepoint_formation::output_type spacepoints_per_event;
        traccc::spacepoint_collection_types::buffer spacepoints_cuda_buffer(0, *mr.host);
        traccc::seed_collection_types::buffer seeds_cuda_buffer(0, *mr.host);

        //Load the current entry tree from data into ientry
        Long64_t ientry = data->LoadTree(jentry);
        if (ientry < 0) break;
        data->fChain->GetEntry(jentry); 
        
        //Create the event data holders (vectors)
        traccc::cell_collection_types::host cells_events_host;
        traccc::cell_module_collection_types::host module_events_host;


        for(int ihit=0; ihit<data->m_Hits_ ; ihit++){
            

            if(data->m_Hits_m_hitType[ihit]!=0) continue;
            if(data->m_Hits_m_detType[ihit]!=1) continue;

            // uint atlasid = data->m_Hits_m_identifierHash[ihit];
            // int ifpixel  = data->m_Hits_m_detType[ihit];
            // //unsigned long int geometry_id = idr->Give50muPixelDetectorID(data->m_Hits_m_identifierHash[ihit],ata->m_Hits_m_detType[ihit]);
            unsigned long int geometry_id = 576461302059245312;
            // int hit_id = ihit;
            // int channel0 = data->m_Hits_m_phiIndex[ihit];
            // int channel1 = data->m_Hits_m_etaIndex[ihit];
            // double timestamp = 0;
            // double value = data->m_Hits_m_ToT[ihit];

            traccc::cell current_cell {data->m_Hits_m_phiIndex[ihit], data->m_Hits_m_etaIndex[ihit], 0, geometry_id};
            traccc::cell_module current_module{};

            // std::cout << "CELL CHANNEL_0 : " << current_cell.channel0 << "\n";
            // std::cout << "CELL CHANNEL_1 : " << current_cell.channel1 << "\n";
            // std::cout << "CELL TIMESTAMP : " << current_cell.time << "\n";
            // std::cout << "CELL GEOM_ID : " << current_cell.module_link << "\n";

            cells_events_host.push_back(current_cell);
            module_events_host.push_back(current_module);

            }

            //Allocate GPU memory and copy the cell and module data from the host type to buffer type
            traccc::cell_collection_types::buffer cells_events_buffer (cells_events_host.size(), mr.main);
            copy(vecmem::get_data(cells_events_host), cells_events_buffer);
            traccc::cell_module_collection_types::buffer module_events_buffer (module_events_host.size(), mr.main);
            copy(vecmem::get_data(module_events_host), module_events_buffer);

            spacepoints_cuda_buffer = ca_cuda(cells_events_buffer, module_events_buffer).first;

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

            //seeds_cuda_buffer = sa_cuda(spacepoints_cuda_buffer);

            //std::cout << vecmem::get_data(cells_events_buffer);


            // Statistics
            n_cells += cells_events_host.size();
            n_spacepoints += spacepoints_per_event.size();
            n_spacepoints_cuda += spacepoints_cuda_buffer.size();

    }


    // // Loop over events
    // for (unsigned int event = common_opts.skip;
    //      event < common_opts.events + common_opts.skip; ++event) {

    //     // Instantiate host containers/collections
    //     traccc::io::cell_reader_output alt_read_out_per_event;
    //     traccc::clusterization_algorithm::output_type measurements_per_event;
    //     traccc::spacepoint_formation::output_type spacepoints_per_event;
    //     traccc::seeding_algorithm::output_type seeds;
    //     traccc::track_params_estimation::output_type params;

    //     // Instantiate cuda containers/collections
    //     traccc::spacepoint_collection_types::buffer spacepoints_cuda_buffer(
    //         0, *mr.host);
    //     traccc::seed_collection_types::buffer seeds_cuda_buffer(0, *mr.host);
    //     traccc::bound_track_parameters_collection_types::buffer
    //         params_cuda_buffer(0, *mr.host);

    //     {
    //         traccc::performance::timer wall_t("Wall time", elapsedTimes);

    //         {
    //             traccc::performance::timer t("File reading  (cpu)",
    //                                          elapsedTimes);
    //             // Read the cells from the relevant event file into host memory.
    //             alt_read_out_per_event = traccc::io::read_cells(
    //                 event, common_opts.input_directory,
    //                 common_opts.input_data_format, &surface_transforms,
    //                 &digi_cfg, &cuda_host_mr);

    //         }  // stop measuring file reading timer

    //         const traccc::cell_collection_types::host& cells_per_event =
    //             alt_read_out_per_event.cells;
    //         const traccc::cell_module_collection_types::host&
    //             modules_per_event = alt_read_out_per_event.modules;

    //         /*-----------------------------
    //             Clusterization and Spacepoint Creation (cuda)
    //         -----------------------------*/
    //         // Create device copy of input collections
    //         traccc::cell_collection_types::buffer cells_buffer(
    //             cells_per_event.size(), mr.main);
    //         copy(vecmem::get_data(cells_per_event), cells_buffer);
    //         traccc::cell_module_collection_types::buffer modules_buffer(
    //             modules_per_event.size(), mr.main);
    //         copy(vecmem::get_data(modules_per_event), modules_buffer);

    //         {
    //             traccc::performance::timer t("Clusterization (cuda)",
    //                                          elapsedTimes);
    //             // Reconstruct it into spacepoints on the device.
    //             spacepoints_cuda_buffer =
    //                 ca_cuda(cells_buffer, modules_buffer).first;
    //             stream.synchronize();
    //         }  // stop measuring clusterization cuda timer

    //         if (run_cpu) {

    //             /*-----------------------------
    //                 Clusterization (cpu)
    //             -----------------------------*/

    //             {
    //                 traccc::performance::timer t("Clusterization  (cpu)",
    //                                              elapsedTimes);
    //                 measurements_per_event =
    //                     ca(cells_per_event, modules_per_event);
    //             }  // stop measuring clusterization cpu timer

    //             /*---------------------------------
    //                 Spacepoint formation (cpu)
    //             ---------------------------------*/

    //             {
    //                 traccc::performance::timer t("Spacepoint formation  (cpu)",
    //                                              elapsedTimes);
    //                 spacepoints_per_event =
    //                     sf(measurements_per_event, modules_per_event);
    //             }  // stop measuring spacepoint formation cpu timer
    //         }

    //         /*----------------------------
    //             Seeding algorithm
    //         ----------------------------*/

    //         // CUDA

    //         {
    //             traccc::performance::timer t("Seeding (cuda)", elapsedTimes);
    //             seeds_cuda_buffer = sa_cuda(spacepoints_cuda_buffer);
    //         }  // stop measuring seeding cuda timer

    //         // CPU

    //         if (run_cpu) {
    //             traccc::performance::timer t("Seeding  (cpu)", elapsedTimes);
    //             seeds = sa(spacepoints_per_event);
    //         }  // stop measuring seeding cpu timer

    //         /*----------------------------
    //         Track params estimation
    //         ----------------------------*/

    //         // CUDA

    //         {
    //             traccc::performance::timer t("Track params (cuda)",
    //                                          elapsedTimes);
    //             params_cuda_buffer =
    //                 tp_cuda(spacepoints_cuda_buffer, seeds_cuda_buffer);
    //         }  // stop measuring track params timer

    //         // CPU

    //         if (run_cpu) {
    //             traccc::performance::timer t("Track params  (cpu)",
    //                                          elapsedTimes);
    //             params = tp(spacepoints_per_event, seeds);
    //         }  // stop measuring track params cpu timer

    //     }  // Stop measuring wall time

    //     /*----------------------------------
    //       compare cpu and cuda result
    //       ----------------------------------*/

    //     traccc::spacepoint_collection_types::host spacepoints_per_event_cuda;
    //     traccc::seed_collection_types::host seeds_cuda;
    //     traccc::bound_track_parameters_collection_types::host params_cuda;
    //     if (run_cpu || i_cfg.check_performance) {
    //         copy(spacepoints_cuda_buffer, spacepoints_per_event_cuda)->wait();
    //         copy(seeds_cuda_buffer, seeds_cuda)->wait();
    //         copy(params_cuda_buffer, params_cuda)->wait();
    //     }

    //     if (run_cpu) {

    //         // Show which event we are currently presenting the results for.
    //         std::cout << "===>>> Event " << event << " <<<===" << std::endl;

    //         // Compare the spacepoints made on the host and on the device.
    //         traccc::collection_comparator<traccc::spacepoint>
    //             compare_spacepoints{"spacepoints"};
    //         compare_spacepoints(vecmem::get_data(spacepoints_per_event),
    //                             vecmem::get_data(spacepoints_per_event_cuda));

    //         // Compare the seeds made on the host and on the device
    //         traccc::collection_comparator<traccc::seed> compare_seeds{
    //             "seeds", traccc::details::comparator_factory<traccc::seed>{
    //                          vecmem::get_data(spacepoints_per_event),
    //                          vecmem::get_data(spacepoints_per_event_cuda)}};
    //         compare_seeds(vecmem::get_data(seeds),
    //                       vecmem::get_data(seeds_cuda));

    //         // Compare the track parameters made on the host and on the device.
    //         traccc::collection_comparator<traccc::bound_track_parameters>
    //             compare_track_parameters{"track parameters"};
    //         compare_track_parameters(vecmem::get_data(params),
    //                                  vecmem::get_data(params_cuda));

    //         /// Statistics
    //         n_modules += alt_read_out_per_event.modules.size();
    //         n_cells += alt_read_out_per_event.cells.size();
    //         n_measurements += measurements_per_event.size();
    //         n_spacepoints += spacepoints_per_event.size();
    //         n_spacepoints_cuda += spacepoints_per_event_cuda.size();
    //         n_seeds_cuda += seeds_cuda.size();
    //         n_seeds += seeds.size();
    //     }
    //     if (i_cfg.check_performance) {

    //         traccc::event_map evt_map(
    //             event, i_cfg.detector_file, i_cfg.digitization_config_file,
    //             common_opts.input_directory, common_opts.input_directory,
    //             common_opts.input_directory, host_mr);
    //         sd_performance_writer.write(
    //             "CUDA", vecmem::get_data(seeds_cuda),
    //             vecmem::get_data(spacepoints_per_event_cuda), evt_map);

    //         if (run_cpu) {
    //             sd_performance_writer.write(
    //                 "CPU", vecmem::get_data(seeds),
    //                 vecmem::get_data(spacepoints_per_event), evt_map);
    //         }
    //     }
    // }

    // if (i_cfg.check_performance) {
    //     sd_performance_writer.finalize();
    // }

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

    std::cout << "- created  (cpu) " << n_seeds << " seeds" << std::endl;
    std::cout << "- created (cuda) " << n_seeds_cuda << " seeds" << std::endl;
    std::cout << "==>Elapsed times...\n" << elapsedTimes << std::endl;


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
