//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 15 11:02:01 2023 by ROOT version 6.24/06
// from TTree HTTEventTree/data
// found on file: singlemu_invPtFlat1_10k_wrap.root
//////////////////////////////////////////////////////////

#ifndef DataList_h
#define DataList_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>

// Header file for the classes stored in the TTree if any.
#include "TObject.h"

using namespace std;

class DataList {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxm_optional_m_OfflineClusters = 1000000;
   static constexpr Int_t kMaxm_optional_m_OfflineTracks = 10000;
   static constexpr Int_t kMaxm_optional_m_TruthTracks = 100000;
   static constexpr Int_t kMaxm_Hits = 2000000;

   // Declaration of leaf types
 //HTTEventInputHeader *HTTEventInputHeader;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   UInt_t          m_event_fUniqueID;
   UInt_t          m_event_fBits;
   ULong_t         m_event_m_run_number;
   ULong_t         m_event_m_event_number;
   Float_t         m_event_m_averageInteractionsPerCrossing;
   Float_t         m_event_m_actualInteractionsPerCrossing;
   Int_t           m_event_m_LB;
   Int_t           m_event_m_BCID;
   UInt_t          m_event_m_extendedLevel1ID;
   UInt_t          m_event_m_level1TriggerType;
   vector<unsigned int> m_event_m_level1TriggerInfo;
   UInt_t          m_optional_fUniqueID;
   UInt_t          m_optional_fBits;
   Int_t           m_optional_m_OfflineClusters_;
   UInt_t          m_optional_m_OfflineClusters_fUniqueID[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   UInt_t          m_optional_m_OfflineClusters_fBits[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
 //vector<HTTHit>  m_optional_m_OfflineClusters_m_hitlist[kMaxm_optional_m_OfflineClusters];
   Int_t           m_optional_m_OfflineClusters_m_clusterEquiv_m_hitType[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   Int_t           m_optional_m_OfflineClusters_m_clusterEquiv_m_detectorZone[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   Int_t           m_optional_m_OfflineClusters_m_clusterEquiv_m_detType[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   UInt_t          m_optional_m_OfflineClusters_m_clusterEquiv_m_identifierHash[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   UInt_t          m_optional_m_OfflineClusters_m_clusterEquiv_m_layer_disk[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   UInt_t          m_optional_m_OfflineClusters_m_clusterEquiv_m_side[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   Int_t           m_optional_m_OfflineClusters_m_clusterEquiv_m_etaModule[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   UInt_t          m_optional_m_OfflineClusters_m_clusterEquiv_m_phiModule[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   UInt_t          m_optional_m_OfflineClusters_m_clusterEquiv_m_etaWidth[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   UInt_t          m_optional_m_OfflineClusters_m_clusterEquiv_m_phiWidth[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   Int_t           m_optional_m_OfflineClusters_m_clusterEquiv_m_layer[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   Int_t           m_optional_m_OfflineClusters_m_clusterEquiv_m_section[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   Int_t           m_optional_m_OfflineClusters_m_clusterEquiv_m_phiIndex[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   Int_t           m_optional_m_OfflineClusters_m_clusterEquiv_m_etaIndex[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   Float_t         m_optional_m_OfflineClusters_m_clusterEquiv_m_x[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   Float_t         m_optional_m_OfflineClusters_m_clusterEquiv_m_y[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   Float_t         m_optional_m_OfflineClusters_m_clusterEquiv_m_z[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   UInt_t          m_optional_m_OfflineClusters_m_clusterEquiv_m_hw_word[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   UInt_t          m_optional_m_OfflineClusters_m_clusterEquiv_m_ToT[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   Long_t          m_optional_m_OfflineClusters_m_clusterEquiv_m_eventindex[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   Long_t          m_optional_m_OfflineClusters_m_clusterEquiv_m_barcode[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   Float_t         m_optional_m_OfflineClusters_m_clusterEquiv_m_barcode_pt[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   ULong_t         m_optional_m_OfflineClusters_m_clusterEquiv_m_parentage_mask[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   UInt_t          m_optional_m_OfflineClusters_m_clusterEquiv_m_truth_fUniqueID[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
   UInt_t          m_optional_m_OfflineClusters_m_clusterEquiv_m_truth_fBits[kMaxm_optional_m_OfflineClusters];   //[m_optional.m_OfflineClusters_]
 //map<pair<long,long>,float> m_optional_m_OfflineClusters_m_clusterEquiv_m_truth_m_truth[kMaxm_optional_m_OfflineClusters];
   Int_t           m_optional_m_OfflineTracks_;
   UInt_t          m_optional_m_OfflineTracks_fUniqueID[kMaxm_optional_m_OfflineTracks];   //[m_optional.m_OfflineTracks_]
   UInt_t          m_optional_m_OfflineTracks_fBits[kMaxm_optional_m_OfflineTracks];   //[m_optional.m_OfflineTracks_]
   Double_t        m_optional_m_OfflineTracks_m_qoverpt[kMaxm_optional_m_OfflineTracks];   //[m_optional.m_OfflineTracks_]
   Double_t        m_optional_m_OfflineTracks_m_eta[kMaxm_optional_m_OfflineTracks];   //[m_optional.m_OfflineTracks_]
   Double_t        m_optional_m_OfflineTracks_m_phi[kMaxm_optional_m_OfflineTracks];   //[m_optional.m_OfflineTracks_]
   Double_t        m_optional_m_OfflineTracks_m_d0[kMaxm_optional_m_OfflineTracks];   //[m_optional.m_OfflineTracks_]
   Double_t        m_optional_m_OfflineTracks_m_z0[kMaxm_optional_m_OfflineTracks];   //[m_optional.m_OfflineTracks_]
   Long_t          m_optional_m_OfflineTracks_m_barcode[kMaxm_optional_m_OfflineTracks];   //[m_optional.m_OfflineTracks_]
   Double_t        m_optional_m_OfflineTracks_m_barcode_frac[kMaxm_optional_m_OfflineTracks];   //[m_optional.m_OfflineTracks_]
 //vector<HTTOfflineHit> m_optional_m_OfflineTracks_m_hits[kMaxm_optional_m_OfflineTracks];
   Int_t           m_optional_m_TruthTracks_;
   UInt_t          m_optional_m_TruthTracks_fUniqueID[kMaxm_optional_m_TruthTracks];   //[m_optional.m_TruthTracks_]
   UInt_t          m_optional_m_TruthTracks_fBits[kMaxm_optional_m_TruthTracks];   //[m_optional.m_TruthTracks_]
   Double_t        m_optional_m_TruthTracks_m_d0[kMaxm_optional_m_TruthTracks];   //[m_optional.m_TruthTracks_]
   Double_t        m_optional_m_TruthTracks_m_z0[kMaxm_optional_m_TruthTracks];   //[m_optional.m_TruthTracks_]
   Double_t        m_optional_m_TruthTracks_m_vtx_x[kMaxm_optional_m_TruthTracks];   //[m_optional.m_TruthTracks_]
   Double_t        m_optional_m_TruthTracks_m_vtx_y[kMaxm_optional_m_TruthTracks];   //[m_optional.m_TruthTracks_]
   Double_t        m_optional_m_TruthTracks_m_vtx_z[kMaxm_optional_m_TruthTracks];   //[m_optional.m_TruthTracks_]
   Double_t        m_optional_m_TruthTracks_m_px[kMaxm_optional_m_TruthTracks];   //[m_optional.m_TruthTracks_]
   Double_t        m_optional_m_TruthTracks_m_py[kMaxm_optional_m_TruthTracks];   //[m_optional.m_TruthTracks_]
   Double_t        m_optional_m_TruthTracks_m_pz[kMaxm_optional_m_TruthTracks];   //[m_optional.m_TruthTracks_]
   Int_t           m_optional_m_TruthTracks_m_q[kMaxm_optional_m_TruthTracks];   //[m_optional.m_TruthTracks_]
   Int_t           m_optional_m_TruthTracks_m_pdgcode[kMaxm_optional_m_TruthTracks];   //[m_optional.m_TruthTracks_]
   Int_t           m_optional_m_TruthTracks_m_barcode[kMaxm_optional_m_TruthTracks];   //[m_optional.m_TruthTracks_]
   Int_t           m_optional_m_TruthTracks_m_evtindex[kMaxm_optional_m_TruthTracks];   //[m_optional.m_TruthTracks_]
   Double_t        m_optional_m_TruthTracks_m_barcode_frac_offline[kMaxm_optional_m_TruthTracks];   //[m_optional.m_TruthTracks_]
   Int_t           m_Hits_;
   Int_t           m_Hits_m_hitType[kMaxm_Hits];   //[m_Hits_]
   Int_t           m_Hits_m_detectorZone[kMaxm_Hits];   //[m_Hits_]
   Int_t           m_Hits_m_detType[kMaxm_Hits];   //[m_Hits_]
   UInt_t          m_Hits_m_identifierHash[kMaxm_Hits];   //[m_Hits_]
   UInt_t          m_Hits_m_layer_disk[kMaxm_Hits];   //[m_Hits_]
   UInt_t          m_Hits_m_side[kMaxm_Hits];   //[m_Hits_]
   Int_t           m_Hits_m_etaModule[kMaxm_Hits];   //[m_Hits_]
   UInt_t          m_Hits_m_phiModule[kMaxm_Hits];   //[m_Hits_]
   UInt_t          m_Hits_m_etaWidth[kMaxm_Hits];   //[m_Hits_]
   UInt_t          m_Hits_m_phiWidth[kMaxm_Hits];   //[m_Hits_]
   Int_t           m_Hits_m_layer[kMaxm_Hits];   //[m_Hits_]
   Int_t           m_Hits_m_section[kMaxm_Hits];   //[m_Hits_]
   Int_t           m_Hits_m_phiIndex[kMaxm_Hits];   //[m_Hits_]
   Int_t           m_Hits_m_etaIndex[kMaxm_Hits];   //[m_Hits_]
   Float_t         m_Hits_m_x[kMaxm_Hits];   //[m_Hits_]
   Float_t         m_Hits_m_y[kMaxm_Hits];   //[m_Hits_]
   Float_t         m_Hits_m_z[kMaxm_Hits];   //[m_Hits_]
   UInt_t          m_Hits_m_hw_word[kMaxm_Hits];   //[m_Hits_]
   UInt_t          m_Hits_m_ToT[kMaxm_Hits];   //[m_Hits_]
   Long_t          m_Hits_m_eventindex[kMaxm_Hits];   //[m_Hits_]
   Long_t          m_Hits_m_barcode[kMaxm_Hits];   //[m_Hits_]
   Float_t         m_Hits_m_barcode_pt[kMaxm_Hits];   //[m_Hits_]
   ULong_t         m_Hits_m_parentage_mask[kMaxm_Hits];   //[m_Hits_]
   UInt_t          m_Hits_m_truth_fUniqueID[kMaxm_Hits];   //[m_Hits_]
   UInt_t          m_Hits_m_truth_fBits[kMaxm_Hits];   //[m_Hits_]
 //map<pair<long,long>,float> m_Hits_m_truth_m_truth[kMaxm_Hits];

   // List of branches
   TBranch        *b_HTTEventInputHeader_fUniqueID;   //!
   TBranch        *b_HTTEventInputHeader_fBits;   //!
   TBranch        *b_HTTEventInputHeader_m_event_fUniqueID;   //!
   TBranch        *b_HTTEventInputHeader_m_event_fBits;   //!
   TBranch        *b_HTTEventInputHeader_m_event_m_run_number;   //!
   TBranch        *b_HTTEventInputHeader_m_event_m_event_number;   //!
   TBranch        *b_HTTEventInputHeader_m_event_m_averageInteractionsPerCrossing;   //!
   TBranch        *b_HTTEventInputHeader_m_event_m_actualInteractionsPerCrossing;   //!
   TBranch        *b_HTTEventInputHeader_m_event_m_LB;   //!
   TBranch        *b_HTTEventInputHeader_m_event_m_BCID;   //!
   TBranch        *b_HTTEventInputHeader_m_event_m_extendedLevel1ID;   //!
   TBranch        *b_HTTEventInputHeader_m_event_m_level1TriggerType;   //!
   TBranch        *b_HTTEventInputHeader_m_event_m_level1TriggerInfo;   //!
   TBranch        *b_HTTEventInputHeader_m_optional_fUniqueID;   //!
   TBranch        *b_HTTEventInputHeader_m_optional_fBits;   //!
   TBranch        *b_HTTEventInputHeader_m_optional_m_OfflineClusters_;   //!
   TBranch        *b_m_optional_m_OfflineClusters_fUniqueID;   //!
   TBranch        *b_m_optional_m_OfflineClusters_fBits;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_hitType;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_detectorZone;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_detType;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_identifierHash;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_layer_disk;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_side;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_etaModule;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_phiModule;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_etaWidth;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_phiWidth;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_layer;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_section;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_phiIndex;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_etaIndex;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_x;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_y;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_z;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_hw_word;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_ToT;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_eventindex;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_barcode;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_barcode_pt;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_parentage_mask;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_truth_fUniqueID;   //!
   TBranch        *b_m_optional_m_OfflineClusters_m_clusterEquiv_m_truth_fBits;   //!
   TBranch        *b_HTTEventInputHeader_m_optional_m_OfflineTracks_;   //!
   TBranch        *b_m_optional_m_OfflineTracks_fUniqueID;   //!
   TBranch        *b_m_optional_m_OfflineTracks_fBits;   //!
   TBranch        *b_m_optional_m_OfflineTracks_m_qoverpt;   //!
   TBranch        *b_m_optional_m_OfflineTracks_m_eta;   //!
   TBranch        *b_m_optional_m_OfflineTracks_m_phi;   //!
   TBranch        *b_m_optional_m_OfflineTracks_m_d0;   //!
   TBranch        *b_m_optional_m_OfflineTracks_m_z0;   //!
   TBranch        *b_m_optional_m_OfflineTracks_m_barcode;   //!
   TBranch        *b_m_optional_m_OfflineTracks_m_barcode_frac;   //!
   TBranch        *b_HTTEventInputHeader_m_optional_m_TruthTracks_;   //!
   TBranch        *b_m_optional_m_TruthTracks_fUniqueID;   //!
   TBranch        *b_m_optional_m_TruthTracks_fBits;   //!
   TBranch        *b_m_optional_m_TruthTracks_m_d0;   //!
   TBranch        *b_m_optional_m_TruthTracks_m_z0;   //!
   TBranch        *b_m_optional_m_TruthTracks_m_vtx_x;   //!
   TBranch        *b_m_optional_m_TruthTracks_m_vtx_y;   //!
   TBranch        *b_m_optional_m_TruthTracks_m_vtx_z;   //!
   TBranch        *b_m_optional_m_TruthTracks_m_px;   //!
   TBranch        *b_m_optional_m_TruthTracks_m_py;   //!
   TBranch        *b_m_optional_m_TruthTracks_m_pz;   //!
   TBranch        *b_m_optional_m_TruthTracks_m_q;   //!
   TBranch        *b_m_optional_m_TruthTracks_m_pdgcode;   //!
   TBranch        *b_m_optional_m_TruthTracks_m_barcode;   //!
   TBranch        *b_m_optional_m_TruthTracks_m_evtindex;   //!
   TBranch        *b_m_optional_m_TruthTracks_m_barcode_frac_offline;   //!
   TBranch        *b_HTTEventInputHeader_m_Hits_;   //!
   TBranch        *b_m_Hits_m_hitType;   //!
   TBranch        *b_m_Hits_m_detectorZone;   //!
   TBranch        *b_m_Hits_m_detType;   //!
   TBranch        *b_m_Hits_m_identifierHash;   //!
   TBranch        *b_m_Hits_m_layer_disk;   //!
   TBranch        *b_m_Hits_m_side;   //!
   TBranch        *b_m_Hits_m_etaModule;   //!
   TBranch        *b_m_Hits_m_phiModule;   //!
   TBranch        *b_m_Hits_m_etaWidth;   //!
   TBranch        *b_m_Hits_m_phiWidth;   //!
   TBranch        *b_m_Hits_m_layer;   //!
   TBranch        *b_m_Hits_m_section;   //!
   TBranch        *b_m_Hits_m_phiIndex;   //!
   TBranch        *b_m_Hits_m_etaIndex;   //!
   TBranch        *b_m_Hits_m_x;   //!
   TBranch        *b_m_Hits_m_y;   //!
   TBranch        *b_m_Hits_m_z;   //!
   TBranch        *b_m_Hits_m_hw_word;   //!
   TBranch        *b_m_Hits_m_ToT;   //!
   TBranch        *b_m_Hits_m_eventindex;   //!
   TBranch        *b_m_Hits_m_barcode;   //!
   TBranch        *b_m_Hits_m_barcode_pt;   //!
   TBranch        *b_m_Hits_m_parentage_mask;   //!
   TBranch        *b_m_Hits_m_truth_fUniqueID;   //!
   TBranch        *b_m_Hits_m_truth_fBits;   //!

   DataList(TTree *tree=0);
   virtual ~DataList();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
  virtual void     WriteCellFiles(Int_t nevents,string dirnm);
  virtual void     WritePixelCellFiles(Int_t nevents,string dirnm);
  void WriteClusters(Int_t nevents, string dirnm);
  void WritePixelClusters(Int_t nevents,string dirnm);
  void DrawVariables();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

  
};

#endif

#ifdef DataList_cxx
DataList::DataList(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("singlemu_invPtFlat1_10k_wrap.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("singlemu_invPtFlat1_10k_wrap.root");
      }
      f->GetObject("HTTEventTree",tree);

   }
   Init(tree);
}

DataList::~DataList()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DataList::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DataList::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void DataList::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_HTTEventInputHeader_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_HTTEventInputHeader_fBits);
   fChain->SetBranchAddress("m_event.fUniqueID", &m_event_fUniqueID, &b_HTTEventInputHeader_m_event_fUniqueID);
   fChain->SetBranchAddress("m_event.fBits", &m_event_fBits, &b_HTTEventInputHeader_m_event_fBits);
   fChain->SetBranchAddress("m_event.m_run_number", &m_event_m_run_number, &b_HTTEventInputHeader_m_event_m_run_number);
   fChain->SetBranchAddress("m_event.m_event_number", &m_event_m_event_number, &b_HTTEventInputHeader_m_event_m_event_number);
   fChain->SetBranchAddress("m_event.m_averageInteractionsPerCrossing", &m_event_m_averageInteractionsPerCrossing, &b_HTTEventInputHeader_m_event_m_averageInteractionsPerCrossing);
   fChain->SetBranchAddress("m_event.m_actualInteractionsPerCrossing", &m_event_m_actualInteractionsPerCrossing, &b_HTTEventInputHeader_m_event_m_actualInteractionsPerCrossing);
   fChain->SetBranchAddress("m_event.m_LB", &m_event_m_LB, &b_HTTEventInputHeader_m_event_m_LB);
   fChain->SetBranchAddress("m_event.m_BCID", &m_event_m_BCID, &b_HTTEventInputHeader_m_event_m_BCID);
   fChain->SetBranchAddress("m_event.m_extendedLevel1ID", &m_event_m_extendedLevel1ID, &b_HTTEventInputHeader_m_event_m_extendedLevel1ID);
   fChain->SetBranchAddress("m_event.m_level1TriggerType", &m_event_m_level1TriggerType, &b_HTTEventInputHeader_m_event_m_level1TriggerType);
   fChain->SetBranchAddress("m_event.m_level1TriggerInfo", &m_event_m_level1TriggerInfo, &b_HTTEventInputHeader_m_event_m_level1TriggerInfo);
   fChain->SetBranchAddress("m_optional.fUniqueID", &m_optional_fUniqueID, &b_HTTEventInputHeader_m_optional_fUniqueID);
   fChain->SetBranchAddress("m_optional.fBits", &m_optional_fBits, &b_HTTEventInputHeader_m_optional_fBits);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters", &m_optional_m_OfflineClusters_, &b_HTTEventInputHeader_m_optional_m_OfflineClusters_);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.fUniqueID", m_optional_m_OfflineClusters_fUniqueID, &b_m_optional_m_OfflineClusters_fUniqueID);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.fBits", m_optional_m_OfflineClusters_fBits, &b_m_optional_m_OfflineClusters_fBits);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_hitType", m_optional_m_OfflineClusters_m_clusterEquiv_m_hitType, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_hitType);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_detectorZone", m_optional_m_OfflineClusters_m_clusterEquiv_m_detectorZone, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_detectorZone);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_detType", m_optional_m_OfflineClusters_m_clusterEquiv_m_detType, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_detType);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_identifierHash", m_optional_m_OfflineClusters_m_clusterEquiv_m_identifierHash, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_identifierHash);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_layer_disk", m_optional_m_OfflineClusters_m_clusterEquiv_m_layer_disk, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_layer_disk);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_side", m_optional_m_OfflineClusters_m_clusterEquiv_m_side, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_side);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_etaModule", m_optional_m_OfflineClusters_m_clusterEquiv_m_etaModule, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_etaModule);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_phiModule", m_optional_m_OfflineClusters_m_clusterEquiv_m_phiModule, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_phiModule);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_etaWidth", m_optional_m_OfflineClusters_m_clusterEquiv_m_etaWidth, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_etaWidth);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_phiWidth", m_optional_m_OfflineClusters_m_clusterEquiv_m_phiWidth, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_phiWidth);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_layer", m_optional_m_OfflineClusters_m_clusterEquiv_m_layer, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_layer);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_section", m_optional_m_OfflineClusters_m_clusterEquiv_m_section, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_section);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_phiIndex", m_optional_m_OfflineClusters_m_clusterEquiv_m_phiIndex, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_phiIndex);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_etaIndex", m_optional_m_OfflineClusters_m_clusterEquiv_m_etaIndex, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_etaIndex);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_x", m_optional_m_OfflineClusters_m_clusterEquiv_m_x, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_x);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_y", m_optional_m_OfflineClusters_m_clusterEquiv_m_y, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_y);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_z", m_optional_m_OfflineClusters_m_clusterEquiv_m_z, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_z);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_hw_word", m_optional_m_OfflineClusters_m_clusterEquiv_m_hw_word, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_hw_word);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_ToT", m_optional_m_OfflineClusters_m_clusterEquiv_m_ToT, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_ToT);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_eventindex", m_optional_m_OfflineClusters_m_clusterEquiv_m_eventindex, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_eventindex);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_barcode", m_optional_m_OfflineClusters_m_clusterEquiv_m_barcode, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_barcode);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_barcode_pt", m_optional_m_OfflineClusters_m_clusterEquiv_m_barcode_pt, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_barcode_pt);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_parentage_mask", m_optional_m_OfflineClusters_m_clusterEquiv_m_parentage_mask, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_parentage_mask);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_truth.fUniqueID", m_optional_m_OfflineClusters_m_clusterEquiv_m_truth_fUniqueID, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_truth_fUniqueID);
   fChain->SetBranchAddress("m_optional.m_OfflineClusters.m_clusterEquiv.m_truth.fBits", m_optional_m_OfflineClusters_m_clusterEquiv_m_truth_fBits, &b_m_optional_m_OfflineClusters_m_clusterEquiv_m_truth_fBits);
   fChain->SetBranchAddress("m_optional.m_OfflineTracks", &m_optional_m_OfflineTracks_, &b_HTTEventInputHeader_m_optional_m_OfflineTracks_);
   fChain->SetBranchAddress("m_optional.m_OfflineTracks.fUniqueID", m_optional_m_OfflineTracks_fUniqueID, &b_m_optional_m_OfflineTracks_fUniqueID);
   fChain->SetBranchAddress("m_optional.m_OfflineTracks.fBits", m_optional_m_OfflineTracks_fBits, &b_m_optional_m_OfflineTracks_fBits);
   fChain->SetBranchAddress("m_optional.m_OfflineTracks.m_qoverpt", m_optional_m_OfflineTracks_m_qoverpt, &b_m_optional_m_OfflineTracks_m_qoverpt);
   fChain->SetBranchAddress("m_optional.m_OfflineTracks.m_eta", m_optional_m_OfflineTracks_m_eta, &b_m_optional_m_OfflineTracks_m_eta);
   fChain->SetBranchAddress("m_optional.m_OfflineTracks.m_phi", m_optional_m_OfflineTracks_m_phi, &b_m_optional_m_OfflineTracks_m_phi);
   fChain->SetBranchAddress("m_optional.m_OfflineTracks.m_d0", m_optional_m_OfflineTracks_m_d0, &b_m_optional_m_OfflineTracks_m_d0);
   fChain->SetBranchAddress("m_optional.m_OfflineTracks.m_z0", m_optional_m_OfflineTracks_m_z0, &b_m_optional_m_OfflineTracks_m_z0);
   fChain->SetBranchAddress("m_optional.m_OfflineTracks.m_barcode", m_optional_m_OfflineTracks_m_barcode, &b_m_optional_m_OfflineTracks_m_barcode);
   fChain->SetBranchAddress("m_optional.m_OfflineTracks.m_barcode_frac", m_optional_m_OfflineTracks_m_barcode_frac, &b_m_optional_m_OfflineTracks_m_barcode_frac);
   fChain->SetBranchAddress("m_optional.m_TruthTracks", &m_optional_m_TruthTracks_, &b_HTTEventInputHeader_m_optional_m_TruthTracks_);
   fChain->SetBranchAddress("m_optional.m_TruthTracks.fUniqueID", m_optional_m_TruthTracks_fUniqueID, &b_m_optional_m_TruthTracks_fUniqueID);
   fChain->SetBranchAddress("m_optional.m_TruthTracks.fBits", m_optional_m_TruthTracks_fBits, &b_m_optional_m_TruthTracks_fBits);
   fChain->SetBranchAddress("m_optional.m_TruthTracks.m_d0", m_optional_m_TruthTracks_m_d0, &b_m_optional_m_TruthTracks_m_d0);
   fChain->SetBranchAddress("m_optional.m_TruthTracks.m_z0", m_optional_m_TruthTracks_m_z0, &b_m_optional_m_TruthTracks_m_z0);
   fChain->SetBranchAddress("m_optional.m_TruthTracks.m_vtx_x", m_optional_m_TruthTracks_m_vtx_x, &b_m_optional_m_TruthTracks_m_vtx_x);
   fChain->SetBranchAddress("m_optional.m_TruthTracks.m_vtx_y", m_optional_m_TruthTracks_m_vtx_y, &b_m_optional_m_TruthTracks_m_vtx_y);
   fChain->SetBranchAddress("m_optional.m_TruthTracks.m_vtx_z", m_optional_m_TruthTracks_m_vtx_z, &b_m_optional_m_TruthTracks_m_vtx_z);
   fChain->SetBranchAddress("m_optional.m_TruthTracks.m_px", m_optional_m_TruthTracks_m_px, &b_m_optional_m_TruthTracks_m_px);
   fChain->SetBranchAddress("m_optional.m_TruthTracks.m_py", m_optional_m_TruthTracks_m_py, &b_m_optional_m_TruthTracks_m_py);
   fChain->SetBranchAddress("m_optional.m_TruthTracks.m_pz", m_optional_m_TruthTracks_m_pz, &b_m_optional_m_TruthTracks_m_pz);
   fChain->SetBranchAddress("m_optional.m_TruthTracks.m_q", m_optional_m_TruthTracks_m_q, &b_m_optional_m_TruthTracks_m_q);
   fChain->SetBranchAddress("m_optional.m_TruthTracks.m_pdgcode", m_optional_m_TruthTracks_m_pdgcode, &b_m_optional_m_TruthTracks_m_pdgcode);
   fChain->SetBranchAddress("m_optional.m_TruthTracks.m_barcode", m_optional_m_TruthTracks_m_barcode, &b_m_optional_m_TruthTracks_m_barcode);
   fChain->SetBranchAddress("m_optional.m_TruthTracks.m_evtindex", m_optional_m_TruthTracks_m_evtindex, &b_m_optional_m_TruthTracks_m_evtindex);
   fChain->SetBranchAddress("m_optional.m_TruthTracks.m_barcode_frac_offline", m_optional_m_TruthTracks_m_barcode_frac_offline, &b_m_optional_m_TruthTracks_m_barcode_frac_offline);
   fChain->SetBranchAddress("m_Hits", &m_Hits_, &b_HTTEventInputHeader_m_Hits_);
   fChain->SetBranchAddress("m_Hits.m_hitType", m_Hits_m_hitType, &b_m_Hits_m_hitType);
   fChain->SetBranchAddress("m_Hits.m_detectorZone", m_Hits_m_detectorZone, &b_m_Hits_m_detectorZone);
   fChain->SetBranchAddress("m_Hits.m_detType", m_Hits_m_detType, &b_m_Hits_m_detType);
   fChain->SetBranchAddress("m_Hits.m_identifierHash", m_Hits_m_identifierHash, &b_m_Hits_m_identifierHash);
   fChain->SetBranchAddress("m_Hits.m_layer_disk", m_Hits_m_layer_disk, &b_m_Hits_m_layer_disk);
   fChain->SetBranchAddress("m_Hits.m_side", m_Hits_m_side, &b_m_Hits_m_side);
   fChain->SetBranchAddress("m_Hits.m_etaModule", m_Hits_m_etaModule, &b_m_Hits_m_etaModule);
   fChain->SetBranchAddress("m_Hits.m_phiModule", m_Hits_m_phiModule, &b_m_Hits_m_phiModule);
   fChain->SetBranchAddress("m_Hits.m_etaWidth", m_Hits_m_etaWidth, &b_m_Hits_m_etaWidth);
   fChain->SetBranchAddress("m_Hits.m_phiWidth", m_Hits_m_phiWidth, &b_m_Hits_m_phiWidth);
   fChain->SetBranchAddress("m_Hits.m_layer", m_Hits_m_layer, &b_m_Hits_m_layer);
   fChain->SetBranchAddress("m_Hits.m_section", m_Hits_m_section, &b_m_Hits_m_section);
   fChain->SetBranchAddress("m_Hits.m_phiIndex", m_Hits_m_phiIndex, &b_m_Hits_m_phiIndex);
   fChain->SetBranchAddress("m_Hits.m_etaIndex", m_Hits_m_etaIndex, &b_m_Hits_m_etaIndex);
   fChain->SetBranchAddress("m_Hits.m_x", m_Hits_m_x, &b_m_Hits_m_x);
   fChain->SetBranchAddress("m_Hits.m_y", m_Hits_m_y, &b_m_Hits_m_y);
   fChain->SetBranchAddress("m_Hits.m_z", m_Hits_m_z, &b_m_Hits_m_z);
   fChain->SetBranchAddress("m_Hits.m_hw_word", m_Hits_m_hw_word, &b_m_Hits_m_hw_word);
   fChain->SetBranchAddress("m_Hits.m_ToT", m_Hits_m_ToT, &b_m_Hits_m_ToT);
   fChain->SetBranchAddress("m_Hits.m_eventindex", m_Hits_m_eventindex, &b_m_Hits_m_eventindex);
   fChain->SetBranchAddress("m_Hits.m_barcode", m_Hits_m_barcode, &b_m_Hits_m_barcode);
   fChain->SetBranchAddress("m_Hits.m_barcode_pt", m_Hits_m_barcode_pt, &b_m_Hits_m_barcode_pt);
   fChain->SetBranchAddress("m_Hits.m_parentage_mask", m_Hits_m_parentage_mask, &b_m_Hits_m_parentage_mask);
   fChain->SetBranchAddress("m_Hits.m_truth.fUniqueID", m_Hits_m_truth_fUniqueID, &b_m_Hits_m_truth_fUniqueID);
   fChain->SetBranchAddress("m_Hits.m_truth.fBits", m_Hits_m_truth_fBits, &b_m_Hits_m_truth_fBits);
   Notify();
}

Bool_t DataList::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DataList::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DataList::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DataList_cxx
