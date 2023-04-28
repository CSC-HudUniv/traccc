#define DataList_cxx
#include "traccc/clusterization/DataList.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "traccc/clusterization/IDreader.h"

void DataList::WriteCellFiles(Int_t nevents, string dirnm=".")
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Int_t nfiles = (nevents < nentries ? nevents : nentries);
   // if(nfiles>5){
   //   nfiles=5;   //for testing purpose
   //   cout<<"Currently, for testing purpose, first 5 events will be read."<<endl;
   // }
     
   cout<<nfiles<<" cell files will be created."<<endl;

   IDreader *idr = new IDreader("trackml-detector.csv");
   //   int nGeoIDs = idr->NIDs();
   
   // Write into csv files
   Long64_t nb=0;
   for (Long64_t jentry=0; jentry<nfiles;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry); 

      Char_t filename[100];
      sprintf(filename,"%s/event0000%05d-cells.csv",dirnm.c_str(),(int)jentry);
      
      ofstream ofs(filename);
      if(!ofs){
	cerr<<"failed to create file"<<endl;
	exit(1);
      }

      ofs<<"geometry_id,hit_id,channel0,channel1,timestamp,value"<<endl;
      cout<<jentry<<" "<<m_Hits_<<endl;
      
      for(int ihit=0; ihit<m_Hits_ ; ihit++){
	if(m_Hits_m_hitType[ihit]!=0) continue;
	uint atlasid = m_Hits_m_identifierHash[ihit];
	int ifpixel  = m_Hits_m_detType[ihit];
	unsigned long int geometry_id
	  = idr->Give50muPixelDetectorID(atlasid,ifpixel);
	int hit_id = ihit;
	int channel0 = m_Hits_m_phiIndex[ihit];
	int channel1 = m_Hits_m_etaIndex[ihit];
	double timestamp = 0;
	double value = m_Hits_m_ToT[ihit];

	ofs<<geometry_id<<","<<hit_id<<","<<channel0<<","<<channel1<<","
	   <<timestamp<<","<<value<<endl;
	// cout<<geometry_id<<","<<hit_id<<","<<channel0<<","<<channel1<<","
	//     <<timestamp<<","<<value<<","
	//     <<m_Hits_m_detType[ihit]<<","<<m_Hits_m_side[ihit]<<endl;
      }

      ofs.close();
      
   }
}

void DataList::WritePixelCellFiles(Int_t nevents, string dirnm=".")
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Int_t nfiles = (nevents < nentries ? nevents : nentries);
     
   cout<<nfiles<<" cell files will be created."<<endl;

   IDreader *idr = new IDreader("trackml-detector.csv");
   
   for (Long64_t jentry=0; jentry<nfiles;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      fChain->GetEntry(jentry); 

      Char_t filename[100];
      sprintf(filename,"%s/pixels/event0000%05d-cells.csv",dirnm.c_str(),(int)jentry);
      
      ofstream ofs(filename);
      if(!ofs){
	cerr<<"failed to create file"<<endl;
	exit(1);
      }

      ofs<<"geometry_id,hit_id,channel0,channel1,timestamp,value"<<endl;
      cout<<jentry<<" "<<m_Hits_<<endl;
      
      for(int ihit=0; ihit<m_Hits_ ; ihit++){
	if(m_Hits_m_hitType[ihit]!=0) continue;
	if(m_Hits_m_detType[ihit]!=1) continue; 
	uint atlasid = m_Hits_m_identifierHash[ihit];
	unsigned long int geometry_id = idr->Give50muPixelDetectorID(atlasid,1);
	int hit_id = ihit;
	int channel0 = m_Hits_m_phiIndex[ihit];
	int channel1 = m_Hits_m_etaIndex[ihit];
	double timestamp = 0;
	double value = m_Hits_m_ToT[ihit];

	ofs<<geometry_id<<","<<hit_id<<","<<channel0<<","<<channel1<<","
	   <<timestamp<<","<<value<<endl;
	// cout<<geometry_id<<","<<hit_id<<","<<channel0<<","<<channel1<<","
	//     <<timestamp<<","<<value<<","
	//     <<m_Hits_m_detType[ihit]<<","<<m_Hits_m_side[ihit]<<endl;
      }

      ofs.close();
      
   }
   cout<<"PixelCellFiles are created"<<endl;
}

void DataList::WriteClusters(Int_t nevents, string dirnm="."){

  if(fChain ==0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Int_t nfiles = (nevents < nentries ? nevents : nentries);
  // if(nfiles>5){
  //   nfiles=5;   //for testing purpose
  //   cout<<"Currently, for testing purpose, first 5 events will be read."<<endl;
  // }
  
  cout<<nfiles<<" cluster files will be created."<<endl;

  IDreader *idr = new IDreader("trackml-detector.csv");

  Char_t filename[100];
  sprintf(filename,"%s/Nclusters.txt",dirnm.c_str());
  ofstream ofs_summary(filename);
 
  for(Long64_t jentry=0; jentry<nfiles; jentry++){
    Long64_t ientry = LoadTree(jentry);
    if(ientry<0) break;
    fChain->GetEntry(jentry);

    sprintf(filename,"%s/event0000%05d-clusters.csv",dirnm.c_str(),(int)jentry);

    ofstream ofs(filename);
    if(!ofs){
      cerr<<"failed to create file"<<endl;
      exit(1);
    }
    
    cout<<m_optional_m_OfflineClusters_<<endl;
    ofs<<"pixel/silicon,atlas_id,geometry_id,side,channel0,channel1,x,y,z"
       <<endl;

    int ncluster=0;
    for(int icls =0; icls < m_optional_m_OfflineClusters_; icls++){

      int ifpixel = m_optional_m_OfflineClusters_m_clusterEquiv_m_detType[icls];
      uint atlasid
	= m_optional_m_OfflineClusters_m_clusterEquiv_m_identifierHash[icls];
      unsigned long int geometry_id =
	idr->Give50muPixelDetectorID(atlasid,ifpixel);
		
      ofs<<m_optional_m_OfflineClusters_m_clusterEquiv_m_detType[icls]<<","
	 <<atlasid<<","
	 <<geometry_id<<","
	 <<m_optional_m_OfflineClusters_m_clusterEquiv_m_side[icls]<<","
	 <<m_optional_m_OfflineClusters_m_clusterEquiv_m_phiIndex[icls]<<","
	 <<m_optional_m_OfflineClusters_m_clusterEquiv_m_etaIndex[icls]<<","
	 <<m_optional_m_OfflineClusters_m_clusterEquiv_m_x[icls]<<","
	 <<m_optional_m_OfflineClusters_m_clusterEquiv_m_y[icls]<<","
	 <<m_optional_m_OfflineClusters_m_clusterEquiv_m_z[icls]
	 <<endl;

      ncluster++;
      
    }
    ofs_summary<<ncluster<<endl;
    ofs.close();
    
  }
  ofs_summary.close();
}
    
void DataList::WritePixelClusters(Int_t nevents, string dirnm="."){

  if(fChain ==0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Int_t nfiles = (nevents < nentries ? nevents : nentries);
  
  cout<<nfiles<<" cluster files will be created."<<endl;

  IDreader *idr = new IDreader("trackml-detector.csv");

  Char_t filename[100];
  sprintf(filename,"%s/pixels/Nclusters.txt",dirnm.c_str());
  ofstream ofs_summary(filename);

  for(Long64_t jentry=0; jentry<nfiles; jentry++){
    Long64_t ientry = LoadTree(jentry);
    if(ientry<0) break;
    fChain->GetEntry(jentry);

    Char_t filename[100];
    sprintf(filename,"%s/pixels/event0000%05d-clusters.csv",dirnm.c_str(),(int)jentry);

    ofstream ofs(filename);
    if(!ofs){
      cerr<<"failed to create file"<<endl;
      exit(1);
    }
    
    cout<<m_optional_m_OfflineClusters_<<endl;
    ofs<<"pixel/silicon,atlas_id,geometry_id,channel0,channel1,x,y,z"<<endl;

    int nclusters=0;
    for(int icls =0; icls < m_optional_m_OfflineClusters_; icls++){
      if(m_optional_m_OfflineClusters_m_clusterEquiv_m_detType[icls]!=1)
	continue;
      
      uint atlasid
	= m_optional_m_OfflineClusters_m_clusterEquiv_m_identifierHash[icls];
      unsigned long int geometry_id = idr->Give50muPixelDetectorID(atlasid,1);
		
      ofs<<m_optional_m_OfflineClusters_m_clusterEquiv_m_detType[icls]<<","
	 <<atlasid<<","
	 <<geometry_id<<","
	 <<m_optional_m_OfflineClusters_m_clusterEquiv_m_phiIndex[icls]<<","
	 <<m_optional_m_OfflineClusters_m_clusterEquiv_m_etaIndex[icls]<<","
	 <<m_optional_m_OfflineClusters_m_clusterEquiv_m_x[icls]<<","
	 <<m_optional_m_OfflineClusters_m_clusterEquiv_m_y[icls]<<","
	 <<m_optional_m_OfflineClusters_m_clusterEquiv_m_z[icls]
	 <<endl;
      nclusters++;

    }
    ofs_summary<<nclusters<<endl;
    ofs.close();
    
  }
  ofs_summary.close();
}
    

void DataList::DrawVariables()
{

  TH2D *h[9];

  TFile *f = new TFile("Variables.root","recreate");
  h[0]= new TH2D("h0","detType vs layer disk",6,0,6,2,0,2);
  h[1]= new TH2D("h1","detectorZone vs layer disk",6,0,6,3,0,3);
  h[2]= new TH2D("h2","detType vs identifierHash",200,0,20000,2,0,2);
  h[3]= new TH2D("h3","detectorZone vs identifierHash",200,0,20000,3,0,3);
  h[4]= new TH2D("h4","layer disk vs identifierHash",200,0,20000,6,0,6);
  h[5]= new TH2D("h5","m_side vs identifierHash",200,0,20000,2,0,2);
  h[6]= new TH2D("h6","etaModule vs identifierHash",200,0,20000,40,-20,20);
  h[7]= new TH2D("h7","phiModule vs identifierHash",200,0,20000,80,0,80);

  TH2D *hpix[3];
  hpix[0] = new TH2D("hpix0","hit r vs z (pixel, Zone=0)",
		     150,-3000,3000,80,0,320);
  hpix[1] = new TH2D("hpix1","hit r vs z (pixel, Zone=1)",
		     150,-3000,3000,80,0,320);
  hpix[2] = new TH2D("hpix2","hit r vs z (pixel, Zone=2)",
		     150,-3000,3000,80,0,320);
    
  Long64_t nentries = fChain->GetEntriesFast();
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      fChain->GetEntry(jentry); 

      for(int ihit=0; ihit<m_Hits_ ; ihit++){
	if(m_Hits_m_hitType[ihit]!=0) continue;

	h[0]->Fill(m_Hits_m_layer_disk[ihit], m_Hits_m_detType[ihit]);
	h[1]->Fill(m_Hits_m_layer_disk[ihit], m_Hits_m_detectorZone[ihit]);
	h[2]->Fill(m_Hits_m_identifierHash[ihit], m_Hits_m_detType[ihit]);
	h[3]->Fill(m_Hits_m_identifierHash[ihit], m_Hits_m_detectorZone[ihit]);
	h[4]->Fill(m_Hits_m_identifierHash[ihit], m_Hits_m_layer_disk[ihit]);
	if(m_Hits_m_detType[ihit]==0)
	  h[5]->Fill(m_Hits_m_identifierHash[ihit], m_Hits_m_side[ihit]);
	h[6]->Fill(m_Hits_m_identifierHash[ihit], m_Hits_m_etaModule[ihit]);
	h[7]->Fill(m_Hits_m_identifierHash[ihit], m_Hits_m_phiModule[ihit]);
		

	if(m_Hits_m_detType[ihit]==1){ //pixel
	  int zoneid = m_Hits_m_detectorZone[ihit];
	  double rr = sqrt(pow(m_Hits_m_x[ihit],2) + pow(m_Hits_m_y[ihit],2));
	  double zz = m_Hits_m_z[ihit];

	  if(zoneid<3 && zoneid>=0) hpix[zoneid]->Fill(zz, rr);
	}
      }
  }

  f->cd();
  TCanvas *c1 = new TCanvas("c1","Canvas",900,600);
  c1->Divide(4,2);
  for(int i=0; i<8; i++){
    c1->cd(i+1);
    h[i]->Draw("colz");
  
  }
  c1->Write();

  TCanvas *c2 = new TCanvas("c2","Pixels",800,400);
  hpix[0]->Draw();
  for(int i=0; i<3; i++){
    hpix[i]->SetMarkerColor(pow(2,i));
    hpix[i]->Draw("same");
  }
  c2->Write(); 
  
  f->Write();
  f->Close();
}
			  
