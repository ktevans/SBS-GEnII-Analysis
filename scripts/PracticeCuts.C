// This script will read in a config file with a list of root files and a simple
// global cut. It will apply the global cut to the events in the root files and
// output a new root file with only events which passed the cut.
// This simple version only considers W2 and ps energy.
// To execute, do the following:
//
// > root
// > .L scripts/SlimGlobal.C
// > SlimGlobal("file")
//
// Note: you must have a directory called "outfiles" created in order to save
// the canvas to a PDF properly.
//
// Kate Evans
// Last modified: July 26 2024

#include <TMath.h>
#include <TF1.h>
#include <TSystem.h>
#include <TChain.h>
#include <TString.h>
#include <TNtuple.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <stack>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TCut.h"
#include "TLatex.h"
#include "TLine.h"

const int maxTracks = 16; //Why 16?

void PracticeCuts(const std::string configfileinput)
{

 TString configfilename = "../config/" + configfileinput + ".cfg";
 cout<<"Reading from configuration file: "<<configfilename<<endl;

 // Add root files to the chain C, can add as many files as you want
 TChain *C = new TChain("T");

 // Set location and name for output file 
 TString outputfilename = "../outfiles/" + configfileinput + ".root";

 // Declare outfile
 TFile *fout = new TFile(outputfilename,"RECREATE");
 cout<<"Writing to file: "<< outputfilename <<endl;

 ifstream configfile(configfilename);
 TString currentline;

 while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endlist") )
 {
    if( !currentline.BeginsWith("#") )
    {
      if(!currentline) cout << "WARNING: No file exists at " << currentline << "." << endl;
      C->Add(currentline);
      cout << "Loaded file at: " << currentline << endl;
    }    
  }

 TCut globalcut = "";
 while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endcut") )
 {
   if( !currentline.BeginsWith("#") )
   {
     globalcut += currentline;
   }    
   cout<< "Global Cut: "<<globalcut<<endl;
 }

 //Defining variables
 double BBtr_vz[maxTracks];
 double BBtr_n, BBtr_x, BBtr_y, BBtr_ph, BBtr_th;
 double BBps_e, BBps_x, BBps_y;
 double BBsh_e, BBsh_x, BBsh_y;
 double BB_etot_over_p;
 double BBgr_size, BBgr_x, BBgr_y, BBgr_tot;
 double ekineW2;
 
 //Set the branch status for variables of interest
  C ->SetBranchStatus("*",0);
 
  C->SetBranchStatus("bb.tr.vz",1);
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.x",1);
  C->SetBranchStatus("bb.tr.y",1);
  C->SetBranchStatus("bb.tr.th",1);
  C->SetBranchStatus("bb.tr.ph",1);

  C->SetBranchStatus("bb.ps.e",1);
  C->SetBranchStatus("bb.ps.x",1);
  C->SetBranchStatus("bb.ps.y",1);

  C->SetBranchStatus("bb.sh.e",1);
  C->SetBranchStatus("bb.sh.x",1);
  C->SetBranchStatus("bb.sh.y",1);

  C->SetBranchStatus("bb.grinch_tdc.clus.tot_mean",1);
  C->SetBranchStatus("bb.grinch_tdc.clus.size",1);
  C->SetBranchStatus("bb.grinch_tdc.clus.x_mean",1);
  C->SetBranchStatus("bb.grinch_tdc.clus.y_mean",1);

  C->SetBranchStatus("bb.etot_over_p",1);
  C->SetBranchStatus("e.kine.W2",1);

 //Define branch addresses
  C->SetBranchAddress("bb.tr.vz", &BBtr_vz);
  C->SetBranchAddress("bb.tr.n", &BBtr_n);
  C->SetBranchAddress("bb.tr.x", &BBtr_x);
  C->SetBranchAddress("bb.tr.y", &BBtr_y);
  C->SetBranchAddress("bb.tr.th", &BBtr_th);
  C->SetBranchAddress("bb.tr.ph", &BBtr_ph);
  
  C->SetBranchAddress("bb.ps.e", &BBps_e);
  C->SetBranchAddress("bb.ps.x", &BBps_x);
  C->SetBranchAddress("bb.ps.y", &BBps_y);

  C->SetBranchAddress("bb.sh.e", &BBsh_e);
  C->SetBranchAddress("bb.sh.x", &BBsh_x);
  C->SetBranchAddress("bb.sh.y", &BBsh_y);

  C->SetBranchAddress("bb.grinch_tdc.clus.tot_mean", &BBgr_tot);
  C->SetBranchAddress("bb.grinch_tdc.clus.size", &BBgr_size);
  C->SetBranchAddress("bb.grinch_tdc.clus.x_mean", &BBgr_x);
  C->SetBranchAddress("bb.grinch_tdc.clus.y_mean", &BBgr_y);

  C->SetBranchAddress("bb.etot_over_p", &BB_etot_over_p);
  C->SetBranchAddress("e.kine.W2", &ekineW2);

 //Make event list for all events which pass the global cut
 C->Draw(">>elist",globalcut);
 TEventList *globalList = (TEventList*)gDirectory->Get("elist");
 globalList->SetReapplyCut(kTRUE);
 C->SetEventList(globalList);

 //Create output tree, P that will contain global cuts and calculated values
 TTree *P= new TTree("P", "Parsed Tree");
    
 // "old" Tree variables we want to keep
 //double BBtr_vz_out[maxTracks];
 double BBtr_n_out, BBtr_x_out, BBtr_y_out, BBtr_ph_out, BBtr_th_out;
 double BBps_e_out, BBps_x_out, BBps_y_out;
 double BBsh_e_out, BBsh_x_out, BBsh_y_out;
 double BB_etot_over_p_out;
 double BBgr_size_out, BBgr_x_out, BBgr_y_out, BBgr_tot_out;
 double ekineW2_out;

 //Create the new branches
 //P->Branch("bb_tr_vz", &BBtr_vz_out, "bb_tr_vz/D");
 P->Branch("bb_tr_n", &BBtr_n_out, "bb_tr_n/D");
 P->Branch("bb_tr_x", &BBtr_x_out, "bb_tr_x/D");
 P->Branch("bb_tr_y", &BBtr_y_out, "bb_tr_y/D");
 P->Branch("bb_tr_th", &BBtr_th_out, "bb_tr_th/D");
 P->Branch("bb_tr_ph", &BBtr_ph_out, "bb_tr_ph/D");

 P->Branch("bb_ps_e", &BBps_e_out, "bb_ps_e/D");
 P->Branch("bb_ps_x", &BBps_x_out, "bb_ps_x/D");
 P->Branch("bb_ps_y", &BBps_y_out, "bb_ps_y/D");

 P->Branch("bb_sh_e", &BBsh_e_out, "bb_sh_e/D");
 P->Branch("bb_sh_x", &BBsh_x_out, "bb_sh_x/D");
 P->Branch("bb_sh_y", &BBsh_y_out, "bb_sh_y/D");

 P->Branch("bb_gr_size", &BBgr_size_out ,"bb_gr_size/D");
 P->Branch("bb_gr_x", &BBgr_x_out ,"bb_gr_x/D");
 P->Branch("bb_gr_y", &BBgr_y_out ,"bb_gr_y/D");
 P->Branch("bb_gr_tot", &BBgr_tot_out ,"bb_gr_tot/D");

 P->Branch("bb_etot_over_p", &BB_etot_over_p_out,"bb_etot_over_p/D");
 P->Branch("e_kine_W2", &ekineW2_out, "e_kine_W2/D");
 

 //Get total number of events in event list
 Long64_t Nevents = globalList->GetN();
 
 //Scan through all the entries in the EventList elist
 //If the rootfiles are empty or don't exist, there will be 0 entries
 //If there are entries, then print out how many
 if(C->GetEntries()==0)
 {
   std::cerr << "\n --- No ROOT file found!! --- \n\n";
   throw;
 }
 else std::cout << "\nFound " << C->GetEntries() << " events. Starting analysis.. \n";

 //Loop over all events in the event list to fill the new tree
 for (size_t iev = 0; iev < Nevents; iev++)
 {

   C->GetEntry(globalList->GetEntry(iev));

   //BBtr_vz_out[maxTracks] = BBtr_vz[maxTracks];
   BBtr_n_out = BBtr_n;
   BBtr_x_out = BBtr_x;
   BBtr_y_out = BBtr_y;
   BBtr_ph_out = BBtr_ph;
   BBtr_th_out = BBtr_th;

   BBps_e_out = BBps_e;
   BBps_x_out = BBps_x;
   BBps_y_out = BBps_y;

   BBsh_e_out = BBsh_e;
   BBsh_x_out = BBsh_x;
   BBsh_y_out = BBsh_y;

   BBgr_size_out = BBgr_size;
   BBgr_x_out = BBgr_x;
   BBgr_y_out = BBgr_y;
   BBgr_tot_out = BBgr_tot;
   
   BB_etot_over_p_out = BB_etot_over_p;
   ekineW2_out = ekineW2;

   P->Fill();

 }//end process over events

 //Write to output rootfile
 fout -> Write();
 P -> Write();

}
