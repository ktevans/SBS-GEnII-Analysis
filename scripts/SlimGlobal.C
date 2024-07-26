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

void SlimGlobal(const std::string configfileinput)
{

 TString configfilename = "../config/" + configfileinput + ".cfg";
 cout<<"Reading from configuration file: "<<configfilename<<endl;

 // Add root files to the chain C, can add as many files as you want
 TChain *C = new TChain("T");

 // Set location and name for output file 
 TString outputfilename = "../outfiles/" + configfileinput + ".root";
 TString outputPDFname = "../outfiles/" + configfileinput + "_plots.pdf";

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
 double BBtr_n, BBps_e;
 double ekineW2;
 
 //Set the branch status for variables of interest
  C ->SetBranchStatus("*",0);
 
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.vz",1);
  C->SetBranchStatus("bb.ps.e",1);
  C->SetBranchStatus("e.kine.W2",1);

 //Define branch addresses
  C->SetBranchAddress("bb.tr.n", &BBtr_n);
  C->SetBranchAddress("bb.tr.vz", BBtr_vz);
  C->SetBranchAddress("bb.ps.e", &BBps_e);
  C->SetBranchAddress("e.kine.W2", &ekineW2);

 //Make event list for all events which pass the global cut
 C->Draw(">>elist",globalcut);
 TEventList *globalList = (TEventList*)gDirectory->Get("elist");
 globalList->SetReapplyCut(kTRUE);
 C->SetEventList(globalList);

 //Create output tree, P that will contain global cuts and calculated values
 TTree *P= new TTree("P", "Parsed Tree");
    
 // "old" Tree variables we want to keep
 double bb_ps_e_out;
 double W2_out;

 //Create the new branches
 P->Branch("bb_ps_e", &bb_ps_e_out, "bb_ps_e/D");
 P->Branch("W2", &W2_out, "W2/D");

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

 //Define 1D histogram to plot W2
 TH1D* h_W = new TH1D("h_W",";W",150,0.0,3.0);
 h_W->GetXaxis()->SetTitle("W^2 [GeV^2]");

 //Define 1D histogram to plot ps_e
 TH1D* h_ps_e = new TH1D("h_ps_e",";ps_e",100,0.0,2.0);
 h_ps_e->GetXaxis()->SetTitle("PreShower Energy [GeV^2]");

 //Loop over all events in the event list to fill the new tree
 for (size_t iev = 0; iev < Nevents; iev++)
 {

   C->GetEntry(globalList->GetEntry(iev));
   
   bb_ps_e_out = BBps_e;
   W2_out = ekineW2;

   P->Fill();

 }//end process over events

 //Loop over all events in the new tree to fill histograms
 for (size_t iev = 0; iev < P->GetEntries(); iev++)
 {

   P->GetEntry(iev);

   h_W->Fill(W2_out);
   h_ps_e->Fill(bb_ps_e_out);

 }//end second process over events

 //Draw canvas, divide into two plots and draw on each one
 TCanvas *c1 = new TCanvas("c1","dxdy",100,100,700,700);
 c1->Divide(1,2);
 c1->cd(1);
 h_W->Draw();
 c1->cd(2);
 h_ps_e->Draw();

 //Save the canvas to a pdf and write to output rootfile
 c1 -> Print(outputPDFname);
 fout -> Write();
 P -> Write();

}
