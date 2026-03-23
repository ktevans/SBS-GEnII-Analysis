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

void QE_comp_fix(const char *kinematic, int kin)
{

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  //TString inputfile = Form("/volatile/halla/sbs/ktevans/KateJackSBSAnalysis/KJ_parsed_GEn_pass2_%s_He3_100.root",kinematic);
  //TString inputfile = "/volatile/halla/sbs/cmjackso/Data/GEN2/QE_Parse_1_GEN2_sbs100p_nucleon_np_model2.root";
  TString inputfile = "/volatile/halla/sbs/vimukthi/outfiles/thesis_final/He3/QE_data_GEN2_sbs100p_nucleon_np_model2_sbstrackingon.root";
  TString outputfile = Form("plots/sean_GEn_pass2_%s_He3_dxdy.pdf",kinematic);
  TString outfile = Form("outfiles/sean_GEn_pass2_%s_He3_dxdy.root",kinematic);
  TFile *fout = new TFile(outfile,"RECREATE");

  TTree *T_data = new TTree("T_data", "Analysis Data Tree");

  double dx_out, dy_out, W2_out;
  int helicity_out;
  T_data->Branch("dx", &dx_out, "dx/D");
  T_data->Branch("dy", &dy_out, "dy/D");
  T_data->Branch("W2", &W2_out, "W2/D");
  T_data->Branch("helicity", &helicity_out, "helicity/I");

  TChain* T = new TChain("Tout");
  T->Add(inputfile);

  double W2;         T->SetBranchAddress("W2", &W2);
  double coin_time;          T->SetBranchAddress("coin_time", &coin_time);
  int helicity;             T->SetBranchAddress("helicity", &helicity);
  double ePS;           T->SetBranchAddress("ePS", &ePS);
  double eSH;           T->SetBranchAddress("eSH", &eSH);
  double trP;           T->SetBranchAddress("trP", &trP);
  double vz;          T->SetBranchAddress("vz", &vz);
  double dx_hcal;           T->SetBranchAddress("dx", &dx_hcal);
  double dy_hcal;           T->SetBranchAddress("dy", &dy_hcal);
  int IHWP;              T->SetBranchAddress("IHWP", &IHWP);

  double coin_mean;
  double coin_sigma;

  int IHWP_flip;

  if(kin==2)
  {
    coin_mean = 0.978;
    coin_sigma = 1.179;
    IHWP_flip = -1;
    std::cout << "\nYou are replaying GEN2!\n";
  }
  else
  {
    coin_mean = 0.0;
    coin_sigma = 400.0;
    IHWP_flip = 1;
    std::cout << "\nWhatever you are replaying does not have the proper cuts applied...\n";
  }

  //Scan through all the entries in the TChain T
  //If the rootfiles are empty or don't exist, there will be 0 entries
  //If there are entries, then print out how many
  if(T->GetEntries()==0)
  {
    std::cerr << "\n --- No ROOT file found!! --- \n\n";
    throw;
  }
  else std::cout << "\nFound " << T->GetEntries() << " events. Starting analysis.. \n";

  TH1D* h_dx = new TH1D("h_dx", ";dx", 70.0, -4.0, 3.0);
  h_dx->GetXaxis()->SetTitle("dx [m]");
  h_dx->SetTitle("dx with Global Cuts and QE Cuts");

  TH1D* h_dy = new TH1D("h_dy", ";dy", 60.0, -3.0, 3.0);
  h_dy->GetXaxis()->SetTitle("dy [m]");
  h_dy->SetTitle("dy with Global Cuts and QE Cuts");

  //Loop over all events to fill the histogram
  for (size_t iev = 0; iev < T->GetEntries(); iev++)
  {
    T->GetEntry(iev);

    if(abs(W2-1)<0.5 && abs(coin_time-coin_mean)<coin_sigma && ePS>0.2 && abs(((ePS+eSH)/trP)-1)<0.2 && abs(vz)<0.27)
    {
      
	  h_dx->Fill(dx_hcal);
      h_dy->Fill(dy_hcal);

      dx_out = dx_hcal;
      dy_out = dy_hcal;
      
      W2_out = W2;
      helicity_out = helicity * IHWP_flip *IHWP;
      T_data->Fill();
    }

  }//end event loop

  TCanvas *c1 = new TCanvas("c1","1D dx and dy Plots",100,100,700,700);
  c1->Divide(1,2);
  c1->cd(1);
  h_dx->Draw();
  c1->cd(2);
  h_dy->Draw();

  printf("You've completed the script!\n");

  //Save the canvas to a pdf
  c1->Print(outputfile);

  fout->Write();
}
