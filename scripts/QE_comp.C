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

void QE_comp(const char *kinematic)
{

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  TString inputfile = Form("/volatile/halla/sbs/ktevans/KateJackSBSAnalysis/KJ_parsed_GEn_pass2_%s_He3_100.root",kinematic);
  TString outputfile = Form("plots/parsed_GEn_pass2_%s_He3_dxdy.pdf",kinematic);
  TString outfile = Form("outfiles/parsed_GEn_pass2_%s_He3_dxdy.root",kinematic);
  TFile *fout = new TFile(outfile,"RECREATE");

  TTree *T_data = new TTree("T_data", "Analysis Data Tree");

  double dx_out, dy_out, W2_out;
  int helicity_out;
  T_data->Branch("dx", &dx_out, "dx/D");
  T_data->Branch("dy", &dy_out, "dy/D");
  T_data->Branch("W2", &W2_out, "W2/D");
  T_data->Branch("helicity", &helicity_out, "helicity/I");

  TChain* T = new TChain("Parse");
  T->Add(inputfile);

  double bb_tr_r_x;         T->SetBranchAddress("bb.tr.r_x", &bb_tr_r_x);
  double bb_tr_r_th;        T->SetBranchAddress("bb.tr.r_th", &bb_tr_r_th);
  double e_kine_W2;         T->SetBranchAddress("e.kine.W2", &e_kine_W2);
  double adc_coin;          T->SetBranchAddress("adc.coin", &adc_coin);
  int helicity;             T->SetBranchAddress("helicity", &helicity);
  double bb_ps_e;           T->SetBranchAddress("bb.ps.e", &bb_ps_e);
  double bb_sh_e;           T->SetBranchAddress("bb.sh.e", &bb_sh_e);
  double bb_tr_p;           T->SetBranchAddress("bb.tr.p", &bb_tr_p);
  double bb_tr_vz;          T->SetBranchAddress("bb.tr.vz", &bb_tr_vz);
  double bb_gr_clus_size;   T->SetBranchAddress("bb.grinch_tdc.clus.size", &bb_gr_clus_size);
  int pass_global;          T->SetBranchAddress("passGlobal", &pass_global);
  double dx_hcal;           T->SetBranchAddress("dx", &dx_hcal);
  double dy_hcal;           T->SetBranchAddress("dy", &dy_hcal);

  double optics_valid_min;
  double optics_valid_max;
  double coin_mean;
  double coin_sigma;

  if(kinematic=="GEN2")
  {
    optics_valid_min = -0.35;
    optics_valid_max = 0.34;
    coin_mean = 129.1;
    coin_sigma = 5.6;
  }
  else if(kinematic=="GEN3")
  {
    optics_valid_min = -0.35;
    optics_valid_max = 0.33;
    coin_mean = 120.3;
    coin_sigma = 6.0;
  }
  else if(kinematic=="GEN4a")
  {
    optics_valid_min = -0.36;
    optics_valid_max = 0.30;
    coin_mean = 121.4;
    coin_sigma = 5.8;
  }
  else if(kinematic=="GEN4b")
  {
    optics_valid_min = -0.37;
    optics_valid_max = 0.32;
    coin_mean = 185.5;
    coin_sigma = 7.0;
  }
  else
  {
    optics_valid_min = -2.0;
    optics_valid_max = 2.0;
    coin_mean = 0.0;
    coin_sigma = 400.0;
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

    if(pass_global==1 && (bb_tr_r_x-0.9*bb_tr_r_th)>optics_valid_min && (bb_tr_r_x-0.9*bb_tr_r_th)<optics_valid_max && abs(adc_coin-coin_mean)<coin_sigma && bb_ps_e>0.2 && abs(((bb_ps_e+bb_sh_e)/bb_tr_p)-1)<0.2 && bb_gr_clus_size>2.0 && abs(bb_tr_vz)<0.27)
    {
      
      if(abs(e_kine_W2-1.0)<0.5)
      {
	h_dx->Fill(dx_hcal);
        h_dy->Fill(dy_hcal);

        dx_out = dx_hcal;
        dy_out = dy_hcal;
      }
      
      W2_out = e_kine_W2;
      helicity_out = helicity;
    }

    T_data->Fill();

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
