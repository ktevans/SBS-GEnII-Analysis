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

void BkgComp()
{

  TString N2File = "outfiles/N2_Corr_SIM_GEn_GEN2_He3_dxdy.root";
  TString InelFile = "outfiles/parsed_SIM_IN_GEn_GEN2_He3_dxdy.root"; //includes QE cuts
  TString DataFile = "/volatile/halla/sbs/ktevans/pass3/QE_data_GEN2_sbs100p_nucleon_np_model2.root"; //Only global cuts
  TString FittedFile = "outfiles/AnalysisResults_GEN2.root"; //scaled and shifted histograms live here

  TFile *T_N2File = TFile::Open(N2File);
  TH1D *hN2dilution_p = nullptr;
  T_N2File->GetObject("hN2dilution_p", hN2dilution_p);

  TFile *T_FittedFile = TFile::Open(FittedFile);
  TH1D *scaled_h_sim_nucleons = nullptr;
  T_FittedFile->GetObject("scaled_h_sim_nucleons", scaled_h_sim_nucleons);

  TH1D *h_N2_scaled = (TH1D*)hN2dilution_p->Clone("h_N2_scaled");
  h_N2_scaled->Multiply(scaled_h_sim_nucleons);
  //h_N2_scaled->Scale(1.0/h_N2_scaled->Integral());
  h_N2_scaled->SetLineColor(kBlue);

  TFile *T_InelFile = TFile::Open(InelFile);
  TH1D *h_dx = nullptr;
  T_InelFile->GetObject("h_dx", h_dx);

  //h_dx->Scale(1.0/h_dx->Integral());
  h_dx->SetLineColor(kGreen);

  TCanvas *c1 = new TCanvas("c1", "Background Shapes", 100,100,800,800);
  c1->cd();
  h_N2_scaled->Draw("HIST");
  h_dx->Draw("HIST SAMES");

}
