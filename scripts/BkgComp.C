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

  TChain* Tout = new TChain("Tout");
  Tout->Add(DataFile);

  Double_t dx;            Tout->SetBranchAddress("dx", &dx);
  Double_t dy;            Tout->SetBranchAddress("dy", &dy);
  int helicity;           Tout->SetBranchAddress("helicity", &helicity);
  double coin;            Tout->SetBranchAddress("adc.coin", &coin);
  double ps_e;            Tout->SetBranchAddress("bb.ps.e", &ps_e);
  double sh_e;            Tout->SetBranchAddress("bb.sh.e", &sh_e);
  double W2;              Tout->SetBranchAddress("e.kine.W2", &W2);
  double Q2;              Tout->SetBranchAddress("e.kine.Q2", &Q2);
  double tr_p;            Tout->SetBranchAddress("bb.tr.p", &tr_p);
  double grinch_track;    Tout->SetBranchAddress("bb.grinch_tdc.clus.trackindex", &grinch_track);
  double grinch_clusSize; Tout->SetBranchAddress("bb.grinch_tdc.clus.size", &grinch_clusSize);

  TH1D* h_dx_acc = new TH1D("h_dx_acc","Accidentals", 100, -3.0, 2.0);
  h_dx_acc->GetXaxis()->SetTitle("dx [m]");
  h_dx_acc->SetLineColor(kRed);

  for (size_t iev = 0; iev < Tout->GetEntries(); iev++)
  {
    Tout->GetEntry(iev);

    if(W2>0.4 && W2<1.6 && abs(coin+0.47385)>6.0 && grinch_track==0.0 && grinch_clusSize>=3.0 && ps_e>0.2 && abs(((ps_e+sh_e)/tr_p)-0.97)<0.2 && abs(dy)<0.88)
    {
      h_dx_acc->Fill(dx);
    }

  }

  h_dx_acc->Scale(0.017);


  TCanvas *c1 = new TCanvas("c1", "Background Shapes", 100,100,800,800);
  c1->cd();
  //scaled_h_sim_nucleons->Draw("HIST");
  h_N2_scaled->Draw("SAMES");
  h_dx->Draw("HIST SAMES");
  h_dx_acc->Draw("HIST SAMES");

}
