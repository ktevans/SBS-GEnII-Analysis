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

void dy_cut()
{

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  TString inputfile = "/volatile/halla/sbs/ktevans/QA/QE_data_GEN2_sbs100p_nucleon_np_model2.root";
  TString outputfile = "plots/dyCut_GEn_pass3_GEN2.pdf";
  TString outfile = "outfiles/dyCut_GEn_pass3_GEN2.root";
  TFile *fout = new TFile(outfile,"RECREATE");

  TTree *T_out = new TTree("T_out", "Analysis Data Tree");

  TChain* T = new TChain("Tout");
  T->Add(inputfile);

  Double_t bb_tr_r_x;             T->SetBranchAddress("bb.tr.r_x", &bb_tr_r_x);
  Double_t bb_tr_r_th;            T->SetBranchAddress("bb.tr.r_th", &bb_tr_r_th);
  Double_t e_kine_W2;             T->SetBranchAddress("e.kine.W2", &e_kine_W2);
  Double_t e_kine_Q2;             T->SetBranchAddress("e.kine.Q2", &e_kine_Q2);
  Double_t p_N;                   T->SetBranchAddress("pN_expect", &p_N);
  Double_t adc_coin;              T->SetBranchAddress("adc.coin", &adc_coin);
  Int_t helicity;                 T->SetBranchAddress("helicity", &helicity);
  Double_t bb_ps_e;               T->SetBranchAddress("bb.ps.e", &bb_ps_e);
  Double_t bb_ps_x;               T->SetBranchAddress("bb.ps.x", &bb_ps_x);
  Double_t bb_ps_y;               T->SetBranchAddress("bb.ps.y", &bb_ps_y);
  Double_t bb_ps_atimeblk;        T->SetBranchAddress("bb.ps.atimeblk", &bb_ps_atimeblk);
  Double_t bb_sh_e;               T->SetBranchAddress("bb.sh.e", &bb_sh_e);
  Double_t bb_sh_x;               T->SetBranchAddress("bb.sh.x", &bb_sh_x);
  Double_t bb_sh_y;               T->SetBranchAddress("bb.sh.y", &bb_sh_y);
  Double_t bb_sh_atimeblk;        T->SetBranchAddress("bb.sh.atimeblk", &bb_sh_atimeblk);
  Double_t bb_tr_p;               T->SetBranchAddress("bb.tr.p", &bb_tr_p);
  Double_t bb_tr_x;               T->SetBranchAddress("bb.tr.x", &bb_tr_x);
  Double_t bb_tr_y;               T->SetBranchAddress("bb.tr.y", &bb_tr_y);
  Double_t bb_tr_vz;              T->SetBranchAddress("bb.tr.vz", &bb_tr_vz);
  Double_t bb_gr_clus_size;       T->SetBranchAddress("bb.grinch_tdc.clus.size", &bb_gr_clus_size);
  Double_t bb_gr_clus_track;      T->SetBranchAddress("bb.grinch_tdc.clus.trackindex", &bb_gr_clus_track);
  Double_t bb_hodotdc_clus_tmean; T->SetBranchAddress("bb.hodotdc.clus.tmean", &bb_hodotdc_clus_tmean);
  Double_t bb_hodotdc_clus_id;    T->SetBranchAddress("bb.hodotdc.clus.id", &bb_hodotdc_clus_id);
  //int pass_global;                T->SetBranchAddress("passGlobal", &pass_global);
  Int_t runnum;                   T->SetBranchAddress("runnum", &runnum);
  Double_t sbs_hcal_e;            T->SetBranchAddress("sbs.hcal.e", &sbs_hcal_e);
  Double_t sbs_hcal_x;            T->SetBranchAddress("sbs.hcal.x", &sbs_hcal_x);
  Double_t sbs_hcal_y;            T->SetBranchAddress("sbs.hcal.y", &sbs_hcal_y);
  Double_t sbs_hcal_atimeblk;     T->SetBranchAddress("sbs.hcal.atimeblk", &sbs_hcal_atimeblk);
  Double_t sbs_hcal_prim_e;       T->SetBranchAddress("sbs.hcal.primary_blk.e", &sbs_hcal_prim_e);
  Double_t sbs_hcal_prim_id;      T->SetBranchAddress("sbs.hcal.primary_blk.id", &sbs_hcal_prim_id);
  Double_t sbs_hcal_sec_e;        T->SetBranchAddress("sbs.hcal.secondary_blk.e", &sbs_hcal_sec_e);
  //Double_t sbs_hcal_clus_blk_id;  T->SetBranchAddress("sbs.hcal.idblk", &sbs_hcal_clus_blk_id);
  Double_t dx_hcal;               T->SetBranchAddress("dx", &dx_hcal);
  Double_t dy_hcal;               T->SetBranchAddress("dy", &dy_hcal);
  Int_t IHWP;                     T->SetBranchAddress("IHWP", &IHWP);
  Double_t gem_tr_ngoodhits;      T->SetBranchAddress("bb.gem.track.ngoodhits", &gem_tr_ngoodhits);
  Double_t gem_tr_chi2ndf;        T->SetBranchAddress("bb.gem.track.chi2ndf", &gem_tr_chi2ndf);

  double optics_valid_min = -0.35;
  double optics_valid_max = 0.34;
  double coin_mean = -0.978;
  double coin_sigma = 3.6;
  double dx_n_mean = -0.121;
  double dx_n_sigma = 0.551;
  double dx_p_mean = -2.752;
  double dx_p_sigma = 0.539;
  double dy_mean = -0.057;
  double dy_sigma = 0.565;
  int firstRun = 2130;
  int lastRun = 2322;
  int totRun = 130;
  double Trp_max = 3.5;
  double Trp_min = 2.0;
  double pN_min = 1.5;
  double pN_max = 3.5;
  double W2_cut_min = 0.0;
  double W2_cut_max = 1.6;

  int runindex = 0;
  int runTrack = 0;

  int IHWP_flip = -1;

  //Scan through all the entries in the TChain T
  //If the rootfiles are empty or don't exist, there will be 0 entries
  //If there are entries, then print out how many
  if(T->GetEntries()==0)
  {
    std::cerr << "\n --- No ROOT file found!! --- \n\n";
    throw;
  }
  else std::cout << "\nFound " << T->GetEntries() << " events. Starting analysis.. \n";

  TH1D* h_dy = new TH1D("h_dy", ";h_dy", 120.0, -8.0, 6.0);
  h_dy->GetXaxis()->SetTitle("dy [m]");

  TH1D* h_dx_anti = new TH1D("h_dx_anti", ";h_dx_anti", 120.0, -8.0, 6.0);
  h_dx_anti->GetXaxis()->SetTitle("dx [m]");

  TH2D* h_dxdy = new TH2D("h_dxdy", ";h_dxdy", 120.0, -8.0, 6.0, 120.0, -8.0, 6.0);
  h_dxdy->GetXaxis()->SetTitle("dy [m]");
  h_dxdy->GetYaxis()->SetTitle("dx [m]");

  //Loop over all events to fill the histogram
  for (size_t iev = 0; iev < T->GetEntries(); iev++)
  {
    T->GetEntry(iev);

    if ((IHWP==-1 || IHWP==1) && (bb_tr_r_x-0.9*bb_tr_r_th)>optics_valid_min && (bb_tr_r_x-0.9*bb_tr_r_th)<optics_valid_max && bb_gr_clus_track==0 && bb_ps_e>0.2 && gem_tr_ngoodhits>=3 && gem_tr_chi2ndf<=15 && abs(adc_coin-coin_mean)<(coin_sigma) && bb_gr_clus_size>2 && abs(((bb_ps_e+bb_sh_e)/bb_tr_p)-0.97)<0.2 && e_kine_W2<2.0)
    {
      //h_dxdy->Fill(dy_hcal,dx_hcal);
      //h_dy->Fill(dy_hcal);

      if (dy_hcal>(dy_mean+(5*dy_sigma)) || dy_hcal<(dy_mean-(5.5*dy_sigma)))
      {
        h_dxdy->Fill(dy_hcal,dx_hcal);
        h_dy->Fill(dy_hcal);
        h_dx_anti->Fill(dx_hcal);
      }
    }

  }//end event loop

  TCanvas *c1 = new TCanvas("c1","dy cut data",100,100,700,700);
  c1->Divide(1,3);
  c1->cd(1);
  h_dxdy->Draw();
  c1->cd(2);
  h_dy->Draw();
  c1->cd(3);
  h_dx_anti->Draw();

  TCanvas *c2 = new TCanvas("c2","dy ant-cut fit",100,100,700,700);
  c2->cd();
  h_dx_anti->Draw("E");

  TF1 *fit_dx = new TF1("fit_dx", "[0] + [1]*x + [2]*TMath::Power(x,2) + [3]*TMath::Power(x,3) + [4]*TMath::Power(x,4)", -8.0, 6.0);
  fit_dx->SetParameters(1,2,3,4,5);
  fit_dx->SetLineColor(kRed);
  h_dx_anti->Fit("fit_dx"); //"W" //Fit using the chi-square method and ignoring the bin uncertainties and skip empty bins.
  fit_dx->Draw("SAMES");

  h_dy->Write();
  h_dx_anti->Write();
  fit_dx->Write();

  fout->Write();
}
