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

  TH1D* h_dy = new TH1D("h_dy", ";h_dy", 60.0, -3.0, 3.0);
  h_dy->GetXaxis()->SetTitle("dy [m]");

  TH1D* h_dy_anti = new TH1D("h_dy_anti", ";h_dy_anti", 60.0, -3.0, 3.0);
  h_dy_anti->GetXaxis()->SetTitle("dy [m]");

  TH2D* h_dxdy = new TH2D("h_dxdy", ";h_dxdy", 90.0, -6.0, 3.0, 60.0, -3.0, 3.0);
  h_dxdy->GetXaxis()->SetTitle("dy [m]");
  h_dxdy->GetYaxis()->SetTitle("dx [m]");

  //Loop over all events to fill the histogram
  for (size_t iev = 0; iev < T->GetEntries(); iev++)
  {
    T->GetEntry(iev);

    if ((IHWP==-1 || IHWP==1) && (bb_tr_r_x-0.9*bb_tr_r_th)>optics_valid_min && (bb_tr_r_x-0.9*bb_tr_r_th)<optics_valid_max && bb_gr_clus_track==0 && bb_ps_e>0.2 && gem_tr_ngoodhits>=3 && gem_tr_chi2ndf<=15 && abs(adc_coin-coin_mean)<(coin_sigma) && bb_gr_clus_size>2 && abs(((bb_ps_e+bb_sh_e)/bb_tr_p)-0.97)<0.2 && abs(e_kine_W2-1.0)<0.5 && (abs(dx_hcal-dx_n_mean)<dx_n_sigma || abs(dx_hcal-dx_p_mean)<dx_p_sigma))
    {
      h_dxdy->Fill(dy_hcal,dx_hcal);
      h_dy->Fill(dy_hcal);

      if (dy_hcal>(dy_mean+dy_sigma)||dy_hcal<(dy_mean-dy_sigma))
      {
        h_dy_anti->Fill(dy_hcal);
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
  h_dy_anti->Draw();

  // ----- come back to this ----
  //c3->Divide(1,2);
  //c3->cd(1);
  //h_prof_pol_n_low->Draw();
  // ----- come back to this ----

  //TF1 *fitn = new TF1("fitn", "[0] + [1]*x + [2]*TMath::Power(x,2) + [3]*TMath::Power(x,3) + [4]*TMath::Power(x,4) + [5]*TMath::Power(x,5) + [6]*TMath::Power(x,6) + [7]*TMath::Power(x,7) + [8]*TMath::Power(x,8) + [9]*TMath::Power(x,9)", -3.5, 3.0);
  //fitn->SetParameters(1,2,3,4,5,6,7,8,9,10);
  //fitn->SetLineColor(kRed);

  //h_prof_pol_n->Fit("fitn");
  //fitn->Draw("SAMES");
  //gStyle->SetOptFit(1111);

  // ----- come back to this ----
  //TF1 *fitn_low = new TF1("fitn_low", "[0]", -4.0, dx_n_cut);
  //fitn_low->SetParameters(1.0);
  //fitn_low->SetLineColor(kBlue);
  //h_prof_pol_n_low->Fit("fitn_low");
  //fitn_low->Draw("SAMES");

  //c3->cd(2);
  //h_prof_pol_n_high->Draw();

  //TF1 *fitn_high = new TF1("fitn_high", "[0]*x*x + [1]*x + [2]", dx_n_cut, 3.0);
  //fitn_high->SetParameters(1.0,2.0,3.0);
  //fitn_high->SetLineColor(kRed);
  //h_prof_pol_n_high->Fit("fitn_high");
  //fitn_high->Draw("SAMES");

  // ----- ^^ come back to this ----


  //*** Testing Fits ***
  //************
  //TF1 *f = new TF1("f",[=](double *x, double */*p*/){return h_prof_pol_p->Interpolate(x[0]);},h_prof_pol_p->GetXaxis()->GetXmin(), h_prof_pol_p->GetXaxis()->GetXmax(), 0);
  //f->Draw("same");
  //f->Write();
  //TSpline3 *spline3 = nullptr;
  //delete spline3;
  //spline3 = new TSpline3(h_prof_pol_p, "b1e1", f->Derivative(-3.5), f->Derivative(3.0));
  // "b1" check first derivative exists for beginning point (f->Derivative(-3.5))
  // "e1" check first derivative exists for ending point (f->Derivative(3.0))
  //spline3->SetLineColor(kGreen);
  //spline3->Draw("same");
  //************
  // ----- come back to this ----
  //c4->Divide(1,2);
  //c4->cd(1);
  //h_prof_pol_p_low->Draw();
  // ----- come back to this ----

  // ----- come back to this ----
  //TF1 *fitp_low = new TF1("fitp_low", "[0] + [1]*cos(x) + [2]*sin(x) + [3]*cos(2*x) + [4]*sin(2*x)+ [5]*cos(3*x) + [6]*sin(3*x)", -4.0, dx_p_cut);
  //fitp_low->SetParameters(1.0,2.0,3.0,4.0,5.0,6.0,7.0);
  //fitp_low->SetLineColor(kBlue);
  //h_prof_pol_p_low->Fit("fitp_low");
  //fitp_low->Draw("SAMES");

  //c4->cd(2);
  //h_prof_pol_p_high->Draw();

  //TF1 *fitp_high = new TF1("fitp_high", "[0]", dx_p_cut, 3.0);
  //fitp_high->SetParameters(0.05);
  //fitp_high->SetLineColor(kRed);
  //h_prof_pol_p_high->Fit("fitp_high");
  //fitp_high->Draw("SAMES");
  // ----- come back to this ----


  //TF1 *fitp = new TF1("fitp", "[0] + [1]*cos(x) + [2]*sin(x) + [3]*cos(2*x) + [4]*sin(2*x) + [5]*cos(3*x) + [6]*sin(3*x) + [7]*cos(4*x) + [8]*sin(4*x)", -3.5, 3.0); // + [5]*cos(3*x) + [6]*sin(3*x)
  //fitp->SetParameters(1,2,3,4,5,6,7,8,9);
  //fitp->SetLineColor(kRed);

  //h_prof_pol_p->Fit("fitp");
  //fitp->Draw("SAMES");
  //gStyle->SetOptFit(1111);

  //fitn_low->Write();
  //fitn_high->Write();
  //fitp_low->Write();
  //fitp_high->Write();

  h_dy->Write();

  fout->Write();
}
