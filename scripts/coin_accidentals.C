#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
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

void coin_accidentals(const char *kinematic)
{

  int kin = 2;

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  TString inputfile = Form("/volatile/halla/sbs/ktevans/QA/QE_data_%s_sbs100p_nucleon_np_model2.root",kinematic);
  TString outputfile = Form("plots/Accidentals_GEn_pass3_%s.pdf",kinematic);
  TString outfile = Form("outfiles/Accidentals_GEn_pass3_%s.root",kinematic);
  TFile *fout = new TFile(outfile,"RECREATE");

  TTree *T_data = new TTree("T_data", "Analysis Data Tree");

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

  double optics_valid_min;
  double optics_valid_max;
  double coin_mean;
  double coin_sigma;
  double dx_n_mean;
  double dx_n_sigma;
  double dx_p_mean;
  double dx_p_sigma;
  double dy_mean;
  double dy_sigma;
  int firstRun;
  int lastRun;
  int totRun;
  double Trp_max;
  double Trp_min;
  double pN_min;
  double pN_max;
  double W2_cut_min;
  double W2_cut_max;

  int runindex = 0;
  int runTrack = 0;

  int IHWP_flip;

  double Mp = 0.93827081;

  if(kin==2)
  {
    optics_valid_min = -0.35;
    optics_valid_max = 0.34;
    coin_mean = -0.978;
    coin_sigma = 2.4;
    IHWP_flip = -1;
    dx_n_mean = -0.121;
    dx_n_sigma = 0.551;
    dx_p_mean = -2.752;
    dx_p_sigma = 0.539;
    dy_mean = -0.066;
    dy_sigma = 0.569;
    firstRun = 2130;
    lastRun = 2322;
    totRun = 130;
    Trp_max = 3.5;
    Trp_min = 2.0;
    pN_min = 1.5;
    pN_max = 3.5;
    W2_cut_min = 0.0;
    W2_cut_max = 1.6;
  }
  else if(kin==3)
  {
    optics_valid_min = -0.35;
    optics_valid_max = 0.33;
    coin_mean = -0.382;
    coin_sigma = 2.2;
    IHWP_flip = 1;
    dx_n_mean = 0.0;
    dx_n_sigma = 0.5;
    dx_p_mean = -1.445;
    dx_p_sigma = 0.315;
    dy_mean = -0.057;
    dy_sigma = 0.374;
    firstRun = 2506;
    lastRun = 3250;
    totRun = 354;
    Trp_max = 3.5;
    Trp_min = 2.0;
    pN_min = 3.5;
    pN_max = 5.5;
    W2_cut_min = 0.0;
    W2_cut_max = 2.0;
  }
  else if(kin==4)
  {
    optics_valid_min = -0.36;
    optics_valid_max = 0.30;
    coin_mean = -0.589;
    coin_sigma = 2.3;
    IHWP_flip = 1;
    dx_n_mean = 0.0;
    dx_n_sigma = 0.5;
    dx_p_mean = -1.108;
    dx_p_sigma = 0.399;
    dy_mean = -0.064;
    dy_sigma = 0.487;
    firstRun = 3510;
    lastRun = 4587;
    totRun = 500;
    Trp_max = 4.0;
    Trp_min = 2.5;
    pN_min = 5.0;
    pN_max = 7.0;
    W2_cut_min = -1.0;
    W2_cut_max = 2.0;
  }
  else if(kin==5)
  {
    optics_valid_min = -0.37;
    optics_valid_max = 0.32;
    coin_mean = -0.568;
    coin_sigma = 2.2;
    IHWP_flip = 1;
    dx_n_mean = 0.0;
    dx_n_sigma = 0.5;
    dx_p_mean = -1.124;
    dx_p_sigma = 0.361;
    dy_mean = -0.071;
    dy_sigma = 0.456;
    firstRun = 5044;
    lastRun = 6083;
    totRun = 480;
    Trp_max = 4.0;
    Trp_min = 2.5;
    pN_min = 5.0;
    pN_max = 7.0;
    W2_cut_min = -1.0;
    W2_cut_max = 2.0;
  }
  else
  {
    optics_valid_min = -2.0;
    optics_valid_max = 2.0;
    coin_mean = 0.0;
    coin_sigma = 400.0;
    IHWP_flip = 1;
    dx_n_mean = 0.0;
    dx_n_sigma = 5.0;
    dx_p_mean = -1.0;
    dx_p_sigma = 5.0;
    dy_mean = 0.0;
    dy_sigma = 5.0;
    firstRun = 0;
    lastRun = 100;
    totRun = 500;
    Trp_max = 3.5;
    Trp_min = 2.0;
    pN_min = 0.0;
    pN_max = 10.0;
    W2_cut_min = -2.0;
    W2_cut_max = 8.0;
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

  TH1D* h_coin = new TH1D("h_coin", "coin", 320.0, -40.0, 40.0);
  h_coin->GetXaxis()->SetTitle("Coincidence Time [ns]");
  h_coin->SetTitle("Coin Time (HCal-BBCal) with Global, Vertex, E/p, PSe, GRINCH, and W2 Cuts");
  int QE_check = 0;
  double hodo_hcal_coin, hodo_sh_coin, hodo_ps_coin, Sf, measE, deltaEfrac, hcalPrimTot, hcalSecPrim;

  //Loop over all events to fill the histogram
  for (size_t iev = 0; iev < T->GetEntries(); iev++)
  {
    T->GetEntry(iev);

    QE_check = 0;

    if (runindex == 0)
    {
      runindex = runnum;
    }

    if (runindex != 0 && runindex !=runnum)
    {
      runindex = runnum;
      runTrack++;
    }

    if((IHWP==-1 || IHWP==1) && (bb_tr_r_x-0.9*bb_tr_r_th)>optics_valid_min && (bb_tr_r_x-0.9*bb_tr_r_th)<optics_valid_max && bb_gr_clus_track==0 && bb_ps_e>0.2 && abs(bb_tr_vz)<0.27 && abs(e_kine_W2-1.0)<0.5 && bb_gr_clus_size>2 && abs(((bb_ps_e+bb_sh_e)/bb_tr_p)-1)<0.2)
    {

      h_coin->Fill(adc_coin);

    }// end global cuts

  }//end event loop

  TCanvas *coin1 = new TCanvas("coin1", "coincidence", 1200, 1000);
  coin1->cd();
  coin1->SetLogy();
  h_coin->Draw();

  coin1->Print("("+outputfile);

  TCanvas *summary = new TCanvas("summary", "summary", 1200, 1000);
  summary->cd();
  TPaveText *pt = new TPaveText(.05,.1,.95,.8);
  pt->AddText("Global Cuts: ");
  pt->AddText(Form("bb.ps.e>0.0 && bb.gem.track.nhits>=3 && bb.gem.track.chi2ndf[0]<=15 && sbs.hcal.e>0.025 && %.1f<e.kine.W2<%.1f",W2_cut_min,W2_cut_max));
  pt->AddText("(IHWP==-1.0 || IHWP==1.0) && optics_valid_min<(bb_tr_r_x-0.9*bb_tr_r_th)<optics_valid_max && bb_gr_clus_track==0");
  pt->AddText("Vertex Cut: abs(bb.tr.vz)<0.27");
  pt->AddText("PreShower Cut: bb.ps.e>0.2");
  pt->AddText("E/p Cut: abs(E/p - 1.0)<0.2");
  pt->AddText("W2 Cut: abs(e.kine.W2 - 1.0)<0.5");
  pt->AddText(Form("Coincidence Cut: abs(adc.coin - %.3f)<%.3f",coin_mean,coin_sigma));
  pt->AddText("GRINCH cut: bb.grinch_tdc.clus.size>2");
  pt->AddText(Form("Proton Spot Cut: %.3f<dx<%.3f && %.3f<dy<%.3f",(dx_p_mean-dx_p_sigma),(dx_p_mean+dx_p_sigma),(dy_mean-dy_sigma),(dy_mean-dy_sigma)));
  pt->AddText(Form("Neutron Spot Cut: %.3f<dx<%.3f && %.3f<dy<%.3f",(dx_n_mean-dx_n_sigma),(dx_n_mean+dx_n_sigma),(dy_mean-dy_sigma),(dy_mean-dy_sigma)));
  pt->Draw();

  summary->Print(outputfile+")");

  T_data->Write("", TObject::kOverwrite);

  h_coin->Write();

  fout->Write();
  //T_data->Delete();
  //T->Delete();
}
