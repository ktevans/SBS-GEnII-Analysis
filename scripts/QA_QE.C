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

void QA_QE(const char *kinematic)
{

  int kin = 3;

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  TString inputfile = Form("/volatile/halla/sbs/ktevans/QA/QE_data_%s_sbs100p_nucleon_np_model2.root",kinematic);
  //TString inputfile = Form("/volatile/halla/sbs/ktevans/KateJackSBSAnalysis/KJ_parsed_GEn_pass2_%s_He3_100.root",kinematic);
  TString outputfile = Form("plots/QA_parsed_GEn_pass3_%s_He3_dxdy.pdf",kinematic);
  TString outfile = Form("outfiles/QA_parsed_GEn_pass3_%s_He3_dxdy.root",kinematic);
  TFile *fout = new TFile(outfile,"RECREATE");

  TTree *T_data = new TTree("T_data", "Analysis Data Tree");

  //double dx_out, dy_out, W2_out;
  //int helicity_out;
  //T_data->Branch("dx", &dx_out, "dx/D");
  //T_data->Branch("dy", &dy_out, "dy/D");
  //T_data->Branch("W2", &W2_out, "W2/D");
  //T_data->Branch("helicity", &helicity_out, "helicity/I");

  //TChain* T = new TChain("Parse");
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
  //Double_t IHWP;                  T->SetBranchAddress("IHWP", &IHWP);
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
    dx_n_mean = -0.147;
    dx_n_sigma = 0.812;
    dx_p_mean = -2.811;
    dx_p_sigma = 0.600;
    dy_mean = -0.096;
    dy_sigma = 0.580;
    firstRun = 2130;
    lastRun = 2322;
    totRun = 130;
    Trp_max = 3.5;
    Trp_min = 2.0;
    pN_min = 0.0;
    pN_max = 6.0;
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
    dy_mean = 0.254;
    dy_sigma = 0.773;
    firstRun = 3510;
    lastRun = 4587;
    totRun = 500;
    Trp_max = 4.0;
    Trp_min = 2.5;
    pN_min = 4.0;
    pN_max = 8.0;
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
    dy_mean = 0.248;
    dy_sigma = 0.769;
    firstRun = 5044;
    lastRun = 6083;
    totRun = 480;
    Trp_max = 4.0;
    Trp_min = 2.5;
    pN_min = 4.0;
    pN_max = 8.0;
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

  // ~~~~~~~~~~~~~~~~~~~~ Basic QE cut plots ~~~~~~~~~~~~~~~~~~~~

  TH1D* h_tr_vz = new TH1D("h_tr_vz", "Track z Vertex", 140.0, -0.8, 0.6);
  h_tr_vz->GetXaxis()->SetTitle("bb.tr.vz [m]");
  h_tr_vz->SetTitle("Track z Vertex with Global Cuts");

  TH1D* h_ps_e_raw = new TH1D("h_ps_e_raw", "PreShower Energy", 100.0, 0.0, 2.0);
  h_ps_e_raw->GetXaxis()->SetTitle("bb.ps.e [GeV]");
  h_ps_e_raw->SetTitle("PreShower Energy with Global and Vertex Cuts");

  TH2D* h2_coin_W2 = new TH2D("h2_coin_W2", "Coin vs W2", 100.0, -2.0, 8.0, 80.0, -20.0, 20.0);
  h2_coin_W2->GetXaxis()->SetTitle("e.kine.W2 [GeV]");
  h2_coin_W2->GetYaxis()->SetTitle("adc.coin [ns]");
  h2_coin_W2->SetTitle("Coincidence Time (HCal-BBCal) vs W2 with Global and Vertex Cuts");

  TH2D* h2_pse_grclus = new TH2D("h2_pse_grclus", "PSe vs GRclus", 100.0, 0.0, 2.0, 20.0, 0.0, 20.0);
  h2_pse_grclus->GetXaxis()->SetTitle("bb.ps.e [GeV]");
  h2_pse_grclus->GetYaxis()->SetTitle("bb.gr.clus.size");
  h2_pse_grclus->SetTitle("PreShower Energy vs GRINCH Cluster Size with Global, Vertex, E/p, and PSe Cuts");

  TH1D* h_W2 = new TH1D("h_W2", "W2", 100.0, 0.0, 2.0);
  h_W2->GetXaxis()->SetTitle("e.kine.W2 [GeV]");
  h_W2->SetTitle("W2 with Global, Vertex, E/p, PSe, Coin, and GRINCH Cuts");

  TH1D* h_W2_raw = new TH1D("h_W2_raw", "W2", 100.0, 0.0, 2.0);
  h_W2_raw->GetXaxis()->SetTitle("e.kine.W2 [GeV]");
  h_W2_raw->SetTitle("W2 with Global Cuts");
  h_W2_raw->SetLineColor(kRed);

  TH1D* h_dx = new TH1D("h_dx", ";dx", 100.0, -6.0, 4.0);
  h_dx->GetXaxis()->SetTitle("dx [m]");
  h_dx->SetTitle("dx with Global, Vertex, E/p, PSe, Coin, GRINCH, and W2 Cuts");

  TH1D* h_dy = new TH1D("h_dy", ";dy", 70.0, -3.0, 4.0);
  h_dy->GetXaxis()->SetTitle("dy [m]");
  h_dy->SetTitle("dy with Global, Vertex, E/p, PSe, Coin, GRINCH, and W2 Cuts");

  TH2D* h2_dxdy = new TH2D("h2_dxdy", "dy vs dx", 70.0, -3.0, 4.0, 100.0, -6.0, 4.0);
  h2_dxdy->GetXaxis()->SetTitle("dy [m]");
  h2_dxdy->GetYaxis()->SetTitle("dx [m]");
  h2_dxdy->SetTitle("dx vs dy with Global, Vertex, E/p, PSe, Coin, GRINCH, and W2 Cuts");

  // ~~~~~~~~~~~~~~~~~~~~ BBCal plots ~~~~~~~~~~~~~~~~~~~~

  TH1D* h_eovp = new TH1D("h_eovp", "E/p", 100.0, 0.0, 2.0);
  h_eovp->GetXaxis()->SetTitle("E/p");
  h_eovp->SetTitle("E/p with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_eovp_runnum = new TH2D("h2_eovp_runnum", "E/p vs runnum", totRun, 0, totRun, 100.0, 0.8, 1.2);
  h2_eovp_runnum->GetYaxis()->SetTitle("E/p");
  h2_eovp_runnum->GetXaxis()->SetTitle("runnum");
  h2_eovp_runnum->SetTitle("E/p vs Run Number with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TProfile* h2_eovp_runnum_prof = new TProfile("h2_eovp_runnum_prof", "Prof E/p", totRun, 0, totRun, 0.8, 1.2);
  h2_eovp_runnum_prof->SetMarkerColor(kRed);
  h2_eovp_runnum_prof->SetMarkerStyle(20);
  h2_eovp_runnum_prof->SetMarkerSize(2);

  TH2D* h2_eovp_trx = new TH2D("h2_eovp_trx", "E/p vs trx", 100, -0.45, 0.6, 100.0, 0.5, 1.5);
  h2_eovp_trx->GetYaxis()->SetTitle("E/p");
  h2_eovp_trx->GetXaxis()->SetTitle("bb.tr.x [m]");
  h2_eovp_trx->SetTitle("E/p vs Track x with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_eovp_try = new TH2D("h2_eovp_try", "E/p vs try", 100, -0.2, 0.15, 100.0, 0.5, 1.5);
  h2_eovp_try->GetYaxis()->SetTitle("E/p");
  h2_eovp_try->GetXaxis()->SetTitle("bb.tr.y [m]");
  h2_eovp_try->SetTitle("E/p vs Track y with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_ps_tot = new TH2D("h2_ps_tot", "ePS vs ePS+eSH", 100, 0.0, 2.0, 100.0, 1.5, 4.5);
  h2_ps_tot->GetXaxis()->SetTitle("bb.ps.e [GeV]");
  h2_ps_tot->GetYaxis()->SetTitle("bb.ps.e+bb.sh.e [GeV]");
  h2_ps_tot->SetTitle("PreShower Energy vs Total BBCal Energy with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH1D* h_ps_e = new TH1D("h_ps_e", "PreShower Energy", 100.0, 0.0, 2.0);
  h_ps_e->GetXaxis()->SetTitle("bb.ps.e [GeV]");
  h_ps_e->SetTitle("PreShower Energy with Global, Vertex, and W2 Cuts");

  TH1D* h_ps_e_anti = new TH1D("h_ps_e_anti", "PreShower Energy (anti-QE cuts)", 100.0, 0.0, 2.0);
  h_ps_e_anti->GetXaxis()->SetTitle("bb.ps.e [GeV]");
  h_ps_e_anti->SetTitle("PreShower Energy with Global, Vertex, and W2 Cuts and Anti - (E/p, Coin, GRINCH and Spot) Cuts");
  h_ps_e_anti->SetLineColor(kRed);

  TH1D* h_ps_e_qe = new TH1D("h_ps_e_qe", "PreShower Energy (QE cuts)", 100.0, 0.0, 2.0);
  h_ps_e_qe->GetXaxis()->SetTitle("bb.ps.e [GeV]");
  h_ps_e_qe->SetTitle("PreShower Energy with Global, Vertex, E/p, Coin, GRINCH, W2, and Spot Cuts");
  h_ps_e_qe->SetLineColor(kBlue);

  TH2D* h2_pse_trx = new TH2D("h2_pse_trx", "PSe vs TrX", 100.0, -0.45, 0.6, 100.0, 0.0, 2.0);
  h2_pse_trx->GetYaxis()->SetTitle("bb.ps.e [GeV]");
  h2_pse_trx->GetXaxis()->SetTitle("bb.tr.x [m]");
  h2_pse_trx->SetTitle("PreShower Energy vs Track x with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_pse_try = new TH2D("h2_pse_try", "PSe vs TrY", 100.0, -0.2, 0.15, 100.0, 0.0, 2.0);
  h2_pse_try->GetYaxis()->SetTitle("bb.ps.e [GeV]");
  h2_pse_try->GetXaxis()->SetTitle("bb.tr.y [m]");
  h2_pse_try->SetTitle("PreShower Energy vs Track y with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_pse_pstime = new TH2D("h2_pse_pstime", "PSe vs PStime", 100.0, 0.0, 2.0, 200.0, -10.0, 25.0);
  h2_pse_pstime->GetXaxis()->SetTitle("bb.ps.e [GeV]");
  h2_pse_pstime->GetYaxis()->SetTitle("bb.ps.atimeblk [ns]");
  h2_pse_pstime->SetTitle("PreShower Energy vs PreShower ADC Time with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_she_trx = new TH2D("h2_she_trx", "SHe vs TrX", 100.0, -0.45, 0.6, 100.0, 0.0, 3.0);
  h2_she_trx->GetYaxis()->SetTitle("bb.sh.e [GeV]");
  h2_she_trx->GetXaxis()->SetTitle("bb.tr.x [m]");
  h2_she_trx->SetTitle("Shower Energy vs Track x with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_she_try = new TH2D("h2_she_try", "SHe vs TrY", 100.0, -0.2, 0.15, 100.0, 0.0, 3.0);
  h2_she_try->GetYaxis()->SetTitle("bb.sh.e [GeV]");
  h2_she_try->GetXaxis()->SetTitle("bb.tr.y [m]");
  h2_she_try->SetTitle("Shower Energy vs Track y with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_she_shtime = new TH2D("h2_she_shtime", "SHe vs SHtime", 100.0, 0.0, 3.0, 200.0, -10.0, 25.0);
  h2_she_shtime->GetXaxis()->SetTitle("bb.sh.e [GeV]");
  h2_she_shtime->GetYaxis()->SetTitle("bb.sh.atimeblk [ns]");
  h2_she_shtime->SetTitle("Shower Energy vs PreShower ADC Time with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH1D* h_sh_e = new TH1D("h_sh_e", "Shower Energy", 100.0, 0.0, 3.0);
  h_sh_e->GetXaxis()->SetTitle("bb.sh.e [GeV]");
  h_sh_e->SetTitle("Shower Energy with Global, Vertex, and W2 Cuts");

  TH1D* h_sh_e_anti = new TH1D("h_sh_e_anti", "Shower Energy (anti-QE cuts)", 100.0, 0.0, 3.0);
  h_sh_e_anti->GetXaxis()->SetTitle("bb.sh.e [GeV]");
  h_sh_e_anti->SetTitle("Shower Energy with Global, Vertex, and W2 Cuts and Anti - (E/p, Coin, GRINCH, and Spot) Cuts");
  h_sh_e_anti->SetLineColor(kRed);

  TH1D* h_sh_e_qe = new TH1D("h_sh_e_qe", "Shower Energy (QE cuts)", 100.0, 0.0, 3.0);
  h_sh_e_qe->GetXaxis()->SetTitle("bb.sh.e [GeV]");
  h_sh_e_qe->SetTitle("Shower Energy with Global, Vertex, E/p, Coin, GRINCH, W2, and Spot Cuts");
  h_sh_e_qe->SetLineColor(kBlue);

  TH2D* h2_ps_xy = new TH2D("h2_ps_xy", "PreShower Position", 2.0, -0.26, 0.26, 26.0, -1.1, 1.1);
  h2_ps_xy->GetXaxis()->SetTitle("bb.ps.y [m]");
  h2_ps_xy->GetYaxis()->SetTitle("bb.ps.x [m]");
  h2_ps_xy->SetTitle("PreShower Position with Global, Vertex, E/p, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_sh_xy = new TH2D("h2_sh_xy", "Shower Position", 7.0, -0.26, 0.26, 27.0, -1.1, 1.1);
  h2_sh_xy->GetXaxis()->SetTitle("bb.sh.y [m]");
  h2_sh_xy->GetYaxis()->SetTitle("bb.sh.x [m]");
  h2_sh_xy->SetTitle("Shower Position with Global, Vertex, E/p, Coin, GRINCH, W2, and Spot Cuts");

  // ~~~~~~~~~~~~~~~~~~~~ HCal plots ~~~~~~~~~~~~~~~~~~~~

  TH2D* h2_hcal_e_x = new TH2D("h2_hcal_e_x", "HCal energy vs HCal x", 24.0, -2.6, 1.1, 200.0, 0.0, 1.5);
  h2_hcal_e_x->GetXaxis()->SetTitle("sbs.hcal.x [m]");
  h2_hcal_e_x->GetYaxis()->SetTitle("sbs.hcal.e [GeV]");
  h2_hcal_e_x->SetTitle("HCal Energy vs HCal x with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_hcal_e_y = new TH2D("h2_hcal_e_y", "HCal energy vs HCal y", 12.0, -0.9, 0.9, 200.0, 0.0, 1.5);
  h2_hcal_e_y->GetXaxis()->SetTitle("sbs.hcal.y [m]");
  h2_hcal_e_y->GetYaxis()->SetTitle("sbs.hcal.e [GeV]");
  h2_hcal_e_y->SetTitle("HCal Energy vs HCal y with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_hcal_time_e = new TH2D("h2_hcal_time_e", "HCal energy vs HCal ADC time", 100.0, -10.0, 10.0, 200.0, 0.0, 1.5);
  h2_hcal_time_e->GetXaxis()->SetTitle("sbs.hcal.atimeblk [ns]");
  h2_hcal_time_e->GetYaxis()->SetTitle("sbs.hcal.e [GeV]");
  h2_hcal_time_e->SetTitle("HCal Energy vs HCal ADC Time with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_hodo_time_hcal_e = new TH2D("h2_hodo_time_hcal_e", "HCal energy vs Hodo Mean time", 100.0, -10.0, 10.0, 200.0, 0.0, 1.5);
  h2_hodo_time_hcal_e->GetXaxis()->SetTitle("bb.hodotdc.clus.tmean [ns]");
  h2_hodo_time_hcal_e->GetYaxis()->SetTitle("sbs.hcal.e [GeV]");
  h2_hodo_time_hcal_e->SetTitle("HCal Energy vs Hodo Mean Time with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH1D* h_pN = new TH1D("h_pN", "Expected HCal Momentum", 100.0, pN_min, pN_max);
  h_pN->GetXaxis()->SetTitle("Expected Scattered Nucleon Momentum [GeV]");
  h_pN->SetTitle("Expected Nucleon Momentum with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH1D* h_Sf = new TH1D("h_Sf", "HCal Sampling Fraction", 100.0, 0.0, 0.5);
  h_Sf->GetXaxis()->SetTitle("Sampling Fraction (sbs.hcal.e*2*Mp/Q2)");
  h_Sf->SetTitle("Sampling Fraction with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH1D* h_measE = new TH1D("h_measE", "HCal Measured Fraction", 100.0, 0.0, 0.5);
  h_measE->GetXaxis()->SetTitle("Measured Energy / Actual Energy (sbs.hcal.e/(sqrt(Mp*Mp+pN*pN)-Mp)");
  h_measE->SetTitle("Measured Energy Fraction with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_delta_Sf_measE = new TH2D("h2_delta_Sf_measE", "Measured Energy Fraction vs Sampling Fraction", 100.0, 0.0, 0.5, 100.0, 0.0, 0.5);
  h2_delta_Sf_measE->GetXaxis()->SetTitle("Sampling Fraction");
  h2_delta_Sf_measE->GetYaxis()->SetTitle("Measured Energy Fraction");
  h2_delta_Sf_measE->SetTitle("Measured Energy Fraction vs Sampling Fraction with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH1D* h_delta_Sf_measE = new TH1D("h_delta_Sf_measE", "Measured Energy Fraction - Sampling Fraction", 100.0, -0.2, 0.2);
  h_delta_Sf_measE->GetXaxis()->SetTitle("Measured Energy Fraction - Sampling Fraction");
  h_delta_Sf_measE->SetTitle("Measured Energy Fraction - Sampling Fraction with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_Sf_x = new TH2D("h2_Sf_x", "Sampling Fraction vs HCal x", 24.0, -2.6, 1.1, 100.0, 0.0, 0.5);
  h2_Sf_x->GetXaxis()->SetTitle("sbs.hcal.x [m]");
  h2_Sf_x->GetYaxis()->SetTitle("Sampling Fraction");
  h2_Sf_x->SetTitle("Sampling Fraction vs HCal x with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_Sf_y = new TH2D("h2_Sf_y", "Sampling Fraction vs HCal y", 12.0, -0.9, 0.9, 100.0, 0.0, 0.5);
  h2_Sf_y->GetXaxis()->SetTitle("sbs.hcal.y [m]");
  h2_Sf_y->GetYaxis()->SetTitle("Sampling Fraction");
  h2_Sf_y->SetTitle("Sampling Fraction vs HCal y with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_hcal_xy = new TH2D("h2_hcal_xy", "Shower Position", 12.0, -0.85, 0.85, 24.0, -2.6, 1.1);
  h2_hcal_xy->GetXaxis()->SetTitle("sbs.hcal.y [m]");
  h2_hcal_xy->GetYaxis()->SetTitle("sbs.hcal.x [m]");
  h2_hcal_xy->SetTitle("HCal Position with Global, Vertex, E/p, Coin, GRINCH, W2, and Spot Cuts");

  TH1D* h_hcal_prim_tot_e = new TH1D("h_hcal_prim_tot_e", "HCal Prim Energy Fraction", 100.0, 0.0, 1.1);
  h_hcal_prim_tot_e->GetXaxis()->SetTitle("sbs.hcal.clusblk.e[0] / sbs.hcal.e");
  h_hcal_prim_tot_e->SetTitle("Energy of Primary HCal Block / Total Cluster Energy with Global, Vertex, E/p, Coin, GRINCH, W2, and Spot Cuts");

  TH1D* h_hcal_sec_prim_e = new TH1D("h_hcal_sec_prim_e", "HCal Prim Sec Energy Fraction", 100.0, 0.0, 1.1);
  h_hcal_sec_prim_e->GetXaxis()->SetTitle("sbs.hcal.clusblk.e[1] / sbs.hcal.clusblk.e[0]");
  h_hcal_sec_prim_e->SetTitle("Energy of Secondary HCal Block / Primary HCal Block with Global, Vertex, E/p, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_hcal_primFrac_blkID = new  TH2D("h2_hcal_primFrac_blkID", "HCal Frac per Block", 288.0, 0.0, 288.0, 100.0, 0.0, 1.0);
  h2_hcal_primFrac_blkID->GetXaxis()->SetTitle("sbs.hcal.clusblk.id[0]");
  h2_hcal_primFrac_blkID->GetYaxis()->SetTitle("sbs.hcal.clusblk.e[0] / sbs.hcal.e");
  h2_hcal_primFrac_blkID->SetTitle("Energy of Primary HCal Block / Total Cluster Energy vs Primary Block ID with Global, Vertex, E/p, Coin, GRINCH, W2, and Spot Cuts");

  TProfile* prof_hcal_primFrac_blkID = new TProfile("prof_hcal_primFrac_blkID", "Prof primFracBLID", 288.0, 0.0, 288.0, 0.0, 1.0);
  prof_hcal_primFrac_blkID->SetMarkerColor(kRed);
  prof_hcal_primFrac_blkID->SetMarkerStyle(20);
  prof_hcal_primFrac_blkID->SetMarkerSize(2);

  TH2D* h2_hcal_primFrac_x = new  TH2D("h2_hcal_primFrac_x", "HCal Frac per X", 24.0, -2.6, 1.1, 100.0, 0.0, 1.0);
  h2_hcal_primFrac_x->GetXaxis()->SetTitle("sbs.hcal.x [m]");
  h2_hcal_primFrac_x->GetYaxis()->SetTitle("sbs.hcal.clusblk.e[0] / sbs.hcal.e");
  h2_hcal_primFrac_x->SetTitle("Energy of Primary HCal Block / Total Cluster Energy vs x with Global, Vertex, E/p, Coin, GRINCH, W2, and Spot Cuts");

  TProfile* prof_hcal_primFrac_x = new TProfile("prof_hcal_primFrac_x", "Prof primFracX", 24.0, -2.6, 1.1, 0.0, 1.0);
  prof_hcal_primFrac_x->SetMarkerColor(kRed);
  prof_hcal_primFrac_x->SetMarkerStyle(20);
  prof_hcal_primFrac_x->SetMarkerSize(2);

  TH2D* h2_hcal_primFrac_y = new  TH2D("h2_hcal_primFrac_y", "HCal Frac per y", 12.0, -0.85, 0.85, 100.0, 0.0, 1.0);
  h2_hcal_primFrac_y->GetXaxis()->SetTitle("sbs.hcal.y [m]");
  h2_hcal_primFrac_y->GetYaxis()->SetTitle("sbs.hcal.clusblk.e[0] / sbs.hcal.e");
  h2_hcal_primFrac_y->SetTitle("Energy of Primary HCal Block / Total Cluster Energy vs y with Global, Vertex, E/p, Coin, GRINCH, W2, and Spot Cuts");

  TProfile* prof_hcal_primFrac_y = new TProfile("prof_hcal_primFrac_y", "Prof primFracY", 12.0, -0.85, 0.85, 0.0, 1.0);
  prof_hcal_primFrac_y->SetMarkerColor(kRed);
  prof_hcal_primFrac_y->SetMarkerStyle(20);
  prof_hcal_primFrac_y->SetMarkerSize(2);

  // ~~~~~~~~~~~~~~~~~~~~ GEM plots ~~~~~~~~~~~~~~~~~~~~

  TH2D* h2_trp_trx = new TH2D("h2_trp_trx", "TrP vs TrX", 100.0, -0.45, 0.6, 100.0, Trp_min, Trp_max);
  h2_trp_trx->GetYaxis()->SetTitle("bb.tr.p [GeV]");
  h2_trp_trx->GetXaxis()->SetTitle("bb.tr.x [m]");
  h2_trp_trx->SetTitle("Track p vs Track x with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_trp_try = new TH2D("h2_trp_try", "TrP vs TrY", 100.0, -0.2, 0.15, 100.0, Trp_min, Trp_max);
  h2_trp_try->GetYaxis()->SetTitle("bb.tr.p [GeV]");
  h2_trp_try->GetXaxis()->SetTitle("bb.tr.y [m]");
  h2_trp_try->SetTitle("Track p vs Track y with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_trp_bbcalE = new TH2D("h2_trp_bbcalE", "TrP vs BBCalE", 100.0, Trp_min, Trp_max, 100.0, 1.5, 4.5);
  h2_trp_bbcalE->GetXaxis()->SetTitle("bb.tr.p [GeV]");
  h2_trp_bbcalE->GetYaxis()->SetTitle("bb.ps.e+bb.sh.e [GeV]");
  h2_trp_bbcalE->SetTitle("Track p vs Track y with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  // ~~~~~~~~~~~~~~~~~~~~ Timing plots ~~~~~~~~~~~~~~~~~~~~

  TH1D* h_coin = new TH1D("h_coin", "coin", 200.0, -10.0, 10.0);
  h_coin->GetXaxis()->SetTitle("Coincidence Time [ns]");
  h_coin->SetTitle("Coin Time (HCal-BBCal) with Global, Vertex, E/p, PSe, GRINCH, W2, and Spot Cuts");

  TH2D* h2_hcal_atime_x = new TH2D("h2_hcal_atime_x", "HCal time vs HCal x", 24.0, -2.6, 1.1, 200.0, -10.0, 25.0);
  h2_hcal_atime_x->GetXaxis()->SetTitle("sbs.hcal.x [m]");
  h2_hcal_atime_x->GetYaxis()->SetTitle("sbs.hcal.atimeblk [ns]");
  h2_hcal_atime_x->SetTitle("HCal ADC Time vs HCal x with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TProfile* prof_hcal_atime_x = new TProfile("prof_hcal_atime_x", "Prof hcalTimeX", 24.0, -2.6, 1.1, -10.0, 25.0);
  prof_hcal_atime_x->SetMarkerColor(kRed);
  prof_hcal_atime_x->SetMarkerStyle(20);
  prof_hcal_atime_x->SetMarkerSize(2);

  TH2D* h2_hcal_atime_y = new TH2D("h2_hcal_atime_y", "HCal time vs HCal y", 12.0, -0.85, 0.85, 200.0, -10.0, 25.0);
  h2_hcal_atime_y->GetXaxis()->SetTitle("sbs.hcal.y [m]");
  h2_hcal_atime_y->GetYaxis()->SetTitle("sbs.hcal.atimeblk [ns]");
  h2_hcal_atime_y->SetTitle("HCal ADC Time vs HCal y with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TProfile* prof_hcal_atime_y = new TProfile("prof_hcal_atime_y", "Prof hcalTimeY", 12.0, -0.85, 0.85, -10.0, 25.0);
  prof_hcal_atime_y->SetMarkerColor(kRed);
  prof_hcal_atime_y->SetMarkerStyle(20);
  prof_hcal_atime_y->SetMarkerSize(2);

  TH2D* h2_ps_atime_x = new TH2D("h2_ps_atime_x", "PS Time vs Tr x", 50.0, -0.45, 0.6, 200.0, -10.0, 25.0);
  h2_ps_atime_x->GetXaxis()->SetTitle("bb.tr.x [m]");
  h2_ps_atime_x->GetYaxis()->SetTitle("bb.ps.atimeblk [ns]");
  h2_ps_atime_x->SetTitle("PS ADC Time vs Track x with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TProfile* prof_ps_atime_x = new TProfile("prof_ps_atime_x", "Prof PSTimeX", 50.0, -0.45, 0.6, -10.0, 25.0);
  prof_ps_atime_x->SetMarkerColor(kRed);
  prof_ps_atime_x->SetMarkerStyle(20);
  prof_ps_atime_x->SetMarkerSize(2);

  TH2D* h2_ps_atime_y = new TH2D("h2_ps_atime_y", "PS Time vs Tr y", 50.0, -0.2, 0.2, 200.0, -10.0, 25.0);
  h2_ps_atime_y->GetXaxis()->SetTitle("bb.tr.y [m]");
  h2_ps_atime_y->GetYaxis()->SetTitle("bb.ps.atimeblk [ns]");
  h2_ps_atime_y->SetTitle("PS ADC Time vs Track y with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TProfile* prof_ps_atime_y = new TProfile("prof_ps_atime_y", "Prof PSTimeY", 50.0, -0.2, 0.2, -10.0, 25.0);
  prof_ps_atime_y->SetMarkerColor(kRed);
  prof_ps_atime_y->SetMarkerStyle(20);
  prof_ps_atime_y->SetMarkerSize(2);

  TH2D* h2_sh_atime_x = new TH2D("h2_sh_atime_x", "SH Time vs Tr x", 50.0, -0.45, 0.6, 200.0, -10.0, 25.0);
  h2_sh_atime_x->GetXaxis()->SetTitle("bb.tr.x [m]");
  h2_sh_atime_x->GetYaxis()->SetTitle("bb.sh.atimeblk [ns]");
  h2_sh_atime_x->SetTitle("SH ADC Time vs Track x with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TProfile* prof_sh_atime_x = new TProfile("prof_sh_atime_x", "Prof SHTimeX", 50.0, -0.45, 0.6, -10.0, 25.0);
  prof_sh_atime_x->SetMarkerColor(kRed);
  prof_sh_atime_x->SetMarkerStyle(20);
  prof_sh_atime_x->SetMarkerSize(2);

  TH2D* h2_sh_atime_y = new TH2D("h2_sh_atime_y", "SH Time vs Tr y", 50.0, -0.2, 0.15, 200.0, -10.0, 25.0);
  h2_sh_atime_y->GetXaxis()->SetTitle("bb.tr.y [m]");
  h2_sh_atime_y->GetYaxis()->SetTitle("bb.sh.atimeblk [ns]");
  h2_sh_atime_y->SetTitle("SH ADC Time vs Track y with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TProfile* prof_sh_atime_y = new TProfile("prof_sh_atime_y", "Prof SHTimeY", 50.0, -0.2, 0.15, -10.0, 25.0);
  prof_sh_atime_y->SetMarkerColor(kRed);
  prof_sh_atime_y->SetMarkerStyle(20);
  prof_sh_atime_y->SetMarkerSize(2);

  TH2D* h2_hcalTH_atime_x = new TH2D("h2_hcalTH_atime_x", "(TH-HCal) Time vs Tr x", 37.0, -2.6, 1.1, 200.0, -20.0, 15.0);
  h2_hcalTH_atime_x->GetXaxis()->SetTitle("bb.tr.x [m]");
  h2_hcalTH_atime_x->GetYaxis()->SetTitle("bb.hodotdc.clus.tmean - sbs.hcal.atimeblk [ns]");
  h2_hcalTH_atime_x->SetTitle("(TH-HCal) Time vs Track x with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_psTH_atime_x = new TH2D("h2_psTH_atime_x", "(TH-PS) Time vs Tr x", 105.0, -0.45, 0.6, 200.0, -20.0, 15.0);
  h2_psTH_atime_x->GetXaxis()->SetTitle("bb.tr.x [m]");
  h2_psTH_atime_x->GetYaxis()->SetTitle("bb.hodotdc.clus.tmean - bb.ps.atimeblk [ns]");
  h2_psTH_atime_x->SetTitle("(TH-PS) Time vs Track x with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_shTH_atime_x = new TH2D("h2_shTH_atime_x", "(TH-SH) Time vs Tr x", 105.0, -0.45, 0.6, 200.0, -20.0, 15.0);
  h2_shTH_atime_x->GetXaxis()->SetTitle("bb.tr.x [m]");
  h2_shTH_atime_x->GetYaxis()->SetTitle("bb.hodotdc.clus.tmean - bb.sh.atimeblk [ns]");
  h2_shTH_atime_x->SetTitle("(TH-SH) Time vs Track x with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_hcalTH_atime_y = new TH2D("h2_hcalTH_atime_y", "(TH-HCal) Time vs Tr y", 18.0, -0.9, 0.9, 200.0, -20.0, 15.0);
  h2_hcalTH_atime_y->GetXaxis()->SetTitle("bb.tr.y [m]");
  h2_hcalTH_atime_y->GetYaxis()->SetTitle("bb.hodotdc.clus.tmean - sbs.hcal.atimeblk [ns]");
  h2_hcalTH_atime_y->SetTitle("(TH-HCal) Time vs Track y with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_psTH_atime_y = new TH2D("h2_psTH_atime_y", "(TH-PS) Time vs Tr y", 35.0, -0.2, 0.15, 200.0, -20.0, 15.0);
  h2_psTH_atime_y->GetXaxis()->SetTitle("bb.tr.y [m]");
  h2_psTH_atime_y->GetYaxis()->SetTitle("bb.hodotdc.clus.tmean - bb.ps.atimeblk [ns]");
  h2_psTH_atime_y->SetTitle("(TH-PS) Time vs Track y with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_shTH_atime_y = new TH2D("h2_shTH_atime_y", "(TH-SH) Time vs Tr y", 35.0, -0.2, 0.15, 200.0, -20.0, 15.0);
  h2_shTH_atime_y->GetXaxis()->SetTitle("bb.tr.y [m]");
  h2_shTH_atime_y->GetYaxis()->SetTitle("bb.hodotdc.clus.tmean - bb.sh.atimeblk [ns]");
  h2_shTH_atime_y->SetTitle("(TH-SH) Time vs Track y with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

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

    // && abs(adc_coin-coin_mean)<coin_sigma && bb_ps_e>0.2 && abs(((bb_ps_e+bb_sh_e)/bb_tr_p)-1)<0.2 && bb_gr_clus_size>2.0 && abs(bb_tr_vz)<0.27
    // pass_global==1 && (IHWP==-1.0 || IHWP==1.0)

    if((IHWP==-1 || IHWP==1) && (bb_tr_r_x-0.9*bb_tr_r_th)>optics_valid_min && (bb_tr_r_x-0.9*bb_tr_r_th)<optics_valid_max && bb_gr_clus_track==0 && bb_ps_e>0.0)
    {

      h_tr_vz->Fill(bb_tr_vz);
      h_W2_raw->Fill(e_kine_W2);

      if (abs(bb_tr_vz)<0.27)
      {

        h_ps_e_raw->Fill(bb_ps_e);

        if (abs(e_kine_W2-1.0)<0.5)
        {
          h_ps_e->Fill(bb_ps_e);
          h_sh_e->Fill(bb_sh_e);

          if (abs(adc_coin-coin_mean)<(coin_sigma) && bb_gr_clus_size>2 && ((abs(dy_hcal-dy_mean)<dy_sigma && abs(dx_hcal-dx_n_mean)<dx_n_sigma) || (abs(dy_hcal-dy_mean)<dy_sigma && abs(dx_hcal-dx_p_mean)<dx_p_sigma)) && abs(((bb_ps_e+bb_sh_e)/bb_tr_p)-1)<0.2)
          {
            h_ps_e_qe->Fill(bb_ps_e);
            h_sh_e_qe->Fill(bb_sh_e);
            QE_check = 1;
          } // end QE cuts

          if (QE_check==0)
          {
            h_ps_e_anti->Fill(bb_ps_e);
            h_sh_e_anti->Fill(bb_sh_e);
          } // en anti-QE check
        }

        if(abs(e_kine_W2-1.0)<0.5 && bb_gr_clus_size>2 && ((abs(dy_hcal-dy_mean)<dy_sigma && abs(dx_hcal-dx_n_mean)<dx_n_sigma) || (abs(dy_hcal-dy_mean)<dy_sigma && abs(dx_hcal-dx_p_mean)<dx_p_sigma)) && abs(((bb_ps_e+bb_sh_e)/bb_tr_p)-1)<0.2)
        {
          h_coin->Fill(adc_coin);
        }

        h2_coin_W2->Fill(e_kine_W2,adc_coin);


        if (bb_ps_e>0.2)
        {

          if (abs(((bb_ps_e+bb_sh_e)/bb_tr_p)-1)<0.2)
          {

            h2_pse_grclus->Fill(bb_ps_e,bb_gr_clus_size);

            if (abs(adc_coin-coin_mean)<(coin_sigma) && bb_gr_clus_size>2)
            {
              h_W2->Fill(e_kine_W2);

              if (abs(e_kine_W2-1.0)<0.5)
              {
                h_dx->Fill(dx_hcal);
                h_dy->Fill(dy_hcal);
                h2_dxdy->Fill(dy_hcal,dx_hcal);

                //dx_out = dx_hcal;
                //dy_out = dy_hcal;

                //T_data->Fill();

                if (abs(dy_hcal-dy_mean)<dy_sigma && ((abs(dx_hcal-dx_n_mean)<dx_n_sigma)||(abs(dx_hcal-dx_p_mean)<dx_p_sigma)))
                {

                  h_eovp->Fill((bb_ps_e+bb_sh_e)/bb_tr_p);
                  h2_eovp_runnum->Fill(runTrack,(bb_ps_e+bb_sh_e)/bb_tr_p);
                  h2_eovp_runnum_prof->Fill(runTrack,(bb_ps_e+bb_sh_e)/bb_tr_p,1);

                  h2_ps_tot->Fill(bb_ps_e,bb_ps_e+bb_sh_e);

                  h2_pse_trx->Fill(bb_tr_x,bb_ps_e);
                  h2_pse_try->Fill(bb_tr_y,bb_ps_e);
                  h2_pse_pstime->Fill(bb_ps_e,bb_ps_atimeblk);

                  h2_she_trx->Fill(bb_tr_x,bb_sh_e);
                  h2_she_try->Fill(bb_tr_y,bb_sh_e);
                  h2_she_shtime->Fill(bb_sh_e,bb_sh_atimeblk);

                  h2_trp_trx->Fill(bb_tr_x,bb_tr_p);
                  h2_trp_try->Fill(bb_tr_y,bb_tr_p);
                  h2_trp_bbcalE->Fill(bb_tr_p,bb_ps_e+bb_sh_e);

                  h2_hcal_e_x->Fill(sbs_hcal_x,sbs_hcal_e);
                  h2_hcal_e_y->Fill(sbs_hcal_y,sbs_hcal_e);
                  h_pN->Fill(p_N);

                  h2_hcal_time_e->Fill(sbs_hcal_atimeblk,sbs_hcal_e);
                  h2_hodo_time_hcal_e->Fill(bb_hodotdc_clus_tmean,sbs_hcal_e);

                  h2_hcal_atime_x->Fill(sbs_hcal_x,sbs_hcal_atimeblk);
                  prof_hcal_atime_x->Fill(sbs_hcal_x,sbs_hcal_atimeblk,1);
                  h2_hcal_atime_y->Fill(sbs_hcal_y,sbs_hcal_atimeblk);
                  prof_hcal_atime_y->Fill(sbs_hcal_y,sbs_hcal_atimeblk,1);

                  h2_ps_atime_x->Fill(bb_tr_x,bb_ps_atimeblk);
                  prof_ps_atime_x->Fill(bb_tr_x,bb_ps_atimeblk,1);
                  h2_ps_atime_y->Fill(bb_tr_y,bb_ps_atimeblk);
                  prof_ps_atime_y->Fill(bb_tr_y,bb_ps_atimeblk,1);

                  h2_sh_atime_x->Fill(bb_tr_x,bb_sh_atimeblk);
                  prof_sh_atime_x->Fill(bb_tr_x,bb_sh_atimeblk,1);
                  h2_sh_atime_y->Fill(bb_tr_y,bb_sh_atimeblk);
                  prof_sh_atime_y->Fill(bb_tr_y,bb_sh_atimeblk,1);

                  hodo_hcal_coin = bb_hodotdc_clus_tmean-sbs_hcal_atimeblk;
                  hodo_ps_coin = bb_hodotdc_clus_tmean-bb_ps_atimeblk;
                  hodo_sh_coin = bb_hodotdc_clus_tmean-bb_sh_atimeblk;

                  h2_hcalTH_atime_x->Fill(sbs_hcal_x,hodo_hcal_coin);
                  h2_psTH_atime_x->Fill(bb_tr_x,hodo_ps_coin);
                  h2_shTH_atime_x->Fill(bb_tr_x,hodo_sh_coin);

                  h2_hcalTH_atime_y->Fill(sbs_hcal_y,hodo_hcal_coin);
                  h2_psTH_atime_y->Fill(bb_tr_y,hodo_ps_coin);
                  h2_shTH_atime_y->Fill(bb_tr_y,hodo_sh_coin);

                  Sf = 2 * sbs_hcal_e * Mp / e_kine_Q2;
                  h_Sf->Fill(Sf);
                  h2_Sf_x->Fill(sbs_hcal_x,Sf);
                  h2_Sf_y->Fill(sbs_hcal_y,Sf);

                  measE = sbs_hcal_e / (TMath::Sqrt((Mp * Mp) + (p_N * p_N)) - Mp);
                  h_measE->Fill(measE);

                  deltaEfrac = measE - Sf;
                  h2_delta_Sf_measE->Fill(measE,Sf);
                  h_delta_Sf_measE->Fill(deltaEfrac);

                  h2_ps_xy->Fill(bb_ps_y,bb_ps_x);
                  h2_sh_xy->Fill(bb_sh_y,bb_sh_x);
                  h2_hcal_xy->Fill(sbs_hcal_y,sbs_hcal_x);

                  hcalPrimTot = sbs_hcal_prim_e / sbs_hcal_e;
                  hcalSecPrim = sbs_hcal_sec_e / sbs_hcal_prim_e;
                  h_hcal_prim_tot_e->Fill(hcalPrimTot);
                  h_hcal_sec_prim_e->Fill(hcalSecPrim);
                  h2_hcal_primFrac_blkID->Fill(sbs_hcal_prim_id,hcalPrimTot);
                  prof_hcal_primFrac_blkID->Fill(sbs_hcal_prim_id,hcalPrimTot,1);
                  h2_hcal_primFrac_x->Fill(sbs_hcal_x,hcalPrimTot);
                  prof_hcal_primFrac_x->Fill(sbs_hcal_x,hcalPrimTot,1);
                  h2_hcal_primFrac_y->Fill(sbs_hcal_y,hcalPrimTot);
                  prof_hcal_primFrac_y->Fill(sbs_hcal_y,hcalPrimTot,1);

                }// end spot cuts

              }// end W2 cut

            }// end coin cut

          }// end E/p cut

        }// end ps cut

      }// end vertex cut

   }// end global cuts

  }//end event loop

  TCanvas *c1 = new TCanvas("c1","1D dx and dy Plots",100,100,1200,1200);
  c1->Divide(1,2);
  c1->cd(1);
  h_dx->Draw();
  c1->cd(2);
  h_dy->Draw();

  //Save the canvas to a pdf
  c1->Print(outputfile+"(");

  TCanvas *c1_2 = new TCanvas("c1_2", "dxdy Plot", 100,100,1200,1200);
  c1_2->cd();
  h2_dxdy->Draw("colz");

  TBox* Bp = new TBox(dy_mean-dy_sigma, dx_p_mean-dx_p_sigma, dy_mean+dy_sigma, dx_p_mean+dx_p_sigma);
  Bp->SetFillStyle(0);
  Bp->SetLineColor(2);
  Bp->SetLineWidth(2);
  Bp->Draw();

  TBox* Bn = new TBox(dy_mean-dy_sigma, dx_n_mean-dx_n_sigma, dy_mean+dy_sigma, dx_n_mean+dx_n_sigma);
  Bn->SetFillStyle(0);
  Bn->SetLineColor(3);
  Bn->SetLineWidth(2);
  Bn->Draw();

  //Save the canvas to a pdf
  c1_2->Print(outputfile);

  TCanvas *c2 = new TCanvas("c2","QE Cuts", 1200, 1000);
  c2->Divide(2,2);
  c2->cd(1);
  h_tr_vz->Draw();
  c2->cd(2);
  h_ps_e_raw->Draw();
  c2->cd(3);
  h2_pse_grclus->Draw("colz");
  c2->cd(4);
  //h_W2_raw->Draw();
  h_W2->Draw();
  //auto legend_W2 = new TLegend(0.1,0.7,0.4,0.9);
  //legend_W2->AddEntry(h_W2_raw, "W2 with Global Cuts", "l");
  //legend_W2->AddEntry(h_W2, "W2 with QE Cuts", "l");
  //legend_W2->Draw();

  //Save the canvas to a pdf
  c2->Print(outputfile);

  TCanvas *c3 = new TCanvas("c3","W2", 1200, 1000);
  c3->Divide(2,2);
  c3->cd(1);
  h2_coin_W2->Draw("colz");
  c3->cd(2);
  h_eovp->Draw();
  c3->cd(3);
  //h2_eovp_runnum->Draw("colz");
  //h2_eovp_runnum_prof->Draw("SAMES");
  c3->cd(4);
  h2_ps_tot->Draw("colz");

  //Save the canvas to a pdf
  c3->Print(outputfile);

  TCanvas* cPos = new TCanvas("cPos", "Cal Positions", 1200, 1000);
  cPos->Divide(3,1);
  cPos->cd(1);
  h2_ps_xy->Draw("colz");
  cPos->cd(2);
  h2_sh_xy->Draw("colz");
  cPos->cd(3);
  h2_hcal_xy->Draw("colz");

  cPos->Print(outputfile);

  TCanvas *cEp = new TCanvas("cEp", "E/p", 1200, 1000);
  cEp->cd();
  h2_eovp_runnum->Draw("colz");
  h2_eovp_runnum_prof->Draw("SAMES");
  cEp->Update();
  TLine *Ep1 = new TLine(0.0, 1.0, totRun, 1.0);
  //Ep1->SetLineColor(kRed);
  Ep1->SetLineWidth(5);
  Ep1->Draw();

  cEp->Print(outputfile);

  TCanvas *c4 = new TCanvas("c4","bbE", 1200, 1000);
  c4->Divide(2,2);
  c4->cd(1);
  h2_pse_trx->Draw("colz");
  c4->cd(2);
  h2_pse_try->Draw("colz");
  c4->cd(3);
  h2_she_trx->Draw("colz");
  c4->cd(4);
  h2_she_try->Draw("colz");

  c4->Print(outputfile);

  TCanvas *c4_2 = new TCanvas("c4_2","bbPS", 1200, 1000);
  c4_2->Divide(1,2);
  c4_2->cd(1);
  h_ps_e->Draw();
  h_ps_e_anti->Draw("same");
  h_ps_e_qe->Draw("same");
  auto legend_ps = new TLegend(0.55,0.7,0.9,0.9);
  legend_ps->AddEntry(h_ps_e, "PS with Global Cuts", "l");
  legend_ps->AddEntry(h_ps_e_qe, "PS with QE Cuts", "l");
  legend_ps->AddEntry(h_ps_e_anti, "PS with Anti-QE Cuts", "l");
  legend_ps->Draw();
  c4_2->cd(2);
  h2_pse_pstime->Draw("colz");

  c4_2->Print(outputfile);

  TCanvas *c4_3 = new TCanvas("c4_3","bbSH", 1200, 1000);
  c4_3->Divide(2,2);
  c4_3->cd(1);
  h_sh_e->Draw();
  h_sh_e_anti->Draw("same");
  h_sh_e_qe->Draw("same");
  auto legend_sh = new TLegend(0.55,0.7,0.9,0.9);
  legend_sh->AddEntry(h_sh_e, "PS with Global Cuts", "l");
  legend_sh->AddEntry(h_sh_e_qe, "PS with QE Cuts", "l");
  legend_sh->AddEntry(h_sh_e_anti, "PS with Anti-QE Cuts", "l");
  legend_sh->Draw();
  c4_3->cd(2);
  h2_she_shtime->Draw("colz");

  c4_3->Print(outputfile);

  TCanvas *cHCal = new TCanvas("cHCal","sbsHCal", 1200, 1000);
  cHCal->Divide(2,2);
  cHCal->cd(1);
  h2_hcal_e_x->Draw("colz");
  cHCal->cd(2);
  h2_hcal_e_y->Draw("colz");
  cHCal->cd(3);
  h2_hcal_time_e->Draw("colz");
  cHCal->cd(4);
  h2_hodo_time_hcal_e->Draw("colz");

  cHCal->Print(outputfile);

  TCanvas *cHCal_2 = new TCanvas("cHCal_2","sbsHCal_2", 1200, 1000);
  cHCal_2->Divide(2,2);
  cHCal_2->cd(1);
  h_Sf->Draw();
  cHCal_2->cd(2);
  h_measE->Draw();
  cHCal_2->cd(3);
  h_delta_Sf_measE->Draw();
  cHCal_2->cd(4);
  h2_delta_Sf_measE->Draw("colz");
  TF1 *f2 = new TF1("f2", "x",0.0,0.4);
  f2->SetLineColor(kRed);
  f2->Draw("SAMES");

  cHCal_2->Print(outputfile);

  TCanvas *cHCal_3 = new TCanvas("cHCal_3","sbsHCal_3", 1200, 1000);
  cHCal_3->Divide(1,2);
  cHCal_3->cd(1);
  h2_Sf_x->Draw("colz");
  cHCal_3->Update();
  TLine *SfxL = new TLine(-2.6, 0.0795, 1.1, 0.0795);
  SfxL->SetLineColor(kRed);
  SfxL->SetLineWidth(3);
  SfxL->Draw();
  cHCal_3->cd(2);
  h2_Sf_y->Draw("colz");
  cHCal_3->Update();
  TLine *SfyL = new TLine(-0.85, 0.0795, 0.85, 0.0795);
  SfyL->SetLineColor(kRed);
  SfyL->SetLineWidth(3);
  SfyL->Draw();

  cHCal_3->Print(outputfile);

  TCanvas* cHCal_4 = new TCanvas("cHCal_4", "sbsHCal_4", 1200, 1000);
  cHCal_4->Divide(2,1);
  cHCal_4->cd(1);
  h_hcal_prim_tot_e->Draw();
  cHCal_4->cd(2);
  h_hcal_sec_prim_e->Draw();

  cHCal_4->Print(outputfile);

  TCanvas* cHCal_5 = new TCanvas("cHCal_5", "sbsHCal_5", 1200, 1000);
  cHCal_5->Divide(1,3);
  cHCal_5->cd(1);
  h2_hcal_primFrac_blkID->Draw("colz");
  prof_hcal_primFrac_blkID->Draw("SAMES")
  cHCal_5->cd(2);
  h2_hcal_primFrac_x->Draw("colz");
  prof_hcal_primFrac_x->Draw("SAMES")
  cHCal_5->cd(3);
  h2_hcal_primFrac_y->Draw("colz");
  prof_hcal_primFrac_y->Draw("SAMES")

  cHCal_5->Print(outputfile);

  TCanvas *c5 = new TCanvas("c5","bbTr", 1200, 1000);
  c5->Divide(2,2);
  c5->cd(1);
  h2_trp_trx->Draw("colz");
  c5->cd(2);
  h2_trp_try->Draw("colz");
  c5->cd(3);
  h2_trp_bbcalE->Draw("colz");
  TF1 *f1 = new TF1("f1", "x",0.0,10.0);
  f1->SetLineColor(kRed);
  f1->Draw("SAMES");
  c5->cd(4);
  h_pN->Draw();

  c5->Print(outputfile);

  TCanvas *coinTime = new TCanvas("coinTime", "Coincidence Time", 1200, 1000);
  coinTime->cd();
  h_coin->Draw("E");
  coinTime->Update();
  TLine *coinMin = new TLine(coin_mean-coin_sigma, 0, coin_mean-coin_sigma, 10000000);
  coinMin->SetLineColor(kRed);
  coinMin->Draw();
  TLine *coinMax = new TLine(coin_mean+coin_sigma, 0, coin_mean+coin_sigma, 10000000);
  coinMax->SetLineColor(kRed);
  coinMax->Draw();
  gPad->Update();

  coinTime->Print(outputfile);

  TCanvas *c6 = new TCanvas("c6","HCalTime", 1200, 1000);
  c6->Divide(1,2);
  c6->cd(1);
  h2_hcal_atime_x->Draw("colz");
  prof_hcal_atime_x->Draw("SAMES");
  c6->cd(2);
  h2_hcal_atime_y->Draw("colz");
  prof_hcal_atime_y->Draw("SAMES");

  c6->Print(outputfile);

  TCanvas *c7 = new TCanvas("c7","PSTime", 1200, 1000);
  c7->Divide(1,2);
  c7->cd(1);
  h2_ps_atime_x->Draw("colz");
  prof_ps_atime_x->Draw("SAMES");
  c7->cd(2);
  h2_ps_atime_y->Draw("colz");
  prof_ps_atime_y->Draw("SAMES");

  c7->Print(outputfile);

  TCanvas *c8 = new TCanvas("c8","SHTime", 1200, 1000);
  c8->Divide(1,2);
  c8->cd(1);
  h2_sh_atime_x->Draw("colz");
  prof_sh_atime_x->Draw("SAMES");
  c8->cd(2);
  h2_sh_atime_y->Draw("colz");
  prof_sh_atime_y->Draw("SAMES");

  c8->Print(outputfile);

  TCanvas *c9 = new TCanvas("c9","TimeX", 1200, 1000);
  c9->Divide(1,3);
  c9->cd(1);
  h2_shTH_atime_x->Draw("colz");
  c9->cd(2);
  h2_psTH_atime_x->Draw("colz");
  c9->cd(3);
  h2_hcalTH_atime_x->Draw("colz");

  c9->Print(outputfile);

  TCanvas *c10 = new TCanvas("c10","TimeY", 1200, 1000);
  c10->Divide(1,3);
  c10->cd(1);
  h2_shTH_atime_y->Draw("colz");
  c10->cd(2);
  h2_psTH_atime_y->Draw("colz");
  c10->cd(3);
  h2_hcalTH_atime_y->Draw("colz");

  c10->Print(outputfile);

  TCanvas *summary = new TCanvas("summary", "summary", 1200, 1000);
  summary->cd();
  TPaveText *pt = new TPaveText(.05,.1,.95,.8);
  pt->AddText("Global Cuts: ");
  pt->AddText("bb.ps.e>0.0 && bb.gem.track.nhits>=3 && bb.gem.track.chi2ndf[0]<=15 && sbs.hcal.e>0.025");
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

  h_dx->Write();
  h_dy->Write();
  //->Write();
  h_tr_vz->Write();
  h_ps_e->Write();
  h2_pse_grclus->Write();
  h_W2->Write();
  h2_coin_W2->Write();
  h_eovp->Write();
  h2_eovp_runnum->Write();
  h2_pse_trx->Write();
  h2_pse_try->Write();
  h2_she_trx->Write();
  h2_she_try->Write();
  h_ps_e_qe->Write();
  h_ps_e_anti->Write();
  h2_pse_pstime->Write();
  h2_trp_trx->Write();
  h2_trp_try->Write();
  h_coin->Write();
  h2_hcal_atime_x->Write();
  h2_hcal_atime_y->Write();
  h2_ps_atime_x->Write();
  h2_ps_atime_y->Write();
  h2_sh_atime_x->Write();
  h2_sh_atime_y->Write();
  h2_shTH_atime_x->Write();
  h2_psTH_atime_x->Write();
  h2_hcalTH_atime_x->Write();
  h2_shTH_atime_y->Write();
  h2_psTH_atime_y->Write();
  h2_hcalTH_atime_y->Write();

  fout->Write();
  T_data->Delete();
  T->Delete();
}
