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

void QA_QE(const char *kinematic)
{

  int kin = 2;

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  TString inputfile = Form("/volatile/halla/sbs/ktevans/KateJackSBSAnalysis/TEST_pass3/KJ_parsed_GEn_pass2_%s_He3_100.root",kinematic);
  TString outputfile = Form("plots/QA_parsed_GEn_pass2_%s_He3_dxdy.pdf",kinematic);
  TString outfile = Form("outfiles/QA_parsed_GEn_pass2_%s_He3_dxdy.root",kinematic);
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

  double bb_tr_r_x;             T->SetBranchAddress("bb.tr.r_x", &bb_tr_r_x);
  double bb_tr_r_th;            T->SetBranchAddress("bb.tr.r_th", &bb_tr_r_th);
  double e_kine_W2;             T->SetBranchAddress("e.kine.W2", &e_kine_W2);
  double adc_coin;              T->SetBranchAddress("adc.coin", &adc_coin);
  int helicity;                 T->SetBranchAddress("helicity", &helicity);
  double bb_ps_e;               T->SetBranchAddress("bb.ps.e", &bb_ps_e);
  double bb_ps_atimeblk;        T->SetBranchAddress("bb.ps.atimeblk", &bb_ps_atimeblk);
  double bb_sh_e;               T->SetBranchAddress("bb.sh.e", &bb_sh_e);
  double bb_sh_atimeblk;        T->SetBranchAddress("bb.sh.atimeblk", &bb_sh_atimeblk);
  double bb_tr_p;               T->SetBranchAddress("bb.tr.p", &bb_tr_p);
  double bb_tr_x;               T->SetBranchAddress("bb.tr.x", &bb_tr_x);
  double bb_tr_y;               T->SetBranchAddress("bb.tr.y", &bb_tr_y);
  double bb_tr_vz;              T->SetBranchAddress("bb.tr.vz", &bb_tr_vz);
  double bb_gr_clus_size;       T->SetBranchAddress("bb.grinch_tdc.clus.size", &bb_gr_clus_size);
  double bb_gr_clus_track;      T->SetBranchAddress("bb.grinch_tdc.clus.trackindex", &bb_gr_clus_track);
  int pass_global;              T->SetBranchAddress("passGlobal", &pass_global);
  int runnum;                   T->SetBranchAddress("runnum", &runnum);
  double sbs_hcal_e;            T->SetBranchAddress("sbs.hcal.e", &sbs_hcal_e);
  double sbs_hcal_x;            T->SetBranchAddress("sbs.hcal.x", &sbs_hcal_x);
  double sbs_hcal_y;            T->SetBranchAddress("sbs.hcal.y", &sbs_hcal_y);
  double sbs_hcal_atimeblk;     T->SetBranchAddress("sbs.hcal.atimeblk", &sbs_hcal_atimeblk);
  int sbs_hcal_clus_blk_id;     T->SetBranchAddress("sbs.hcal.clus_blk.id", &sbs_hcal_clus_blk_id);
  double dx_hcal;               T->SetBranchAddress("dx", &dx_hcal);
  double dy_hcal;               T->SetBranchAddress("dy", &dy_hcal);
  double IHWP;                  T->SetBranchAddress("IHWP", &IHWP);

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

  int IHWP_flip;

  if(kin==2)
  {
    optics_valid_min = -0.35;
    optics_valid_max = 0.34;
    //coin_mean = 129.1;
    coin_mean = -0.08;
    //coin_sigma = 5.6;
    coin_sigma = 1.883;
    IHWP_flip = -1;
    dx_n_mean = -0.147;
    dx_n_sigma = 0.812;
    dx_p_mean = -2.811;
    dx_p_sigma = 0.600;
    dy_mean = 0.431;
    dy_sigma = 1.265;
    firstRun = 2130;
    lastRun = 2322;
  }
  else if(kin==3)
  {
    optics_valid_min = -0.35;
    optics_valid_max = 0.33;
    //coin_mean = 120.3;
    coin_mean = 0.4239;
    //coin_sigma = 6.0;
    coin_sigma = 2.728;
    IHWP_flip = 1;
    dx_n_mean = 0.0;
    dx_n_sigma = 1.0;
    dx_p_mean = -1.541;
    dx_p_sigma = 0.365;
    dy_mean = 0.336;
    dy_sigma = 0.913;
    firstRun = 2506;
    lastRun = 3250;
  }
  else if(kin==4)
  {
    optics_valid_min = -0.36;
    optics_valid_max = 0.30;
    //coin_mean = 121.4;
    coin_mean = 0.2289;
    //coin_sigma = 5.8;
    coin_sigma = 2.017;
    IHWP_flip = 1;
    dx_n_mean = 0.0;
    dx_n_sigma = 0.5;
    dx_p_mean = -1.124;
    dx_p_sigma = 0.464;
    dy_mean = 0.254;
    dy_sigma = 0.773;
    firstRun = 3510;
    lastRun = 4587;
  }
  else if(kin==5)
  {
    optics_valid_min = -0.37;
    optics_valid_max = 0.32;
    //coin_mean = 185.5;
    coin_mean = 0.2546;
    //coin_sigma = 7.0;
    coin_sigma = 2.695;
    IHWP_flip = 1;
    dx_n_mean = 0.0;
    dx_n_sigma = 0.5;
    dx_p_mean = -1.124;
    dx_p_sigma = 0.361;
    dy_mean = 0.248;
    dy_sigma = 0.769;
    firstRun = 5044;
    lastRun = 6083;
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

  TH1D* h_tr_vz = new TH1D("h_tr_vz", "Track z Vertex", 140.0, -0.8, 0.6);
  h_tr_vz->GetXaxis()->SetTitle("bb.tr.vz [m]");
  h_tr_vz->SetTitle("Track z Vertex with Global Cuts");

  TH1D* h_ps_e = new TH1D("h_ps_e", "PreShower Energy", 100.0, 0.0, 2.0);
  h_ps_e->GetXaxis()->SetTitle("bb.ps.e [GeV]");
  h_ps_e->SetTitle("PreShower Energy with Global and Vertex Cuts");

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

  TH1D* h_eovp = new TH1D("h_eovp", "E/p", 100.0, 0.5, 1.5);
  h_eovp->GetXaxis()->SetTitle("E/p");
  h_eovp->SetTitle("E/p with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_pse_trx = new TH2D("h2_pse_trx", "PSe vs TrX", 100.0, -0.45, 0.6, 100.0, 0.0, 3.0);
  h2_pse_trx->GetYaxis()->SetTitle("bb.ps.e [GeV]");
  h2_pse_trx->GetXaxis()->SetTitle("bb.tr.x [m]");
  h2_pse_trx->SetTitle("PreShower Energy vs Track x with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_pse_try = new TH2D("h2_pse_try", "PSe vs TrY", 100.0, -0.2, 0.15, 100.0, 0.0, 3.0);
  h2_pse_try->GetYaxis()->SetTitle("bb.ps.e [GeV]");
  h2_pse_try->GetXaxis()->SetTitle("bb.tr.y [m]");
  h2_pse_try->SetTitle("PreShower Energy vs Track y with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_she_trx = new TH2D("h2_she_trx", "SHe vs TrX", 100.0, -0.45, 0.6, 100.0, 0.0, 3.0);
  h2_she_trx->GetYaxis()->SetTitle("bb.sh.e [GeV]");
  h2_she_trx->GetXaxis()->SetTitle("bb.tr.x [m]");
  h2_she_trx->SetTitle("Shower Energy vs Track x with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_she_try = new TH2D("h2_she_try", "SHe vs TrY", 100.0, -0.2, 0.15, 100.0, 0.0, 3.0);
  h2_she_try->GetYaxis()->SetTitle("bb.sh.e [GeV]");
  h2_she_try->GetXaxis()->SetTitle("bb.tr.y [m]");
  h2_she_try->SetTitle("Shower Energy vs Track y with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_trp_trx = new TH2D("h2_trp_trx", "TrP vs TrX", 100.0, -0.45, 0.6, 100.0, 2.0, 3.5);
  h2_trp_trx->GetYaxis()->SetTitle("bb.tr.p [GeV]");
  h2_trp_trx->GetXaxis()->SetTitle("bb.tr.x [m]");
  h2_trp_trx->SetTitle("Track p vs Track x with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_trp_try = new TH2D("h2_trp_try", "TrP vs TrY", 100.0, -0.2, 0.15, 100.0, 2.0, 3.5);
  h2_trp_try->GetYaxis()->SetTitle("bb.tr.p [GeV]");
  h2_trp_try->GetXaxis()->SetTitle("bb.tr.y [m]");
  h2_trp_try->SetTitle("Track p vs Track y with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH1D* h_coin = new TH1D("h_coin", "coin", 100.0, -5.0, 5.0);
  h_coin->GetXaxis()->SetTitle("Coincidence Time [ns]");
  h_coin->SetTitle("Coin Time (HCal-BBCal) with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_eovp_runnum = new TH2D("h2_eovp_runnum", "E/p vs runnum", lastRun-firstRun, firstRun, lastRun, 100.0, 0.5, 1.5);
  h2_eovp_runnum->GetYaxis()->SetTitle("E/p");
  h2_eovp_runnum->GetXaxis()->SetTitle("runnum");
  h2_eovp_runnum->SetTitle("E/p vs Run Number with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_hcal_atime_x = new TH2D("h2_hcal_atime_x", "HCal time vs HCal x", 37.0, -2.6, 1.1, 100.0, -10.0, 15.0);
  h2_hcal_atime_x->GetXaxis()->SetTitle("sbs.hcal.x [m]");
  h2_hcal_atime_x->GetYaxis()->SetTitle("sbs.hcal.atimeblk [ns]");
  h2_hcal_atime_x->SetTitle("HCal ADC Time vs HCal x with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_hcal_atime_y = new TH2D("h2_hcal_atime_y", "HCal time vs HCal y", 18.0, -0.9, 0.9, 100.0, -10.0, 15.0);
  h2_hcal_atime_y->GetXaxis()->SetTitle("sbs.hcal.y [m]");
  h2_hcal_atime_y->GetYaxis()->SetTitle("sbs.hcal.atimeblk [ns]");
  h2_hcal_atime_y->SetTitle("HCal ADC Time vs HCal y with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_ps_atime_x = new TH2D("h2_ps_atime_x", "PS Time vs Tr x", 105.0, -0.45, 0.6, 100.0, -10.0, 15.0);
  h2_ps_atime_x->GetXaxis()->SetTitle("bb.tr.x [m]");
  h2_ps_atime_x->GetYaxis()->SetTitle("bb.ps.atimeblk [ns]");
  h2_ps_atime_x->SetTitle("PS ADC Time vs Track x with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_ps_atime_y = new TH2D("h2_ps_atime_y", "PS Time vs Tr y", 35.0, -0.2, 0.15, 100.0, -10.0, 15.0);
  h2_ps_atime_y->GetXaxis()->SetTitle("bb.tr.y [m]");
  h2_ps_atime_y->GetYaxis()->SetTitle("bb.ps.atimeblk [ns]");
  h2_ps_atime_y->SetTitle("PS ADC Time vs Track y with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_sh_atime_x = new TH2D("h2_sh_atime_x", "SH Time vs Tr x", 105.0, -0.45, 0.6, 100.0, -10.0, 15.0);
  h2_sh_atime_x->GetXaxis()->SetTitle("bb.tr.x [m]");
  h2_sh_atime_x->GetYaxis()->SetTitle("bb.sh.atimeblk [ns]");
  h2_sh_atime_x->SetTitle("SH ADC Time vs Track x with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  TH2D* h2_sh_atime_y = new TH2D("h2_sh_atime_y", "SH Time vs Tr y", 35.0, -0.2, 0.15, 100.0, -10.0, 15.0);
  h2_sh_atime_y->GetXaxis()->SetTitle("bb.tr.y [m]");
  h2_sh_atime_y->GetYaxis()->SetTitle("bb.sh.atimeblk [ns]");
  h2_sh_atime_y->SetTitle("SH ADC Time vs Track y with Global, Vertex, E/p, PSe, Coin, GRINCH, W2, and Spot Cuts");

  //Loop over all events to fill the histogram
  for (size_t iev = 0; iev < T->GetEntries(); iev++)
  {
    T->GetEntry(iev);

    // && abs(adc_coin-coin_mean)<coin_sigma && bb_ps_e>0.2 && abs(((bb_ps_e+bb_sh_e)/bb_tr_p)-1)<0.2 && bb_gr_clus_size>2.0 && abs(bb_tr_vz)<0.27

    if(pass_global==1 && (IHWP==-1.0 || IHWP==1.0) && (bb_tr_r_x-0.9*bb_tr_r_th)>optics_valid_min && (bb_tr_r_x-0.9*bb_tr_r_th)<optics_valid_max && bb_gr_clus_track==0)
    {

      h_tr_vz->Fill(bb_tr_vz);

      if (abs(bb_tr_vz)<0.27)
      {
        h_ps_e->Fill(bb_ps_e);

        h2_coin_W2->Fill(e_kine_W2,adc_coin);

        if (bb_ps_e>0.2)
        {

          if (abs(((bb_ps_e+bb_sh_e)/bb_tr_p)-1)<0.2)
          {

            h2_pse_grclus->Fill(bb_ps_e,bb_gr_clus_size);

            if (abs(adc_coin-coin_mean)<(2*coin_sigma) && bb_gr_clus_size>2)
            {
              h_W2->Fill(e_kine_W2);

              if (abs(e_kine_W2-1.0)<0.5)
              {
                h_dx->Fill(dx_hcal);
                h_dy->Fill(dy_hcal);
                h2_dxdy->Fill(dy_hcal,dx_hcal);

                dx_out = dx_hcal;
                dy_out = dy_hcal;

                if (abs(dy_hcal-dy_mean)<dy_sigma && (abs(dx_hcal-dx_n_mean)<dx_n_sigma)||(abs(dx_hcal-dx_p_mean)<dx_p_sigma))
                {

                  h_eovp->Fill((bb_ps_e+bb_sh_e)/bb_tr_p);
                  h2_eovp_runnum->Fill(runnum,(bb_ps_e+bb_sh_e)/bb_tr_p);

                  h2_pse_trx->Fill(bb_tr_x,bb_ps_e);
                  h2_pse_try->Fill(bb_tr_y,bb_ps_e);

                  h2_she_trx->Fill(bb_tr_x,bb_sh_e);
                  h2_she_try->Fill(bb_tr_y,bb_sh_e);

                  h2_trp_trx->Fill(bb_tr_x,bb_tr_p);
                  h2_trp_try->Fill(bb_tr_y,bb_tr_p);

                  h_coin->Fill(adc_coin);

                  h2_hcal_atime_x->Fill(sbs_hcal_x,sbs_hcal_atimeblk);
                  h2_hcal_atime_y->Fill(sbs_hcal_y,sbs_hcal_atimeblk);

                  h2_ps_atime_x->Fill(bb_tr_x,bb_ps_atimeblk);
                  h2_ps_atime_y->Fill(bb_tr_y,bb_ps_atimeblk);

                  h2_sh_atime_x->Fill(bb_tr_x,bb_sh_atimeblk);
                  h2_sh_atime_y->Fill(bb_tr_y,bb_sh_atimeblk);

                }// end spot cuts

              }// end W2 cut

            }// end coin cut

          }// end E/p cut

        }// end ps cut

      }// end vertex cut

      W2_out = e_kine_W2;
      helicity_out = helicity * (IHWP * IHWP_flip); //* IHWP;//-1* IHWP * helicity * IHWP_flip;
   }// end global cuts

    T_data->Fill();

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
  h_ps_e->Draw();
  c2->cd(3);
  h2_pse_grclus->Draw("colz");
  c2->cd(4);
  h_W2->Draw();

  //Save the canvas to a pdf
  c2->Print(outputfile);

  TCanvas *c3 = new TCanvas("c3","W2", 1200, 1000);
  c3->Divide(1,3);
  c3->cd(1);
  h2_coin_W2->Draw("colz");
  c3->cd(2);
  h_eovp->Draw();
  c3->cd(3);
  h2_eovp_runnum->Draw("colz");

  //Save the canvas to a pdf
  c3->Print(outputfile);

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

  TCanvas *c5 = new TCanvas("c5","bbTr", 1200, 1000);
  c5->Divide(1,2);
  c5->cd(1);
  h2_trp_trx->Draw("colz");
  c5->cd(2);
  h2_trp_try->Draw("colz");

  c5->Print(outputfile);

  TCanvas *coinTime = new TCanvas("coinTime", "Coincidence Time", 1200, 1000);
  coinTime->cd();
  h_coin->Draw();
  coinTime->Update();
  TLine *coinMin = new TLine(coin_mean-2*coin_sigma, 0, coin_mean-2*coin_sigma, 10000000);
  coinMin->SetLineColor(kRed);
  coinMin->Draw();
  TLine *coinMax = new TLine(coin_mean+2*coin_sigma, 0, coin_mean+2*coin_sigma, 10000000);
  coinMax->SetLineColor(kRed);
  coinMax->Draw();
  gPad->Update();

  coinTime->Print(outputfile);

  TCanvas *c6 = new TCanvas("c6","HCalTime", 1200, 1000);
  c6->Divide(1,2);
  c6->cd(1);
  h2_hcal_atime_x->Draw("colz");
  c6->cd(2);
  h2_hcal_atime_y->Draw("colz");

  c6->Print(outputfile);

  TCanvas *c7 = new TCanvas("c7","PSTime", 1200, 1000);
  c7->Divide(1,2);
  c7->cd(1);
  h2_ps_atime_x->Draw("colz");
  c7->cd(2);
  h2_ps_atime_y->Draw("colz");

  c7->Print(outputfile);

  TCanvas *c8 = new TCanvas("c8","SHTime", 1200, 1000);
  c8->Divide(1,2);
  c8->cd(1);
  h2_sh_atime_x->Draw("colz");
  c8->cd(2);
  h2_sh_atime_y->Draw("colz");

  c8->Print(outputfile);

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
  pt->AddText(Form("Coincidence Cut: abs(adc.coin - %.3f)<(2 * %.3f)",coin_mean,coin_sigma));
  pt->AddText("GRINCH cut: bb.grinch_tdc.clus.size>2");
  pt->AddText(Form("Proton Spot Cut: %.3f<dx<%.3f && %.3f<dy<%.3f",(dx_p_mean-dx_p_sigma),(dx_p_mean+dx_p_sigma),(dy_mean-dy_sigma),(dy_mean-dy_sigma)));
  pt->AddText(Form("Neutron Spot Cut: %.3f<dx<%.3f && %.3f<dy<%.3f",(dx_n_mean-dx_p_sigma),(dx_n_mean+dx_p_sigma),(dy_mean-dy_sigma),(dy_mean-dy_sigma)));
  pt->Draw();

  summary->Print(outputfile+")");

  fout->Write();
}
