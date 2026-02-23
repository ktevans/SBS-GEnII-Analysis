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

  int IHWP_flip;

  if(kinematic=="GEN2")
  {
    optics_valid_min = -0.35;
    optics_valid_max = 0.34;
    //coin_mean = 129.1;
    coin_mean = -0.08;
    //coin_sigma = 5.6;
    coin_sigma = 1.883;
    IHWP_flip = -1;
  }
  else if(kinematic=="GEN3")
  {
    optics_valid_min = -0.35;
    optics_valid_max = 0.33;
    //coin_mean = 120.3;
    coin_mean = 0.4239;
    //coin_sigma = 6.0;
    coin_sigma = 2.728;
    IHWP_flip = 1;
  }
  else if(kinematic=="GEN4a")
  {
    optics_valid_min = -0.36;
    optics_valid_max = 0.30;
    //coin_mean = 121.4;
    coin_mean = 0.2289;
    //coin_sigma = 5.8;
    coin_sigma = 2.017;
    IHWP_flip = 1;
  }
  else if(kinematic=="GEN4b")
  {
    optics_valid_min = -0.37;
    optics_valid_max = 0.32;
    //coin_mean = 185.5;
    coin_mean = 0.2546;
    //coin_sigma = 7.0;
    coin_sigma = 2.695;
    IHWP_flip = 1;
  }
  else
  {
    optics_valid_min = -2.0;
    optics_valid_max = 2.0;
    coin_mean = 0.0;
    coin_sigma = 400.0;
    IHWP_flip = 1;
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
  h2_coin_W2->SetTitle("Coincidence Time (HCal-BBCal) vs W2 with Global, Vertex, and PSe Cuts");

  TH2D* h2_pse_grclus = new TH2D("h2_pse_grclus", "PSe vs GRclus", 100.0, 0.0, 2.0, 100.0, 0.0, 20.0);
  h2_pse_grclus->GetXaxis()->SetTitle("bb.ps.e [GeV]");
  h2_pse_grclus->GetYaxis()->SetTitle("bb.gr.clus.size");
  h2_pse_grclus->SetTitle("PreShower Energy vs GRINCH Cluster Size with Global, Vertex, and PSe Cuts");

  TH1D* h_W2 = new TH1D("h_W2", "W2", 100.0, 0.0, 2.0);
  h_W2->GetXaxis()->SetTitle("e.kine.W2 [GeV]");
  h_W2->SetTitle("W2 with Global, Vertex, and PSe Cuts");

  TH1D* h_dx = new TH1D("h_dx", ";dx", 70.0, -4.0, 3.0);
  h_dx->GetXaxis()->SetTitle("dx [m]");
  h_dx->SetTitle("dx with Global, Vertex, PSe, and W2 Cuts");

  TH1D* h_dy = new TH1D("h_dy", ";dy", 60.0, -3.0, 3.0);
  h_dy->GetXaxis()->SetTitle("dy [m]");
  h_dy->SetTitle("dy with Global, Vertex, PSe, and W2 Cuts");

  //Loop over all events to fill the histogram
  for (size_t iev = 0; iev < T->GetEntries(); iev++)
  {
    T->GetEntry(iev);

    // && abs(adc_coin-coin_mean)<coin_sigma && bb_ps_e>0.2 && abs(((bb_ps_e+bb_sh_e)/bb_tr_p)-1)<0.2 && bb_gr_clus_size>2.0 && abs(bb_tr_vz)<0.27

    if(pass_global==1 && (IHWP==-1.0 || IHWP==1.0) && (bb_tr_r_x-0.9*bb_tr_r_th)>optics_valid_min && (bb_tr_r_x-0.9*bb_tr_r_th)<optics_valid_max&&bb_gr_clus_track==0)
    {

      h_tr_vz->Fill(bb_tr_vz);

      if (abs(bb_tr_vz)<0.27)
      {
        h_ps_e->Fill(bb_ps_e);

        if (bb_ps_e>0.2)
        {

          if (abs(((bb_ps_e+bb_sh_e)/bb_tr_p)-1)<0.2)
          {
            h2_coin_W2->Fill(adc_coin,e_kine_W2);

            h2_pse_grclus->Fill(bb_ps_e,bb_gr_clus_size);

            if (abs(adc_coin-coin_mean)<(2*coin_sigma))
            {
              h_W2->Fill(e_kine_W2);

              if(abs(e_kine_W2-1.0)<0.5)
              {
                h_dx->Fill(dx_hcal);
                h_dy->Fill(dy_hcal);

                dx_out = dx_hcal;
                dy_out = dy_hcal;
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

  TCanvas *c1 = new TCanvas("c1","1D dx and dy Plots",100,100,700,700);
  c1->Divide(1,2);
  c1->cd(1);
  h_dx->Draw();
  c1->cd(2);
  h_dy->Draw();

  //Save the canvas to a pdf
  c1->Print(outputfile+"(");

  TCanvas *c2 = new TCanvas("c2","QE Cuts",100,100,700,700);
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

  TCanvas *c3 = new TCanvas("c3","W2",100,100,700,700);
  c3->cd();
  h2_coin_W2->Draw("colz");

  //Save the canvas to a pdf
  c3->Print(outputfile+")");

  fout->Write();
}
