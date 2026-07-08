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

void N2_Corr(const char *kinematic)
{

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  TString inputfile_He3 = Form("/volatile/halla/sbs/ktevans/QE_sim/QE_sim_%s_sbs100p_nucleon_np_model2_elastic.root",kinematic);
  TString inputfile_N2 = Form("/volatile/halla/sbs/ktevans/2026SIM/QE_sim_%s_sbs100p_nucleon_np_model2_N2_elastic.root",kinematic);
  TString outputfile = Form("plots/N2_Corr_SIM_GEn_%s_He3_dxdy.pdf",kinematic);
  TString outfile = Form("outfiles/N2_Corr_SIM_GEn_%s_He3_dxdy.root",kinematic);
  TFile *fout = new TFile(outfile,"RECREATE");

  TTree *T_sim = new TTree("T_sim", "Analysis Data Tree");

  double dx_out_He3, dy_out_He3, fnucl_out_He3, weight_out_He3, dx_out_N2, dy_out_N2, fnucl_out_N2, weight_out_N2;
  T_sim->Branch("dx.He3", &dx_out_He3, "dx.He3/D");
  T_sim->Branch("dy.He3", &dy_out_He3, "dy.He3/D");
  T_sim->Branch("fnucl.He3", &fnucl_out_He3, "fnucl.He3/D");
  T_sim->Branch("weight.He3", &weight_out_He3, "weight.He3/D");
  T_sim->Branch("dx.N2", &dx_out_N2, "dx.N2/D");
  T_sim->Branch("dy.N2", &dy_out_N2, "dy.N2/D");
  T_sim->Branch("fnucl.N2", &fnucl_out_N2, "fnucl.N2/D");
  T_sim->Branch("weight.N2", &weight_out_N2, "weight.N2/D");

  TH1D* h_N2_tot_dx = new TH1D("h_N2_tot_dx","Total N2 Events", 130, -6.0, 7.0);
  h_N2_tot_dx->GetXaxis()->SetTitle("dx [m]");

  TH1D* h_N2_QE_dx = new TH1D("h_N2_QE_dx","QE N2 Events", 130, -6.0, 7.0);
  h_N2_QE_dx->GetXaxis()->SetTitle("dx [m]");

  TH1D* h_He3_tot_dx = new TH1D("h_He3_tot_dx","Total He3 Events", 130, -6.0, 7.0);
  h_He3_tot_dx->GetXaxis()->SetTitle("dx [m]");

  TH1D* h_He3_QE_dx = new TH1D("h_He3_QE_dx","QE He3 Events", 130, -6.0, 7.0);
  h_He3_QE_dx->GetXaxis()->SetTitle("dx [m]");

  TChain* T_He3 = new TChain("Tout");
  T_He3->Add(inputfile_He3);

  double e_kine_W2_He3;         T_He3->SetBranchAddress("W2", &e_kine_W2_He3);
  double bb_ps_e_He3;           T_He3->SetBranchAddress("ePS", &bb_ps_e_He3);
  double bb_sh_e_He3;           T_He3->SetBranchAddress("eSH", &bb_sh_e_He3);
  double bb_tr_p_He3;           T_He3->SetBranchAddress("trP", &bb_tr_p_He3);
  double bb_tr_vz_He3;          T_He3->SetBranchAddress("vz", &bb_tr_vz_He3);
  double sbs_hcal_e_He3;        T_He3->SetBranchAddress("eHCAL", &sbs_hcal_e_He3);
  double dx_hcal_He3;           T_He3->SetBranchAddress("dx", &dx_hcal_He3);
  double dy_hcal_He3;           T_He3->SetBranchAddress("dy", &dy_hcal_He3);
  double fnucl_hcal_He3;        T_He3->SetBranchAddress("fnucl", &fnucl_hcal_He3);
  double weight_hcal_He3;       T_He3->SetBranchAddress("weight", &weight_hcal_He3);

  //Scan through all the entries in the TChain T
  //If the rootfiles are empty or don't exist, there will be 0 entries
  //If there are entries, then print out how many
  if(T_He3->GetEntries()==0)
  {
    std::cerr << "\n --- No ROOT file found!! --- \n\n";
    throw;
  }
  else std::cout << "\nFound " << T_He3->GetEntries() << " events. Starting analysis.. \n";

  //Loop over all events to fill the histogram
  for (size_t iev = 0; iev < T_He3->GetEntries(); iev++)
  {
    T_He3->GetEntry(iev);

    h_He3_tot_dx->Fill(dx_hcal_He3,weight_hcal_He3);

    if(abs(e_kine_W2_He3-1.0)<0.5 && bb_ps_e_He3>0.2 && abs(((bb_ps_e_He3+bb_sh_e_He3)/bb_tr_p_He3)-0.97)<0.2 && abs(bb_tr_vz_He3)<0.27 && sbs_hcal_e_He3>0.025)
    {
      h_He3_QE_dx->Fill(dx_hcal_He3,weight_hcal_He3);

      dx_out_He3 = dx_hcal_He3;
      dy_out_He3 = dy_hcal_He3;
      fnucl_out_He3 = fnucl_hcal_He3;
      weight_out_He3 = weight_hcal_He3;

      T_sim->Fill();
    }
  }//end event loop

  TChain* T_N2 = new TChain("Tout");
  T_N2->Add(inputfile_N2);

  double e_kine_W2_N2;         T_N2->SetBranchAddress("W2", &e_kine_W2_N2);
  double bb_ps_e_N2;           T_N2->SetBranchAddress("ePS", &bb_ps_e_N2);
  double bb_sh_e_N2;           T_N2->SetBranchAddress("eSH", &bb_sh_e_N2);
  double bb_tr_p_N2;           T_N2->SetBranchAddress("trP", &bb_tr_p_N2);
  double bb_tr_vz_N2;          T_N2->SetBranchAddress("vz", &bb_tr_vz_N2);
  double sbs_hcal_e_N2;        T_N2->SetBranchAddress("eHCAL", &sbs_hcal_e_N2);
  double dx_hcal_N2;           T_N2->SetBranchAddress("dx", &dx_hcal_N2);
  double dy_hcal_N2;           T_N2->SetBranchAddress("dy", &dy_hcal_N2);
  double fnucl_hcal_N2;        T_N2->SetBranchAddress("fnucl", &fnucl_hcal_N2);
  double weight_hcal_N2;       T_N2->SetBranchAddress("weight", &weight_hcal_N2);

  //Scan through all the entries in the TChain T
  //If the rootfiles are empty or don't exist, there will be 0 entries
  //If there are entries, then print out how many
  if(T_N2->GetEntries()==0)
  {
    std::cerr << "\n --- No ROOT file found!! --- \n\n";
    throw;
  }
  else std::cout << "\nFound " << T_N2->GetEntries() << " events. Starting analysis.. \n";

  //Loop over all events to fill the histogram
  for (size_t iev = 0; iev < T_N2->GetEntries(); iev++)
  {
    T_N2->GetEntry(iev);

    h_N2_tot_dx->Fill(dx_hcal_N2,weight_hcal_N2);

    if(abs(e_kine_W2_N2-1.0)<0.5 && bb_ps_e_N2>0.2 && abs(((bb_ps_e_N2+bb_sh_e_N2)/bb_tr_p_N2)-0.97)<0.2 && abs(bb_tr_vz_N2)<0.27 && sbs_hcal_e_N2>0.025)
    {
      h_N2_QE_dx->Fill(dx_hcal_N2,weight_hcal_N2);

      dx_out_N2 = dx_hcal_N2;
      dy_out_N2 = dy_hcal_N2;
      fnucl_out_N2 = fnucl_hcal_N2;
      weight_out_N2 = weight_hcal_N2;

      T_sim->Fill();
    }
  }//end event loop

  TH1D* hN2frac = (TH1D*) h_N2_QE_dx->Clone("hN2frac");
  hN2frac->Divide(h_N2_tot_dx);

  TH1D* hHe3frac = (TH1D*) h_He3_QE_dx->Clone("hHe3frac");
  hHe3frac->Divide(h_He3_tot_dx);
  hHe3frac->Scale(14.0);

  TH1D* hN2dilution = (TH1D*) hN2frac->Clone("hN2dilution");
  hN2dilution->Divide(hHe3frac);
  hN2dilution->SetTitle("N2 Dilution Simulation Scaling");

  TCanvas *c1 = new TCanvas("c1","N2 Dilution as a Function of dx",100,100,700,700);
  c1->cd();
  hN2dilution->Draw();

  printf("You've completed the script!\n");

  //Save the canvas to a pdf
  c1->Print(outputfile);

  hN2frac->Write();
  hHe3frac->Write();
  hN2dilution->Write();
  fout->Write();
}
