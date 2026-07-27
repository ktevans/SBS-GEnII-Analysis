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

void Inel_Corr(const char *kinematic)
{

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  double numberBins = 96;
  double dx_min_d = -2.8;
  double dx_max_d = 2.0;

  TString inputfile_He3 = Form("/volatile/halla/sbs/ktevans/2026SIM/kte_QE_sim_%s_sbs100p_nucleon_np_model2_elastic.root",kinematic);
  TString inputfile_inel = Form("/volatile/halla/sbs/ktevans/2026SIM/QE_sim_%s_sbs100p_nucleon_np_model2_inelastic.root",kinematic);
  TString outputfile = Form("plots/Inel_Corr_SIM_GEn_%s_He3_dxdy.pdf",kinematic);
  TString outfile = Form("outfiles/Inel_Corr_SIM_GEn_%s_He3_dxdy.root",kinematic);
  TFile *fout = new TFile(outfile,"RECREATE");

  TTree *T_sim = new TTree("T_sim", "Analysis Data Tree");

  double dx_out_p_He3, dx_out_n_He3, dy_out_He3, fnucl_out_He3, weight_out_He3, dx_out_p_inel, dx_out_n_inel, dy_out_inel, fnucl_out_inel, weight_out_inel;
  T_sim->Branch("dx.p.He3", &dx_out_p_He3, "dx.p.He3/D");
  T_sim->Branch("dx.n.He3", &dx_out_n_He3, "dx.n.He3/D");
  T_sim->Branch("dy.He3", &dy_out_He3, "dy.He3/D");
  T_sim->Branch("fnucl.He3", &fnucl_out_He3, "fnucl.He3/D");
  T_sim->Branch("weight.He3", &weight_out_He3, "weight.He3/D");
  T_sim->Branch("dx.p.inel", &dx_out_p_inel, "dx.p.inel/D");
  T_sim->Branch("dx.n.inel", &dx_out_n_inel, "dx.n.inel/D");
  T_sim->Branch("dy.inel", &dy_out_inel, "dy.inel/D");
  T_sim->Branch("fnucl.inel", &fnucl_out_inel, "fnucl.inel/D");
  T_sim->Branch("weight.inel", &weight_out_inel, "weight.inel/D");

  TH1D* h_inel_tot_dx_p = new TH1D("h_inel_tot_dx_p","Total inel Proton Events", numberBins, dx_min_d, dx_max_d);
  h_inel_tot_dx_p->GetXaxis()->SetTitle("dx [m]");

  TH1D* h_inel_QE_dx_p = new TH1D("h_inel_QE_dx_p","QE inel Proton Events", numberBins, dx_min_d, dx_max_d);
  h_inel_QE_dx_p->GetXaxis()->SetTitle("dx [m]");

  TH1D* h_He3_tot_dx_p = new TH1D("h_He3_tot_dx_p","Total He3 Proton Events", numberBins, dx_min_d, dx_max_d);
  h_He3_tot_dx_p->GetXaxis()->SetTitle("dx [m]");

  TH1D* h_He3_QE_dx_p = new TH1D("h_He3_QE_dx_p","QE He3 Proton Events", numberBins, dx_min_d, dx_max_d);
  h_He3_QE_dx_p->GetXaxis()->SetTitle("dx [m]");

  TH1D* h_inel_tot_dx_n = new TH1D("h_inel_tot_dx_n","Total inel Neutron Events", numberBins, dx_min_d, dx_max_d);
  h_inel_tot_dx_n->GetXaxis()->SetTitle("dx [m]");

  TH1D* h_inel_QE_dx_n = new TH1D("h_inel_QE_dx_n","QE inel Neutron Events", numberBins, dx_min_d, dx_max_d);
  h_inel_QE_dx_n->GetXaxis()->SetTitle("dx [m]");

  TH1D* h_He3_tot_dx_n = new TH1D("h_He3_tot_dx_n","Total He3 Neutron Events", numberBins, dx_min_d, dx_max_d);
  h_He3_tot_dx_n->GetXaxis()->SetTitle("dx [m]");

  TH1D* h_He3_QE_dx_n = new TH1D("h_He3_QE_dx_n","QE He3 Neutron Events", numberBins, dx_min_d, dx_max_d);
  h_He3_QE_dx_n->GetXaxis()->SetTitle("dx [m]");


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

    if(fnucl_hcal_He3==1.0)
    {
      h_He3_tot_dx_p->Fill(dx_hcal_He3,weight_hcal_He3);
    }

    if(fnucl_hcal_He3==0.0)
    {
      h_He3_tot_dx_n->Fill(dx_hcal_He3,weight_hcal_He3);
    }

    if(abs(e_kine_W2_He3-1.0)<0.5 && bb_ps_e_He3>0.2 && abs(((bb_ps_e_He3+bb_sh_e_He3)/bb_tr_p_He3)-0.97)<0.2 && abs(bb_tr_vz_He3)<0.27 && sbs_hcal_e_He3>0.025 && abs(dy_hcal_He3)<0.5)
    {

      if(fnucl_hcal_He3==1.0)
      {
        h_He3_QE_dx_p->Fill(dx_hcal_He3,weight_hcal_He3);
        dx_out_p_He3 = dx_hcal_He3;
      }

      if(fnucl_hcal_He3==0.0)
      {
        h_He3_QE_dx_n->Fill(dx_hcal_He3,weight_hcal_He3);
        dx_out_n_He3 = dx_hcal_He3;
      }

      dy_out_He3 = dy_hcal_He3;
      fnucl_out_He3 = fnucl_hcal_He3;
      weight_out_He3 = weight_hcal_He3;

      T_sim->Fill();
    }
  }//end event loop

  TChain* T_inel = new TChain("Tout");
  T_inel->Add(inputfile_inel);

  double e_kine_W2_inel;         T_inel->SetBranchAddress("W2", &e_kine_W2_inel);
  double bb_ps_e_inel;           T_inel->SetBranchAddress("ePS", &bb_ps_e_inel);
  double bb_sh_e_inel;           T_inel->SetBranchAddress("eSH", &bb_sh_e_inel);
  double bb_tr_p_inel;           T_inel->SetBranchAddress("trP", &bb_tr_p_inel);
  double bb_tr_vz_inel;          T_inel->SetBranchAddress("vz", &bb_tr_vz_inel);
  double sbs_hcal_e_inel;        T_inel->SetBranchAddress("eHCAL", &sbs_hcal_e_inel);
  double dx_hcal_inel;           T_inel->SetBranchAddress("dx", &dx_hcal_inel);
  double dy_hcal_inel;           T_inel->SetBranchAddress("dy", &dy_hcal_inel);
  double fnucl_hcal_inel;        T_inel->SetBranchAddress("fnucl", &fnucl_hcal_inel);
  double weight_hcal_inel;       T_inel->SetBranchAddress("weight", &weight_hcal_inel);

  //Scan through all the entries in the TChain T
  //If the rootfiles are empty or don't exist, there will be 0 entries
  //If there are entries, then print out how many
  if(T_inel->GetEntries()==0)
  {
    std::cerr << "\n --- No ROOT file found!! --- \n\n";
    throw;
  }
  else std::cout << "\nFound " << T_inel->GetEntries() << " events. Starting analysis.. \n";

  //Loop over all events to fill the histogram
  for (size_t iev = 0; iev < T_inel->GetEntries(); iev++)
  {
    T_inel->GetEntry(iev);

    if(fnucl_hcal_inel==1.0)
    {
      h_inel_tot_dx_p->Fill(dx_hcal_inel,weight_hcal_inel);
    }

    if(fnucl_hcal_inel==0.0)
    {
      h_inel_tot_dx_n->Fill(dx_hcal_inel,weight_hcal_inel);
    }

    if(abs(e_kine_W2_inel-1.0)<0.5 && bb_ps_e_inel>0.2 && abs(((bb_ps_e_inel+bb_sh_e_inel)/bb_tr_p_inel)-0.97)<0.2 && abs(bb_tr_vz_inel)<0.27 && sbs_hcal_e_inel>0.025 && abs(dy_hcal_inel)<0.99)
    {

      if(fnucl_hcal_inel==1.0)
      {
        h_inel_QE_dx_p->Fill(dx_hcal_inel,weight_hcal_inel);
        dx_out_p_inel = dx_hcal_inel;
      }

      if(fnucl_hcal_inel==0.0)
      {
        h_inel_QE_dx_n->Fill(dx_hcal_inel,weight_hcal_inel);
        dx_out_n_inel = dx_hcal_inel;
      }

      dy_out_inel = dy_hcal_inel;
      fnucl_out_inel = fnucl_hcal_inel;
      weight_out_inel = weight_hcal_inel;

      T_sim->Fill();
    }
  }//end event loop

  TH1D* hinelfrac_n = (TH1D*) h_inel_QE_dx_n->Clone("hinelfrac_n");
  hinelfrac_n->Divide(h_inel_tot_dx_n);

  TH1D* hHe3frac_n = (TH1D*) h_He3_QE_dx_n->Clone("hHe3frac_n");
  hHe3frac_n->Divide(h_He3_tot_dx_n);
  //hHe3frac_n->Scale(14.0);

  TH1D* hineldilution_n = (TH1D*) hinelfrac_n->Clone("hineldilution_n");
  hineldilution_n->Divide(hHe3frac_n);
  hineldilution_n->SetTitle("inel Neutron Dilution Simulation Scaling");

  TH1D* hinelfrac_p = (TH1D*) h_inel_QE_dx_p->Clone("hinelfrac_p");
  hinelfrac_p->Divide(h_inel_tot_dx_p);

  TH1D* hHe3frac_p = (TH1D*) h_He3_QE_dx_p->Clone("hHe3frac_p");
  hHe3frac_p->Divide(h_He3_tot_dx_p);
  //hHe3frac_p->Scale(7.0);

  TH1D* hineldilution_p = (TH1D*) hinelfrac_p->Clone("hineldilution_p");
  hineldilution_p->Divide(hHe3frac_p);
  hineldilution_p->SetTitle("inel Proton Dilution Simulation Scaling");

  TCanvas *c1 = new TCanvas("c1","inel Dilution as a Function of dx",100,100,700,700);
  c1->cd();
  //hineldilution_p->Scale(0.0957);
  //hineldilution_n->Scale(0.1747);
  gPad->Update();
  hineldilution_p->Add(hineldilution_n);
  hineldilution_p->Draw();
  //hineldilution_n->Draw("SAMES");
  //gPad->Update();

  printf("You've completed the script!\n");

  //Save the canvas to a pdf
  c1->Print(outputfile);

  hinelfrac_p->Write();
  hHe3frac_p->Write();
  hineldilution_p->Write();
  hinelfrac_n->Write();
  hHe3frac_n->Write();
  hineldilution_n->Write();
  fout->Write();
}
