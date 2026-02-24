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

void MissingMom(const char *kinematic)
{

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  TString inputfile = Form("outfiles/parsed_SIM_GEn_%s_He3_dxdy.root",kinematic);
  TString outputfile = Form("plots/parsed_GEn_pass2_%s_simulation.pdf",kinematic);
  TString outfile = Form("outfiles/parsed_GEn_pass2_%s_simulation.root",kinematic);
  TFile *fout = new TFile(outfile,"RECREATE");

  TTree *T_out = new TTree("Tout", "Simulation Data Tree");

  double dx_out, dy_out, W2_out, ePS_out, epx_out, epy_out, epz_out, fnucl_out, npx_out, npy_out, npz_out, weight_out;
  T_out->Branch("dx",     &dx_out,     "dx/D");
  T_out->Branch("dy",     &dy_out,     "dy/D");
  T_out->Branch("W2",     &W2_out,     "W2/D");

  double missing_mom, pol_p, pol_n;
  T_out->Branch("missing_mom",   &missing_mom,   "missing_mom/D");
  T_out->Branch("pol_p",         &pol_p,         "pol_p/D");
  T_out->Branch("pol_n",         &pol_n,         "pol_n/D");

  TChain* T = new TChain("T_sim");
  T->Add(inputfile);

  double dx;         T->SetBranchAddress("dx", &dx);
  double dy;         T->SetBranchAddress("dy", &dy);
  double W2;         T->SetBranchAddress("W2", &W2);
  double ePS;        T->SetBranchAddress("ePS", &ePS);
  double epx;        T->SetBranchAddress("epx", &epx);
  double epy;        T->SetBranchAddress("epy", &epy);
  double epz;        T->SetBranchAddress("epz", &epz);
  double fnucl;      T->SetBranchAddress("fnucl", &fnucl);
  double npx;        T->SetBranchAddress("npx", &npx);
  double npy;        T->SetBranchAddress("npy", &npy);
  double npz;        T->SetBranchAddress("npz", &npz);
  double weight;     T->SetBranchAddress("weight", &weight);

  //double dx_p_shift = 0.375;
  //double dx_p_shift = 0.5;
  double dx_p_shift = 0.0;

  //double dx_n_shift = 0.7;
  //double dx_n_shift = 0.7;
  double dx_n_shift = 0.7;

  //double beam_e = 4.291;
  //double beam_e = 6.373;
  double beam_e = 8.448;

  //double n_min = -0.95;
  //double n_min = -0.5;
  double n_min = -0.3;

  //double n_max = 0.95;
  //double n_max = 0.5;
  double n_max = 0.3;

  //Scan through all the entries in the TChain T
  //If the rootfiles are empty or don't exist, there will be 0 entries
  //If there are entries, then print out how many
  if(T->GetEntries()==0)
  {
    std::cerr << "\n --- No ROOT file found!! --- \n\n";
    throw;
  }
  else std::cout << "\nFound " << T->GetEntries() << " events. Starting analysis.. \n";

  TH1D* h_dx_p = new TH1D("h_dx_p", ";dx_p", 140.0, -4.0, 3.0);
  h_dx_p->GetXaxis()->SetTitle("dx [m]");
  h_dx_p->SetTitle("dx for Protons");

  TH1D* h_dx_n = new TH1D("h_dx_n", ";dx_n", 140.0, -4.0, 3.0);
  h_dx_n->GetXaxis()->SetTitle("dx [m]");
  h_dx_n->SetTitle("dx for Neutrons");

  TH2D* h_dx_missing_mom_p = new TH2D("h_dx_missing_mom_p", ";h_dx_missing_mom_p", 140.0, -4.0, 3.0, 80.0, 0.0, 0.4);
  h_dx_missing_mom_p->GetXaxis()->SetTitle("dx [m]");
  h_dx_missing_mom_p->GetYaxis()->SetTitle("Missing Momentum [GeV]");
  h_dx_missing_mom_p->SetTitle("dx for Protons");

  TH2D* h_dx_missing_mom_n = new TH2D("h_dx_missing_mom_n", ";h_dx_missing_mom_n", 140.0, -4.0, 3.0, 80.0, 0.0, 0.4);
  h_dx_missing_mom_n->GetXaxis()->SetTitle("dx [m]");
  h_dx_missing_mom_n->GetYaxis()->SetTitle("Missing Momentum [GeV]");
  h_dx_missing_mom_n->SetTitle("dx for Neutrons");

  TH2D* h_dx_pol_p = new TH2D("h_dx_pol_p", ";h_dx_pol_p", 140.0, -4.0, 3.0, 80.0, -0.5, 1.5);
  h_dx_pol_p->GetXaxis()->SetTitle("dx [m]");
  h_dx_pol_p->GetYaxis()->SetTitle("Nucleon Effective Polarization");
  h_dx_pol_p->SetTitle("dx for Protons");

  TProfile* h_prof_pol_p = new TProfile("h_prof_pol_p", "pol_prof_p", 70.0, -4.0, 3.0, -0.5, 1.5);
  TProfile* h_prof_pol_n = new TProfile("h_prof_pol_n", "pol_prof_n", 70.0, -4.0, 3.0, -0.5, 1.5);

  TH2D* h_dx_pol_n = new TH2D("h_dx_pol_n", ";h_dx_pol_n", 140.0, -4.0, 3.0, 80.0, -0.5, 1.5);
  h_dx_pol_n->GetXaxis()->SetTitle("dx [m]");
  h_dx_pol_n->GetYaxis()->SetTitle("Nucleon Effective Polarization");
  h_dx_pol_n->SetTitle("dx for Neutrons");

  TH2D* h_dx_pol_p_inWindow = new TH2D("h_dx_pol_p_inWindow", ";h_dx_pol_p_inWindow", 200.0, n_min-0.05, n_max+0.05, 200.0, -0.5, 1.5);
  h_dx_pol_p_inWindow->GetXaxis()->SetTitle("dx [m]");
  h_dx_pol_p_inWindow->GetYaxis()->SetTitle("Nucleon Effective Polarization");
  h_dx_pol_p_inWindow->SetTitle("Proton Effective Polarization in the Neutron Window");

  TH2D* h_dx_pol_n_inWindow = new TH2D("h_dx_pol_n_inWindow", ";h_dx_pol_n_inWindow", 200.0, n_min-0.05, n_max+0.05, 200.0, -0.5, 1.5);
  h_dx_pol_n_inWindow->GetXaxis()->SetTitle("dx [m]");
  h_dx_pol_n_inWindow->GetYaxis()->SetTitle("Nucleon Effective Polarization");
  h_dx_pol_n_inWindow->SetTitle("Neutron Effective Polarization in the Neutron Window");

  //Loop over all events to fill the histogram
  for (size_t iev = 0; iev < T->GetEntries(); iev++)
  {
    T->GetEntry(iev);

    missing_mom = TMath::Sqrt(TMath::Power((-beam_e)+(npz)+(epz),2)+TMath::Power((npx)+(epx),2)+TMath::Power((npy)+(epy),2));

      if (fnucl==1.0)
      {
        h_dx_p->Fill(dx-dx_p_shift,weight);
        h_dx_missing_mom_p->Fill(dx-dx_p_shift,missing_mom,weight);
        if (missing_mom<0.2)
        {
          pol_p = (191.16 * TMath::Power(missing_mom,4)) - (111.53 * TMath::Power(missing_mom,3)) + (18.734 * TMath::Power(missing_mom,2)) - (0.1466 * missing_mom) - 0.0934;
        }
        if (missing_mom>=0.2)
        {
          pol_p = (585.36 * TMath::Power(missing_mom,4)) - (617.45 * TMath::Power(missing_mom,3)) + (228.85 * TMath::Power(missing_mom,2)) - (36.301 * missing_mom) - 2.1489;
        }
        h_dx_pol_p->Fill(dx-dx_p_shift,pol_p,weight);
        h_prof_pol_p->Fill(dx-dx_p_shift,pol_p,weight);
        if (dx-dx_p_shift<n_max&&dx-dx_p_shift>n_min)
        {
          h_dx_pol_p_inWindow->Fill(dx-dx_p_shift,pol_p,weight);
        }
      }

      if (fnucl==0.0)
      {
        h_dx_n->Fill(dx-dx_n_shift,weight);
        h_dx_missing_mom_n->Fill(dx-dx_n_shift,missing_mom,weight);
        if (missing_mom<0.2)
        {
          pol_n = (-149.54 * TMath::Power(missing_mom,4)) + (26.742 * TMath::Power(missing_mom,3)) - (3.1007 * TMath::Power(missing_mom,2)) + (0.064 * missing_mom) + 0.9999;
        }
        if (missing_mom>=0.2)
        {
          pol_n = (-409.23 * TMath::Power(missing_mom,4)) + (660.3 * TMath::Power(missing_mom,3)) - (364.59 * TMath::Power(missing_mom,2)) + (77.325 * missing_mom) - 4.6319;
        }
        h_dx_pol_n->Fill(dx-dx_n_shift,pol_n,weight);
        h_prof_pol_n->Fill(dx-dx_n_shift,pol_n,weight);
        if (dx-dx_n_shift<n_max&&dx-dx_n_shift>n_min)
        {
          h_dx_pol_n_inWindow->Fill(dx-dx_n_shift,pol_n,weight);
        }
      }

      dx_out = dx-dx_n_shift;
      dy_out = dy;
      W2_out = W2;

    T_out->Fill();

  }//end event loop

  double mean_p = 0.0;
  double mean_n = 0.0;

  mean_p = h_dx_pol_p_inWindow->GetMean(2);
  mean_n = h_dx_pol_n_inWindow->GetMean(2);

  TCanvas *c1 = new TCanvas("c1","Missing Momentum and Polarization 2D Plots",100,100,700,700);
  c1->Divide(1,2);
  c1->cd(1);
  h_dx_missing_mom_p->Draw();
  h_dx_missing_mom_n->Draw("SAMES");
  c1->cd(2);
  h_dx_pol_p->Draw();
  h_dx_pol_n->Draw("SAMES");

  TCanvas *c2 = new TCanvas("c2","Effective Polarization in Neutron Window",100,100,700,700);
  c2->Divide(1,2);
  c2->cd(1);
  h_dx_pol_p_inWindow->Draw();
  c2->cd(2);
  h_dx_pol_n_inWindow->Draw();

  printf("Mean Proton Effective Polarization in Neutron Window = %f \n", mean_p);
  printf("Mean Neutron Effective Polarization in Neutron Window = %f \n", mean_n);
  printf("You've completed the script!\n");

  //Save the canvas to a pdf
  //c1->Print(outputfile);
  c2->Print(outputfile);

  fout->Write();
}
