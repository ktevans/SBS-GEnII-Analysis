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

  TString inputfile = Form("/volatile/halla/sbs/ktevans/QE_sim/QE_sim_%s_sbs100p_nucleon_np_model2_elastic.root",kinematic);
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

  TChain* T = new TChain("Tout");
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

  double dx_p_shift = 0.0;
  double dx_n_shift = 0.0;
  double beam_e = 4.291;

  if (kinematic == "GEN2")
  {
    beam_e = 4.291;
  }
  if (kinematic == "GEN3")
  {
    beam_e = 6.373;
  }
  if (kinematic == "GEN4")
  {
    beam_e = 8.448;
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

  TH1D* h_dx_p = new TH1D("h_dx_p", ";dx_p", 70.0, -4.0, 3.0);
  h_dx_p->GetXaxis()->SetTitle("dx [m]");
  h_dx_p->SetTitle("dx for Protons");

  TH1D* h_dx_n = new TH1D("h_dx_n", ";dx_n", 70.0, -4.0, 3.0);
  h_dx_n->GetXaxis()->SetTitle("dx [m]");
  h_dx_n->SetTitle("dx for Neutrons");

  TH2D* h_dx_missing_mom_p = new TH2D("h_dx_missing_mom_p", ";h_dx_missing_mom_p", 70.0, -4.0, 3.0, 40.0, 0.0, 0.4);
  h_dx_missing_mom_p->GetXaxis()->SetTitle("dx [m]");
  h_dx_missing_mom_p->GetYaxis()->SetTitle("Missing Momentum [GeV]");
  h_dx_missing_mom_p->SetTitle("dx for Protons");

  TH2D* h_dx_missing_mom_n = new TH2D("h_dx_missing_mom_n", ";h_dx_missing_mom_n", 70.0, -4.0, 3.0, 40.0, 0.0, 0.4);
  h_dx_missing_mom_n->GetXaxis()->SetTitle("dx [m]");
  h_dx_missing_mom_n->GetYaxis()->SetTitle("Missing Momentum [GeV]");
  h_dx_missing_mom_n->SetTitle("dx for Neutrons");

  TH2D* h_dx_pol_p = new TH2D("h_dx_pol_p", ";h_dx_pol_p", 70.0, -4.0, 3.0, 20.0, -0.5, 1.5);
  h_dx_pol_p->GetXaxis()->SetTitle("dx [m]");
  h_dx_pol_p->GetYaxis()->SetTitle("Nucleon Effective Polarization");
  h_dx_pol_p->SetTitle("dx for Protons");

  TH2D* h_dx_pol_n = new TH2D("h_dx_pol_n", ";h_dx_pol_n", 70.0, -4.0, 3.0, 20.0, -0.5, 1.5);
  h_dx_pol_n->GetXaxis()->SetTitle("dx [m]");
  h_dx_pol_n->GetYaxis()->SetTitle("Nucleon Effective Polarization");
  h_dx_pol_n->SetTitle("dx for Neutrons");

  //Loop over all events to fill the histogram
  for (size_t iev = 0; iev < T->GetEntries(); iev++)
  {
    T->GetEntry(iev);

    missing_mom = TMath::Sqrt(TMath::Power((-beam_e)+(npz)+(epz),2)+TMath::Power((npx)+(epx),2)+TMath::Power((npy)+(epy),2));

    if(ePS>0.2&&W2<1.5)
    {

      if (fnucl==1.0)
      {
        h_dx_p->Fill(dx,weight);
        h_dx_missing_mom_p->Fill(dx,missing_mom,weight);
        if (missing_mom<0.2)
        {
          pol_p = (191.16 * TMath::Power(missing_mom,4)) - (111.53 * TMath::Power(missing_mom,3)) + (18.734 * TMath::Power(missing_mom,2)) - (0.1466 * missing_mom) - 0.0934;
        }
        if (missing_mom>=0.2)
        {
          pol_p = (585.36 * TMath::Power(missing_mom,4)) - (617.45 * TMath::Power(missing_mom,3)) + (228.85 * TMath::Power(missing_mom,2)) - (36.301 * missing_mom) - 2.1489;
        }
        h_dx_pol_p->Fill(dx,pol_p,weight);
      }

      if (fnucl==0.0)
      {
        h_dx_n->Fill(dx,weight);
        h_dx_missing_mom_n->Fill(dx,missing_mom,weight);
        if (missing_mom<0.2)
        {
          pol_n = (-149.54 * TMath::Power(missing_mom,4)) + (26.742 * TMath::Power(missing_mom,3)) - (3.1007 * TMath::Power(missing_mom,2)) + (0.064 * missing_mom) + 0.9999;
        }
        if (missing_mom>=0.2)
        {
          pol_n = (-409.23 * TMath::Power(missing_mom,4)) + (660.3 * TMath::Power(missing_mom,3)) - (364.59 * TMath::Power(missing_mom,2)) + (77.325 * missing_mom) - 4.6319;
        }
        h_dx_pol_n->Fill(dx,pol_n,weight);
      }

      dx_out = dx;
      dy_out = dy;
      W2_out = W2;

    }

    T_out->Fill();

  }//end event loop

  TCanvas *c1 = new TCanvas("c1","1D dx and dy Plots",100,100,700,700);
  c1->Divide(1,2);
  c1->cd(1);
  h_dx_missing_mom_p->Draw();
  h_dx_missing_mom_n->Draw("SAMES");
  c1->cd(2);
  h_dx_pol_p->Draw();
  h_dx_pol_n->Draw("SAMES");

  printf("You've completed the script!\n");

  //Save the canvas to a pdf
  c1->Print(outputfile);

  fout->Write();
}
