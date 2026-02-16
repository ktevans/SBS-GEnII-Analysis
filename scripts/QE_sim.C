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

void QE_sim(const char *kinematic)
{

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  TString inputfile = Form("/volatile/halla/sbs/ktevans/QE_sim/QE_sim_%s_sbs100p_nucleon_np_model2.root",kinematic);
  TString outputfile = Form("plots/parsed_SIM_GEn_%s_He3_dxdy.pdf",kinematic);
  TString outfile = Form("outfiles/parsed_SIM_GEn_%s_He3_dxdy.root",kinematic);
  TFile *fout = new TFile(outfile,"RECREATE");

  TTree *T_sim = new TTree("T_sim", "Analysis Data Tree");

  double dx_out, dy_out, fnucl_out, weight_out, epx_out, epy_out, epz_out, npx_out, npy_out, npz_out;
  T_sim->Branch("dx", &dx_out, "dx/D");
  T_sim->Branch("dy", &dy_out, "dy/D");
  T_sim->Branch("fnucl", &fnucl_out, "fnucl/D");
  T_sim->Branch("weight", &weight_out, "weight/D");
  T_sim->Branch("epx", &epx_out, "epx/D");
  T_sim->Branch("epy", &epy_out, "epy/D");
  T_sim->Branch("epz", &epz_out, "epz/D");
  T_sim->Branch("npx", &npx_out, "npx/D");
  T_sim->Branch("npy", &npy_out, "npy/D");
  T_sim->Branch("npz", &npz_out, "npz/D");

  TChain* T = new TChain("Tout");
  T->Add(inputfile);

  double e_kine_W2;         T->SetBranchAddress("W2", &e_kine_W2);
  double bb_ps_e;           T->SetBranchAddress("ePS", &bb_ps_e);
  double bb_sh_e;           T->SetBranchAddress("eSH", &bb_sh_e);
  double bb_tr_p;           T->SetBranchAddress("trP", &bb_tr_p);
  double bb_tr_vz;          T->SetBranchAddress("vz", &bb_tr_vz);
  double sbs_hcal_e;        T->SetBranchAddress("eHCAL", &sbs_hcal_e);
  double bb_gr_clus_size;   T->SetBranchAddress("grinch_clus_size", &bb_gr_clus_size);
  double bb_gr_clus_track;  T->SetBranchAddress("grinch_track", &bb_gr_clus_track);
  double dx_hcal;           T->SetBranchAddress("dx", &dx_hcal);
  double dy_hcal;           T->SetBranchAddress("dy", &dy_hcal);
  double fnucl_hcal;        T->SetBranchAddress("fnucl", &fnucl_hcal);
  double weight_hcal;       T->SetBranchAddress("weight", &weight_hcal);
  double epx_in;            T->SetBranchAddress("epx", &epx_in);
  double epy_in;            T->SetBranchAddress("epy", &epy_in);
  double epz_in;            T->SetBranchAddress("epz", &epz_in);
  double npx_in;            T->SetBranchAddress("npx", &npx_in);
  double npy_in;            T->SetBranchAddress("npy", &npy_in);
  double npz_in;            T->SetBranchAddress("npz", &npz_in);

  //Scan through all the entries in the TChain T
  //If the rootfiles are empty or don't exist, there will be 0 entries
  //If there are entries, then print out how many
  if(T->GetEntries()==0)
  {
    std::cerr << "\n --- No ROOT file found!! --- \n\n";
    throw;
  }
  else std::cout << "\nFound " << T->GetEntries() << " events. Starting analysis.. \n";

  TH1D* h_dx = new TH1D("h_dx", ";dx", 70.0, -4.0, 3.0);
  h_dx->GetXaxis()->SetTitle("dx [m]");
  h_dx->SetTitle("Simulated dx with QE Cuts");

  TH1D* h_dy = new TH1D("h_dy", ";dy", 60.0, -3.0, 3.0);
  h_dy->GetXaxis()->SetTitle("dy [m]");
  h_dy->SetTitle("Simulated dy with QE Cuts");

  //Loop over all events to fill the histogram
  for (size_t iev = 0; iev < T->GetEntries(); iev++)
  {
    T->GetEntry(iev);

    if(abs(e_kine_W2-1.0)<0.7 && bb_ps_e>0.2 && abs(((bb_ps_e+bb_sh_e)/bb_tr_p)-1)<0.2 && abs(bb_tr_vz)<0.27 && sbs_hcal_e>0.025)
    {
      h_dx->Fill(dx_hcal,weight_hcal);
      h_dy->Fill(dy_hcal,weight_hcal);

      dx_out = dx_hcal;
      dy_out = dy_hcal;
      fnucl_out = fnucl_hcal;
      weight_out = weight_hcal;
      epx_out = epx_in;
      epy_out = epy_in;
      epz_out = epz_in;
      npx_out = npx_in;
      npy_out = npy_in;
      npz_out = npz_in;
    }

    T_sim->Fill();

  }//end event loop

  TCanvas *c1 = new TCanvas("c1","1D dx and dy Plots",100,100,700,700);
  c1->Divide(1,2);
  c1->cd(1);
  h_dx->Draw();
  c1->cd(2);
  h_dy->Draw();

  printf("You've completed the script!\n");

  //Save the canvas to a pdf
  c1->Print(outputfile);

  fout->Write();
}
