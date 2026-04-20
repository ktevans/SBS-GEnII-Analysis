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

void MissingMom(const char *kinematic, int kin)
{

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  TString inputfile = Form("outfiles/parsed_SIM_GEn_%s_He3_dxdy.root",kinematic);
  TString outputfile = Form("plots/parsed_GEn_pass2_%s_simulation.pdf",kinematic);
  TString outfile = Form("outfiles/parsed_GEn_pass2_%s_simulation.root",kinematic);
  TFile *fout = new TFile(outfile,"RECREATE");

  TTree *T_out = new TTree("Tout", "Simulation Data Tree");

  double dx_out, dy_out, epx_out, epy_out, epz_out, fnucl_out, npx_out, npy_out, npz_out, weight_out;
  T_out->Branch("dx",     &dx_out,     "dx/D");
  T_out->Branch("dy",     &dy_out,     "dy/D");
  //T_out->Branch("W2",     &W2_out,     "W2/D");

  double missing_mom, pol_p, pol_n, pol_p_w, pol_n_w;
  T_out->Branch("missing_mom",   &missing_mom,   "missing_mom/D");
  T_out->Branch("pol_p",         &pol_p,         "pol_p/D");
  T_out->Branch("pol_n",         &pol_n,         "pol_n/D");
  T_out->Branch("pol_p_w",       &pol_p_w,       "pol_p_w/D");
  T_out->Branch("pol_n_w",       &pol_n_w,       "pol_n_w/D");

  TChain* T = new TChain("T_sim");
  T->Add(inputfile);

  double dx;         T->SetBranchAddress("dx", &dx);
  double dy;         T->SetBranchAddress("dy", &dy);
  //double W2;         T->SetBranchAddress("W2", &W2);
  //double ePS;        T->SetBranchAddress("ePS", &ePS);
  double epx;        T->SetBranchAddress("epx", &epx);
  double epy;        T->SetBranchAddress("epy", &epy);
  double epz;        T->SetBranchAddress("epz", &epz);
  double fnucl;      T->SetBranchAddress("fnucl", &fnucl);
  double npx;        T->SetBranchAddress("npx", &npx);
  double npy;        T->SetBranchAddress("npy", &npy);
  double npz;        T->SetBranchAddress("npz", &npz);
  double weight;     T->SetBranchAddress("weight", &weight);

  double dx_p_shift;
  double dx_n_shift;
  double beam_e;
  double n_min;
  double n_max;
  double dx_p_cut;
  double dx_n_cut;

  if (kin == 2)
  {
    dx_p_shift = 0.438;
    dx_n_shift = 0.118;
    beam_e = 4.291;
    n_min = -0.95;
    n_max = 0.95;
    dx_p_cut = 0.0;
    dx_n_cut = -1.0;
  }

  if (kin == 3)
  {
    dx_p_shift = 0.5;
    dx_n_shift = 0.7;
    beam_e = 6.373;
    n_min = -0.5;
    n_max = 0.5;
    dx_p_cut = 0.0;
    dx_n_cut = -1.0;
  }

  if (kin == 4)
  {
    dx_p_shift = 0.0;
    dx_n_shift = 0.7;
    beam_e = 8.448;
    n_min = -0.3;
    n_max = 0.3;
    dx_p_cut = 0.0;
    dx_n_cut = -1.0;
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

  TH2D* h_dx_pol_p_low = new TH2D("h_dx_pol_p_low", ";h_dx_pol_p_low", 140.0, -4.0, 3.0, 80.0, -0.5, 1.5);
  TH2D* h_dx_pol_p_high = new TH2D("h_dx_pol_p_high", ";h_dx_pol_p_high", 140.0, -4.0, 3.0, 80.0, -0.5, 1.5);
  TProfile* h_prof_pol_p_low = new TProfile("h_prof_pol_p_low", "pol_prof_p_low", 70.0, -3.5, 3.0, -0.5, 1.5);
  TProfile* h_prof_pol_p_high = new TProfile("h_prof_pol_p_high", "pol_prof_p_high", 70.0, -3.5, 3.0, -0.5, 1.5);

  TH2D* h_dx_pol_n = new TH2D("h_dx_pol_n", ";h_dx_pol_n", 140.0, -4.0, 3.0, 80.0, -0.5, 1.5);
  h_dx_pol_n->GetXaxis()->SetTitle("dx [m]");
  h_dx_pol_n->GetYaxis()->SetTitle("Nucleon Effective Polarization");
  h_dx_pol_n->SetTitle("dx for Neutrons");

  TH2D* h_dx_pol_n_low = new TH2D("h_dx_pol_n_low", ";h_dx_pol_n_low", 140.0, -4.0, 3.0, 80.0, -0.5, 1.5);
  TH2D* h_dx_pol_n_high = new TH2D("h_dx_pol_n_high", ";h_dx_pol_n_high", 140.0, -4.0, 3.0, 80.0, -0.5, 1.5);
  TProfile* h_prof_pol_n_low = new TProfile("h_prof_pol_n_low", "pol_prof_n_low", 70.0, -3.5, 3.0, -0.5, 1.5);
  TProfile* h_prof_pol_n_high = new TProfile("h_prof_pol_n_high", "pol_prof_n_high", 70.0, -3.5, 3.0, -0.5, 1.5);

  TProfile* h_prof_pol_p = new TProfile("h_prof_pol_p", "pol_prof_p", 70.0, -3.5, 3.0, -0.5, 1.5);
  TProfile* h_prof_pol_n = new TProfile("h_prof_pol_n", "pol_prof_n", 70.0, -3.5, 3.0, -0.5, 1.5);

  TH2D* h_dx_pol_p_inWindow = new TH2D("h_dx_pol_p_inWindow", ";h_dx_pol_p_inWindow", 200.0, n_min-0.05, n_max+0.05, 200.0, -0.5, 1.5);
  h_dx_pol_p_inWindow->GetXaxis()->SetTitle("dx [m]");
  h_dx_pol_p_inWindow->GetYaxis()->SetTitle("Nucleon Effective Polarization");
  h_dx_pol_p_inWindow->SetTitle("Proton Effective Polarization in the Neutron Window");

  TH2D* h_dx_pol_n_inWindow = new TH2D("h_dx_pol_n_inWindow", ";h_dx_pol_n_inWindow", 200.0, n_min-0.05, n_max+0.05, 200.0, -0.5, 1.5);
  h_dx_pol_n_inWindow->GetXaxis()->SetTitle("dx [m]");
  h_dx_pol_n_inWindow->GetYaxis()->SetTitle("Nucleon Effective Polarization");
  h_dx_pol_n_inWindow->SetTitle("Neutron Effective Polarization in the Neutron Window");

  double a, a_err, b, b_err, c, c_err, d, d_err, e, e_err, sigma;

  //Loop over all events to fill the histogram
  for (size_t iev = 0; iev < T->GetEntries(); iev++)
  {
    T->GetEntry(iev);

    missing_mom = TMath::Sqrt(TMath::Power((-beam_e)+(npz)+(epz),2)+TMath::Power((npx)+(epx),2)+TMath::Power((npy)+(epy),2));
    sigma = 1.0/weight;

      if (fnucl==1.0)
      {
        h_dx_p->Fill(dx-dx_p_shift,weight);
        h_dx_missing_mom_p->Fill(dx-dx_p_shift,missing_mom,weight);
        if (missing_mom<0.2)
        {
          a =  198.53995712204684 ;
          a_err =  TMath::Power(24.594485711656095,2) ;
          b =  -112.1798011587985 ;
          b_err =  TMath::Power(10.642529332523253,2) ;
          c =  18.35846646718483 ;
          c_err =  TMath::Power(1.6103041867447094,2) ;
          d =  -0.09501747239412166 ;
          d_err =  TMath::Power(0.1018222797267054,2) ;
          e =  -0.0946833878842563 ;
          e_err =  TMath::Power(0.002402729508830831,2) ;
        }
        if (missing_mom>=0.2)
        {
          a =  408.547770052618 ;
          a_err =  TMath::Power(36.453926608120966,2) ;
          b =  -426.06473793939705 ;
          b_err =  TMath::Power(37.76476287606027,2) ;
          c =  152.68913846154405 ;
          c_err =  TMath::Power(14.507465096982482,2) ;
          d =  -23.07351286601611 ;
          d_err =  TMath::Power(2.450405880413825,2) ;
          e =  1.3017802399082927 ;
          e_err =  TMath::Power(0.15359726205698873,2) ;
        }
        pol_p = (a * TMath::Power(missing_mom,4)) + (b * TMath::Power(missing_mom,3)) + (c * TMath::Power(missing_mom,2)) + (d * missing_mom) + e;

        pol_p_w = 1.0/((TMath::Power(missing_mom,8)*a_err) + (TMath::Power(missing_mom,6)*b_err) + (TMath::Power(missing_mom,4)*c_err) + (TMath::Power(missing_mom,2)*d_err) + e_err + (16*a*a*TMath::Power(missing_mom,6)*sigma) + (24*a*b*TMath::Power(missing_mom,5)*sigma) + ((16*a*c + 9*b*b)*TMath::Power(missing_mom,4)*sigma) + ((8*a*d + 12*b*c)*TMath::Power(missing_mom,3)*sigma) + ((6*b*d + 4*c*c)*TMath::Power(missing_mom,2)*sigma) + (4*c*d*missing_mom*sigma) + (d*d*sigma));

        // ^^ This weight is wrong. I need to incorperate the cross terms. However, the fit is good enough that the fit error is much smaller than the error simply due to the simulation, so I am going to ignore the fit error and just use the simulation weight.

        h_dx_pol_p->Fill(dx-dx_p_shift,pol_p,weight);
        h_prof_pol_p->Fill(dx-dx_p_shift,pol_p,weight);
        if (dx-dx_p_shift<n_max&&dx-dx_p_shift>n_min)
        {
          h_dx_pol_p_inWindow->Fill(dx-dx_p_shift,pol_p,weight);
        }
        if (dx-dx_p_shift<dx_p_cut)
        {
          h_dx_pol_p_low->Fill(dx-dx_p_shift,pol_p,weight);
          h_prof_pol_p_low->Fill(dx-dx_p_shift,pol_p,weight);
        }
        if (dx-dx_p_shift>=dx_p_cut)
        {
          h_dx_pol_p_high->Fill(dx-dx_p_shift,pol_p,weight);
          h_prof_pol_p_high->Fill(dx-dx_p_shift,pol_p,weight);
        }
      }

      if (fnucl==0.0)
      {
        h_dx_n->Fill(dx-dx_n_shift,weight);
        h_dx_missing_mom_n->Fill(dx-dx_n_shift,missing_mom,weight);
        if (missing_mom<0.2)
        {
          a =  -206.6702353176982 ;
          a_err =  TMath::Power(5.515059269519871,2) ;
          b =  54.385593370006745 ;
          b_err =  TMath::Power(2.922097091242897,2) ;
          c =  -7.751073896982325 ;
          c_err =  TMath::Power(0.5568231435537595,2) ;
          d =  0.37912363572302266 ;
          d_err =  TMath::Power(0.04496132325371255,2) ;
          e =  0.9928502121808195 ;
          e_err =  TMath::Power(0.0012929231472295886,2) ;
        }
        if (missing_mom>=0.2)
        {
          a =  883.3727552389017 ;
          a_err =  TMath::Power(130.37671654356538,2) ;
          b =  -792.0583539603548 ;
          b_err =  TMath::Power(151.20956516982548,2) ;
          c =  230.03074404947085 ;
          c_err =  TMath::Power(64.81973223961523,2) ;
          d =  -27.750236381289398 ;
          d_err =  TMath::Power(12.155125256774213,2) ;
          e =  2.1343364619401908 ;
          e_err =  TMath::Power(0.8406721498772978,2) ;
        }
        pol_n = (a * TMath::Power(missing_mom,4)) + (b * TMath::Power(missing_mom,3)) + (c * TMath::Power(missing_mom,2)) + (d * missing_mom) + e;

        pol_n_w = 1.0/((TMath::Power(missing_mom,8)*a_err) + (TMath::Power(missing_mom,6)*b_err) + (TMath::Power(missing_mom,4)*c_err) + (TMath::Power(missing_mom,2)*d_err) + e_err + (16*a*a*TMath::Power(missing_mom,6)*sigma) + (24*a*b*TMath::Power(missing_mom,5)*sigma) + ((16*a*c + 9*b*b)*TMath::Power(missing_mom,4)*sigma) + ((8*a*d + 12*b*c)*TMath::Power(missing_mom,3)*sigma) + ((6*b*d + 4*c*c)*TMath::Power(missing_mom,2)*sigma) + (4*c*d*missing_mom*sigma) + (d*d*sigma));

        h_dx_pol_n->Fill(dx-dx_n_shift,pol_n,weight);
        h_prof_pol_n->Fill(dx-dx_n_shift,pol_n,weight);
        if (dx-dx_n_shift<n_max&&dx-dx_n_shift>n_min)
        {
          h_dx_pol_n_inWindow->Fill(dx-dx_n_shift,pol_n,weight);
        }
        if (dx-dx_n_shift<dx_n_cut)
        {
          h_dx_pol_n_low->Fill(dx-dx_n_shift,pol_n,weight);
          h_prof_pol_n_low->Fill(dx-dx_n_shift,pol_n,weight);
        }
        if (dx-dx_n_shift>=dx_n_cut)
        {
          h_dx_pol_n_high->Fill(dx-dx_n_shift,pol_n,weight);
          h_prof_pol_n_high->Fill(dx-dx_n_shift,pol_n,weight);
        }
      }

      dx_out = dx-dx_n_shift;
      dy_out = dy;
      T_out->Fill();

  }//end event loop

  for (size_t iev = 0; iev < T->GetEntries(); iev++)
  {
    T->GetEntry(iev);

    missing_mom = TMath::Sqrt(TMath::Power((-beam_e)+(npz)+(epz),2)+TMath::Power((npx)+(epx),2)+TMath::Power((npy)+(epy),2));



  }// end second event loop

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

  c2->Print(outputfile);

  TCanvas *c3 = new TCanvas("c3", "Neutron Profile Fitting", 100,100,700,700);
  c3->Divide(1,2);
  c3->cd(1);
  h_prof_pol_n_low->Draw();

  //TF1 *fitn = new TF1("fitn", "[0] + [1]*x + [2]*TMath::Power(x,2) + [3]*TMath::Power(x,3) + [4]*TMath::Power(x,4) + [5]*TMath::Power(x,5) + [6]*TMath::Power(x,6) + [7]*TMath::Power(x,7) + [8]*TMath::Power(x,8) + [9]*TMath::Power(x,9)", -3.5, 3.0);
  //fitn->SetParameters(1,2,3,4,5,6,7,8,9,10);
  //fitn->SetLineColor(kRed);

  //h_prof_pol_n->Fit("fitn");
  //fitn->Draw("SAMES");
  //gStyle->SetOptFit(1111);

  TF1 *fitn_low = new TF1("fitn_low", "[0]", -3.5, dx_n_cut);
  fitn_low->SetParameters(1.0);
  fitn_low->SetLineColor(kBlue);
  h_prof_pol_n_low->Fit("fitn_low");
  fitn_low->Draw("SAMES");

  c3->cd(2);
  h_prof_pol_n_high->Draw();

  TF1 *fitn_high = new TF1("fitn_high", "[0]*x*x + [1]*x + [2]", dx_n_cut, 3.0);
  fitn_high->SetParameters(1.0,2.0,3.0);
  fitn_high->SetLineColor(kRed);
  h_prof_pol_n_high->Fit("fitn_high");
  fitn_high->Draw("SAMES");

  //h_prof_pol_n->Draw("SAMES");

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

  TCanvas *c4 = new TCanvas("c4", "Proton Profile Fitting", 100,100,700,700);
  c4->Divide(1,2);
  c4->cd(1);
  h_prof_pol_p_low->Draw();

  TF1 *fitp_low = new TF1("fitp_low", "[0] + [1]*cos(x) + [2]*sin(x) + [3]*cos(2*x) + [4]*sin(2*x) + [5]*cos(3*x) + [6]*sin(3*x)", -3.5, dx_p_cut);
  fitp_low->SetParameters(1.0,2.0,3.0,4.0,5.0,6.0,7.0);
  fitp_low->SetLineColor(kBlue);
  h_prof_pol_p_low->Fit("fitp_low");
  fitp_low->Draw("SAMES");

  c4->cd(2);
  h_prof_pol_p_high->Draw();

  TF1 *fitp_high = new TF1("fitp_high", "[0]", dx_p_cut, 3.0);
  fitp_high->SetParameters(0.05);
  fitp_high->SetLineColor(kRed);
  h_prof_pol_p_high->Fit("fitp_high");
  fitp_high->Draw("SAMES");

  //TF1 *fitp = new TF1("fitp", "[0] + [1]*cos(x) + [2]*sin(x) + [3]*cos(2*x) + [4]*sin(2*x) + [5]*cos(3*x) + [6]*sin(3*x) + [7]*cos(4*x) + [8]*sin(4*x)", -3.5, 3.0); // + [5]*cos(3*x) + [6]*sin(3*x)
  //fitp->SetParameters(1,2,3,4,5,6,7,8,9);
  //fitp->SetLineColor(kRed);

  //h_prof_pol_p->Fit("fitp");
  //fitp->Draw("SAMES");
  //gStyle->SetOptFit(1111);

  fout->Write();
}
