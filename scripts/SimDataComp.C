#include <TF1.h>
#include <TAxis.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TError.h>

#include <string>
#include <cmath>
#include <iostream>
#include <TMath.h>

//To run script:
// Go to SBS-GEnII-Analysis directory and type "root"
// Type ".L scripts/SimDataComp.C"
//   This will compile the code and show any errors but will not run the code.
// To run the code, type one of the following commands:
//  for GEN2:  "SimDataComp(2)"
//  for GEN3:  "SimDataComp(3)"
//  for GEN4a: "SimDataComp(4)"
//  for GEN4b: "SimDataComp(5)"

//TH1::SetDefaultSumw2();

TH1D *h_data_dx;
TH1D *h_simIN_dx;
TH1D *h_sim_proton_dx;
TH1D *h_sim_neutron_dx;
TH1D *h_prob_proton_dx;
TH1D *h_prob_neutron_dx;
TH1D *h_prob_bckgrnd_dx;
TH1D *h_prob_proton_dx_polW;
TH1D *h_prob_neutron_dx_polW;
//TF1 *fitn_low_in;
//TF1 *fitn_high_in;
//TF1 *fitp_low_in;
//TF1 *fitp_high_in;

double fitsim( double *x, double *par)
{

  const double dx   = x[0];
  const double Rp   = par[0];
  const double shp  = par[1];
  const double Rn   = par[2];
  const double shn  = par[3];
  const double Rbg  = par[4];
  const double shbg = par[5];

  //condition ? expression_if_true : expression_if_false
  const double p  = (h_sim_proton_dx  ? h_sim_proton_dx->Interpolate(dx - shp)  :0.0);
  const double n  = (h_sim_neutron_dx ? h_sim_neutron_dx->Interpolate(dx - shn) :0.0);
  const double bg = (h_simIN_dx       ? h_simIN_dx->Interpolate(dx - shbg)      :0.0);

  return Rp * p + Rn * n + Rbg * bg;

}//end fitsim

double fitAsym(double *xA, double *parA)
{

  const double dxA  = xA[0];
  const double Ap   = parA[0];
  const double An   = parA[1];
  const double Abg  = parA[2];

  //condition ? expression_if_true : expression_if_false
  const double Prob_p  = (h_prob_proton_dx  ? h_prob_proton_dx->Interpolate(dxA)  :0.0);
  const double Prob_n  = (h_prob_neutron_dx ? h_prob_neutron_dx->Interpolate(dxA) :0.0);
  const double Prob_bg = (h_prob_bckgrnd_dx ? h_prob_bckgrnd_dx->Interpolate(dxA) :0.0);

  return (Ap * Prob_p) + (An * Prob_n) + (Abg * Prob_bg);

}//end fitAsym

void SimDataComp(int kin)
{

  TString data_file;
  TString nucleon_sim_file;
  TString inel_sim_file;
  TString pol_func_file;
  TString N2_dilution_file;
  TString title_words;

  int npar = 6;
  double dx_min_d, dx_max_d;
  int dx_min_i, dx_max_i;
  double pol_beam, pol_beam_err, pol_targ, pol_targ_err;

  if(kin == 2)
  {

    data_file = "outfiles/parsed_GEn_pass3_GEN2_He3_dxdy.root";
    nucleon_sim_file = "outfiles/parsed_SIM_GEn_GEN2_He3_dxdy.root";
    inel_sim_file = "outfiles/parsed_SIM_IN_GEn_GEN2_He3_dxdy.root";
    pol_func_file = "outfiles/parsed_GEn_pass2_GEN2_simulation.root";
    N2_dilution_file = "outfiles/N2_Corr_SIM_GEn_GEN2_He3_dxdy.root";
    title_words = "GEN2";
    dx_min_d = -2.8;
    dx_min_i = -3;
    dx_max_d = 2.0;
    dx_max_i = 2;
    pol_beam = 0.8409;
    pol_beam_err = 0.0018;
    pol_targ = 0.3557;
    pol_targ_err = 0.0179;
    cout<<"\nHi! You are analyzing GEN2!\n";

  }

  if(kin == 3)
  {

    data_file = "outfiles/parsed_GEn_pass3_GEN3_He3_dxdy.root";
    nucleon_sim_file = "outfiles/parsed_SIM_GEn_GEN3_He3_dxdy.root";
    inel_sim_file = "outfiles/parsed_SIM_IN_GEn_GEN3_He3_dxdy.root";
    pol_func_file = "outfiles/parsed_GEn_pass2_GEN3_simulation.root";
    N2_dilution_file = "outfiles/N2_Corr_SIM_GEn_GEN3_He3_dxdy.root";
    title_words = "GEN3";
    dx_min_d = -2.5;
    dx_min_i = -3;
    dx_max_d = 1.7;
    dx_max_i = 2;
    pol_beam = 0.8649;
    pol_beam_err = 0.0008;
    pol_targ = 0.4212;
    pol_targ_err = 0.011;
    cout<<"\nHi! You are analyzing GEN3!\n";

  }

  if(kin == 4)
  {

    data_file = "outfiles/parsed_GEn_pass2_GEN4a_He3_dxdy.root";
    nucleon_sim_file = "outfiles/parsed_SIM_GEn_GEN4_He3_dxdy.root";
    inel_sim_file = "outfiles/parsed_SIM_IN_GEn_GEN4_He3_dxdy.root";
    pol_func_file = "outfiles/parsed_GEn_pass2_GEN4a_simulation.root";
    N2_dilution_file = "outfiles/N2_Corr_SIM_GEn_GEN4_He3_dxdy.root";
    title_words = "GEN4a";
    dx_min_d = -3.0;
    dx_min_i = -3;
    dx_max_d = 2.0;
    dx_max_i = 2;
    pol_beam = 0.8293;
    pol_beam_err = 0.0010;
    pol_targ = 0.4908;
    pol_targ_err = 0.0202;
    cout<<"\nHi! You are analyzing GEN4a!\n";

  }

  if(kin == 5)
  {

    data_file = "outfiles/parsed_GEn_pass2_GEN4b_He3_dxdy.root";
    nucleon_sim_file = "outfiles/parsed_SIM_GEn_GEN4_He3_dxdy.root";
    inel_sim_file = "outfiles/parsed_SIM_IN_GEn_GEN4_He3_dxdy.root";
    pol_func_file = "outfiles/parsed_GEn_pass2_GEN4b_simulation.root";
    N2_dilution_file = "outfiles/N2_Corr_SIM_GEn_GEN4_He3_dxdy.root";
    title_words = "GEN4b";
    dx_min_d = -3.0;
    dx_min_i = -3;
    dx_max_d = 2.0;
    dx_max_i = 2;
    pol_beam = 0.8293;
    pol_beam_err = 0.0010;
    pol_targ = 0.4908;
    pol_targ_err = 0.0202;
    cout<<"\nHi! You are analyzing GEN4b!\n";

  }

  //else
  //{
    //data_file = "null";
    //nucleon_sim_file = "null";
    //inel_sim_file = "null";
    //title_words = "null";
    //dx_min_d = 0.0;
    //dx_min_i = 0;
    //dx_max_d = 0.0;
    //dx_max_i = 0.0;
  //}

  //gErrorIgnoreLevel = kError;

  int numberBins = 185;
  int asymBinning = 192;

  TChain* T_data = new TChain("T_data");
  T_data->Add(data_file);

  Double_t dx;      T_data->SetBranchAddress("dx", &dx);
  Double_t dy;      T_data->SetBranchAddress("dy", &dy);
  int helicity;     T_data->SetBranchAddress("helicity", &helicity);

  if(T_data->GetEntries()==0)
  {
    std::cerr << "\n --- No ROOT file found!! --- \n\n";
    throw;
  }
  else std::cout << "\nFound " << T_data->GetEntries() << " events. \n";

  h_data_dx = new TH1D("h_data_dx", ";dx", numberBins, dx_min_d, dx_max_d);
  h_data_dx->GetXaxis()->SetTitle("dx [m]");
  h_data_dx->Sumw2();

  TH1D* h_neg_hel_dx = new TH1D("h_neg_hel_dx",";-hel", asymBinning, dx_min_d, dx_max_d);
  h_neg_hel_dx->GetXaxis()->SetTitle("dx [m]");
  h_neg_hel_dx->SetLineColor(kRed);
  h_neg_hel_dx->Sumw2();

  TH1D* h_pos_hel_dx = new TH1D("h_pos_hel_dx",";+hel", asymBinning, dx_min_d, dx_max_d);
  h_pos_hel_dx->GetXaxis()->SetTitle("dx [m]");
  h_pos_hel_dx->SetLineColor(kBlue);
  h_pos_hel_dx->Sumw2();

  double Npos = 0.0;
  double Nneg = 0.0;

  double dxEllipse;

  for (size_t iev = 0; iev < T_data->GetEntries(); iev++)
  {
    T_data->GetEntry(iev);

    h_data_dx->Fill(dx);
    dxEllipse = (((dx-0.0)*(dx-0.0))/(0.99*0.99)) + (((dy-0.0)*(dy-0.0))/(0.99*0.99));

    if(helicity==1)
    {
      h_pos_hel_dx->Fill(dx);

      if(dxEllipse<=1.0)
      {
        Npos++;
      }

    }

    if(helicity==-1)
    {
      h_neg_hel_dx->Fill(dx);

      if(dxEllipse<=1.0)
      {
        Nneg++;
      }
    }

  }//end loop over events

  std::cout << "\nFound " << Npos << " positive helicity events in the neutron window. \n";
  std::cout << "\nFound " << Nneg << " negative helicity events in the neutron window. \n";

  TChain* T_sim = new TChain("T_sim");
  T_sim->Add(nucleon_sim_file);

  Double_t dx_sim;     T_sim->SetBranchAddress("dx", &dx_sim);
  Double_t fnucl;      T_sim->SetBranchAddress("fnucl", &fnucl);
  Double_t weight;     T_sim->SetBranchAddress("weight", &weight);

  if(T_sim->GetEntries()==0)
  {
    std::cerr << "\n --- No ROOT file found!! --- \n\n";
    throw;
  }
  else std::cout << "\nFound " << T_sim->GetEntries() << " events. \n";

  h_sim_proton_dx = new TH1D("h_sim_proton_dx", ";dx_sim_p", numberBins, dx_min_d, dx_max_d);
  h_sim_proton_dx->GetXaxis()->SetTitle("dx [m]");
  h_sim_proton_dx->SetLineColor(kGreen);

  h_sim_neutron_dx = new TH1D("h_sim_neutron_dx", ";dx_sim_n", numberBins, dx_min_d, dx_max_d);
  h_sim_neutron_dx->GetXaxis()->SetTitle("dx [m]");
  h_sim_neutron_dx->SetLineColor(kRed);

  for (size_t iev = 0; iev < T_sim->GetEntries(); iev++)
  {
    T_sim->GetEntry(iev);

    if(fnucl == 1.0)
    {
      h_sim_proton_dx->Fill(dx_sim,weight);
    }

    if(fnucl == 0.0)
    {
      h_sim_neutron_dx->Fill(dx_sim,weight);
    }

  }//end loop over events

  TChain* T_simIN = new TChain("T_sim");
  T_simIN->Add(inel_sim_file);

  Double_t dx_simIN;      T_simIN->SetBranchAddress("dx", &dx_simIN);
  Double_t fnucl_IN;      T_simIN->SetBranchAddress("fnucl", &fnucl_IN);
  Double_t weight_IN;     T_simIN->SetBranchAddress("weight", &weight_IN);

  if(T_simIN->GetEntries()==0)
  {
    std::cerr << "\n --- No ROOT file found!! --- \n\n";
    throw;
  }
  else std::cout << "\nFound " << T_simIN->GetEntries() << " events. \n";

  h_simIN_dx = new TH1D("h_simIN_dx", ";dxIN", numberBins, dx_min_d, dx_max_d);
  h_simIN_dx->GetXaxis()->SetTitle("dx [m]");
  h_simIN_dx->SetLineColor(kBlue);

  for (size_t iev = 0; iev < T_simIN->GetEntries(); iev++)
  {
    T_simIN->GetEntry(iev);

    h_simIN_dx->GetSumw2();
    weight_IN = weight_IN * 1e31;

    h_simIN_dx->Fill(dx_simIN,weight_IN);

  }//end loop over events

  TFile *polFile = TFile::Open(pol_func_file);
  TF1 *fitp = nullptr;
  TF1 *fitn = nullptr;
  polFile->GetObject("fitp", fitp);
  polFile->GetObject("fitn", fitn);

  TFile *N2File = TFile::Open(N2_dilution_file);
  TH1D *hN2dilution_p = nullptr;
  N2File->GetObject("hN2dilution_p", hN2dilution_p);

  //std::cout << "sqrt(n) error " << h_simIN_dx->GetBinError(250) << std::endl;
  //std::cout << "Get Sumw2 Error: " << h_simIN_dx->GetSumw2() << std::endl;

  //******************************************************

  //Normalize each distribution using the sum of the dataset
  double scale = h_data_dx->Integral();
  //h_data_dx->Scale(1.0/h_data_dx->Integral());
  h_data_dx->Scale(1.0/scale);
  h_pos_hel_dx->Scale(1.0/h_pos_hel_dx->Integral());
  h_neg_hel_dx->Scale(1.0/h_neg_hel_dx->Integral());
  h_sim_proton_dx->Scale(1.0/h_sim_proton_dx->Integral());
  h_sim_neutron_dx->Scale(1.0/h_sim_neutron_dx->Integral());
  h_simIN_dx->Scale(1.0/h_simIN_dx->Integral());

  int nbins = h_data_dx->GetNbinsX();
  double xmin = h_data_dx->GetXaxis()->GetBinLowEdge(1);
  double xmax = h_data_dx->GetXaxis()->GetBinUpEdge(nbins);

  TF1 *FitFunc = new TF1("FitFunc",&fitsim,dx_min_i,dx_max_i,6); //-6,4,6

  FitFunc->SetNpx(numberBins);

  //----- GEN2 -----
  double startpar[] = {1.0,-0.5,0.5,-0.7,0.1,-1.0};
  FitFunc->SetParameters(startpar);
  FitFunc->SetParLimits(0,0.0,100);   // proton scale
  FitFunc->SetParLimits(1,-4.0,4.0);  // proton shift
  FitFunc->SetParLimits(2,0.0,100);   // neutron scale
  FitFunc->SetParLimits(3,-4.0,4.0);  // neutron shift
  FitFunc->SetParLimits(4,0.0,100);   // background scale
  FitFunc->SetParLimits(5,-4.0,4.0);  // background shift

  //----- GEN3 -----
  //double startpar[] = {0.2,-0.5,0.1,-0.7,0.4,-1.0};
  //FitFunc->SetParameters(startpar);
  //FitFunc->SetParLimits(0,0.1,5.0);   // proton scale
  //FitFunc->SetParLimits(1,-1.0,1.0);  // proton shift
  //FitFunc->SetParLimits(2,0.0,5.0);   // neutron scale
  //FitFunc->SetParLimits(3,-1.0,0.0);  // neutron shift
  //FitFunc->SetParLimits(4,0.0,0.5);   // background scale
  //FitFunc->SetParLimits(5,-1.0,2.0);   // background shift

  h_data_dx->Fit(FitFunc,"0","",xmin,xmax);

  std::cout << "Proton Shift: " << FitFunc->GetParameter(0)*FitFunc->GetParameter(1) << std::endl;
  std::cout << "Neutron Shift: " << FitFunc->GetParameter(2)*FitFunc->GetParameter(3) << std::endl;
  std::cout << "Background Shift: " << FitFunc->GetParameter(4)*FitFunc->GetParameter(5) << std::endl;

  TH1D* shifted_h_sim_proton_dx = (TH1D*)h_sim_proton_dx->Clone("shifted_h_sim_proton_dx");
  shifted_h_sim_proton_dx->Reset("ICESM");
  //shifted_h_sim_proton_dx->Sumw2();

  TH1D* shifted_h_sim_neutron_dx = (TH1D*)h_sim_neutron_dx->Clone("shifted_h_sim_neutron_dx");
  shifted_h_sim_neutron_dx->Reset("ICESM");
  //shifted_h_sim_neutron_dx->Sumw2();

  TH1D* shifted_h_simIN_dx = (TH1D*)h_simIN_dx->Clone("shifted_h_simIN_dx");
  shifted_h_simIN_dx->Reset("ICESM");
  //shifted_h_simIN_dx->Sumw2();

  for(int i =1; i <= nbins; i++)
  {

    const double b_p  = shifted_h_sim_proton_dx->GetBinCenter(i);
    const double b_n  = shifted_h_sim_neutron_dx->GetBinCenter(i);
    const double b_in = shifted_h_simIN_dx->GetBinCenter(i);

    const double v_p         = h_sim_proton_dx->Interpolate(b_p  - FitFunc->GetParameter(1));
    const double v_p_error   = h_sim_proton_dx->GetBinError(b_p  - FitFunc->GetParameter(1));
    const double v_n         = h_sim_neutron_dx->Interpolate(b_n - FitFunc->GetParameter(3));
    const double v_n_error   = h_sim_proton_dx->GetBinError(b_n  - FitFunc->GetParameter(3));
    const double v_in        = h_simIN_dx->Interpolate(b_in      - FitFunc->GetParameter(5));
    const double v_in_error  = h_simIN_dx->GetBinError(b_in      - FitFunc->GetParameter(5));

    shifted_h_sim_proton_dx->SetBinContent(i, v_p);
    shifted_h_sim_proton_dx->SetBinError(i, v_p_error);
    shifted_h_sim_neutron_dx->SetBinContent(i, v_n);
    shifted_h_sim_neutron_dx->SetBinError(i, v_n_error);
    shifted_h_simIN_dx->SetBinContent(i, v_in);
    shifted_h_simIN_dx->SetBinError(i, v_in_error);

  }//end loop over bins

  shifted_h_sim_proton_dx->Scale(FitFunc->GetParameter(0));
  shifted_h_sim_neutron_dx->Scale(FitFunc->GetParameter(2));
  shifted_h_simIN_dx->Scale(FitFunc->GetParameter(4));

  TH1F* h_total_dx = new TH1F("h_total_dx","",nbins,xmin,xmax);

  for(int ibin = 0; ibin < nbins; ibin++)
  {
    h_total_dx->SetBinContent(ibin, shifted_h_sim_neutron_dx->GetBinContent(ibin) + shifted_h_sim_proton_dx->GetBinContent(ibin) + shifted_h_simIN_dx->GetBinContent(ibin));
  }

  h_data_dx                ->Scale(scale);
  h_pos_hel_dx             ->Scale(scale);
  h_neg_hel_dx             ->Scale(scale);
  shifted_h_sim_proton_dx  ->Scale(scale);
  shifted_h_sim_neutron_dx ->Scale(scale);
  shifted_h_simIN_dx       ->Scale(scale);
  h_total_dx               ->Scale(scale);

  //Account for diluitions

  TH1D* scaled_h_sim_nucleons = (TH1D*)shifted_h_sim_proton_dx->Clone("scaled_h_sim_nucleons");
  scaled_h_sim_nucleons->Add(shifted_h_sim_neutron_dx);

  TCanvas *step1 = new TCanvas("step1","testing",100,100,1500,500);
  step1->cd();
  scaled_h_sim_nucleons->Draw("HIST");

  TH1D* scaled_hN2dilution = (TH1D*)hN2dilution_p->Clone("scaled_hN2dilution");
  scaled_hN2dilution->SetLineColor(kRed); //line, marker, fill
  scaled_hN2dilution->Reset("ICESM");

  scaled_hN2dilution->Multiply(scaled_h_sim_nucleons);

  //Create residual and probability plots

  double totalentries = h_data_dx->GetEntries();
  int h_Nbins = h_data_dx->GetNbinsX();
  double h_minX = h_data_dx->GetXaxis()->GetXmin();
  double h_maxX = h_data_dx->GetXaxis()->GetXmax();

  TH1D* h_resid = new TH1D("h_resid",";Residuals",h_Nbins,h_minX,h_maxX);
  h_resid->GetXaxis()->SetTitle("Residuals (Data - Fit) [m]");

  h_prob_proton_dx = new TH1D("h_prob_proton_dx:",";prob dx proton",h_Nbins,h_minX,h_maxX);
  h_prob_proton_dx->GetXaxis()->SetTitle("dx [m]");
  h_prob_proton_dx->GetYaxis()->SetTitle("Proton Probability");
  h_prob_proton_dx->GetYaxis()->SetRangeUser(0.0,1.0);

  h_prob_neutron_dx = new TH1D("h_prob_neutron_dx:",";prob dx neutron",h_Nbins,h_minX,h_maxX);
  h_prob_neutron_dx->GetXaxis()->SetTitle("dx [m]");
  h_prob_neutron_dx->GetYaxis()->SetTitle("Neutron Probability");
  h_prob_neutron_dx->GetYaxis()->SetRangeUser(0.0,1.0);

  h_prob_bckgrnd_dx = new TH1D("h_prob_bckgrnd_dx:",";prob dx background",h_Nbins,h_minX,h_maxX);
  h_prob_bckgrnd_dx->GetXaxis()->SetTitle("dx [m]");
  h_prob_bckgrnd_dx->GetYaxis()->SetTitle("Background Probability");
  h_prob_bckgrnd_dx->GetYaxis()->SetRangeUser(0.0,1.0);

  TH1D* h_asym = new TH1D("h_asym",";Raw Asymmetry",h_Nbins,h_minX,h_maxX);
  h_asym->GetXaxis()->SetTitle("dx [m]");
  h_asym->SetTitle("Raw Asymmetry (N+ - N-)/(N+ + N-)");

  TH1D* hAsymDiff = (TH1D*) h_pos_hel_dx->Clone("hAsymDiff");
  hAsymDiff->Sumw2();
  hAsymDiff->Add(h_neg_hel_dx, -1.0);

  TH1D* hAsymSum = (TH1D*) h_pos_hel_dx->Clone("hAsymSum");
  hAsymSum->Sumw2();
  hAsymSum->Add(h_neg_hel_dx);

  TH1D* hAsym = (TH1D*) hAsymDiff->Clone("hAsym");
  hAsym->Divide(hAsymSum);
  hAsym->Sumw2();
  //hAsym->Rebin();

  TH1D* h_fullProb = new TH1D("h_fullProb","100 Percent",h_Nbins,h_minX,h_maxX);
  h_fullProb->GetYaxis()->SetRangeUser(0.0,1.0);
  h_fullProb->Sumw2();

  for (int bin = 0; bin < h_Nbins; bin++)
  {

    double hist_val    = h_data_dx->GetBinContent(bin);
    double hist_error  = h_data_dx->GetBinError(bin);
    double fit_val     = h_total_dx->GetBinContent(bin);

    double new_val = (hist_val - fit_val);
    h_resid->SetBinContent(bin,new_val);
    h_resid->SetBinError(bin,hist_error);

    double c_p       = shifted_h_sim_proton_dx->GetBinContent(bin);
    double c_p_err   = shifted_h_sim_proton_dx->GetBinError(bin);
    double c_n       = shifted_h_sim_neutron_dx->GetBinContent(bin);
    double c_n_err   = shifted_h_sim_neutron_dx->GetBinError(bin);
    double c_bg      = shifted_h_simIN_dx->GetBinContent(bin); //Backgrounds with asymmetry contributions
    double c_bg_err  = shifted_h_simIN_dx->GetBinError(bin);
    double c_dil     = scaled_hN2dilution->GetBinContent(bin); //Backgrounds with no asymmetry contribution
    double c_dil_err = scaled_hN2dilution->GetBinError(bin);

    double P_tot = c_p + c_n + c_bg + c_dil;

    //For asymmetry fit, we only care about if background contributes to asymmetry, so there is no need to track the probability of being a dilution particle. The total on the denominator should include dilution considerations.

    double P_p       = 0.0;
    double P_p_err   = 1.0;
    double P_n       = 0.0;
    double P_n_err   = 1.0;
    double P_bg      = 0.0;
    double P_bg_err  = 1.0;

    if (P_tot != 0.0)
    {

      P_p  = c_p / P_tot;
      P_n  = c_n / P_tot;
      P_bg = c_bg / P_tot;

      P_p_err  = TMath::Sqrt( (c_n+c_bg+c_dil)*(c_n+c_bg+c_dil)*c_p_err*c_p_err + c_p*c_p*c_n_err*c_n_err + c_p*c_p*c_bg_err*c_bg_err + c_p*c_p*c_dil_err*c_dil_err ) / ( (c_p+c_n+c_bg+c_dil)*(c_p+c_n+c_bg+c_dil) );
      P_n_err  = TMath::Sqrt( (c_bg+c_p+c_dil)*(c_bg+c_p+c_dil)*c_n_err*c_n_err + c_n*c_n*c_p_err*c_p_err + c_n*c_n*c_bg_err*c_bg_err + c_n*c_n*c_dil_err*c_dil_err ) / ( (c_p+c_n+c_bg+c_dil)*(c_p+c_n+c_bg+c_dil) );
      P_bg_err = TMath::Sqrt( (c_n+c_p+c_dil)*(c_n+c_p+c_dil)*c_bg_err*c_bg_err + c_bg*c_bg*c_p_err*c_p_err + c_bg*c_bg*c_p_err*c_p_err + c_bg*c_bg*c_dil_err*c_dil_err ) / ( (c_p+c_n+c_bg+c_dil)*(c_p+c_n+c_bg+c_dil) );

    }//end if total probability is nonzero

    //std::cout << "Proton probability: " << P_p << std::endl;
    h_prob_proton_dx->SetBinContent(bin,P_p);
    h_prob_proton_dx->SetBinError(bin,P_p_err);
    h_prob_neutron_dx->SetBinContent(bin,P_n);
    h_prob_neutron_dx->SetBinError(bin,P_n_err);
    h_prob_bckgrnd_dx->SetBinContent(bin,P_bg);
    h_prob_bckgrnd_dx->SetBinError(bin,P_bg_err);
    h_fullProb->SetBinContent(bin,1.0);

  }//end loop over bins

  h_resid           -> SetEntries(totalentries);
  h_prob_proton_dx  -> SetEntries(totalentries);
  h_prob_neutron_dx -> SetEntries(totalentries);
  h_prob_bckgrnd_dx -> SetEntries(totalentries);
  h_fullProb        -> SetEntries(totalentries);

  double pol_combo = pol_beam * pol_targ;
  h_fullProb->Add(hN2dilution_p,-1.0);
  h_fullProb->Add(h_prob_bckgrnd_dx,-1.0);

  TCanvas *step3 = new TCanvas("step3","1 - (N2 Dilution) - (Inel Prob)",100,100,1500,500);
  step3->cd();
  h_fullProb->Draw("HIST");

  //Apply effective nucleon polarization and beam+target polarization
  h_prob_proton_dx  -> Multiply(fitp, pol_combo);
  h_prob_neutron_dx -> Multiply(fitn, pol_combo);

  //Subtract off dilution
  //h_prob_proton_dx  -> Divide(h_fullProb);
  //h_prob_neutron_dx -> Divide(h_fullProb);

  TCanvas *step5 = new TCanvas("step5","Nucleon Probability Scaled by Polarization and Dilution",100,100,1500,500);
  step5->Divide(1,2);
  step5->cd(1);
  h_prob_proton_dx->Draw("E");
  step5->cd(2);
  h_prob_neutron_dx->Draw("E");

  TF1 *AsymFitFunc = new TF1("AsymFitFunc",&fitAsym,dx_min_i,dx_max_i,3); //-6,3,3

  AsymFitFunc->SetNpx(numberBins);
  double Asymstartpar[] = {0.0,0.04,0.0};
  AsymFitFunc->SetParameters(Asymstartpar);
  AsymFitFunc->SetParLimits(0,-0.1,2.0); // proton asymmetry
  AsymFitFunc->SetParLimits(1,0.0,0.4);  // neutron asymmetry
  AsymFitFunc->SetParLimits(2,-0.1,0.1); // background asymmetry

  hAsym->Fit(AsymFitFunc,"0","",h_minX,h_maxX);

  TH1D* scaled_proton_prob  = (TH1D*)h_prob_proton_dx->Clone("scaled_proton_prob");
  TH1D* scaled_neutron_prob = (TH1D*)h_prob_neutron_dx->Clone("scaled_neutron_prob");
  TH1D* scaled_bckgrnd_prob = (TH1D*)h_prob_bckgrnd_dx->Clone("scaled_bckgrnd_prob");

  scaled_proton_prob  -> Scale(AsymFitFunc->GetParameter(0));
  scaled_neutron_prob -> Scale(AsymFitFunc->GetParameter(1));
  scaled_bckgrnd_prob -> Scale(AsymFitFunc->GetParameter(2));

  std::cout << "Proton Asymmetry: "     << AsymFitFunc->GetParameter(0) << std::endl;
  std::cout << "Neutron Asymmetry: "    << AsymFitFunc->GetParameter(1) << std::endl;
  std::cout << "Background Asymmetry: " << AsymFitFunc->GetParameter(2) << std::endl;

  //Create canvas and draw all the histograms

  TCanvas *c1 = new TCanvas("c1", "dx", 100,100,1600,800);
  gStyle->SetOptStat(0);

  TPad *pad1 = new TPad("pad1","Pad with the plot",0.0,0.3,1.0,1.0);
  TPad *pad2 = new TPad("pad2","Pad with resid",0.0,0.0,1.0,0.3);

  pad1->Draw();
  pad2->Draw();
  pad1->SetBottomMargin(0.1); // Upper and lower plot are not joined
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);

  pad1->cd();
  pad1->SetTitle(title_words);
  h_total_dx->Draw("HIST");
  shifted_h_simIN_dx->Draw("HIST SAMES");
  shifted_h_sim_proton_dx->Draw("HIST SAMES");
  shifted_h_sim_neutron_dx->Draw("HIST SAMES");
  h_data_dx->Draw("E SAMES");
  h_data_dx->GetXaxis()->SetTitle("dx [m]");

  auto legend = new TLegend(0.55,0.70,0.99,0.99);
  legend->SetTextSize(0.03);
  legend->SetHeader("Fitting","C");
  legend->AddEntry(h_data_dx,"Data","lep");
  legend->AddEntry(shifted_h_sim_proton_dx,Form("Rp=%.3g, sh_p=%.3g m",FitFunc->GetParameter(0),FitFunc->GetParameter(0)*FitFunc->GetParameter(1)),"le");
  legend->AddEntry(shifted_h_sim_neutron_dx,Form("Rn=%.3g, sh_n=%.3g m",FitFunc->GetParameter(2),FitFunc->GetParameter(2)*FitFunc->GetParameter(3)),"le");
  legend->AddEntry(shifted_h_simIN_dx,Form("Rbg=%.3g, sh_bg=%.3g m",FitFunc->GetParameter(4),FitFunc->GetParameter(4)*FitFunc->GetParameter(5)),"le");
  legend->AddEntry(h_total_dx,"Total Fit: (Rp*h_p+sh_p)+(Rn*h_n+sh_n)+(Rbg*h_bg+sh_bg)","l");
  legend->Draw();

  pad2->cd();
  h_resid->Draw("E");
  h_resid->GetXaxis()->SetTitle("dx [m]");
  h_resid->GetYaxis()->SetTitle("Residuals (Data - Fit)");

  c1->cd();  // c1 is the TCanvas
  TPad *pad5 = new TPad("all","all",0,0,1,1);
  pad5->SetFillStyle(4000);  // transparent
  pad5->Draw();
  pad5->cd();

  TLatex *lat = new TLatex();
  lat->DrawLatexNDC(.4,.95,title_words);

  TCanvas *c2 = new TCanvas("c2","Probabilities",100,100,1500,500);
  c2->Divide(3,1);

  c2->cd(1);
  h_prob_proton_dx->Draw();
  c2->cd(2);
  h_prob_neutron_dx->Draw();
  c2->cd(3);
  h_prob_bckgrnd_dx->Draw();

  TCanvas *c3 = new TCanvas("c3","Asymmetry Fit",100,100,800,500);
  //c3->Divide(3,1);
  c3->cd();
  gStyle->SetOptStat(0);
  c3->Update();

  //h_pos_hel_dx->Draw("E");
  //h_neg_hel_dx->Draw("E SAMES");
  //auto legendHEL = new TLegend(0.7,0.8,0.99,0.99);
  //legendHEL->SetTextSize(0.03);
  //legendHEL->AddEntry(h_pos_hel_dx,"Positive Helicity","lep");
  //legendHEL->AddEntry(h_neg_hel_dx,"Negative Helicity","lep");
  //legendHEL->Draw();

  hAsym->Draw("E");
  AsymFitFunc->Draw("SAMES");

  auto legendA = new TLegend(0.1,0.8,0.9,0.9);
  legendA->SetTextSize(0.02);
  legendA->AddEntry(hAsym,"Data","lep");
  legendA->AddEntry(AsymFitFunc,Form("A(dx) = %.4gPp(dx) + %.4gPn(dx) + %.4gPbg(dx)",AsymFitFunc->GetParameter(0),AsymFitFunc->GetParameter(1),AsymFitFunc->GetParameter(2)),"l");
  legendA->Draw();

  //c3->cd(3);
  //scaled_proton_prob->SetTitle("Scaled Probabilities");
  //gPad->DrawFrame(dx_min_d, -0.05, dx_max_d, 0.1);
  //scaled_neutron_prob->SetLineColor(kBlue);
  //scaled_neutron_prob->Draw("E SAME");
  //scaled_proton_prob->SetLineColor(kRed);
  //scaled_proton_prob->Draw("E SAME");
  //scaled_bckgrnd_prob->SetLineColor(kGreen);
  //scaled_bckgrnd_prob->Draw("E SAME");
  //auto legendProb = new TLegend(0.7,0.8,0.9,0.9);
  //legendProb->SetTextSize(0.02);
  //legendProb->AddEntry(scaled_proton_prob, "Proton", "l");
  //legendProb->AddEntry(scaled_neutron_prob, "Neutron", "l");
  //legendProb->AddEntry(scaled_bckgrnd_prob, "Background", "l");
  //legendProb->Draw();

  TCanvas *c4 = new TCanvas("c4","Asymmetry",100,100,1000,1000);
  c4->Divide(3,1);
  c4->cd(1);
  shifted_h_sim_proton_dx->Draw("E");
  gStyle->SetOptStat(11001111);
  gPad->Update();
  c4->cd(2);
  shifted_h_sim_neutron_dx->Draw("E");
  gStyle->SetOptStat(11001111);
  gPad->Update();
  c4->cd(3);
  shifted_h_simIN_dx->Draw("E");
  gStyle->SetOptStat(11001111);
  gPad->Update();

  TCanvas *c5 = new TCanvas("c5","AsymmetryHistOperator",100,100,1000,1000);
  c5->cd();
  hAsym->Draw("E");
  c5->SetGrid();
  c5->Update();

  //TCanvas *c6 = new TCanvas("c6","testing",100,100,1500,500);
  //c6->cd();
  //h_data_asym_raw_early->Draw();

  delete T_data;
  delete T_sim;
  delete T_simIN;

}//end main
