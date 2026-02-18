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

//TH1::SetDefaultSumw2();

TH1D *h_data_dx;
TH1D *h_simIN_dx;
TH1D *h_sim_proton_dx;
TH1D *h_sim_neutron_dx;
TH1D *h_prob_proton_dx;
TH1D *h_prob_neutron_dx;
TH1D *h_prob_bckgrnd_dx;

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

void SimDataComp()
{

  int kin = 2;
  auto data_file = "null";
  auto nucleon_sim_file = "null";
  auto inel_sim_file = "null";
  auto title_words = "null";

  int npar = 6;
  double dx_min_d, dx_max_d;
  int dx_min_i, dx_max_i;

  if(kin == 2)
  {

    data_file = "outfiles/parsed_GEn_pass2_GEN2_He3_dxdy.root";
    nucleon_sim_file = "outfiles/parsed_SIM_GEn_GEN2_He3_dxdy.root";
    inel_sim_file = "outfiles/parsed_SIM_IN_GEn_GEN2_He3_dxdy.root";
    title_words = "GEN2";
    dx_min_d = -5.0;
    dx_min_i = -5;
    dx_max_d = 2.0;
    dx_max_i = 2;

  }

  if(kin == 3)
  {

    data_file = "outfiles/parsed_GEn_pass2_GEN3_He3_dxdy.root";
    nucleon_sim_file = "outfiles/parsed_SIM_GEn_GEN3_He3_dxdy.root";
    inel_sim_file = "outfiles/parsed_SIM_IN_GEn_GEN3_He3_dxdy.root";
    title_words = "GEN3";
    dx_min_d = -6.0;
    dx_min_i = -6;
    dx_max_d = 4.0;
    dx_max_i = 4;

  }

  if(kin == 4)
  {

    data_file = "outfiles/parsed_GEn_pass2_GEN4a_He3_dxdy.root";
    nucleon_sim_file = "outfiles/parsed_SIM_GEn_GEN4a_He3_dxdy.root";
    inel_sim_file = "outfiles/parsed_SIM_IN_GEn_GEN4_He3_dxdy.root";
    title_words = "GEN4a";
    dx_min_d = -3.0;
    dx_min_i = -3;
    dx_max_d = 2.0;
    dx_max_i = 2;

  }

  if(kin == 5)
  {

    data_file = "outfiles/parsed_GEn_pass2_GEN4b_He3_dxdy.root";
    nucleon_sim_file = "outfiles/parsed_SIM_GEn_GEN4b_He3_dxdy.root";
    inel_sim_file = "outfiles/parsed_SIM_IN_GEn_GEN4_He3_dxdy.root";
    title_words = "GEN4b";
    dx_min_d = -3.0;
    dx_min_i = -3;
    dx_max_d = 2.0;
    dx_max_i = 2;

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

  gErrorIgnoreLevel = kError;

  int numberBins = 250;

  TChain* T_data = new TChain("T_data");
  T_data->Add(data_file);

  Double_t dx;      T_data->SetBranchAddress("dx", &dx);
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

  TH1D* h_neg_hel_dx = new TH1D("h_neg_hel_dx",";-hel", numberBins, dx_min_d, dx_max_d);
  h_neg_hel_dx->GetXaxis()->SetTitle("dx [m]");
  h_neg_hel_dx->SetLineColor(kRed);
  //h_neg_hel_dx->Sumw2();

  TH1D* h_pos_hel_dx = new TH1D("h_pos_hel_dx",";+hel", numberBins, dx_min_d, dx_max_d);
  h_pos_hel_dx->GetXaxis()->SetTitle("dx [m]");
  h_pos_hel_dx->SetLineColor(kBlue);
  //h_pos_hel_dx->Sumw2();

  int helPosArray[numberBins];
  int helNegArray[numberBins];

  double binSize = 10.0 / numberBins;

  for(int j = 0; j < numberBins; j++)
  {
    helPosArray[j] = 0;
    helNegArray[j] = 0;
  }

  for (size_t iev = 0; iev < T_data->GetEntries(); iev++)
  {
    T_data->GetEntry(iev);

    h_data_dx->Fill(dx);

    int binAt = (int) ((dx + 5.0) / binSize); // what???

    if(helicity==1)  //h_pos_hel_dx->Fill(dx);
    {
      h_pos_hel_dx->Fill(dx);
      helPosArray[binAt]++;
    }

    if(helicity==-1)  //h_neg_hel_dx->Fill(dx);
    {
      h_neg_hel_dx->Fill(dx);
      helNegArray[binAt]++;
    }

  }//end loop over events

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

    h_simIN_dx->Fill(dx_simIN);//,weight_IN);

  }//end loop over events

  //std::cout << "sqrt(n) error " << h_simIN_dx->GetBinError(250) << std::endl;
  //std::cout << "Get Sumw2 Error: " << h_simIN_dx->GetSumw2() << std::endl;

  //******************************************************

  double scale = h_data_dx->Integral();
  //h_data_dx->Scale(1.0/h_data_dx->Integral());
  h_data_dx->Scale(1.0/scale);
  h_sim_proton_dx->Scale(1.0/h_sim_proton_dx->Integral());
  h_sim_neutron_dx->Scale(1.0/h_sim_neutron_dx->Integral());
  h_simIN_dx->Scale(1.0/h_simIN_dx->Integral());

  int nbins = h_data_dx->GetNbinsX();
  double xmin = h_data_dx->GetXaxis()->GetBinLowEdge(1);
  double xmax = h_data_dx->GetXaxis()->GetBinUpEdge(nbins);

  TF1 *FitFunc = new TF1("FitFunc",&fitsim,dx_min_i,dx_max_i,6); //-6,4,6

  FitFunc->SetNpx(numberBins);
  
  //----- GEN2 -----
  double startpar[] = {1.0,-0.5,0.3,-0.7,1.0,-1.0};
  FitFunc->SetParameters(startpar);
  FitFunc->SetParLimits(0,0.1,100);   // proton scale
  FitFunc->SetParLimits(1,-1.0,1.0);  // proton shift
  FitFunc->SetParLimits(2,0.0,100);   // neutron scale
  FitFunc->SetParLimits(3,-1.0,0.0);  // neutron shift
  FitFunc->SetParLimits(4,0.0,0.2);   // background scale
  FitFunc->SetParLimits(5,-1.0,1.0);  // background shift

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
  shifted_h_sim_proton_dx  ->Scale(scale);
  shifted_h_sim_neutron_dx ->Scale(scale);
  shifted_h_simIN_dx       ->Scale(scale);
  h_total_dx               ->Scale(scale);

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
  //h_asym->Sumw2();

  double A_array[numberBins];
  double A_err_array[numberBins];

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
    double c_bg      = shifted_h_simIN_dx->GetBinContent(bin);
    double c_bg_err  = shifted_h_simIN_dx->GetBinError(bin);

    double P_tot = c_p + c_n + c_bg;

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

      P_p_err  = TMath::Sqrt( (c_n+c_bg)*(c_n+c_bg)*c_p_err*c_p_err + c_p*c_p*c_n_err*c_n_err + c_p*c_p*c_bg_err*c_bg_err ) / ( (c_p+c_n+c_bg)*(c_p+c_n+c_bg) );
      P_n_err  = TMath::Sqrt( (c_bg+c_p)*(c_bg+c_p)*c_n_err*c_n_err + c_n*c_n*c_p_err*c_p_err + c_n*c_n*c_bg_err*c_bg_err ) / ( (c_p+c_n+c_bg)*(c_p+c_n+c_bg) );
      P_bg_err = TMath::Sqrt( (c_n+c_p)*(c_n+c_p)*c_bg_err*c_bg_err + c_bg*c_bg*c_p_err*c_p_err + c_bg*c_bg*c_p_err*c_p_err ) / ( (c_p+c_n+c_bg)*(c_p+c_n+c_bg) );

    }//end if total probability is nonzero

    //std::cout << "Proton probability: " << P_p << std::endl;
    h_prob_proton_dx->SetBinContent(bin,P_p);
    h_prob_proton_dx->SetBinError(bin,P_p_err);
    h_prob_neutron_dx->SetBinContent(bin,P_n);
    h_prob_neutron_dx->SetBinError(bin,P_n_err);
    h_prob_bckgrnd_dx->SetBinContent(bin,P_bg);
    h_prob_bckgrnd_dx->SetBinError(bin,P_bg_err);

    double c_pos      = h_pos_hel_dx->GetBinContent(bin);
    double c_pos_err  = h_pos_hel_dx->GetBinError(bin);
    double c_neg      = h_neg_hel_dx->GetBinContent(bin);
    double c_neg_err  = h_neg_hel_dx->GetBinError(bin);

    A_array[bin] = (helNegArray[bin] - helPosArray[bin])*1.0 / (helPosArray[bin] + helNegArray[bin]);
    A_err_array[bin] = std::sqrt(std::max(0.0,(4.0*helPosArray[bin]*helNegArray[bin])/std::pow((helPosArray[bin] + helNegArray[bin]),3)));
    //TMath::Sqrt((1 - A_array[bin]*A_array[bin])/(helPosArray[bin]+helNegArray[bin]));

    double Asym_tot = c_pos + c_neg;
    double Asym_diff = c_pos - c_neg;

    double Asym_raw = 0.0;
    double Asym_raw_err = 1.0;
    //if counts add up to zero then the stat error should be infinity, but in order to avoid the code exploding, i'll just make the error really large if we dont have any counts in that bin. this is just for safety. every bin should have counts.

    if (Asym_tot != 0.0)
    {
      Asym_raw = Asym_diff / Asym_tot;
      Asym_raw_err = 2 * TMath::Sqrt( c_neg*c_neg*c_pos_err*c_pos_err + c_pos*c_pos*c_neg_err*c_neg_err ) / ( (c_pos+c_neg)*(c_pos+c_neg) );
      //std::sqrt(std::max(0.0,(4.0*c_pos*c_neg)/std::pow(Asym_tot,3)));
      //^^method from JAck's code
    }//end if the total asym is nonzero

    h_asym->SetBinContent(bin, Asym_raw);
    //h_asym->SetBinContent(bin, A_array[bin]);
    h_asym->SetBinError(bin, Asym_raw_err);
    //h_asym->SetBinError(bin, A_err_array[bin]);

  }//end loop over bins

  h_resid           -> SetEntries(totalentries);
  h_prob_proton_dx  -> SetEntries(totalentries);
  h_prob_neutron_dx -> SetEntries(totalentries);
  h_prob_bckgrnd_dx -> SetEntries(totalentries);
  h_asym            -> SetEntries(totalentries);

  TF1 *AsymFitFunc = new TF1("AsymFitFunc",&fitAsym,dx_min_i,dx_max_i,3); //-6,3,3

  AsymFitFunc->SetNpx(numberBins);
  double Asymstartpar[] = {0.0,0.04,0.0};
  AsymFitFunc->SetParameters(Asymstartpar);
  AsymFitFunc->SetParLimits(0,-0.1,0.1); // proton asymmetry
  AsymFitFunc->SetParLimits(1,0.0,0.1);  // neutron asymmetry
  AsymFitFunc->SetParLimits(2,-0.1,0.1); // background asymmetry

  h_asym->Fit(AsymFitFunc,"0","",h_minX,h_maxX);

  TH1D* scaled_proton_prob  = (TH1D*)h_prob_proton_dx->Clone("scaled_proton_prob");
  TH1D* scaled_neutron_prob = (TH1D*)h_prob_neutron_dx->Clone("scaled_neutron_prob");
  TH1D* scaled_bckgrnd_prob = (TH1D*)h_prob_bckgrnd_dx->Clone("scaled_bckgrnd_prob");

  scaled_proton_prob->Scale(AsymFitFunc->GetParameter(0));
  scaled_neutron_prob->Scale(AsymFitFunc->GetParameter(1));
  scaled_bckgrnd_prob->Scale(AsymFitFunc->GetParameter(2));

  std::cout << "Proton Asymmetry: " << AsymFitFunc->GetParameter(0) << std::endl;
  std::cout << "Neutron Asymmetry: " << AsymFitFunc->GetParameter(1) << std::endl;
  std::cout << "Background Asymmetry: " << AsymFitFunc->GetParameter(2) << std::endl;

  //TH1F* h_Asym_dx = new TH1F("h_Asym_dx","",h_Nbins,h_minX,h_maxX);

  //for(int ibin = 0; ibin < h_Nbins; ibin++)
  //{
    //h_Asym_dx->SetBinContent(ibin, scaled_proton_prob->GetBinContent(ibin) + scaled_neutron_prob->GetBinContent(ibin) + scaled_bckgrnd_prob->GetBinContent(ibin));
  //}

  //h_Asym_dx->Scale(scale);

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

  auto legend = new TLegend(0.55,0.75,0.99,0.99);
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

  TCanvas *c3 = new TCanvas("c3","Asymmetry",100,100,1500,500);
  c3->Divide(3,1);

  gStyle->SetOptStat(0);

  c3->cd(1);
  h_pos_hel_dx->Draw("E");
  h_neg_hel_dx->Draw("E SAMES");

  auto legendHEL = new TLegend(0.7,0.8,0.99,0.99);
  legendHEL->SetTextSize(0.03);
  legendHEL->AddEntry(h_pos_hel_dx,"Positive Helicity","lep");
  legendHEL->AddEntry(h_neg_hel_dx,"Negative Helicity","lep");
  legendHEL->Draw();

  c3->cd(2);
  h_asym->Draw("E");
  //h_Asym_dx->Draw("SAMES");
  AsymFitFunc->Draw("SAMES");

  auto legendA = new TLegend(0.1,0.8,0.9,0.9);
  legendA->SetTextSize(0.02);
  legendA->AddEntry(h_asym,"Data","lep");
  //legendA->AddEntry(h_Asym_dx,"Scaled Fit","l");
  //legendA->AddEntry(AsymFitFunc,"Fit","l");
  legendA->AddEntry(AsymFitFunc,Form("A(dx) = %.4gPp(dx) + %.4gPn(dx) + %.4gPbg(dx)",AsymFitFunc->GetParameter(0),AsymFitFunc->GetParameter(1),AsymFitFunc->GetParameter(2)),"l");
  legendA->Draw();

  c3->cd(3);
  scaled_proton_prob->SetTitle("Scaled Probabilities");

  scaled_proton_prob->SetLineColor(kRed);
  scaled_proton_prob->Draw("E");

  scaled_neutron_prob->SetLineColor(kBlue);
  scaled_neutron_prob->Draw("E SAMES");

  scaled_bckgrnd_prob->SetLineColor(kGreen);
  scaled_bckgrnd_prob->Draw("E SAMES");

  auto legendProb = new TLegend(0.7,0.8,0.9,0.9);
  legendProb->SetTextSize(0.02);
  legendProb->AddEntry(scaled_proton_prob, "Proton", "l");
  legendProb->AddEntry(scaled_neutron_prob, "Neutron", "l");
  legendProb->AddEntry(scaled_bckgrnd_prob, "Background", "l");
  legendProb->Draw();

  TCanvas *c4 = new TCanvas("c4","Asymmetry",100,100,1000,1000);
  c4->cd();
  //h_asym->GetYaxis()->SetRangeUser(-0.4,0.6);
  h_asym->Draw("E");
  c4->SetGrid();
  c4->Update();

}//end main
