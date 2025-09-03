#include <TF1.h>

TH1D *h_data_dx;
TH1D *h_simIN_dx;
TH1D *h_sim_proton_dx;
TH1D *h_sim_neutron_dx;

double fitsim( double *x, double *par)
  {
    double dx = x[0];

    double Norm_overall = par[0];
    double R_pn = par[1];
    double Bg_norm = par[2];
    double bg = h_simIN_dx->Interpolate(dx);
  
    double simu = Norm_overall * (h_sim_proton_dx->Interpolate(dx) + R_pn * h_sim_neutron_dx->Interpolate(dx) + Bg_norm * bg);
    return simu;   
  }//end fitsim

void SimDataComp()
{

  TChain* T_data = new TChain("T_data");
  T_data->Add("outfiles/parsed_GEn_pass2_GEN2_He3_dxdy.root");

  Double_t dx;       T_data->SetBranchAddress("dx", &dx);

  if(T_data->GetEntries()==0)
  {
    std::cerr << "\n --- No ROOT file found!! --- \n\n";
    throw;
  }
  else std::cout << "\nFound " << T_data->GetEntries() << " events. \n";

  h_data_dx = new TH1D("h_data_dx", ";dx", 500, -6.0, 4.0);
  h_data_dx->GetXaxis()->SetTitle("dx [m]");

  for (size_t iev = 0; iev < T_data->GetEntries(); iev++)
  {
    T_data->GetEntry(iev);

    h_data_dx->Fill(dx);
  }//end loop over events

  TChain* T_sim = new TChain("T_sim");
  T_sim->Add("outfiles/parsed_SIM_GEn_GEN2_He3_dxdy.root");

  Double_t dx_sim;     T_sim->SetBranchAddress("dx", &dx_sim);
  Double_t fnucl;      T_sim->SetBranchAddress("fnucl", &fnucl);
  Double_t weight;     T_sim->SetBranchAddress("weight", &weight);

  if(T_sim->GetEntries()==0)
  {
    std::cerr << "\n --- No ROOT file found!! --- \n\n";
    throw;
  }
  else std::cout << "\nFound " << T_sim->GetEntries() << " events. \n";

  h_sim_proton_dx = new TH1D("h_sim_proton_dx", ";dx_sim_p", 500, -6.0, 4.0);
  h_sim_proton_dx->GetXaxis()->SetTitle("dx [m]");
  h_sim_proton_dx->SetLineColor(kGreen);

  h_sim_neutron_dx = new TH1D("h_sim_neutron_dx", ";dx_sim_n", 500, -6.0, 4.0);
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
  T_simIN->Add("outfiles/parsed_SIM_IN_GEn_GEN2_He3_dxdy.root");

  Double_t dx_simIN;      T_simIN->SetBranchAddress("dx", &dx_simIN);
  Double_t fnucl_IN;      T_simIN->SetBranchAddress("fnucl", &fnucl_IN);
  Double_t weight_IN;     T_simIN->SetBranchAddress("weight", &weight_IN);

  if(T_simIN->GetEntries()==0)
  {
    std::cerr << "\n --- No ROOT file found!! --- \n\n";
    throw;
  }
  else std::cout << "\nFound " << T_simIN->GetEntries() << " events. \n";

  h_simIN_dx = new TH1D("h_simIN_dx", ";dxIN", 500, -6.0, 4.0);
  h_simIN_dx->GetXaxis()->SetTitle("dx [m]");
  h_simIN_dx->SetLineColor(kBlue);

  for (size_t iev = 0; iev < T_simIN->GetEntries(); iev++)
  {
    T_simIN->GetEntry(iev);

    h_simIN_dx->Fill(dx_simIN,weight_IN);
  }//end loop over events

  //******************************************************

  double scale = h_data_dx->Integral();
  h_data_dx->Scale(1.0/h_data_dx->Integral());
  h_sim_proton_dx->Scale(1.0/h_sim_proton_dx->Integral());
  h_sim_neutron_dx->Scale(1.0/h_sim_neutron_dx->Integral());
  h_simIN_dx->Scale(1.0/h_simIN_dx->Integral());

  int nbins = h_data_dx->GetNbinsX();
  double xmin = h_data_dx->GetXaxis()->GetBinLowEdge(1);
  double xmax = h_data_dx->GetXaxis()->GetBinUpEdge(nbins);

  int npar = 3;

  TF1 *FitFunc = new TF1("FitFunc",&fitsim,0,10,3);

  FitFunc->SetNpx(1000);
  double startpar[] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  FitFunc->SetParameters(startpar);
  FitFunc->SetParLimits(0,0.1,100);
  FitFunc->SetParLimits(1,0.1,100);
  FitFunc->SetParLimits(2,0,100);

  h_data_dx->Fit(FitFunc,"0","",xmin,xmax);

  h_sim_proton_dx->Scale(FitFunc->GetParameter(0));
  h_sim_neutron_dx->Scale(FitFunc->GetParameter(0)*FitFunc->GetParameter(1));
  h_simIN_dx->Scale(FitFunc->GetParameter(0)*FitFunc->GetParameter(2));

  TH1F* h_total_dx = new TH1F("h_total_dx","",nbins,xmin,xmax);
  h_total_dx->SetMarkerStyle('.');
  for(int ibin = 0; ibin < nbins; ibin++)
  {
    h_total_dx->SetBinContent(ibin, h_sim_neutron_dx->GetBinContent(ibin) + h_sim_proton_dx->GetBinContent(ibin) + h_simIN_dx->GetBinContent(ibin));
  }

  h_data_dx->Scale(scale);
  h_sim_proton_dx->Scale(scale);
  h_sim_neutron_dx->Scale(scale);
  h_simIN_dx->Scale(scale);
  h_total_dx->Scale(scale);

  TCanvas *c1 = new TCanvas("c1", "dx", 100,100,1600,800);
  c1->cd();
  h_data_dx->Draw();
  h_simIN_dx->Draw("SAMES");
  h_sim_proton_dx->Draw("SAMES");
  h_sim_neutron_dx->Draw("SAMES");
  h_total_dx->Draw("SAMES");
}
