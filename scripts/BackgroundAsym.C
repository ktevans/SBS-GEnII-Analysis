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
//  for GEN2:  "BackgroundAsym(2)"
//  for GEN3:  "BackgroundAsym(3)"
//  for GEN4a: "BackgroundAsym(4)"
//  for GEN4b: "BackgroundAsym(5)"

void BackgroundAsym(int kin)
{

  TString data_file;
  TString title_words;

  int npar = 6;
  double dx_min_d, dx_max_d;
  int dx_min_i, dx_max_i;

  if(kin == 2)
  {

    data_file = "/volatile/halla/sbs/ktevans/pass3/QE_data_GEN2_sbs100p_nucleon_np_model2.root";
    title_words = "GEN2";
    dx_min_d = -2.8;
    dx_min_i = -3;
    dx_max_d = 2.0;
    dx_max_i = 2;
    cout<<"\nHi! You are analyzing GEN2!\n";

  }

  if(kin == 3)
  {

    data_file = "/volatile/halla/sbs/ktevans/pass3/QE_data_GEN3_sbs100p_nucleon_np_model2.root";
    title_words = "GEN3";
    dx_min_d = -2.5;
    dx_min_i = -3;
    dx_max_d = 1.7;
    dx_max_i = 2;
    cout<<"\nHi! You are analyzing GEN3!\n";

  }

  if(kin == 4)
  {

    data_file = "/volatile/halla/sbs/ktevans/pass3/QE_data_GEN4a_sbs100p_nucleon_np_model2.root";
    title_words = "GEN4a";
    dx_min_d = -3.0;
    dx_min_i = -3;
    dx_max_d = 2.0;
    dx_max_i = 2;
    cout<<"\nHi! You are analyzing GEN4a!\n";

  }

  if(kin == 5)
  {

    data_file = "/volatile/halla/sbs/ktevans/pass3/QE_data_GEN4b_sbs100p_nucleon_np_model2.root";
    title_words = "GEN4b";
    dx_min_d = -3.0;
    dx_min_i = -3;
    dx_max_d = 2.0;
    dx_max_i = 2;
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

  int numberBins = 96;

  TChain* Tout = new TChain("Tout");
  Tout->Add(data_file);

  Double_t dx;            Tout->SetBranchAddress("dx", &dx);
  Double_t dy;            Tout->SetBranchAddress("dy", &dy);
  int helicity;           Tout->SetBranchAddress("helicity", &helicity);
  double coin;            Tout->SetBranchAddress("adc.coin", &coin);
  double ps_e;            Tout->SetBranchAddress("bb.ps.e", &ps_e);
  double sh_e;            Tout->SetBranchAddress("bb.sh.e", &sh_e);
  double W2;              Tout->SetBranchAddress("e.kine.W2", &W2);
  double Q2;              Tout->SetBranchAddress("e.kine.Q2", &Q2);
  double tr_p;            Tout->SetBranchAddress("bb.tr.p", &tr_p);
  double grinch_track;    Tout->SetBranchAddress("bb.grinch_tdc.clus.trackindex", &grinch_track);
  double grinch_clusSize; Tout->SetBranchAddress("bb.grinch_tdc.clus.size", &grinch_clusSize);

  if(Tout->GetEntries()==0)
  {
    std::cerr << "\n --- No ROOT file found!! --- \n\n";
    throw;
  }
  else std::cout << "\nFound " << Tout->GetEntries() << " events. \n";

  TH1D* h_neg_hel_dx_inel = new TH1D("h_neg_hel_dx_inel",";-hel", numberBins, dx_min_d, dx_max_d);
  h_neg_hel_dx_inel->GetXaxis()->SetTitle("dx [m]");
  h_neg_hel_dx_inel->Sumw2();

  TH1D* h_pos_hel_dx_inel = new TH1D("h_pos_hel_dx_inel",";+hel", numberBins, dx_min_d, dx_max_d);
  h_pos_hel_dx_inel->GetXaxis()->SetTitle("dx [m]");
  h_pos_hel_dx_inel->Sumw2();

  for (size_t iev = 0; iev < Tout->GetEntries(); iev++)
  {
    Tout->GetEntry(iev);

    if(W2>2.0 && abs(coin+0.47385)<3.6 && grinch_track==0.0 && grinch_clusSize>=3.0 && abs(((ps_e+sh_e)/tr_p)-0.97)<0.2 && abs(dy)<0.99)
    {
      if(helicity==-1)
      {
        h_neg_hel_dx_inel->Fill(dx);
      }

      if(helicity==1)
      {
        h_pos_hel_dx_inel->Fill(dx);
      }
    }

  }//end loop over events

  TH1D* hAsymDiff_inel = (TH1D*) h_pos_hel_dx_inel->Clone("hAsymDiff_inel");
  hAsymDiff_inel->Sumw2();
  hAsymDiff_inel->Add(h_neg_hel_dx_inel, -1.0);

  TH1D* hAsymSum_inel = (TH1D*) h_pos_hel_dx_inel->Clone("hAsymSum_inel");
  hAsymSum_inel->Sumw2();
  hAsymSum_inel->Add(h_neg_hel_dx_inel);

  TH1D* hAsym_inel = (TH1D*) hAsymDiff_inel->Clone("hAsym_inel");
  hAsym_inel->Divide(hAsymSum_inel);
  hAsym_inel->Sumw2();

  TCanvas *c1 = new TCanvas("c1","Inelastic Asymmetry",100,100,1500,500);
  c1->cd();
  hAsym_inel->Draw();

}//end main
