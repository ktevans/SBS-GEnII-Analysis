#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TLegend.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include <TStopwatch.h>
#include <list>

void ElectronArmCuts(const std::string cutType)
{
  TFile *inputFile = new TFile("../outfiles/GEN2_He3_run2321.root","READ");
  TTree *C = (TTree*)inputFile->Get("P");
  
  //PreShower Histograms
  TH1D* h_ps_e_1 = new TH1D("h_ps_e_1", "PreShower Energy", 100, 0.0, 2.0);
  h_ps_e_1 -> GetXaxis() -> SetTitle("PreShower Energy [GeV]");
  TH1D* h_ps_e_2 = new TH1D("h_ps_e_2", "PreShower Energy", 100, 0.0, 2.0);
  h_ps_e_2 -> GetXaxis() -> SetTitle("PreShower Energy [GeV]");
  TH1D* h_ps_e_3 = new TH1D("h_ps_e_3", "PreShower Energy", 100, 0.0, 2.0);
  h_ps_e_3 -> GetXaxis() -> SetTitle("PreShower Energy [GeV]");
  TH1D* h_ps_e_4 = new TH1D("h_ps_e_4", "PreShower Energy", 100, 0.0, 2.0);
  h_ps_e_4 -> GetXaxis() -> SetTitle("PreShower Energy [GeV]");
  TH1D* h_ps_e_5 = new TH1D("h_ps_e_5", "PreShower Energy", 100, 0.0, 2.0);
  h_ps_e_5 -> GetXaxis() -> SetTitle("PreShower Energy [GeV]");
  h_ps_e_1->SetLineColor(2);
  h_ps_e_2->SetLineColor(3);
  h_ps_e_3->SetLineColor(4);
  h_ps_e_4->SetLineColor(7);
  h_ps_e_5->SetLineColor(6);

  //Shower Histograms
  TH1D* h_sh_e_1 = new TH1D("h_sh_e_1", "Shower Energy", 250, 0.0, 5.0);
  h_sh_e_1 -> GetXaxis() -> SetTitle("Shower Energy [GeV]");
  TH1D* h_sh_e_2 = new TH1D("h_sh_e_2", "Shower Energy", 250, 0.0, 5.0);
  h_sh_e_2 -> GetXaxis() -> SetTitle("Shower Energy [GeV]");
  TH1D* h_sh_e_3 = new TH1D("h_sh_e_3", "Shower Energy", 250, 0.0, 5.0);
  h_sh_e_3 -> GetXaxis() -> SetTitle("Shower Energy [GeV]");
  TH1D* h_sh_e_4 = new TH1D("h_sh_e_4", "Shower Energy", 250, 0.0, 5.0);
  h_sh_e_4 -> GetXaxis() -> SetTitle("Shower Energy [GeV]");
  TH1D* h_sh_e_5 = new TH1D("h_sh_e_5", "Shower Energy", 250, 0.0, 5.0);
  h_sh_e_5 -> GetXaxis() -> SetTitle("Shower Energy [GeV]");
  h_sh_e_1->SetLineColor(2);
  h_sh_e_2->SetLineColor(3);
  h_sh_e_3->SetLineColor(4);
  h_sh_e_4->SetLineColor(7);
  h_sh_e_5->SetLineColor(6);

  //GRINCH Time over Threshold Histograms
  TH1D* h_gr_tot_1 = new TH1D("h_gr_tot_1", "GRINCH Cluster Mean Time over Threshold", 100, 0.0, 100.0);
  h_gr_tot_1 -> GetXaxis() -> SetTitle("GRINCH Cluster Mean Time over Threshold [ns]");
  TH1D* h_gr_tot_2 = new TH1D("h_gr_tot_2", "GRINCH Cluster Mean Time over Threshold", 100, 0.0, 100.0);
  h_gr_tot_2 -> GetXaxis() -> SetTitle("GRINCH Cluster Mean Time over Threshold [ns]");
  TH1D* h_gr_tot_3 = new TH1D("h_gr_tot_3", "GRINCH Cluster Mean Time over Threshold", 100, 0.0, 100.0);
  h_gr_tot_3 -> GetXaxis() -> SetTitle("GRINCH Cluster Mean Time over Threshold [ns]");
  TH1D* h_gr_tot_4 = new TH1D("h_gr_tot_4", "GRINCH Cluster Mean Time over Threshold", 100, 0.0, 100.0);
  h_gr_tot_4 -> GetXaxis() -> SetTitle("GRINCH Cluster Mean Time over Threshold [ns]");
  TH1D* h_gr_tot_5 = new TH1D("h_gr_tot_5", "GRINCH Cluster Mean Time over Threshold", 100, 0.0, 100.0);
  h_gr_tot_5 -> GetXaxis() -> SetTitle("GRINCH Cluster Mean Time over Threshold [ns]");
  h_gr_tot_1->SetLineColor(2);
  h_gr_tot_2->SetLineColor(3);
  h_gr_tot_3->SetLineColor(4);
  h_gr_tot_4->SetLineColor(7);
  h_gr_tot_5->SetLineColor(6);

  //Total Energy over Momentum Histograms
  TH1D* h_EoP_1 = new TH1D("h_EoP_1", "E over P", 150, 0.0, 3.0);
  h_EoP_1 -> GetXaxis() -> SetTitle("PS + SH Energy over Momentum");
  TH1D* h_EoP_2 = new TH1D("h_EoP_2", "E over P", 150, 0.0, 3.0);
  h_EoP_2 -> GetXaxis() -> SetTitle("PS + SH Energy over Momentum");
  TH1D* h_EoP_3 = new TH1D("h_EoP_3", "E over P", 150, 0.0, 3.0);
  h_EoP_3 -> GetXaxis() -> SetTitle("PS + SH Energy over Momentum");
  TH1D* h_EoP_4 = new TH1D("h_EoP_4", "E over P", 150, 0.0, 3.0);
  h_EoP_4 -> GetXaxis() -> SetTitle("PS + SH Energy over Momentum");
  TH1D* h_EoP_5 = new TH1D("h_EoP_5", "E over P", 150, 0.0, 3.0);
  h_EoP_5 -> GetXaxis() -> SetTitle("PS + SH Energy over Momentum");
  h_EoP_1->SetLineColor(2);
  h_EoP_2->SetLineColor(3);
  h_EoP_3->SetLineColor(4);
  h_EoP_4->SetLineColor(7);
  h_EoP_5->SetLineColor(6);

  //Define your cuts based on user input
  TCut cut1 = "";
  TCut cut2 = "";
  TCut cut3 = "";
  TCut cut4 = "";
  TCut cut5 = "";
  TString cutLabel1 = "";
  TString cutLabel2 = "";
  TString cutLabel3 = "";
  TString cutLabel4 = "";
  TString cutLabel5 = "";
  TString var = "";
  
  if(cutType=="ps")
  {
    //PreShower Energy Cuts
    cut1 = "bb_ps_e>0.0";
    cut2 = "bb_ps_e>0.1";
    cut3 = "bb_ps_e>0.2";
    cut4 = "bb_ps_e>0.5";
    cut5 = "bb_ps_e>1.0";
    cutLabel1 = "bb_ps_e>0.0";
    cutLabel2 = "bb_ps_e>0.1";
    cutLabel3 = "bb_ps_e>0.2";
    cutLabel4 = "bb_ps_e>0.5";
    cutLabel5 = "bb_ps_e>1.0";
    var = "PreShower";
  }
  if(cutType=="sh")
  {
    //Shower Energy Cuts
    cut1 = "bb_sh_e>0.0&&bb_sh_e<4.0";
    cut2 = "bb_sh_e>0.5&&bb_sh_e<3.0";
    cut3 = "bb_sh_e>0.6&&bb_sh_e<2.8";
    cut4 = "bb_sh_e>0.8&&bb_sh_e<2.5";
    cut5 = "bb_sh_e>1.0&&bb_sh_e<2.0";
    cutLabel1 = "bb_sh_e>0.0&&bb_sh_e<4.0";
    cutLabel2 = "bb_sh_e>0.5&&bb_sh_e<3.0";
    cutLabel3 = "bb_sh_e>0.6&&bb_sh_e<2.8";
    cutLabel4 = "bb_sh_e>0.8&&bb_sh_e<2.5";
    cutLabel5 = "bb_sh_e>1.0&&bb_sh_e<2.0";
    var = "Shower";
  }
  if(cutType=="gr")
  {
    //GRINCH Cluster Size Cuts
    cut1 = "bb_gr_size>0.0";
    cut2 = "bb_gr_size>1.0";
    cut3 = "bb_gr_size>2.0";
    cut4 = "bb_gr_size>3.0";
    cut5 = "bb_gr_size>4.0";
    cutLabel1 = "bb_gr_size>0.0";
    cutLabel2 = "bb_gr_size>1.0";
    cutLabel3 = "bb_gr_size>2.0";
    cutLabel4 = "bb_gr_size>3.0";
    cutLabel5 = "bb_gr_size>4.0";
    var = "GRINCH_ClusSize";
  }
  if(cutType=="EoP")
  {
    //E/p Cuts
    cut1 = "bb_etot_over_p>0.0&&bb_etot_over_p<2.5";
    cut2 = "bb_etot_over_p>0.5&&bb_etot_over_p<1.5";
    cut3 = "bb_etot_over_p>0.7&&bb_etot_over_p<1.3";
    cut4 = "bb_etot_over_p>0.8&&bb_etot_over_p<1.2";
    cut5 = "bb_etot_over_p>0.9&&bb_etot_over_p<1.1";
    cutLabel1 = "bb_etot_over_p>0.0&&bb_etot_over_p<2.5";
    cutLabel2 = "bb_etot_over_p>0.5&&bb_etot_over_p<1.5";
    cutLabel3 = "bb_etot_over_p>0.7&&bb_etot_over_p<1.3";
    cutLabel4 = "bb_etot_over_p>0.8&&bb_etot_over_p<1.2";
    cutLabel5 = "bb_etot_over_p>0.9&&bb_etot_over_p<1.1";
    var = "EoverP";
  }

  //Draw all of the histograms with the 5 different cuts
  C->Draw(Form("%s>>h_ps_e_1","bb_ps_e"),cut1);
  C->Draw(Form("%s>>h_ps_e_2","bb_ps_e"),cut2,"SAMES");
  C->Draw(Form("%s>>h_ps_e_3","bb_ps_e"),cut3,"SAMES");
  C->Draw(Form("%s>>h_ps_e_4","bb_ps_e"),cut4,"SAMES");
  C->Draw(Form("%s>>h_ps_e_5","bb_ps_e"),cut5,"SAMES");

  C->Draw(Form("%s>>h_sh_e_1","bb_sh_e"),cut1);
  C->Draw(Form("%s>>h_sh_e_2","bb_sh_e"),cut2,"SAMES");
  C->Draw(Form("%s>>h_sh_e_3","bb_sh_e"),cut3,"SAMES");
  C->Draw(Form("%s>>h_sh_e_4","bb_sh_e"),cut4,"SAMES");
  C->Draw(Form("%s>>h_sh_e_5","bb_sh_e"),cut5,"SAMES");

  C->Draw(Form("%s>>h_gr_tot_1","bb_gr_tot"),cut1);
  C->Draw(Form("%s>>h_gr_tot_2","bb_gr_tot"),cut2,"SAMES");
  C->Draw(Form("%s>>h_gr_tot_3","bb_gr_tot"),cut3,"SAMES");
  C->Draw(Form("%s>>h_gr_tot_4","bb_gr_tot"),cut4,"SAMES");
  C->Draw(Form("%s>>h_gr_tot_5","bb_gr_tot"),cut5,"SAMES");

  C->Draw(Form("%s>>h_EoP_1","bb_etot_over_p"),cut1);
  C->Draw(Form("%s>>h_EoP_2","bb_etot_over_p"),cut2,"SAMES");
  C->Draw(Form("%s>>h_EoP_3","bb_etot_over_p"),cut3,"SAMES");
  C->Draw(Form("%s>>h_EoP_4","bb_etot_over_p"),cut4,"SAMES");
  C->Draw(Form("%s>>h_EoP_5","bb_etot_over_p"),cut5,"SAMES");

  TCanvas *c1 = new TCanvas("c1","ElectronArm",1000,1000,1000,1000);

  auto cutLegend = new TPaveLabel(0.0,0.0,1.0,0.1,"test","brNDC");
  cutLegend->Draw();
  gPad->Update();

  TPad *pad1 = new TPad("pad1","title",0.0,0.85,1.0,1.0);
  TPad *pad2 = new TPad("pad2","plots",0.0,0.0,1.0,0.82);
  pad1->Draw();
  pad2->Draw();

  pad2->cd();
  pad2->Divide(2,2);

  gStyle->SetOptStat(0);

  pad2->cd(1);
  h_ps_e_1->Draw();
  h_ps_e_2->Draw("SAMES");
  h_ps_e_3->Draw("SAMES");
  h_ps_e_4->Draw("SAMES");
  h_ps_e_5->Draw("SAMES");
  //gPad->BuildLegend();

  pad2->cd(2);
  h_sh_e_1->Draw();
  h_sh_e_2->Draw("SAMES");
  h_sh_e_3->Draw("SAMES");
  h_sh_e_4->Draw("SAMES");
  h_sh_e_5->Draw("SAMES");
  //gPad->BuildLegend();

  pad2->cd(3);
  h_gr_tot_1->Draw();
  h_gr_tot_2->Draw("SAMES");
  h_gr_tot_3->Draw("SAMES");
  h_gr_tot_4->Draw("SAMES");
  h_gr_tot_5->Draw("SAMES");
  //gPad->BuildLegend();

  pad2->cd(4);
  h_EoP_1->Draw();
  h_EoP_2->Draw("SAMES");
  h_EoP_3->Draw("SAMES");
  h_EoP_4->Draw("SAMES");
  h_EoP_5->Draw("SAMES");
  //gPad->BuildLegend();

  pad1->cd();

  TText *t1 = new TText(0.5,0.9,cutLabel1);
  t1->SetTextColor(2);
  t1->SetTextSize(0.2);
  t1->SetTextAlign(22);
  t1->Draw();

  TText *t2 = new TText(0.5,0.7,cutLabel2);
  t2->SetTextColor(3);
  t2->SetTextSize(0.2);
  t2->SetTextAlign(22);
  t2->Draw();

  TText *t3 = new TText(0.5,0.5,cutLabel3);
  t3->SetTextColor(4);
  t3->SetTextSize(0.2);
  t3->SetTextAlign(22);
  t3->Draw();

  TText *t4 = new TText(0.5,0.3,cutLabel4);
  t4->SetTextColor(7);
  t4->SetTextSize(0.2);
  t4->SetTextAlign(22);
  t4->Draw();

  TText *t5 = new TText(0.5,0.1,cutLabel5);
  t5->SetTextColor(6);
  t5->SetTextSize(0.2);
  t5->SetTextAlign(22);
  t5->Draw();
  
  
  c1->Update();
  c1->Print("../outfiles/ElectronArm_"+var+".pdf");
  
}  
