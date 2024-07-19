#include <TSystem.h>
#include <TChain.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
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
#include<math.h>
#include <stack>
#include "TLorentzVector.h"
#include "TCut.h"
#include "TLatex.h"
#include "TLine.h"
using namespace std;


//global variables
const int maxTracks = 16; 
const double PI = TMath::Pi();
const double M_e = 0.0051;
const double M_p = 0.938272;
const double M_n = 0.939565;

//const double hcalheight = -0.2897; // Need for Pass0 or Pass1 data
//const double hcalheight = 0; // should put in config file maybe


////  Detectors
  // HCAL 
  // Dimensions
  static const Int_t hcalchan = 288;
  static const Int_t hcalcol = 12;
  static const Int_t hcalrow = 24;
  static const Double_t hcal_hrange = 1.85928;    //m, total range in horizontal direction of HCal (end-to-end)
  static const Double_t hcal_vrange = 3.81;       //m, total range in vertical direction of HCal (end-to-end)
  static const Double_t hcaladc_binw = 4.;        //ns, width of each ADC bin
  static const Double_t hcalblk_w_p0 = 0.15;      //m, width of a HCAL block, pass0/1
  static const Double_t hcalblk_h_p0 = 0.15;      //m, height of a HCAL block, pass0/1 
  static const Double_t hcalblk_w = 0.1524;       //m, width of a HCAL block
  static const Double_t hcalblk_h = 0.1524;       //m, height of a HCAL block
  static const Double_t hcalblk_div_h = 0.15494;  //m, horizontal center-to-center dist.
  static const Double_t hcalblk_div_v = 0.15875;  //m, vertical center-to-center dist.
  static const Double_t hcalblk_div_hyp = 0.15875;//m, division corner-to-cornter dist.
  static const Double_t hcalblk_gap_h = 0.00254;  //m, horiz. gap bet. two blocks
  static const Double_t hcalblk_gap_v = 0.00635;  //m, vert. gap bet. two blocks
  // Positions (mc)
  static const Double_t hcalposXi_mc = -2.655;    //m, distance from beam center to top of HCal w/75cm offset
  static const Double_t hcalposXf_mc = 1.155;     //m, distance from beam center to bottom of HCal w/75cm offset
  static const Double_t hcalposYi_mc = -0.92964;  //m, distance from beam center to opposite-beam side of HCal
  static const Double_t hcalposYf_mc = 0.92964;   //m, distance from beam center to beam side of HCal
  // Pass0/1 (no block spacing)
  static const Double_t hcalposXi_p0 = -2.16014;  //m, distance from beam center to top of HCal w/75cm offset
  static const Double_t hcalposXf_p0 = 1.43826;   //m, distance from beam center to bottom of HCal w/75cm offset
  static const Double_t hcalposYi_p0 = -0.9;      //m, distance from beam center to opposite-beam side of HCal
  static const Double_t hcalposYf_p0 = 0.9;  
//m, distance from beam center to beam side of HCal
  // Positions (data fits)
  static const Double_t hcalposXi = -2.268095;    //m, distance from beam center to top of HCal (obsolete)
  static const Double_t hcalposXf = 1.538095;     //m, distance from beam center to bottom of HCal (obsolete)
  static const Double_t hcalposYi = -0.931545;    //m, distance from beam center to opposite-beam side of HCal (obsolete)
  static const Double_t hcalposYf = 0.931545;     //m, distance from beam center to beam side of HCal (obsolete)
  // Global
  //static const Double_t hcalvoff = -0.2897;       //m, height of the center of hcal above beam (m)
  //static const Double_t hcalvoff = -0.3735;       //m, height of the center of hcal above beam (m) (sbs8)
  //static const Double_t hcalvoff = -0.75;       //m, height of the center of hcal above beam (m)
  //static const Double_t hcalvoff = 0.0;         //m, height of the center of hcal above beam (m) after pass2 corr

 // target
  static const Double_t l_tgt = 0.15;  //m 
  static const Double_t celldiameter = 1.6*2.54; //cm, 
// shielding
  static const Double_t rho_al = 2.7; //g/cc
  static const Double_t aldEdx = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV

  // LH2
  static const Double_t lh2tarrho = 0.0723;     //g/cc, target density
  static const Double_t lh2cthick = 0.02;       //cm, target cell thickness
  static const Double_t lh2uwallthick = 0.0145; //cm, upstream wall thickness
  static const Double_t lh2dwallthick = 0.015;  //cm, downstream wall thickness
  static const Double_t lh2dEdx = 0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy. On the other hand, target group -> 0.00480
  // LD2
  static const Double_t ld2tarrho = 0.169;      //g/cc, target density
  static const Double_t ld2dEdx = 0.00581;      //According to https://open.library.ubc.ca/media/stream/pdf/831/1.0085416/1, pick up a factor of 1.012. On the other hand, target group -> 0.00240
  static const Double_t ld2uwallthick = 0.0145; //cm, assume same as hydrogen for now
  static const Double_t ld2dwallthick = 0.015;  //cm, assume same as hydrogen for now



const Double_t atime_diff_mean = 0;
const Double_t atime_diff_sigma = 3;
const Double_t Nsigma_atime_diff =5;

const Double_t HCAL_atime_mean = 0;
const Double_t HCAL_atime_sigma = 5;
const Double_t Nsigma_HCAL_atime =5;

//function declarations
std::vector<Double_t> hcalaa_data (int exblkN_x=1, int exblkN_y=1);
std::vector<Double_t> hcalaa_mc (int exblkN_x=1, int exblkN_y=1);
bool hcalaaON (Double_t hcalx, Double_t hcaly, std::vector<Double_t> hcalaa);
std::vector<Double_t> hcalfid (Double_t dxsig_p, Double_t dxsig_n, Double_t dysig, std::vector<Double_t> hcalaa);
bool hcalfidIN (Double_t hcalx_exp, Double_t hcaly_exp, Double_t dx_pn, vector<Double_t> fid);
Double_t calculate_nsigma_fid_y(Double_t hcaly_exp,Double_t dysig, std::vector<Double_t> hcalaa);
Double_t calculate_nsigma_fid_x(Double_t hcalx_exp,Double_t dxsig_n, Double_t dxsig_p, Double_t dx_pn, std::vector<Double_t> hcalaa);


void dxdy_data_v2(int entries_input = -1, TString configfileinput="sbs4_30p" ){//main
  
  gStyle->SetNumberContours(255); 
  
  double M_N = M_p;// mass of the nucleon. We'll determine if it's deuterium or hydrogen based on the config file input;

  TChain *C = new TChain("T");
 
  //SBS 4
  double E_e = 3.728;
  double BB_d = 1.7988;
  double BB_th = 36.0* TMath::DegToRad();
  double HCal_d =11.0;
  double HCal_th = 31.9 * TMath::DegToRad();
  double W2_mean = 1.00;
  double W2_sig = 0.24;
  //const double hcalheight = -0.2897; // Need for Pass0 or Pass1 data
  double hcalheight = 0; // 
  double dxsig_p = 0.167; // sbs4 30% data 
  double dxsig_n = 0.18; // sbs4 30% data 
  double dysig = 0.27; // sbs4 30% data
  double dx_pn = 1.18; // 0.639 for sbs4 30% data
  double nsigma_p = 1;
  double nsigma_n =1;
  double nsigma_y =1;
  TString target = "deuterium";
  int is_deuterium = 0;
  int is_hydrogen = 0; 
  double mag_field =30;

  //cout<<"kine: "<<kine<<endl;
  //cout<<"config file"<<configfileinput<<endl;

  //// set location to config file
   //string configfilename = Form("../config/sbs%d.cfg",kine);
   TString configfilename = "../config/" + configfileinput + ".cfg";
  // string configfilename = "/w/halla-scshelf2102/sbs/msatnik/GMn/config/sbs4_30p.cfg";
   cout<<"reading from config file: "<<configfilename<<endl;

   // set location and name for output file 
   TString outputfilename = "../output/" + configfileinput + ".root";
   //TString outputfilename = "../output/test.root";

   // Declare outfile
   // TFile *fout = new TFile( Form("../output/sbs%d.root",kine), "RECREATE" );
   //TFile *fout = new TFile("../output/sbs4_30p.root", "RECREATE" );
  TFile *fout = new TFile(outputfilename,"RECREATE");
  cout<<"writing to file: "<< outputfilename <<endl;

    //cout<<"testing reading in from a config file"<<endl;

  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      if(!currentline) cout << "WARNING: No file exists at " << currentline << "." << endl;
      C->Add(currentline);
      cout << "Loaded file at: " << currentline << endl;
    }    
  }
  TCut globalcut = "";
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
    }    
    cout<< "Global Cut: "<<globalcut<<endl;
  }
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ){
    TObjArray *tokens = currentline.Tokenize(" ");
    Int_t ntokens = tokens->GetEntries();
    if( ntokens>1 ){
      TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
      if( skey == "E_e" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	E_e = sval.Atof();
	cout << "Loading beam energy: " << E_e << endl;
      }
      if( skey == "HCal_d" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_d = sval.Atof();
	cout << "Loading HCal distance: " << HCal_d << endl;
      }
      if( skey == "HCal_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_th = sval.Atof();	
	cout << "Loading HCal angle: " << HCal_th << endl;
	HCal_th = HCal_th * TMath::DegToRad();
      }
      if( skey == "hcalheight" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	hcalheight = sval.Atof();	
	cout << "Loading HCal 'height': " << hcalheight<< endl;
      }
      if( skey == "BB_d" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	BB_d = sval.Atof();
	cout << "Loading BB distance: " << BB_d << endl;
      }
      if( skey == "BB_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	BB_th = sval.Atof();	
	cout << "Loading BB angle: " << BB_th << endl;
	BB_th = BB_th* TMath::DegToRad();
      }
      if( skey == "W2_mean" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W2_mean = sval.Atof();
	cout << "Loading W2 mean cut: " << W2_mean << endl;
      }
      if( skey == "W2_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W2_sig = sval.Atof();
	cout << "Loading W2 sigma cut: " << W2_sig << endl;
      }
      if( skey == "dxsig_p" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        dxsig_p = sval.Atof();
	cout << "Loading dxsig_p: " << dxsig_p << endl;
      }
      if( skey == "dxsig_n" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        dxsig_n = sval.Atof();
	cout << "Loading dxsig_n: " << dxsig_n << endl;
      }
      if( skey == "dysig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        dysig = sval.Atof();
	cout << "Loading dysig: " << dysig << endl;
      }
      if( skey == "dx_pn" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        dx_pn = sval.Atof();
	cout << "Loading dx_pn: " << dx_pn << endl;
      }
      if( skey == "nsigma_p" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        nsigma_p = sval.Atof();
	cout << "Loading nsigma_p: " << nsigma_p << endl;
      }
      if( skey == "nsigma_n" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        nsigma_n = sval.Atof();
	cout << "Loading nsigma_n: " << nsigma_n<< endl;
      }
      if( skey == "nsigma_y" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        nsigma_y = sval.Atof();
	cout << "Loading nsigma_y: " << nsigma_y<< endl;
      }
      if( skey == "target" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        target = sval;
	cout << "Loading target type: " << target << endl;
      }
    if( skey == "mag_field" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        mag_field = sval.Atof();
	cout << "Loading magnetic field: " << mag_field<< endl;
      }
    

    }
    delete tokens;
  }

if (target == "deuterium")
    {
      std::cout<<"doing deuterium"<<std::endl;
      M_N = (M_p + M_n)/2;
      is_deuterium =1;
      is_hydrogen = 0;
    }
  else   if (target == "hydrogen")
    {
      std::cout<<"doing hydrogen"<<std::endl;
      M_N = M_p;
      is_deuterium =0;
      is_hydrogen =1;
    }
  else {
    std::cout<<"ERROR: invalid target type"<<std::endl;
    return;
  }


  //TCut globalcut = "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&sbs.hcal.e>0.025&&bb.ps.e+bb.sh.e>1.7";

  // cout<<"Applying Global Cut to rootfiles. May take a few minutes...."<<endl;
  // TEventList *elist = new TEventList("elist","Elastic Event List");
  // C->Draw(">>elist",globalcut);


 double BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks],BBgem_track_nhits[maxTracks],
    BBgem_track_ngoodhits[maxTracks], BBgem_track_chi2ndf[maxTracks];
  double BBtr_vz[maxTracks];
  double BBtgt_x[maxTracks], BBtgt_y[maxTracks], BBtgt_th[maxTracks], BBtgt_ph[maxTracks];
  double BBtr_r_th[maxTracks], BBtr_r_ph[maxTracks],BBtr_r_x[maxTracks],BBtr_r_y[maxTracks], BBtr_th[maxTracks],BBtr_ph[maxTracks];
  double BBtr_n, BBps_x, BBps_y, BBps_e, BBps_rowblk, BBps_colblk, BBsh_x, BBsh_y, BBsh_e,BBsh_rowblk, BBsh_colblk,BBsh_atimeblk;
  double HCALx, HCALy, HCALe, ekineW2;
  double HCALx_default, HCALy_default,HCALe_default,HCAL_atime_default;

  double HCAL_clus_e[1000],HCAL_clus_x[1000],HCAL_clus_y[1000],HCAL_clus_row[1000], HCAL_clus_col[1000],HCAL_clus_nblk[1000],HCAL_clus_id[1000],HCAL_clus_tdctime[1000],HCAL_clus_atime[1000];
  int Ndata_HCAL_clus_id;

  double BBgr_clus_size,BBgr_clus_adc,BBgr_clus_trackindex;

  C ->SetBranchStatus("*",0);

  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus("sbs.hcal.atimeblk",1);
  C->SetBranchStatus("sbs.hcal.clus.id",1);
  C->SetBranchStatus("sbs.hcal.clus.e",1);
  C->SetBranchStatus("sbs.hcal.clus.x",1);
  C->SetBranchStatus("sbs.hcal.clus.y",1);
  C->SetBranchStatus("sbs.hcal.clus.row",1);
  C->SetBranchStatus("sbs.hcal.clus.col",1);
  C->SetBranchStatus("sbs.hcal.clus.nblk",1);
  C->SetBranchStatus("sbs.hcal.clus.tdctime",1);
  C->SetBranchStatus("sbs.hcal.clus.atime",1);
  C->SetBranchStatus("Ndata.sbs.hcal.clus.id",1);

  C->SetBranchStatus("bb.grinch_tdc.clus.size",1);
  C->SetBranchStatus("bb.grinch_tdc.clus.adc",1);
  C->SetBranchStatus("bb.grinch_tdc.clus.trackindex",1);

  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.p",1);
  C->SetBranchStatus("bb.tr.vz",1);
  C->SetBranchStatus("bb.tr.th",1);
  C->SetBranchStatus("bb.tr.ph",1);
  C->SetBranchStatus("bb.tr.r_th",1);
  C->SetBranchStatus("bb.tr.r_ph",1);
  C->SetBranchStatus("bb.tr.r_x",1);
  C->SetBranchStatus("bb.tr.r_y",1);
  C->SetBranchStatus("bb.gem.track.nhits",1);
  C->SetBranchStatus("bb.gem.track.ngoodhits",1);
  C->SetBranchStatus("bb.gem.track.chi2ndf",1);
 

  C->SetBranchStatus("bb.ps.e",1);
  C->SetBranchStatus("bb.ps.x",1);
  C->SetBranchStatus("bb.ps.y",1);
  C->SetBranchStatus("bb.ps.rowblk",1);
  C->SetBranchStatus("bb.ps.colblk",1);
  C->SetBranchStatus("bb.sh.e",1);
  C->SetBranchStatus("bb.sh.x",1);
  C->SetBranchStatus("bb.sh.y",1);
  C->SetBranchStatus("bb.sh.atimeblk",1);
  C->SetBranchStatus("bb.sh.rowblk",1);
  C->SetBranchStatus("bb.sh.colblk",1);

  C->SetBranchStatus("e.kine.W2",1);

  C->SetBranchAddress("sbs.hcal.x", &HCALx_default);
  C->SetBranchAddress("sbs.hcal.y", &HCALy_default);
  C->SetBranchAddress("sbs.hcal.e", &HCALe_default);
  C->SetBranchAddress("sbs.hcal.atimeblk",&HCAL_atime_default);
  C->SetBranchAddress("sbs.hcal.clus.id",&HCAL_clus_id);
  C->SetBranchAddress("sbs.hcal.clus.e",&HCAL_clus_e);
  C->SetBranchAddress("sbs.hcal.clus.x",&HCAL_clus_x);
  C->SetBranchAddress("sbs.hcal.clus.y",&HCAL_clus_y);
  C->SetBranchAddress("sbs.hcal.clus.row",&HCAL_clus_row);
  C->SetBranchAddress("sbs.hcal.clus.col",&HCAL_clus_col);
  C->SetBranchAddress("sbs.hcal.clus.nblk",&HCAL_clus_nblk);
  C->SetBranchAddress("sbs.hcal.clus.tdctime",&HCAL_clus_tdctime);
  C->SetBranchAddress("sbs.hcal.clus.atime",&HCAL_clus_atime);
  C->SetBranchAddress("Ndata.sbs.hcal.clus.id",&Ndata_HCAL_clus_id);

  C->SetBranchAddress("bb.grinch_tdc.clus.size",&BBgr_clus_size);
  C->SetBranchAddress("bb.grinch_tdc.clus.adc",&BBgr_clus_adc);
  C->SetBranchAddress("bb.grinch_tdc.clus.trackindex",&BBgr_clus_trackindex);
  
  C->SetBranchAddress("bb.tr.n", &BBtr_n);
  C->SetBranchAddress("bb.tr.vz", BBtr_vz);
  C->SetBranchAddress("bb.tr.px", BBtr_px);
  C->SetBranchAddress("bb.tr.py", BBtr_py);
  C->SetBranchAddress("bb.tr.pz", BBtr_pz);
  C->SetBranchAddress("bb.tr.p", BBtr_p);
  C->SetBranchAddress("bb.tr.th",BBtr_th);
  C->SetBranchAddress("bb.tr.ph",BBtr_ph);
  C->SetBranchAddress("bb.tr.r_th",BBtr_r_th);
  C->SetBranchAddress("bb.tr.r_ph",BBtr_r_ph);
  C->SetBranchAddress("bb.tr.r_y",BBtr_r_y);
  C->SetBranchAddress("bb.tr.r_x",BBtr_r_x);
  C->SetBranchAddress("bb.gem.track.nhits",BBgem_track_nhits);
  C->SetBranchAddress("bb.gem.track.ngoodhits",BBgem_track_ngoodhits);
  C->SetBranchAddress("bb.gem.track.chi2ndf",BBgem_track_chi2ndf);
  

  C->SetBranchAddress("bb.ps.e", &BBps_e);
  C->SetBranchAddress("bb.ps.x", &BBps_x);
  C->SetBranchAddress("bb.ps.y", &BBps_y);
  C->SetBranchAddress("bb.ps.rowblk", &BBps_rowblk);
  C->SetBranchAddress("bb.ps.colblk", &BBps_colblk);
  C->SetBranchAddress("bb.sh.e", &BBsh_e);
  C->SetBranchAddress("bb.sh.x", &BBsh_x);
  C->SetBranchAddress("bb.sh.y", &BBsh_y);
  C->SetBranchAddress("bb.sh.rowblk", &BBsh_rowblk);
  C->SetBranchAddress("bb.sh.colblk", &BBsh_colblk);
  C->SetBranchAddress("e.kine.W2", &ekineW2);
  C->SetBranchAddress("bb.sh.colblk", &BBsh_colblk);
  C->SetBranchAddress("bb.sh.atimeblk",&BBsh_atimeblk);


  // Create output tree
    TTree *P= new TTree("P", "Parsed Tree");
    
    // "old" Tree variables we want to keep
    double bb_tr_vz_out;
    double bb_tr_n_out;
    double bb_tr_p_out;
    double bb_tr_th_out;
    double bb_tr_ph_out;
    double bb_tr_r_th_out;
    double bb_tr_r_ph_out;
    double bb_tr_r_x_out;
    double bb_tr_r_y_out;
    double bb_ps_e_out;
    double bb_ps_rowblk_out;
    double bb_ps_colblk_out;
    double bb_sh_e_out;
    double bb_sh_rowblk_out;
    double bb_sh_colblk_out;
    //double bb_hodotdc_clus_tmean_out;
    //double bb_grinch_tdc_clus_size_out;
    //double bb_grinch_tdc_clus_trackindex_out;
    double bb_gem_track_nhits_out;
    double bb_gem_track_ngoodhits_out;
    double bb_gem_track_chi2ndf_out;
    double bb_etot_over_p_out;
    double bb_grinch_clus_size_out;
    double bb_grinch_clus_adc_out;
    double bb_grinch_clus_trackindex_out;

    // hcal best cluster 
    double hcal_x, hcal_y; // 
    double hcal_e;

    // hcal default cluster
    double hcal_x_def, hcal_y_def;
    double hcal_e_def;

    // calculated variables 
    double hcal_dx, hcal_dy; // 
    double hcal_x_exp, hcal_y_exp; // 
    double Q2_out, W2_out, nu_out, tau_out, epsilon_out, precon_out;

     // mc info. These will all be filled with zeros. I want the branches to match for ease of use later. 
    int is_proton_out=0;
    int is_neutron_out=0;
    double corrected_weight_out=1; // when weight ==1, it's as if it was unweighted. 
    double luminosity_out=0;
    double generation_volume_out=0;
    double max_weight_out=0;
    double Ntried_sum_out=0;

    // which target it is. 
    int is_deuterium_out;
    int is_hydrogen_out;
   
    // fiducial and active area params
    double dxsig_p_out;
    double dxsig_n_out;
    double dysig_out;
    double dx_pn_out; 

    double nsigx_fid_out, nsigy_fid_out;

    // other
    double e_over_p_out;

    // HCal - Shower timing best cluster
    double hcal_clus_atime_out;
    double bb_sh_atimeblk_out;
    double hcal_sh_atime_diff_out;
    int passed_atime_cuts_out;
   

    
    
    


  // Create new branches for the output tree
    P->Branch("bb_tr_n", &bb_tr_n_out, "bb_tr_n/D");
    P->Branch("bb_tr_vz", &bb_tr_vz_out, "bb_tr_vz/D");
    P->Branch("bb_tr_p", &bb_tr_p_out, "bb_tr_p/D");
    P->Branch("bb_tr_th", &bb_tr_th_out, "bb_tr_th/D");
    P->Branch("bb_tr_ph", &bb_tr_ph_out, "bb_tr_ph/D");
    P->Branch("bb_tr_r_x", &bb_tr_r_x_out, "bb_tr_r_x/D");
    P->Branch("bb_tr_r_y", &bb_tr_r_y_out, "bb_tr_r_y/D");
    P->Branch("bb_tr_r_th", &bb_tr_r_th_out, "bb_tr_r_th/D"); // point where straightline tracks begin 
    P->Branch("bb_tr_r_ph", &bb_tr_r_ph_out, "bb_tr_r_ph/D"); // point where straightline tracks begin 
    P->Branch("bb_ps_e", &bb_ps_e_out, "bb_ps_e/D");
    P->Branch("bb_ps_rowblk", &bb_ps_rowblk_out, "bb_ps_rowblk/D");
    P->Branch("bb_ps_colblk", &bb_ps_colblk_out, "bb_ps_colblk/D");
    P->Branch("bb_sh_e", &bb_sh_e_out, "bb_sh_e/D");
    P->Branch("bb_sh_rowblk", &bb_sh_rowblk_out, "bb_sh_rowblk/D");
    P->Branch("bb_sh_colblk", &bb_sh_colblk_out, "bb_sh_colblk/D");
    //P->Branch("bb_hodotdc_clus_tmean", &bb_hodotdc_clus_tmean_out, "bb_hodotdc_clus_tmean/D");
    //P->Branch("bb_grinch_tdc_clus_size", &bb_grinch_tdc_clus_size_out, "bb_grinch_tdc_clus_size/D");
    //P->Branch("bb_grinch_tdc_clus_trackindex", &bb_grinch_tdc_clus_trackindex_out, "bb_grinch_tdc_clus_trackindex/D");
    P->Branch("bb_gem_track_nhits", &bb_gem_track_nhits_out, "bb_gem_track_nhits/D");
    P->Branch("bb_gem_track_ngoodhits", &bb_gem_track_ngoodhits_out, "bb_gem_track_ngoodhits/D");
    P->Branch("bb_gem_track_chi2ndf", &bb_gem_track_chi2ndf_out, "bb_gem_track_chi2ndf/D");
    //P->Branch("bb_etot_over_p", &bb_etot_over_p_out, "bb_etot_over_p/D");
    P->Branch("bb_grinch_clus_size",&bb_grinch_clus_size_out,"bb_grinch_clus_size/D");
    P->Branch("bb_grinch_clus_adc",&bb_grinch_clus_adc_out,"bb_grinch_clus_adc/D");
    P->Branch("bb_grinch_clus_trackindex",&bb_grinch_clus_trackindex_out,"bb_grinch_clus_trackindex/D");

    // calculated variables
    P->Branch("hcal_x", &hcal_x, "hcal_x/D");
    P->Branch("hcal_y", &hcal_y, "hcal_y/D");
    P->Branch("hcal_dx", &hcal_dx, "hcal_dx/D");
    P->Branch("hcal_dy", &hcal_dy, "hcal_dy/D");
    P->Branch("hcal_x_exp", &hcal_x_exp, "hcal_x_exp/D");
    P->Branch("hcal_y_exp", &hcal_y_exp, "hcal_y_exp/D");

    P->Branch("Q2", &Q2_out, "Q2/D");
    P->Branch("W2", &W2_out, "W2/D");
    P->Branch("nu", &nu_out, "nu/D");
    P->Branch("tau", &tau_out, "tau/D");
    P->Branch("epsilon", &epsilon_out, "epsilon/D");
    P->Branch("precon", &precon_out, "precon/D");

    // hcal best cluster
    P->Branch("hcal_e", &hcal_e,"hcal_e/D");

    // Target for Data
    P->Branch("is_deuterium", &is_deuterium_out, "is_deuterium/I");
    P->Branch("is_hydrogen", &is_hydrogen_out, "is_hydrogen/I");
    

    // fiducial and active area params
    P->Branch("dxsig_p", &dxsig_p_out, "dxsig_p/D");
    P->Branch("dxsig_n", &dxsig_n_out, "dxsig_n/D");
    P->Branch("dysig", &dysig_out, "dysig/D");
    P->Branch("dx_pn", &dx_pn, "dx_pn/D");
    P->Branch("nsigx_fid", &nsigx_fid_out, "nsigx_fid/D");
    P->Branch("nsigy_fid", &nsigy_fid_out, "nsigy_fid/D");

    // others
    P->Branch("e_over_p", &e_over_p_out, "e_over_p/D");
    P->Branch("mag_field", &mag_field, "mag_field/D");

    // HCal - Shower timing
    P->Branch("hcal_clus_atime", &hcal_clus_atime_out,"hcal_clus_atime/D");
    P->Branch("bb_sh_atimeblk", &bb_sh_atimeblk_out,"bb_sh_atimeblk/D");
    P->Branch("hcal_sh_atime_diff", &hcal_sh_atime_diff_out,"hcal_sh_atime_diff/D");
    P->Branch("passed_atime_cuts",&passed_atime_cuts_out,"passed_atime_cuts/I");

    // mc info. These will all be filled with zeros. I want the branches to match between data and MC for ease later. 
    P->Branch("is_proton", &is_proton_out, "is_proton/I");
    P->Branch("is_neutron", &is_neutron_out, "is_neutron/I");
    P->Branch("corrected_weight", &corrected_weight_out, "corrected_weight/D");
    P->Branch("luminosity", &luminosity_out, "luminosity/D");
    P->Branch("generation_volume", &generation_volume_out, "generation_volume/D");
    P->Branch("max_weight", &max_weight_out, "max_weight/D");
    P->Branch("Ntried_sum", &Ntried_sum_out, "Ntried_sum/D");


  TH1D *h_dx_HCAL = new TH1D("h_dx_HCAL", " ; x_{HCAL} - x_{exp} (m)  ", 200,-4,4);
  TH1D *h_dy_HCAL = new TH1D("h_dy_HCAL", " ; y_{HCAL} - y_{exp} (m)  ", 200,-2,2);
  TH2D *h_dxdy_HCAL = new TH2D("h_dxdy_HCAL","HCAL_{actual} - HCAL_{expected} positions  ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 100, -2, 2, 200, -4, 4);
  TH1D* h_W2 = new TH1D("h_W2", "W2",100,0,2);
  TH1D* h_W2_test = new TH1D("h_W2_test", "W2",100,0,2);

   TH1D *h_dx_HCAL_failed_ps = new TH1D("h_dx_HCAL_failed_ps", " ; x_{HCAL} - x_{exp} (m)  ", 200,-4,4);
   TH1D* h_W2_failed_ps = new TH1D("h_W2_failed_ps", "W2",100,0,2);



  TH1D *h_dx_HCAL_noFid = new TH1D("h_dx_HCAL_noFid", " ; x_{HCAL} - x_{exp} (m)  ", 200,-4,4);
  TH1D *h_dy_HCAL_noFid = new TH1D("h_dy_HCAL_noFid", " ; y_{HCAL} - y_{exp} (m)  ", 100,-2,2);
  TH2D *h_dxdy_HCAL_noFid = new TH2D("h_dxdy_HCAL_noFid","HCAL_{actual} - HCAL_{expected} positions  ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 100, -2, 2, 200, -4, 4 );
  TH1D* h_W2_noFid = new TH1D("h_W2_noFid", "W2",100,0,2);

  TH1D *h_dx_HCAL_noAA = new TH1D("h_dx_HCAL_noAA", " ; x_{HCAL} - x_{exp} (m)  ", 200,-4,4);
  TH1D *h_dy_HCAL_noAA = new TH1D("h_dy_HCAL_noAA", " ; y_{HCAL} - y_{exp} (m)  ", 100,-2,2);
  TH2D *h_dxdy_HCAL_noAA = new TH2D("h_dxdy_HCAL_noAA","HCAL_{actual} - HCAL_{expected} positions  ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 100, -2, 2, 200, -4, 4 );
  TH1D *h_dx_HCAL_failedFid = new TH1D("h_dx_HCAL_failedFid", " ; x_{HCAL} - x_{exp} (m)  ", 200,-4,4);
  TH1D *h_dy_HCAL_failedFid = new TH1D("h_dy_HCAL_failedFid", " ; y_{HCAL} - y_{exp} (m)  ", 100,-2,2);
  TH1D* h_W2_noAA = new TH1D("h_W2_noAA", "W2",100,0,2);

  TH2D* h_xy_HCAL = new TH2D("h_xy_HCAL", "Actual HCAL Positions ;HCAL Y_{actual} ; HCAL X_{actual}",100,-1,1,225,-3,1.5);
  TH2D* h_xy_HCAL_exp = new TH2D("h_xy_HCAL_exp", "Expected HCAL positions;HCAL Y_{expected} ; HCAL X_{expected}",100,-1,1,225,-3,1.5);
  TH2D* h_xy_HCAL_exp_failedFid = new TH2D("h_xy_HCAL_exp_failedFid", "Expected: Failed Fiducial;HCAL Y_{expected} ; HCAL X_{expected}",100,-1,1,225,-3,1.5);
  TH2D* h_xy_HCAL_exp_passedFid = new TH2D("h_xy_HCAL_exp_passedFid", "Expected: Passed Fiducial;HCAL Y_{expected} ; HCAL X_{expected}",100,-1,1,225,-3,1.5);
  TH2D* h_xy_HCAL_passedFid_n = new TH2D("h_xy_HCAL_passedFid_n", "Expected: Passed Fiducial for Neutron hyp.;HCAL Y_{expected} ; HCAL X_{expected}",100,-1,1,225,-3,2);
  TH2D* h_xy_HCAL_passedFid_p = new TH2D("h_xy_HCAL_passedFid_p", "Expected: Passed Fiducial for Proton hyp.;HCAL Y_{expected} ; HCAL X_{expected}",100,-1,1,225,-3,2);

   TH1D* h_nsigx_fid = new TH1D("h_nsigx_fid", "nsigx_fid",200,-20,20);
   TH1D* h_nsigy_fid = new TH1D("h_nsigy_fid", "nsigy_fid",200,-20,20);
   TH1D *h_dx_HCAL_nsig = new TH1D("h_dx_HCAL_nsig ", " ; x_{HCAL} - x_{exp} (m)  ", 200,-4,4);


  TH2D* h_W2_dx =  new TH2D("h_W2_dx",";HCAL dx; W2",500,-4,3,100,0,4);
  TH2D* h_p_dx = new TH2D("h_p_dx",";HCAL dx;tr.px;",350,-4,3,100,0,2);

  TH1D* h_atime_diff = new TH1D("h_atime_diff","HCAL_clus_atime[i] - BBsh_atimeblk",300,-100,100);
  TH1D* h_atime_diff_failed = new TH1D("h_atime_diff_failed","HCAL_clus_atime[i] - BBsh_atimeblk",300,-100,100);
  TH1D* h_atime_diff_passed = new TH1D("h_atime_diff_passed","HCAL_clus_atime[i] - BBsh_atimeblk",300,-100,100);
  TH1D* h_HCAL_clus_atime = new TH1D("h_HCAL_clus_atime","sbs.hcal.clus.atime", 300,-100,200);
  TH1D* h_BBsh_atime =  new TH1D("h_BBsh_atime","bb.sh.atimeblk",400,-50,50);
  TH1D* h_HCAL_max_clus_e =  new TH1D("h_HCAL_max_clus_e","max energy for HCAL cluster that also passed adc coin cut ",100,0,0.5);
  TH1D* h_max_e_index = new TH1D("h_max_e_index","index of best hcal cluster",21,-1,20);
  TH1D* h_atime_HCAL_failed = new TH1D("h_atime_HCAL_failed","HCAL_clus_atime[i] - BBsh_atimeblk",300,-100,100);
  TH1D* h_atime_HCAL_passed = new TH1D("h_atime_HCAL_passed","HCAL_clus_atime[i] - BBsh_atimeblk",300,-100,100);

  //// define hcal active area
  std::vector<Double_t> hcalaa = hcalaa_mc(1, 1);//
  cout<< "HCal acitve area cut params: "<< hcalaa[0] << ", " << hcalaa[1] << ", " << hcalaa[2] << ", " << hcalaa[3] <<endl;

  

  //// define hcal fiducial cut 
  std::vector<Double_t> hcalfid_data =hcalfid(nsigma_p*dxsig_p, nsigma_n*dxsig_n, nsigma_y*dysig, hcalaa);
  cout<< "HCal fiducial cut params: "<< hcalfid_data[0] << ", " << hcalfid_data[1] << ", " << hcalfid_data[2] << ", " << hcalfid_data[3] <<endl;
  cout<< "dxsig_p =  "<< dxsig_p <<", dxsig_n = "<<dxsig_n<< ", dysig = " <<dysig<<", dx_pn = "<<dx_pn<<endl;





  //Long64_t Nevents = elist ->GetN();

 // Set long int to keep track of total entries
  cout<<"Loading branches. Hold tight! "<<endl;
  Long64_t Nevents = C->GetEntries();
  UInt_t run_number = 0;
 
  cout<<"Entries: "<<Nevents<<endl;


  Int_t max = 0;
    if (entries_input == -1){
      max = Nevents;
    }
    else{
      max = entries_input; 
    }
    if (max > Nevents){ max = Nevents;}
    cout<<"max = "<<max<<endl;


  cout<< "# Events = "<<Nevents<<endl;
  for( Long64_t nevent = 1; nevent <max; nevent++){

    if (nevent%50000==0) cout << " Entry = " << nevent << endl;
    
    C->GetEntry(nevent); 

      // doing a global cut manually here for now. 
    if (BBtr_n ==1 && abs(BBtr_vz[0])<0.08 && BBps_e > 0.1 && BBgem_track_nhits[0]>2){ 
  
 double Eloss_outgoing = celldiameter/2.0/sin(BB_th) * ld2tarrho * ld2dEdx;

      double ebeam_c = E_e - ( (BBtr_vz[0]+l_tgt/2.0) * lh2tarrho * lh2dEdx + lh2uwallthick * rho_al * aldEdx );//correct for energy loss of beam electron as it passes into the target material towards the vertex position.

      TVector3 vertex(0,0,BBtr_vz[0]);

      //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
      Double_t precon = BBtr_p[0] + Eloss_outgoing; //correct the momentum of e' after the vertex as it passes out of the target material.
      
      TLorentzVector Pbeam(0,0,ebeam_c,ebeam_c); //beam momentum
      //TLorentzVector pe( px[0], py[0], pz[0], p[0] ); //e' momentum
      TLorentzVector pe( precon*BBtr_px[0]/BBtr_p[0], precon*BBtr_py[0]/BBtr_p[0], precon*BBtr_pz[0]/BBtr_p[0], precon ); //e' recon plvect
      TLorentzVector Ptarg(0, 0, 0, M_N);
      TLorentzVector q = Pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
      TLorentzVector pN; //N' momentum
      //	TVector3 pNhat; //Unit N' 3-vector

      double etheta =  acos( pe.Pz() / pe.E() );      
      double ephi = atan2( pe.Py(), pe.Px() ); 
      double pcent =E_e/( 1. + ( E_e/M_N )*( 1.0 - cos(etheta) ) );
      double phNexp = ephi + PI;
      double Q2, W2, nu, thNexp, pNexp, ebeam_o, tau, epsilon;
      ebeam_o =  ebeam_c - celldiameter/2.0/sin(etheta) * ld2tarrho * ld2dEdx; //Second energy correction accounting for energy loss leaving target. (Although S doesn't seem to use this)
       
      //reconstruct track with angles (better resolution with GEMs)
      nu = Pbeam.E() - pcent;
      pNexp = sqrt( pow(nu, 2.) + 2. * M_N * nu );
      thNexp = acos( (Pbeam.E()-pcent*cos(etheta)) / pNexp );
      TVector3 pNhat( sin(thNexp) * cos(phNexp), sin(thNexp) * sin(phNexp), cos(thNexp) );  //Unit N' 3-vector
      pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+Ptarg.E() );
      Q2 = 2.0 * Pbeam.E()* pe.E() * ( 1.0 - cos(etheta) );
      W2 = pow( M_N, 2.0 ) + 2.0*M_N * (Pbeam.E()-pe.E()) - Q2;
      tau =  Q2 / ( 4.0 * pow( M_N, 2.0 ) );
      epsilon = 1.0 / (1.0 + 2.0 * (1.0 + tau) * pow(tan(etheta / 2.0), 2));

      // HCAL 
      TVector3 hcal_zaxis(sin(-HCal_th),0,cos(-HCal_th)); // Clock-wise rotation about Y axis
      TVector3 hcal_xaxis( 0, -1, 0 );  // -Y axis of Hall CoS = X axis of hcal CoS
      TVector3 hcal_yaxis = hcal_zaxis.Cross(hcal_xaxis).Unit();
      TVector3 hcal_origin = HCal_d * hcal_zaxis + hcalheight * hcal_xaxis;
      double sintersect = ( hcal_origin - vertex).Dot( hcal_zaxis) / (pNhat.Dot( hcal_zaxis));    // intersection of a ray with a plane in hcal coordinates
      TVector3 hcal_intersect = vertex + sintersect * pNhat; 
      double yexpect_HCAL = (hcal_intersect - hcal_origin).Dot( hcal_yaxis);
      double xexpect_HCAL =  (hcal_intersect - hcal_origin).Dot( hcal_xaxis);

      double e_over_p = ( BBps_e+ BBsh_e ) / BBtr_p[0];

      Double_t nsigy_fid = calculate_nsigma_fid_y( yexpect_HCAL,  dysig, hcalaa);
      //cout<<"n_sigy = "<<nsigy_fid<<endl;
      Double_t nsigx_fid = calculate_nsigma_fid_x(xexpect_HCAL, dxsig_n, dxsig_p,  dx_pn,  hcalaa);
      //cout<<"n_sigx = "<<nsigx_fid<<endl;

      if (W2>3) continue; 

      

      // /// what i had before 
      // //Phyiscs calcuations
    
      // double  ebeam_c = E_e - ( (BBtr_vz[0]+l_tgt/2.0) * ld2tarrho * ld2dEdx + ld2uwallthick * rho_al * aldEdx ); // where do we get the ebeam from? Because I'm just reading it in from a file right now. 


      // double etheta = acos( BBtr_pz[0]/BBtr_p[0]);
      // double ephi = atan2(BBtr_py[0],BBtr_px[0]);
    
      // TVector3 vertex(0,0,BBtr_vz[0]);
      // TLorentzVector Pbeam(0,0,E_e,E_e);
      // TLorentzVector kprime(BBtr_px[0], BBtr_py[0], BBtr_pz[0], BBtr_p[0]);
      // TLorentzVector Ptarg(0, 0, 0, M_p);

      // TLorentzVector q = Pbeam - kprime;
      // TLorentzVector PgammaN = Ptarg +q; //should go through and write this out. Momentum of virtual photon

      // double pel =( E_e/ (1. +E_e/M_p*(1.-cos(etheta))) + E_e/ (1. +E_e/M_n*(1.-cos(etheta))) ) /2 ;
      // //centural momentum of elastically scattered electron. We don't know if it's a proton or a neutron so we are going to take the average.  
      // double nu = E_e -BBtr_p[0]; //kinetic energy of the elasticlly scattered electron 
      // double pp = sqrt(pow(nu,2)+2 *( (M_p+M_n)/2 )*nu); 
      // double phinucleon = ephi + PI; //coplanar 
      // //double thetanucleon = acos((E_e - pel*cos(etheta)) / pp); //
      // double thetanucleon = acos((E_e - BBtr_pz[0])/pp); //
      // // cout<<thetanucleon<<" "<<thetanucleon2<<endl;

      // TVector3 pNhat( sin(thetanucleon) *cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));
    
      // TVector3 HCAL_zaxis(sin(-HCal_th),0,cos(-HCal_th));
      // TVector3 HCAL_xaxis(0,-1,0);
      // TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();

      // TVector3 HCAL_origin = HCal_d * HCAL_zaxis + hcalheight * HCAL_xaxis;

      // double sintersect = ( HCAL_origin - vertex).Dot( HCAL_zaxis) / (pNhat.Dot( HCAL_zaxis));
    
      // TVector3 HCAL_intersect = vertex + sintersect * pNhat; 

      // double yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis);
      // double xexpect_HCAL =  (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis);

      // double Q2 = 2*E_e*BBtr_p[0]*( 1-BBtr_pz[0]/BBtr_p[0]);
    
      // double W = PgammaN.M();
      // double W2 = ekineW2;


      h_BBsh_atime ->Fill(BBsh_atimeblk);

      ///////////////////////////////////////////////
      // HCAL best cluster analysis
      /////////////////////////////////////////////
      std::vector<bool> passed_atime_diff;
      std::vector<bool> passed_HCAL_atime;
      // Loop over the clusters and see if they are within the expected HCAL time and the coinicidence time with the shower. 
      for (int i = 0; i< Ndata_HCAL_clus_id; i++)
	{
	  double atime_diff = HCAL_clus_atime[i] - BBsh_atimeblk;
	  h_atime_diff ->Fill(atime_diff); 
	  h_HCAL_clus_atime ->Fill(HCAL_clus_atime[i]);

	  // see if the HCAL_clus_atime for this cluster is within a few sigma of the mean of the distribution 
	  if (abs(HCAL_clus_atime[i] -HCAL_atime_mean)<=Nsigma_HCAL_atime*HCAL_atime_sigma){
	    passed_HCAL_atime.push_back(true); 
	    h_atime_HCAL_passed ->Fill(HCAL_clus_atime[i]);
	  } 
	  else {
	    passed_HCAL_atime.push_back(false);
	    h_atime_HCAL_failed ->Fill(HCAL_clus_atime[i]);
	  }

	  // see if the difference in adc time between HCAL cluster and the shower cluster is within a few sigma of the mean of the distribution
	  if (abs(atime_diff - atime_diff_mean)<=Nsigma_atime_diff*atime_diff_sigma){
	    passed_atime_diff.push_back(true); 
	    h_atime_diff_passed ->Fill(atime_diff);
	  } 
	  else {
	    passed_atime_diff.push_back(false);
	    h_atime_diff_failed ->Fill(atime_diff);
	  }
	}

      // find the cluster with the max energy that passes the adc time checks we made above
      Double_t max_e = 0;
      Int_t max_e_index = -1;
      for (int i = 0; i< Ndata_HCAL_clus_id; i++)
	{
	  if (passed_atime_diff[i] && passed_HCAL_atime[i]){
	    if (HCAL_clus_e[i] > max_e){
	      max_e = HCAL_clus_e[i];
	      max_e_index = i;
	    }
	  }	
	}

      h_max_e_index ->Fill(max_e_index);

      
      int passed_atime_cuts = 1; // flag to say if the a cluster passed the in-time and coin-time cuts
      
      if(max_e_index == -1){// no clusters passed the adc time checks
	max_e_index = 0; // let it fill with the default highest energy values. We'll have to cut on hcal time and coin time later.
	passed_atime_cuts = 0; // put up a flag that it failed one or both of the atime cuts and is using the default cluster. 
      }
    
     
      
      //if (abs(ekineW2- W2_mean)>W2_sig ) continue; // I probably should remove this when I turn this into a parsing script. Or just make it > 1.5 or something. 
      // if (ekineW2>2 ) continue; // I probably should remove this when I turn this into a parsing script. Or just make it > 1.5 or something. 


      h_HCAL_max_clus_e ->Fill(HCAL_clus_e[max_e_index]);

      HCALy = HCAL_clus_y[max_e_index];
      HCALx = HCAL_clus_x[max_e_index];
      HCALe = HCAL_clus_e[max_e_index];

      if( HCALe < 0.01) continue;// may also want to eleminate this when making a parsing script? Or at least lower it a bit. // 0.025 is what others were using. 

      // HCALy = HCALy_default;
      // HCALx = HCALx_default;
      // HCALe = HCALe_default;

      double dx = HCALx - xexpect_HCAL;
      double dy = HCALy - yexpect_HCAL;

      h_W2_dx ->Fill(dx, ekineW2); 

      if(BBps_e<0.2)
	{
	  h_dx_HCAL_failed_ps ->Fill(dx);
	  h_W2_failed_ps->Fill(W2);
	}

      h_xy_HCAL ->Fill(HCALy, HCALx);
      h_xy_HCAL_exp ->Fill(yexpect_HCAL, xexpect_HCAL);

      h_dx_HCAL_noAA ->Fill(dx);
      h_dy_HCAL_noAA ->Fill(dy);
      h_dxdy_HCAL_noAA ->Fill(dy,dx);
      h_W2_noAA ->Fill(ekineW2);



      //// check if it's on the hcal active area
      if ( hcalaaON( HCALx , HCALy , hcalaa))// Uses HCal actual positions. 
	{
	  h_dx_HCAL_noFid ->Fill(dx);
	  h_dy_HCAL_noFid ->Fill(dy);
	  h_dxdy_HCAL_noFid ->Fill(dy,dx);
	  h_W2_noFid ->Fill(ekineW2);
	  //// check if it's within the fiducial cuts 
	  if(  hcalfidIN(xexpect_HCAL, yexpect_HCAL, dx_pn, hcalfid_data ) ) 
	    {
	      h_dx_HCAL ->Fill(dx);
	      h_dy_HCAL ->Fill(dy);
	      h_dxdy_HCAL ->Fill(dy,dx);
	      h_W2 ->Fill(ekineW2);
	      h_xy_HCAL_exp_passedFid ->Fill(yexpect_HCAL, xexpect_HCAL);
	      h_p_dx ->Fill(dx, BBtr_px[0]);
	      h_xy_HCAL_passedFid_n ->Fill(yexpect_HCAL, xexpect_HCAL);
	      h_xy_HCAL_passedFid_p ->Fill(yexpect_HCAL, xexpect_HCAL - dx_pn);
	    }
	}
      if (!hcalaaON( xexpect_HCAL , yexpect_HCAL , hcalaa) || !hcalfidIN(xexpect_HCAL, yexpect_HCAL, dx_pn, hcalfid_data ))
	{
	  h_dx_HCAL_failedFid ->Fill(dx);
	  h_dy_HCAL_failedFid ->Fill(dy);
	  h_xy_HCAL_exp_failedFid ->Fill(yexpect_HCAL, xexpect_HCAL);
	}

      // testing my nsigx and nsigy 
      h_nsigx_fid ->Fill(nsigx_fid);
      h_nsigy_fid ->Fill(nsigy_fid);
      if (nsigx_fid > 1 && nsigy_fid >1)
	{
	  h_dx_HCAL_nsig->Fill(dx);
	}

 // Update and fill the output branches 
      //Fill old output tree
      bb_tr_p_out = BBtr_p[0];
      bb_tr_n_out = BBtr_n;
      bb_tr_vz_out = BBtr_vz[0];
      bb_tr_th_out = BBtr_th[0];
      bb_tr_ph_out = BBtr_ph[0];
      bb_tr_r_x_out = BBtr_r_x[0];
      bb_tr_r_y_out = BBtr_r_y[0];
      bb_tr_r_th_out = BBtr_r_th[0];
      bb_tr_r_ph_out = BBtr_r_ph[0];
      bb_ps_e_out = BBps_e;
      bb_ps_rowblk_out = BBps_rowblk;
      bb_ps_colblk_out = BBps_colblk;
      bb_sh_e_out = BBsh_e;
      bb_sh_rowblk_out = BBsh_rowblk;
      bb_sh_colblk_out = BBsh_colblk;
      //bb_hodotdc_clus_tmean_out = hodotmean[0];
      //bb_grinch_tdc_clus_size_out = grinchClusSize;
      //bb_grinch_tdc_clus_trackindex_out = grinchClusTrIndex;
      bb_gem_track_nhits_out = BBgem_track_nhits[0];
      bb_gem_track_ngoodhits_out = BBgem_track_ngoodhits[0];
      bb_gem_track_chi2ndf_out = BBgem_track_chi2ndf[0];
      //bb_etot_over_p_out = eop;
      bb_grinch_clus_size_out = BBgr_clus_size;
      bb_grinch_clus_adc_out = BBgr_clus_adc;
      bb_grinch_clus_trackindex_out= BBgr_clus_trackindex;

      // HCAL projections and variables 
      hcal_x = HCALx;
      hcal_y = HCALy;
      hcal_dx = dx;
      hcal_dy = dy;
      hcal_x_exp = xexpect_HCAL;
      hcal_y_exp = yexpect_HCAL;
      Q2_out = Q2;
      W2_out = W2;
      nu_out = nu;
      tau_out = tau;
      epsilon_out = epsilon;
      precon_out = precon;

      // HCAL best cluster 
      hcal_e = HCAL_clus_e[max_e_index];

      // mc info
      is_hydrogen_out = is_hydrogen;
      is_deuterium_out  = is_deuterium;

      // fiducial and active area params
      dxsig_p_out =  dxsig_p;
      dxsig_n_out =  dxsig_n;
      dysig_out =  dysig;
      dx_pn_out = dx_pn; 
      nsigx_fid_out = nsigx_fid;
      nsigy_fid_out = nsigy_fid;

      // other
      e_over_p_out = e_over_p;

      // Hcal - Shower timing

      hcal_clus_atime_out = HCAL_clus_atime[max_e_index];
      bb_sh_atimeblk_out = BBsh_atimeblk;
      hcal_sh_atime_diff_out = HCAL_clus_atime[max_e_index] - BBsh_atimeblk;
      passed_atime_cuts_out = passed_atime_cuts;
   
      P -> Fill();
    }// end global cut  
 

  }// end event loop

  //// Drawing the active area cuts and the fiducial cuts to visulaize them. 

  double aaXi = hcalaa[0];
  double aaXf = hcalaa[1];
  double aaYi = hcalaa[2];
  double aaYf = hcalaa[3]; 

  double fidXi = hcalfid_data[0];
  double fidXf = hcalfid_data[1];
  double fidYi = hcalfid_data[2];
  double fidYf = hcalfid_data[3];
  

 
  TLine *LineXi = new TLine(aaYi,aaXi, aaYf, aaXi); //horizontal line at aaXi from Yi to Yf.  geometry: (x1, y1, x2, y2)
  LineXi->SetLineWidth(2);
  // LineXi->SetLineStyle(2);
  LineXi ->SetLineColor(kRed);
  TLine *LineXf =  new TLine(aaYi,aaXf, aaYf, aaXf); //horizontal line at aaXf from Yi to Yf.
  LineXf->SetLineWidth(2);
  //LineXf->SetLineStyle(2);
  LineXf ->SetLineColor(kRed);
  TLine *LineYi = new TLine(aaYi,aaXi, aaYi, aaXf); //vert line at aaYi from Xi to Xf.
  LineYi->SetLineWidth(2);
  //LineYi->SetLineStyle(2);
  LineYi ->SetLineColor(kRed);
  TLine *LineYf = new TLine(aaYf,aaXi, aaYf, aaXf); //vert line at aaYf from Xi to Xf.
  LineYf->SetLineWidth(2);
  //LineYf->SetLineStyle(2);
  LineYf ->SetLineColor(kRed);

  TLine *LineFidXi = new TLine(fidYi,fidXi, fidYf, fidXi); //horizontal line at fidXi from Yi to Yf.  geometry: (x1, y1, x2, y2)
  LineFidXi->SetLineWidth(2);
  // LineFidXi->SetLineStyle(2);
  LineFidXi ->SetLineColor(kMagenta);
  TLine *LineFidXf =  new TLine(fidYi,fidXf, fidYf, fidXf); //horizontal line at fidXf from Yi to Yf.
  LineFidXf->SetLineWidth(2);
  //LineFidXf->SetLineStyle(2);
  LineFidXf ->SetLineColor(kMagenta);
  TLine *LineFidXproton =  new TLine(fidYi, fidXi + dx_pn, fidYf, fidXi + dx_pn); //horizontal line at fidXi+dx_pn from Yi to Yf.
  LineFidXproton->SetLineWidth(2);
  LineFidXproton->SetLineStyle(2);
  LineFidXproton ->SetLineColor(kMagenta);
  TLine *LineFidYi = new TLine(fidYi,fidXi, fidYi, fidXf); //vert line at fidYi from Xi to Xf.
  LineFidYi->SetLineWidth(2);
  // LineFidYi->SetLineStyle(2);
  LineFidYi ->SetLineColor(kMagenta);
  TLine *LineFidYf = new TLine(fidYf,fidXi, fidYf, fidXf); //vert line at fidYf from Xi to Xf.
  LineFidYf->SetLineWidth(2);
  //LineFidYf->SetLineStyle(2);
  LineFidYf ->SetLineColor(kMagenta);

  TLine *LinePosXi = new TLine(hcalposYi_mc,hcalposXi_mc, hcalposYf_mc, hcalposXi_mc); //horizontal line at hcalposXi from Yi to Yf.  geometry: (x1, y1, x2, y2)
  LinePosXi->SetLineWidth(2);
  // LinePosXi->SetLineStyle(2);
  LinePosXi ->SetLineColor(kGreen);
  TLine *LinePosXf =  new TLine(hcalposYi_mc,hcalposXf_mc,hcalposYf_mc, hcalposXf_mc); //horizontal line at hcalposXf from Yi to Yf.
  LinePosXf->SetLineWidth(2);
  //LinePosXf->SetLineStyle(2);
  LinePosXf ->SetLineColor(kGreen);
  TLine *LinePosYi = new TLine(hcalposYi_mc,hcalposXi_mc,hcalposYi_mc, hcalposXf_mc); //vert line at hcalposYi from Xi to Xf.
  LinePosYi->SetLineWidth(2);
  // LinePosYi->SetLineStyle(2);
  LinePosYi ->SetLineColor(kGreen);
  TLine *LinePosYf = new TLine(hcalposYf_mc,hcalposXi_mc,hcalposYf_mc,hcalposXf_mc); //vert line at hcalposYf from Xi to Xf.
  LinePosYf->SetLineWidth(2);
  //LinePosYf->SetLineStyle(2);
  LinePosYf ->SetLineColor(kGreen);
  



  TCanvas* canvas = new TCanvas("canvas", "Canvas");
  canvas->Divide(3,1);
  canvas ->cd(1);
  h_xy_HCAL ->Draw("colz");
  LineXi ->Draw("same");
  LineXf ->Draw("same");
  LineYi ->Draw("same");
  LineYf->Draw("same");
  LineFidXi ->Draw("same");
  LineFidXf ->Draw("same");
  LineFidYi ->Draw("same");
  LineFidYf->Draw("same");
  LinePosXi ->Draw("same");
  LinePosXf ->Draw("same");
  LinePosYi ->Draw("same");
  LinePosYf->Draw("same");
  canvas ->cd(2);
  h_xy_HCAL_exp_passedFid->Draw("colz");
  LineXi ->Draw("same");
  LineXf ->Draw("same");
  LineYi ->Draw("same");
  LineYf->Draw("same");
  LineFidXi ->Draw("same");
  LineFidXf ->Draw("same");
  LineFidYi ->Draw("same");
  LineFidYf->Draw("same");
  LineFidXproton->Draw("same");
  LinePosXi ->Draw("same");
  LinePosXf ->Draw("same");
  LinePosYi ->Draw("same");
  LinePosYf->Draw("same");
  canvas ->cd(3);
  h_xy_HCAL_exp_failedFid->Draw("colz");
  LineXi ->Draw("same");
  LineXf ->Draw("same");
  LineYi ->Draw("same");
  LineYf->Draw("same");
  LineFidXi ->Draw("same");
  LineFidXf ->Draw("same");
  LineFidYi ->Draw("same");
  LineFidYf->Draw("same");
  LineFidXproton->Draw("same");
  LinePosXi ->Draw("same");
  LinePosXf ->Draw("same");
  LinePosYi ->Draw("same");
  LinePosYf->Draw("same");

  TCanvas* canvas2 = new TCanvas("canvas2", "Canvas");
  canvas2->Divide(2,1);
  canvas2 ->cd(1);
  h_xy_HCAL_passedFid_n->Draw("colz");
  LineFidXi ->Draw("same");
  LineFidXf ->Draw("same");
  LineFidYi ->Draw("same");
  LineFidYf->Draw("same");
  LineFidXproton->Draw("same");
  canvas2 ->cd(2);
  h_xy_HCAL_passedFid_p->Draw("colz");
  LineFidXi ->Draw("same");
  LineFidXf ->Draw("same");
  LineFidYi ->Draw("same");
  LineFidYf->Draw("same");
  LineFidXproton->Draw("same");

  //// end visualizig cuts 

  fout->Write();
  P->Write();

}// end main

// BEGIN FUNCTIONS 

 // Establish hcal active area excluding N blks from edge, Pass0/1 DB
  std::vector<Double_t> hcalaa_data(int exblkN_x=1, int exblkN_y=1) {
    std::vector<Double_t> hcalaa;
    Double_t hcalaaXi = hcalposXi_p0 + exblkN_x*hcalblk_w_p0;
    Double_t hcalaaXf = hcalposXf_p0 - exblkN_x*hcalblk_w_p0;
    Double_t hcalaaYi = hcalposYi_p0 + exblkN_y*hcalblk_h_p0;
    Double_t hcalaaYf = hcalposYf_p0 - exblkN_y*hcalblk_h_p0;
    hcalaa.push_back( hcalaaXi ); 
    hcalaa.push_back( hcalaaXf );
    hcalaa.push_back( hcalaaYi );
    hcalaa.push_back( hcalaaYf );
    return hcalaa;
  }

// Establish hcal active area excluding N blks from edge, MC DB
std::vector<Double_t> hcalaa_mc (int exblkN_x=1, int exblkN_y=1) {
  std::vector<Double_t> hcalaa;
  Double_t hcalaaXi = hcalposXi_mc + exblkN_x*hcalblk_div_v;
  Double_t hcalaaXf = hcalposXf_mc - exblkN_x*hcalblk_div_v;
  Double_t hcalaaYi = hcalposYi_mc + exblkN_y*hcalblk_div_h;
  Double_t hcalaaYf = hcalposYf_mc - exblkN_y*hcalblk_div_h;
  hcalaa.push_back( hcalaaXi ); 
  hcalaa.push_back( hcalaaXf );
  hcalaa.push_back( hcalaaYi );
  hcalaa.push_back( hcalaaYf );
  return hcalaa;
}

 // Check position per event against hcal active area (TRUE if detection on active area)
  bool hcalaaON (Double_t hcalx, Double_t hcaly, std::vector<Double_t> hcalaa) {
    bool on = false;
    // active area dimensions
    Double_t hcalx_t = hcalaa[0];
    Double_t hcalx_b = hcalaa[1];
    Double_t hcaly_r = hcalaa[2];
    Double_t hcaly_l = hcalaa[3];
    on = hcaly>hcaly_r && hcaly<hcaly_l && hcalx>hcalx_t && hcalx<hcalx_b;
    return on;
  } 


// Overload for most cases where dy margin negligable and dx top/bottom configurable by p/n
  std::vector<Double_t> hcalfid (Double_t dxsig_p, Double_t dxsig_n, Double_t dysig, std::vector<Double_t> hcalaa) {
    std::vector<Double_t> fid;
    Double_t hcalx_t = hcalaa[0] + dxsig_p;  // top margin (relevant for proton)
    Double_t hcalx_b = hcalaa[1] - dxsig_n;  // bottom margin (relevant for neutron)
    Double_t hcaly_r = hcalaa[2] + dysig;    // right margin
    Double_t hcaly_l = hcalaa[3] - dysig;    // left margin
    fid.push_back( hcalx_t ); 
    fid.push_back( hcalx_b );
    fid.push_back( hcaly_r );
    fid.push_back( hcaly_l );
    return fid;
  }


  // Check position by event and verify in hcal fiducial area
  bool hcalfidIN (Double_t hcalx_exp, Double_t hcaly_exp, Double_t dx_pn, vector<Double_t> fid) {
    Double_t hcalx_t = fid[0];
    Double_t hcalx_b = fid[1];
    Double_t hcaly_r = fid[2];
    Double_t hcaly_l = fid[3];

    Double_t hcalx_exp_p = hcalx_exp - dx_pn;      //define the exp pos of a proton from obs dx peak diff

    bool infid = hcaly_exp>hcaly_r && hcaly_exp<hcaly_l &&      //dy same for protons and neutrons
	         hcalx_exp>hcalx_t && hcalx_exp<hcalx_b &&      //dx for neutrons
	         hcalx_exp_p>hcalx_t && hcalx_exp_p<hcalx_b;	//dx for protons							     
    return infid;
  } 

Double_t calculate_nsigma_fid_y(Double_t hcaly_exp,Double_t dysig, std::vector<Double_t> hcalaa)
{
  Double_t hcalaa_t = hcalaa[0];
  Double_t hcalaa_b = hcalaa[1] ;
  Double_t hcalyaa_neg = hcalaa[2]; // I'm a bit confused as for if the "right" and "left" are from hcal's perspective or the transport coord perspective, so I'll just call them "neg" and "pos" for now. hcalyaa_neg = -0.77 and hcalyaa_pos = 0.77 with 1 block
  Double_t hcalyaa_pos = hcalaa[3];
  
  // goal: how many sigma away from each active area boundary is the expected value?
  // Also want to incorporate information on if it is on the wrong side of the boundary: ie allow for negative values. 
  // Example: say hcaly_exp = 0.1
  Double_t distance_to_negative_y_aa = hcaly_exp - hcalyaa_neg; // 0.1 - (-0.77) = 0.87. So it is 0.87 away from the negative active area boundary.  
  Double_t nsigma_neg = distance_to_negative_y_aa  / dysig ;  // express in terms of how many dy sigma that is. 
  // example: if it was outside the boundary, say hcaly_exp = -0.8:  -0.8 - (-0.77) = -0.03. The sign shows it's outside the boundary. 
  Double_t distance_to_positive_y_aa = hcalyaa_pos - hcaly_exp; // 0.77 - 0.1 = 0.67 
  // example: if it was outside the boundary, say hcaly_exp = 0.8:  0.77 - 0.8 = -0.03. The sign shows it's outside the boundary. 
  Double_t nsigma_pos =  distance_to_positive_y_aa  / dysig ;// express in terms of how many dy sigma that is. 

  // if(nsigma_neg <0 ||nsigma_pos<0){
  //   cout<< "hcalyaa_neg = "<< hcalyaa_neg << " , hcalyaa_pos = "<< hcalyaa_pos << " , hcaly_exp = " <<hcaly_exp<<endl;
  //   cout << "nsigma_neg = "<<nsigma_neg<<" , nsigma_pos = "<<nsigma_pos<<endl;
  // }
  return std::min(nsigma_pos,nsigma_neg);
}

Double_t calculate_nsigma_fid_x(Double_t hcalx_exp,Double_t dxsig_n, Double_t dxsig_p, Double_t dx_pn, std::vector<Double_t> hcalaa)
{
  Double_t hcalaa_t = hcalaa[0]; // about -2.5 for 1 block aa cut
  Double_t hcalaa_b = hcalaa[1] ;// about 1.0 for 1 block aa cut 
  
   Double_t hcalx_exp_p = hcalx_exp - dx_pn;      //define the exp pos of a proton from obs dx peak diff

// goal: how many sigma away from each active area boundary is the expected value?
  // Also want to incorporate information on if it is on the wrong side of the boundary: ie allow for negative values. 
   // example: say hcalx_exp = 0.1
   Double_t distance_to_top_x_aa_n = hcalx_exp - hcalaa_t; // 0.1 - (-2.5) = 2.6. So it is 2.6 away from the top active area boundary.
   Double_t distance_to_top_x_aa_p = hcalx_exp_p - hcalaa_t; // proton hypothesis 
   Double_t nsigma_top_n = distance_to_top_x_aa_n / dxsig_p; // express in terms of how many dxsig_p that is. 
   Double_t nsigma_top_p = distance_to_top_x_aa_p / dxsig_p;// proton hypotheis. 
   /// what if outside the aa boundary? say hcalx_exp = -2.6:  -2.6 - (-2.5) = -0.1. The sign shows it's outside the boundary. 

   Double_t distance_to_bottom_x_aa_n = hcalaa_b - hcalx_exp; // 1.0 - 0.1 = 0.9. So it is 0.9 away from the bottom active area boundary.
   Double_t distance_to_bottom_x_aa_p = hcalaa_b - hcalx_exp_p; // proton hypothesis.
   Double_t nsigma_bottom_n = distance_to_bottom_x_aa_n / dxsig_n; // express in terms of how many dxsig_n that is. 
   Double_t nsigma_bottom_p = distance_to_bottom_x_aa_p / dxsig_n;// proton hypotheis. 
   /// what if outside the aa boundary? say hcalx_exp = 2.0:  1.0 - 2.0  = -1. The sign shows it's outside the boundary. 


   //cout<< "hcalaa_t = "<< hcalaa_t << " , hcalyaa_b = "<< hcalaa_b<< " , hcalx_exp = " <<hcalx_exp<<" , hcalx_exp_p = " <<hcalx_exp_p<<endl;
   //cout << "nsigma_bottom_n= "<<nsigma_bottom_n<<" nsigma_bottom_p, = "<<nsigma_bottom_p << " nsigma_top_n = "<<nsigma_top_n<<" , nsigma_top_p = "<<nsigma_top_p<<endl;

  return std::min({nsigma_bottom_n, nsigma_bottom_p ,nsigma_top_n,nsigma_top_p});
}



// END FUNCTIONS
