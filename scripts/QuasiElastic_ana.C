//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Provakar Datta
//   Modified by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified July 7, 2023
//
//
//   The purpose of this script is to take a configuraiton
//   file for some SBS experiment and to produced some
//   analyzed output with cuts on good elastic events.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
#include <vector>
#include <iostream>

#include "TCut.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TChain.h"
#include "TVector3.h"
#include "TStopwatch.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"

#include "../include/gen-ana.h"


DBparse::DBInfo DBInfo;

// Load database files
void getDB(TString cfg){

  cout<<"Attempting to load DB File"<<endl;
  cout<<"---------------------------------------------------------------"<<endl;

   vector<DBparse::DBrequest> request = {
     {"He3 Polarization","He3 target polarization", 1},
     {"Beam Polarization","Beam Polarization values",1},
     {"Helicity Quality","Helicity readback good? (0/1 = bad/good)",1},
     {"Moller Quality","Moller measurements known? (0/1 = no/yes)",1},
     {"Field Measurement","Magnetic field direction",1}
  };

  DBInfo.cfg = cfg;
  DBInfo.var_req = request;

  DB_load(DBInfo);

  cout<<"---------------------------------------------------------------"<<endl;

}


int QuasiElastic_ana(const std::string configfilename, std::string filebase="../outfiles/QE_data")
{

  string configdir = "../config/";

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading input config file ---------------------------------------
  Utilities::KinConf kin_info = Utilities::LoadKinConfig(configdir + configfilename,1);

  getDB(kin_info.conf);

  // parsing trees
  TChain *C = LoadRawRootFiles(kin_info, 1);
  TChain *E = LoadRawRootFiles_E(kin_info, 1);

  // seting up the desired SBS configuration
  TString conf = kin_info.conf;
  int sbsmag = kin_info.sbsmag;
  SBSconfig sbsconf(conf, sbsmag);
  sbsconf.Print();

  // Choosing the model of calculation
  // model 0 => uses reconstructed p as independent variable
  // model 1 => uses reconstructed angles as independent variable
  // model 2 => uses 4-vector calculation
  int model = kin_info.model;
  if (model == 0) std::cout << "Using model 0 [recon. p as indep. var.] for analysis.." << std::endl;
  else if (model == 1) std::cout << "Using model 1 [recon. angle as indep. var.] for analysis.." << std::endl;
  else if (model == 2) std::cout << "Using model 2 [4-vector calculation] for analysis.." << std::endl;
  else { std::cerr << "Enter a valid model number! **!**" << std::endl; throw; }

  // choosing nucleon type
  TString Ntype = kin_info.Ntype;

  double GEMpitch = 10.0*TMath::DegToRad();
  TVector3 GEMzaxis(-sin(GEMpitch),0,cos(GEMpitch));
  TVector3 GEMyaxis(0,1,0);
  TVector3 GEMxaxis = (GEMyaxis.Cross(GEMzaxis)).Unit();

  // setting up ROOT tree branch addresses ---------------------------------------

  // setting up global cuts
  TCut globalcut = kin_info.globalcut;
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);

  int maxNtr=1000;
  C->SetBranchStatus("*",0);
  // beam energy - Probably we should take an average over 100 events
  // double HALLA_p;
  // setrootvar::setbranch(C, "HALLA_p", "", &HALLA_p);

  // taking beam energy from the epics
  Long64_t evnum;
  double HALLA_p;

  E->SetBranchAddress("evnum",&evnum);
  E->SetBranchAddress("HALLA_p",&HALLA_p);

  std::vector<double> epics_evnum;
  std::vector<double> epics_HALLA_p;
  
  epics_evnum.reserve(E->GetEntries());
  epics_HALLA_p.reserve(E->GetEntries());

  for(int i =0;i<E->GetEntries();++i)
  {
    E->GetEntry(i);
    epics_evnum.push_back(static_cast<double>(evnum));
    epics_HALLA_p.push_back(HALLA_p);
  }

  // Some CODA event information
  double evtime;   setrootvar::setbranch(C,"g","evtime",&evtime);
  double evnum_T;  setrootvar::setbranch(C,"g","evnum",&evnum_T);

  // bbcal sh clus var
  double eSH,xSH,ySH,atimeSH;
  int idSH;
  std::vector<std::string> bbcalclvar = {"e","x","y","atimeblk","idblk"};
  std::vector<void*> bbcalclvar_mem = {&eSH,&xSH,&ySH,&atimeSH,&idSH};
  setrootvar::setbranch(C,"bb.sh",bbcalclvar,bbcalclvar_mem);

  // bbcal ps clus var
  double ePS,xPS,yPS,atimePS;
  int idPS;
  std::vector<std::string> bbcalpsclvar = {"e","x","y","atimeblk","idblk"};
  std::vector<void*> bbcalpsclvar_mem = {&ePS,&xPS,&yPS,&atimePS,&idPS};
  setrootvar::setbranch(C,"bb.ps",bbcalpsclvar,bbcalpsclvar_mem);

  int maxhcal = 1000;

  // hcal clus var
  double eHCAL[maxhcal], xHCAL[maxhcal], yHCAL[maxhcal], rblkHCAL[maxhcal], cblkHCAL[maxhcal], idblkHCAL[maxhcal],tdctimeHCAL[maxhcal],atimeHCAL[maxhcal];
  std::vector<std::string> hcalclvar = {"e","x","y","rowblk","colblk","idblk","tdctimeblk","atimeblk"};
  std::vector<void*> hcalclvar_mem = {&eHCAL,&xHCAL,&yHCAL,&rblkHCAL,&cblkHCAL,&idblkHCAL,&tdctimeHCAL,&atimeHCAL};
  setrootvar::setbranch(C, "sbs.hcal", hcalclvar, hcalclvar_mem);
  double hcal_blk_e[maxhcal], hcal_blk_id[maxhcal];
  setrootvar::setbranch(C, "sbs.hcal.clus_blk", "e", &hcal_blk_e);
  setrootvar::setbranch(C, "sbs.hcal.clus_blk", "id", &hcal_blk_id);

  // grinch var
  double grinch_track[maxhcal], grinch_clus_size[maxhcal];
  std::vector<std::string> grinchvar = {"trackindex","size"};
  std::vector<void*> grinchvar_mem = {&grinch_track,&grinch_clus_size};
  setrootvar::setbranch(C, "bb.grinch_tdc.clus", grinchvar, grinchvar_mem);
  double grinch_pmt[maxhcal], grinch_time[maxhcal], grinch_time_corr[maxhcal];
  setrootvar::setbranch(C, "bb.grinch_tdc.hit", "pmtnum", &grinch_pmt);
  setrootvar::setbranch(C, "bb.grinch_tdc.hit", "time", &grinch_time);
  setrootvar::setbranch(C, "bb.grinch_tdc.hit", "time_corr", &grinch_time_corr);

  // hodoscope
  const int maxClus = 1000;
  double hodo_time[maxClus], hodo_id[maxClus];
  int nhodo_clus;
  setrootvar::setbranch(C,"bb.hodotdc.clus.bar.tdc","meantime",&hodo_time);
  setrootvar::setbranch(C,"bb.hodotdc.clus.bar.tdc","id",&hodo_id);
  setrootvar::setbranch(C,"bb.hodotdc.clus","tfinal",&hodo_time);
  setrootvar::setbranch(C,"Ndata.bb.hodotdc.clus.bar.tdc","meantime",&nhodo_clus);

  const int maxTracks = 1000;

  // GEMs
  double gem_nhits[maxTracks], gem_ngood[maxTracks], gem_chi2ndf[maxTracks];
  setrootvar::setbranch(C,"bb.gem.track","nhits",&gem_nhits);
  setrootvar::setbranch(C,"bb.gem.track","ngoodhits",&gem_ngood);
  setrootvar::setbranch(C,"bb.gem.track","chi2ndf",&gem_chi2ndf);

  double sbs_gem_nhits[maxTracks];
  int nsbs_gem_nhits;

  double bb_gem_nhits[maxTracks];
  double bb_gem_track_x[maxTracks];
  double bb_gem_track_y[maxTracks];
  double bb_gem_track_xp[maxTracks];
  double bb_gem_track_yp[maxTracks];

  int nbb_gem_nhits;

  setrootvar::setbranch(C,"sbs.gem.track","nhits",&sbs_gem_nhits);
  setrootvar::setbranch(C,"bb.gem.track","nhits",&bb_gem_nhits);  
  setrootvar::setbranch(C,"bb.gem.track","x",&bb_gem_track_x);  
  setrootvar::setbranch(C,"bb.gem.track","y",&bb_gem_track_y);  
  setrootvar::setbranch(C,"bb.gem.track","xp",&bb_gem_track_xp);  
  setrootvar::setbranch(C,"bb.gem.track","yp",&bb_gem_track_yp);  
  //setrootvar::setbranch(C,"Ndata.sbs.gem.track","nhits",&nsbs_gem_nhits);

  // track var
  double ntrack, p[maxNtr],px[maxNtr],py[maxNtr],pz[maxNtr],xTr[maxNtr],yTr[maxNtr],thTr[maxNtr],phTr[maxNtr];
  double vx[maxNtr],vy[maxNtr],vz[maxNtr];
  double xtgt[maxNtr],ytgt[maxNtr],thtgt[maxNtr],phtgt[maxNtr];
  double xfp[maxNtr],yfp[maxNtr],thfp[maxNtr],phfp[maxNtr];
  std::vector<std::string> trvar = {"n","p","px","py","pz","vx","vy","vz","tg_x","tg_y","tg_th","tg_ph","r_x","r_y","r_th","r_ph","x","y","ph","th"};
  std::vector<void*> trvar_mem = {&ntrack,&p,&px,&py,&pz,&vx,&vy,&vz,&xtgt,&ytgt,&thtgt,&phtgt,&xfp,&yfp,&thfp,&phfp,&xTr,&yTr,&phTr,&thTr};
  setrootvar::setbranch(C,"bb.tr",trvar,trvar_mem);

  //sbs track var
  double ntrack_sbs, chi2_sbs[maxNtr], p_sbs[maxNtr],px_sbs[maxNtr],py_sbs[maxNtr],pz_sbs[maxNtr],xTr_sbs[maxNtr],yTr_sbs[maxNtr],thTr_sbs[maxNtr],phTr_sbs[maxNtr];
  double vx_sbs[maxNtr],vy_sbs[maxNtr],vz_sbs[maxNtr];
  double ytgt_sbs[maxNtr],thtgt_sbs[maxNtr],phtgt_sbs[maxNtr];
  double xfp_sbs[maxNtr],yfp_sbs[maxNtr],thfp_sbs[maxNtr],phfp_sbs[maxNtr];
  std::vector<std::string> sbstrvar ={"n","chi2","p","px","py","pz","vx","vy","vz","tg_y","tg_th","tg_ph","r_x","r_y","r_th","r_ph","x","y","th","ph"};
  std::vector<void*> sbstrvar_mem = {&ntrack_sbs,&chi2_sbs,&p_sbs,&px_sbs,&py_sbs,&pz_sbs,&vx_sbs,&vy_sbs,&vz_sbs,&ytgt_sbs,&thtgt_sbs,&phtgt_sbs,&xfp_sbs,&yfp_sbs,&thfp_sbs,&phfp_sbs,&xTr_sbs,&yTr_sbs,&thTr_sbs,&phTr_sbs};
  //setrootvar::setbranch(C,"sbs.tr",sbstrvar,sbstrvar_mem);

  // tdctrig variable (N/A for simulation)
  int tdcElemN;
  int Ndata_bb_tdctrig_tdc;
  double tdcTrig[maxNtr], tdcElem[maxNtr];

  std::vector<std::string> tdcvar = {"tdcelemID","tdcelemID","tdc"};
  std::vector<void*> tdcvar_mem = {&tdcElem,&tdcElemN,&tdcTrig};
  setrootvar::setbranch(C,"bb.tdctrig",tdcvar,tdcvar_mem,1);
  setrootvar::setbranch(C,"Ndata.bb.tdctrig","tdc",&Ndata_bb_tdctrig_tdc);

  //Beam helicity variables
  double helicity;
  setrootvar::setbranch(C,"scalhel","hel",&helicity);

  //IHWP State
  double IHWP;
  setrootvar::setbranch(C,"IGL1I00OD16_16","",&IHWP);

  //trigbits
  double trigbits;
  setrootvar::setbranch(C,"g","trigbits",&trigbits);

  //BPMA
  double BPMAx, BPMAy;
  std::vector<std::string> BPMAvar = {"x","y"};
  std::vector<void*> BPMAvar_mem = {&BPMAx,&BPMAy};
  setrootvar::setbranch(C,"Lrb.BPMA",BPMAvar,BPMAvar_mem);

  //Raster
  double rawcurx, rawcury;
  std::vector<std::string> Rastervar = {"x","y"};
  std::vector<void*> Rastervar_mem = {&rawcurx,&rawcury};
  setrootvar::setbranch(C,"Lrb.Raster.rawcur",Rastervar,Rastervar_mem);

  //Raster 2
  double rawcur2x, rawcur2y;
  std::vector<std::string> Raster2var = {"x","y"};
  std::vector<void*> Raster2var_mem = {&rawcur2x,&rawcur2y};
  setrootvar::setbranch(C,"Lrb.Raster2.rawcur",Raster2var,Raster2var_mem);

  // turning on the remaining branches we use for the globalcut
  C->SetBranchStatus("bb.gem.track.nhits", 1);
  C->SetBranchStatus("bb.etot_over_p", 1);
  C->SetBranchStatus("sbs.hcal.nclus", 1);

  // defining the outputfile
  TString outFile = Form("%s_" + sbsconf.GetSBSconf() + "_sbs%dp_nucleon_%s_model%d.root",
			 filebase.c_str(),  sbsconf.GetSBSmag(), Ntype.Data(), model);
  TFile *fout = new TFile(outFile.Data(), "RECREATE");

  // defining histograms
  TH1F *h_W = Utilities::TH1FhW("h_W");
  TH1F *h_W_cut = Utilities::TH1FhW("h_W_cut");
  TH1F *h_W_acut = Utilities::TH1FhW("h_W_acut");
  TH1D *h_dpel = new TH1D("h_dpel",";p/p_{elastic}(#theta)-1;",100,-0.3,0.3);

  TH1F *h_Q2 = Utilities::TH1FhQ2("h_Q2", conf);
  vector<double> hdx_lim = kin_info.hdx_lim;
  vector<double> hdy_lim = kin_info.hdy_lim;
  TH1F *h_dxHCAL = new TH1F("h_dxHCAL","; x_{HCAL} - x_{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);
  TH1F *h_dyHCAL = new TH1F("h_dyHCAL","; y_{HCAL} - y_{exp} (m);",int(hdy_lim[0]),hdy_lim[1],hdy_lim[2]);
  TH1F *h_coin_time = new TH1F("h_coin_time", "Coincidence time (ns)", 200, 380, 660);

  TH2F *h2_rcHCAL = Utilities::TH2FHCALface_rc("h2_rcHCAL");
  TH2F *h2_dxdyHCAL = Utilities::TH2FdxdyHCAL("h2_dxdyHCAL");
  TH2F *h2_xyHCAL_p = Utilities::TH2FHCALface_xy_data("h2_xyHCAL_p");
  TH2F *h2_xyHCAL_n = Utilities::TH2FHCALface_xy_data("h2_xyHCAL_n");

  // Defining interesting ROOT tree branches
  TTree *Tout = new TTree("Tout", "");
  int T_runnum;         Tout->Branch("runnum", &T_runnum, "runnum/I");
  TDatime T_datetime;   Tout->Branch("datetime", "TDatime", &T_datetime);
  //cuts
  bool WCut;            Tout->Branch("WCut", &WCut, "WCut/B");
  bool pCut;            Tout->Branch("pCut", &pCut, "pCut/B");
  bool nCut;            Tout->Branch("nCut", &nCut, "nCut/B");
  bool fiduCut;         Tout->Branch("fiduCut", &fiduCut, "fiduCut/B");
  bool coinCut;         Tout->Branch("coinCut", &coinCut, "coinCut/B");
  //kine
  double T_ebeam;       Tout->Branch("ebeam", &T_ebeam, "ebeam/D");
  double T_nu;          Tout->Branch("e.kine.nu", &T_nu, "e.kine.nu/D");
  double T_Q2;          Tout->Branch("e.kine.Q2", &T_Q2, "e.kine.Q2/D");
  double T_W2;          Tout->Branch("e.kine.W2", &T_W2, "e.kine.W2/D");
  double T_dpel;        Tout->Branch("dpel", &T_dpel, "dpel/D");
  double T_ephi;        Tout->Branch("ephi", &T_ephi, "ephi/D");
  double T_etheta;      Tout->Branch("etheta", &T_etheta, "etheta/D");
  double T_pcentral;    Tout->Branch("pcentral", &T_pcentral, "pcentral/D");
  double T_pN_expect; 	Tout->Branch("pN_expect", &T_pN_expect, "pN_expect/D");
  double T_ptheta_cal;	Tout->Branch("ptheta_cal",&T_ptheta_cal,"ptheta_cal/D");
  double T_ptheta; 	    Tout->Branch("ptheta",&T_ptheta,"ptheta/D");
  double T_pphi_cal;    Tout->Branch("pphi_cal",&T_pphi_cal,"pphi_cal/D");
  double T_pphi;        Tout->Branch("pphi",&T_pphi,"pphi/D");
  //track
  double T_vz;          Tout->Branch("bb.tr.vz", &T_vz, "bb.tr.vz/D");
  double T_vx;          Tout->Branch("bb.tr.vx", &T_vx, "bb.tr.vx/D");
  double T_vy;          Tout->Branch("bb.tr.vy", &T_vy, "bb.tr.vy/D");
  double T_xTr;         Tout->Branch("bb.tr.x", &T_xTr, "bb.tr.x/D");
  double T_yTr;         Tout->Branch("bb.tr.y", &T_yTr, "bb.tr.y/D");
  double T_thTr;        Tout->Branch("bb.tr.th", &T_thTr, "bb.tr.th/D");
  double T_phTr;        Tout->Branch("bb.tr.ph", &T_phTr, "bb.tr.ph/D");
  double T_xtgt;        Tout->Branch("xtgt", &T_xtgt, "xtgt/D");
  double T_ytgt;        Tout->Branch("ytgt", &T_ytgt, "ytgt/D");
  double T_thtgt;       Tout->Branch("thtgt", &T_thtgt, "thtgt/D");
  double T_phtgt;       Tout->Branch("phtgt", &T_phtgt, "phtgt/D");
  double T_thetabend;   Tout->Branch("thetabend", &T_thetabend, "thetabend/D");
  double T_xfp;         Tout->Branch("bb.tr.r_x", &T_xfp, "bb.tr.r_x/D");
  double T_yfp;         Tout->Branch("bb.tr.r_y", &T_yfp, "bb.tr.r_y/D");
  double T_thfp;        Tout->Branch("bb.tr.r_th", &T_thfp, "bb.tr.r_th/D");
  double T_phfp;        Tout->Branch("bb.tr.r_ph", &T_phfp, "bb.tr.r_ph/D");
  double T_trP;         Tout->Branch("bb.tr.p", &T_trP, "bb.tr.p/D");
  double T_trPx;        Tout->Branch("bb.tr.px", &T_trPx, "bb.tr.px/D");
  double T_trPy;        Tout->Branch("bb.tr.py", &T_trPy, "bb.tr.py/D");
  double T_trPz;        Tout->Branch("bb.tr.pz", &T_trPz, "bb.tr.pz/D");
  double T_trN;         Tout->Branch("bb.tr.n", &T_trN, "bb.tr.n/D");
  //GEMs
  double T_gem_nhits;   Tout->Branch("bb.gem.track.nhits", &T_gem_nhits, "bb.gem.track.nhits/D");
  double T_gem_ngood;   Tout->Branch("bb.gem.track.ngoodhits", &T_gem_ngood, "bb.gem.track.ngoodhits/D");
  double T_gem_chi2ndf; Tout->Branch("bb.gem.track.chi2ndf", &T_gem_chi2ndf, "bb.gem.track.chi2ndf/D");
  double T_gem_x;       Tout->Branch("bb.gem.track.x", &T_gem_x, "bb.gem.track.x/D");
  double T_gem_y;       Tout->Branch("bb.gem.track.y", &T_gem_y, "bb.gem.track.y/D");
  double T_gem_xp;      Tout->Branch("bb.gem.track.xp", &T_gem_xp, "bb.gem.track.xp/D");
  double T_gem_yp;      Tout->Branch("bb.gem.track.yp", &T_gem_yp, "bb.gem.track.yp/D");
  //BBCAL
  double T_ePS;         Tout->Branch("bb.ps.e", &T_ePS, "bb.ps.e/D");
  double T_xPS;         Tout->Branch("bb.ps.x", &T_xPS, "bb.ps.x/D");
  double T_yPS;         Tout->Branch("bb.ps.y", &T_yPS, "bb.ps.y/D");
  double T_eSH;         Tout->Branch("bb.sh.e", &T_eSH, "bb.sh.e/D");
  double T_xSH;         Tout->Branch("bb.sh.x", &T_xSH, "bb.sh.x/D");
  double T_ySH;         Tout->Branch("bb.sh.y", &T_ySH, "bb.sh.y/D");
  //HCAL
  double T_eHCAL;           Tout->Branch("sbs.hcal.e", &T_eHCAL, "sbs.hcal.e/D");
  double T_xHCAL;           Tout->Branch("sbs.hcal.x", &T_xHCAL, "sbs.hcal.x/D");
  double T_yHCAL;           Tout->Branch("sbs.hcal.y", &T_yHCAL, "sbs.hcal.y/D");
  double T_idblkHCAL;       Tout->Branch("sbs.hcal.idblk", &T_idblkHCAL, "sbs.hcal.idblk/D");
  double T_hcal_primblk_e;  Tout->Branch("sbs.hcal.primary_blk.e", &T_hcal_primblk_e, "sbs.hcal.primary_blk.e/D");
  double T_hcal_primblk_id; Tout->Branch("sbs.hcal.primary_blk.id", &T_hcal_primblk_id, "sbs.hcal.primary_blk.id/D");
  double T_hcal_secblk_e;   Tout->Branch("sbs.hcal.secondary_blk.e", &T_hcal_secblk_e, "sbs.hcal.secondary_blk.e/D");
  double T_hcal_secblk_id;  Tout->Branch("sbs.hcal.secondary_blk.id", &T_hcal_secblk_id, "sbs.hcal.secondary_blk.id/D");
  double T_xHCAL_exp;       Tout->Branch("xHCAL_exp", &T_xHCAL_exp, "xHCAL_exp/D");
  double T_yHCAL_exp;       Tout->Branch("yHCAL_exp", &T_yHCAL_exp, "yHCAL_exp/D");
  double T_dx;              Tout->Branch("dx", &T_dx, "dx/D");
  double T_dy;              Tout->Branch("dy", &T_dy, "dy/D");
  double T_theta_pq;        Tout->Branch("theta_pq", &T_theta_pq, "theta_pq/D");
  //GRINCH
  double T_grinch_track;        Tout->Branch("bb.grinch_tdc.clus.trackindex", &T_grinch_track, "bb.grinch_tdc.clus.trackindex/D");
  double T_grinch_clus_size;    Tout->Branch("bb.grinch_tdc.clus.size", &T_grinch_clus_size, "bb.grinch_tdc.clus.size/D");
  double T_grinch_pmt;          Tout->Branch("bb.grinch_tdc.hit.pmtnum", &T_grinch_pmt, "bb.grinch_tdc.hit.pmtnum/D");
  double T_grinch_time;         Tout->Branch("bb.grinch_tdc.hit.time", &T_grinch_time, "bb.grinch_tdc.hit.time/D");
  double T_grinch_time_corr;    Tout->Branch("bb.grinch_tdc.hit.time_corr", &T_grinch_time_corr, "bb.grinch_tdc.hit.time_corr/D");

  //Timing Information
  double T_coin_time;            Tout->Branch("adc.coin", &T_coin_time, "adc.coin/D");
  double T_hcal_time;            Tout->Branch("sbs.hcal.atimeblk", &T_hcal_time, "sbs.hcal.atimeblk/D");
  double T_bbcal_time;           Tout->Branch("bb.sh.atimeblk", &T_bbcal_time, "bb.sh.atimeblk/D");
  double T_ps_time;              Tout->Branch("bb.ps.atimeblk", &T_ps_time, "bb.ps.atimeblk/D");
  int T_nhodo_clus;              Tout->Branch("nhodo_clus", &T_nhodo_clus, "nhodo_clus/I");
  double T_hodo_id[maxClus];     Tout->Branch("bb.hodotdc.clus.id", &T_hodo_id, "bb.hodotdc.clus.id[nhodo_clus]/D");
  double T_hodo_time[maxClus];   Tout->Branch("bb.hodotdc.clus.tmean", &T_hodo_time, "bb.hodotdc.clus.tmean[nhodo_clus]/D");

  //BPM and Raster information
  double T_BPMAx;     Tout->Branch("BPMAx", &T_BPMAx, "BPMAx/D");
  double T_BPMAy;     Tout->Branch("BPMAy", &T_BPMAy, "BPMAy/D");
  double T_rawcurx;  Tout->Branch("Rasterx", &T_rawcurx, "Rasterx/D");
  double T_rawcury;  Tout->Branch("Rastery", &T_rawcury, "Rastery/D");
  double T_rawcur2x;  Tout->Branch("Raster2x", &T_rawcur2x, "Raster2x/D");
  double T_rawcur2y;  Tout->Branch("Raster2y", &T_rawcur2y, "Raster2y/D");

  //trigbits
  double T_trigbits; Tout->Branch("trigbits", &T_trigbits, "trigbits/D");

  //Beam/Target information
  int T_helicity;        Tout->Branch("helicity", &T_helicity, "helicity/I");
  int T_IHWP;            Tout->Branch("IHWP", &T_IHWP, "IHWP/I");
  double T_He3Pol;       Tout->Branch("He3Pol", &T_He3Pol, "He3Pol/D");
  double T_err_He3Pol;   Tout->Branch("err_He3Pol", &T_err_He3Pol, "err_He3Pol/D");

  // Do the energy loss calculation here ...........

  // HCAL cut definitions
  double sbs_kick = kin_info.sbs_kick;
  vector<double> dx_p = kin_info.dx_p;
  vector<double> dy_p = kin_info.dy_p;
  double Nsigma_cut_dx_p = kin_info.Nsigma_dx_p;
  double Nsigma_cut_dy_p = kin_info.Nsigma_dy_p;
  vector<double> dx_n = kin_info.dx_n;
  vector<double> dy_n = kin_info.dy_n;
  double Nsigma_cut_dx_n = kin_info.Nsigma_dx_n;
  double Nsigma_cut_dy_n = kin_info.Nsigma_dy_n;
  double coin_min = -2.0; //fix later
  double coin_max = 2.0; //fix later
  vector<double> hcal_active_area = cut::hcal_active_area_data(); // Exc. 1 blk from all 4 sides
  vector<double> hcal_safety_margin = cut::hcal_safety_margin(dx_p[1], dx_n[1], dy_p[1], hcal_active_area);

  // elastic cut limits
  double W2min = kin_info.W2min;
  double W2max = kin_info.W2max;

  // costruct axes of HCAL CoS in Hall CoS
  double hcal_voffset = kin_info.hcal_voffset;
  double hcal_hoffset = kin_info.hcal_hoffset;
  vector<TVector3> HCAL_axes; kine::SetHCALaxes(sbsconf.GetSBStheta_rad(), HCAL_axes);
  TVector3 HCAL_origin = sbsconf.GetHCALdist()*HCAL_axes[2] + hcal_voffset*HCAL_axes[0] + hcal_hoffset*HCAL_axes[1];

  // looping through the tree ---------------------------------------
  std::cout << std::endl;
  long nevent = 0, nevents = C->GetEntries();
  int treenum = 0, currenttreenum = 0, currentrunnum = 0;
  int IHWP_run = -100;
  time_t run_time_unix;

  cout<<"Processing "<<nevents<<" events"<<endl;

  while (C->GetEntry(nevent++)) {

    // print progress
    if( nevent % 1000 == 0 ) std::cout << nevent*100.0/nevents << "% \r";
    std::cout.flush();

    // apply global cuts efficiently (AJRP method)
    currenttreenum = C->GetTreeNumber();
    if (nevent == 1 || currenttreenum != treenum) {
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();

      //Get the run number
      string s = C->GetFile()->GetName();
      int start = s.find("_stream0");
      start -= 4;
      int end = start + 4;
      T_runnum = stoi(s.substr(start,end - start));

      //Get the time this run started
      auto* Run_Data = C->GetFile()->Get<THaRunBase>("Run_Data");
      TDatime run_time = Run_Data->GetDate();
      run_time.Set(run_time.GetYear(),run_time.GetMonth(),run_time.GetDay(),run_time.GetHour(),run_time.GetMinute(),0);
      run_time_unix = run_time.Convert();

      // We must loop over a small subset to get the correct IHWP state for the entire run
      //First we do some stuff to max sure we dont reach the file limit
      int file_nevents = C->GetTree()->GetEntries();
      int max_events = 8000;
      if(file_nevents < max_events) max_events = file_nevents - 100;

      int start_event = nevent;
      while (C->GetEntry(start_event++) && start_event < nevent + max_events){
	    if(IHWP == 1 || IHWP == -1) IHWP_run = IHWP;
      }
      C->GetEntry(nevent - 1);

    }

    bool passedgCut = GlobalCut->EvalInstance(0) != 0;

    if (!passedgCut) continue;

    double bbcal_trig_time=0., hcal_trig_time=0.;
    for(int ihit=0; ihit<tdcElemN; ihit++){
      if(tdcElem[ihit]==5) bbcal_trig_time=tdcTrig[ihit];
      if(tdcElem[ihit]==0) hcal_trig_time=tdcTrig[ihit];
    }

    //Calculate absolute time of this event
    double time_interval = 4;  //in ns
    int time_rel = evtime*time_interval*1e-9 / 60; //in min, rounded

    TDatime time_abs(run_time_unix + time_rel * 60); //Add the relative minutes

    auto it = DBInfo.He3Pol.find(time_abs); // Get Polarization from the table
    if(it == DBInfo.He3Pol.end())
    {
      T_He3Pol = -1;
      T_err_He3Pol = -1;
    }
    else
    {
      T_He3Pol = it->second.first;
      T_err_He3Pol = it->second.second;
    }

    T_datetime = time_abs;

    //Timing Information
    T_hcal_time = atimeHCAL[0];
    T_bbcal_time = atimeSH;
    T_ps_time = atimePS;

    T_nhodo_clus = nhodo_clus;
    for(int iclus = 0; iclus < nhodo_clus; iclus++)
    {
        T_hodo_time[iclus] = hodo_time[iclus];
        T_hodo_id[iclus] = hodo_id[iclus];
    }

    double coin_time = atimeHCAL[0] - atimeSH;
    T_coin_time = coin_time;
    h_coin_time->Fill(coin_time);

    coinCut = coin_time > coin_min && coin_time < coin_max;

    //Grinch info
    T_grinch_track = grinch_track[0];
    T_grinch_clus_size = grinch_clus_size[0];
    T_grinch_pmt = grinch_pmt[0];
    T_grinch_time = grinch_time[0];
    T_grinch_time_corr = grinch_time_corr[0];

    //trigbits
    T_trigbits = trigbits;

    //Beam helicity information
    T_helicity = helicity;
    T_IHWP = IHWP_run;

    //BPMs and Rasters
    T_BPMAx = BPMAx;
    T_BPMAy = BPMAy;
    T_rawcurx = rawcurx;
    T_rawcury = rawcury;
    T_rawcur2x = rawcur2x;
    T_rawcur2y = rawcur2y;

     // kinematic parameters
     //double ebeam = sbsconf.GetEbeam();       // Expected beam energy (GeV) [Get it from EPICS, eventually]

     size_t j   = evnum_T; //up(evnum_T);        // first EPICS evnum > current event
     size_t idx = (j==0) ? 0 : j-1;    // take the latest EPICS sample <= event evnum
    
     double pMeV = 0.0;

     if (epics_HALLA_p.size()>0)
     {
       if(idx >= epics_HALLA_p.size()) idx = epics_HALLA_p.size()-1;
       pMeV = epics_HALLA_p[idx];
     }
     else
     {
       pMeV= sbsconf.GetEbeam() * 1000; 
     }

     double ebeam = pMeV/1000.0;

     double ebeam_corr = ebeam; //- MeanEloss;
     double precon = p[0]; //+ MeanEloss_outgoing

    // constructing the 4 vectors
    // Reaction    : e + e' -> N + N'
    // Conservation: Pe + Peprime = PN + PNprime
    TVector3 vertex(0, 0, vz[0]);
    TLorentzVector Pe(0,0,ebeam_corr,ebeam_corr);   // incoming e-
    TLorentzVector Peprime(px[0] * (precon/p[0]),   // scattered e-
			   py[0] * (precon/p[0]),
			   pz[0] * (precon/p[0]),
			   precon);
    TLorentzVector PN;                              // target nucleon [Ntype ??]
    kine::SetPN(Ntype, PN);

    double etheta = kine::etheta(Peprime);
    double ephi = kine::ephi(Peprime);
    double pcentral = kine::pcentral(ebeam_corr, etheta, Ntype.Data());

    double nu = 0.;                   // energy of the virtual photon
    double pN_expect = 0.;            // expected recoil nucleon momentum
    double thetaN_expect = 0.;        // expected recoil nucleon theta
    double phiN_expect = ephi + constant::pi;
    // Different modes of calculation. Goal is to achieve the best resolution
    // model 0 = uses reconstructed p as independent variable
    // model 1 = uses reconstructed angles as independent variable
    // model 2 = uses 4-vector calculation
    TVector3 pNhat;                   // 3-momentum of the recoil nucleon (Unit)
    double Q2recon, W2recon;
    if (model == 0) {
      nu = Pe.E() - Peprime.E();
      pN_expect = kine::pN_expect(nu, Ntype.Data());
      thetaN_expect = acos((Pe.E() - Peprime.Pz()) / pN_expect);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);
      Q2recon = kine::Q2(Pe.E(), Peprime.E(), etheta);
      W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype.Data());
    } else if (model == 1) {
      nu = Pe.E() - pcentral;
      pN_expect = kine::pN_expect(nu, Ntype.Data());
      thetaN_expect = acos((Pe.E() - pcentral*cos(etheta)) / pN_expect);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);
      Q2recon = kine::Q2(Pe.E(), Peprime.E(), etheta);
      W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype.Data());
    } else if (model == 2) {
      TLorentzVector q = Pe - Peprime; // 4-momentum of virtual photon
      nu = q.E();
      pNhat = q.Vect().Unit();
      pN_expect = q.P();//spatial magnitude
      Q2recon = -q.M2();
      W2recon = (PN + q).M2();
      thetaN_expect = q.Theta();
    }
    h_Q2->Fill(Q2recon);
    double Wrecon = sqrt(max(0., W2recon));
    double dpel = Peprime.E()/pcentral - 1.0;
    h_dpel->Fill(dpel);

    TVector3 enhat_tgt( thtgt[0], phtgt[0], 1.0 );
    enhat_tgt = enhat_tgt.Unit();

    TVector3 enhat_fp( thfp[0], phfp[0], 1.0 );
    enhat_fp = enhat_fp.Unit();

    TVector3 enhat_fp_rot = enhat_fp.X() * GEMxaxis + enhat_fp.Y() * GEMyaxis + enhat_fp.Z() * GEMzaxis;

    double thetabend = acos( enhat_fp_rot.Dot( enhat_tgt ) );
    double pphi = atan(py_sbs[0]/px_sbs[0]) ;
    double ptheta = acos(pz_sbs[0]/p_sbs[0]);

    //T_pN_expect = pN_expect;

    T_ebeam = Pe.E();

    T_nu = nu;
    T_Q2 = Q2recon;
    T_W2 = W2recon;
    T_dpel = dpel;
    T_ephi = ephi;
    T_etheta = etheta;
    T_pcentral = pcentral;
    T_pN_expect = pN_expect;
    T_ptheta_cal = thetaN_expect;
    T_pphi_cal = phiN_expect;
    T_ptheta = ptheta;
    T_pphi = pphi;

    T_vz = vz[0];
    T_vx = vx[0];
    T_vy = vy[0];
    T_xtgt = xtgt[0];
    T_ytgt = ytgt[0];
    T_thtgt = thtgt[0];
    T_phtgt = phtgt[0];
    T_xTr = xTr[0];
    T_yTr = yTr[0];
    T_thTr = thTr[0];
    T_phTr = phTr[0];
    T_thetabend = thetabend;
    T_xfp = xfp[0];
    T_yfp = yfp[0];
    T_thfp = thfp[0];
    T_phfp = phfp[0];
    T_trP = p[0];
    T_trPx = px[0];
    T_trPy = py[0];
    T_trPz = pz[0];
    T_trN = ntrack;

    T_gem_nhits = gem_nhits[0];
    T_gem_ngood = gem_ngood[0];
    T_gem_chi2ndf = gem_chi2ndf[0];
    T_gem_x = bb_gem_track_x[0];
    T_gem_y = bb_gem_track_y[0];
    T_gem_xp = bb_gem_track_xp[0];
    T_gem_yp = bb_gem_track_yp[0];

    T_ePS = ePS;
    T_xPS = xPS;
    T_yPS = yPS;
    T_eSH = eSH;
    T_xSH = xSH;
    T_ySH = ySH;

    T_eHCAL = eHCAL[0];
    T_xHCAL = xHCAL[0];
    T_yHCAL = yHCAL[0];
    T_idblkHCAL = idblkHCAL[0];
    T_hcal_primblk_e = hcal_blk_e[0];
    T_hcal_primblk_id = hcal_blk_id[0];
    T_hcal_secblk_e = hcal_blk_e[1];
    T_hcal_secblk_id = hcal_blk_id[1];

    // Expected position of the q vector at HCAL
    vector<double> xyHCAL_exp; // xyHCAL_exp[0] = xHCAL_exp & xyHCAL_exp[1] = yHCAL_exp
    double theta_pq = kine::theta_pq(vertex, HCAL_origin, pNhat, xHCAL[0], yHCAL[0]);
    kine::GetxyHCALexpect(vertex, pNhat, HCAL_origin, HCAL_axes, xyHCAL_exp);
    double dx = xHCAL[0] - xyHCAL_exp[0];
    double dy = yHCAL[0] - xyHCAL_exp[1];

    T_theta_pq = theta_pq;
    T_xHCAL_exp = xyHCAL_exp[0];
    T_yHCAL_exp = xyHCAL_exp[1];
    T_dx = dx;
    T_dy = dy;

    // HCAL active area and safety margin cuts [Fiducial region]
    bool AR_cut = cut::inHCAL_activeA(xHCAL[0], yHCAL[0], hcal_active_area);
    bool FR_cut = cut::inHCAL_fiducial(xyHCAL_exp[0], xyHCAL_exp[1], sbs_kick, hcal_safety_margin);

    // HCAL cuts
    pCut = pow((dx-dx_p[0]) / (dx_p[1]*Nsigma_cut_dx_p), 2) + pow((dy-dy_p[0]) / (dy_p[1]*Nsigma_cut_dy_p), 2) <= 1.;
    nCut = pow((dx-dx_n[0]) / (dx_n[1]*Nsigma_cut_dx_n), 2) + pow((dy-dy_n[0]) / (dy_n[1]*Nsigma_cut_dy_n), 2) <= 1.;

    fiduCut = AR_cut && FR_cut;
    WCut = W2recon > W2min && W2recon < W2max;

    // W cut and coincidence cut
    if (WCut && coinCut) {
      // fiducial cut
      //if (fiduCut) { //No fiducial cut for GEN?
	h_dxHCAL->Fill(dx);
	h_dyHCAL->Fill(dy);
	h2_rcHCAL->Fill(cblkHCAL[0], rblkHCAL[0]);
	h2_dxdyHCAL->Fill(dy, dx);

	if (pCut) h2_xyHCAL_p->Fill(xyHCAL_exp[1], xyHCAL_exp[0] - sbs_kick);
	if (nCut) h2_xyHCAL_n->Fill(xyHCAL_exp[1], xyHCAL_exp[0]);
	//}
    }
    h_W->Fill(Wrecon);
    // fiducial cut but no W cut
    //if (fiduCut) { //No fiducial cut for GEN
    //h_W_cut->Fill(Wrecon);
      if (pCut || nCut) {
	h_W_cut->Fill(Wrecon);
      } else {
	h_W_acut->Fill(Wrecon);
      }
      //}


      Tout->Fill();
  } // event loop
  std::cout << std::endl << std::endl;

  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1000);
  c1->Divide(2,2);

  c1->cd(1); h2_dxdyHCAL->Draw("colz");
  TEllipse Ep_p;
  Ep_p.SetFillStyle(0); Ep_p.SetLineColor(2); Ep_p.SetLineWidth(2);
  Ep_p.DrawEllipse(dy_p[0], dx_p[0], Nsigma_cut_dy_p*dy_p[1], Nsigma_cut_dx_p*dx_p[1], 0,360,0);
  TEllipse Ep_n;
  Ep_n.SetFillStyle(0); Ep_n.SetLineColor(3); Ep_n.SetLineWidth(2);
  Ep_n.DrawEllipse(dy_n[0], dx_n[0], Nsigma_cut_dy_n*dy_n[1], Nsigma_cut_dx_n*dx_n[1], 0,360,0);

  c1->cd(2);
  h_W->Draw(); h_W->SetLineColor(1);
  h_W_cut->Draw("same"); h_W_cut->SetLineColor(2);
  h_W_acut->Draw("same");

  c1->cd(3);
  h2_xyHCAL_p->Draw("colz");
  Utilities::DrawArea(hcal_active_area);
  Utilities::DrawArea(hcal_safety_margin,4);

  c1->cd(4);
  h2_xyHCAL_n->Draw("colz");
  Utilities::DrawArea(hcal_active_area);
  Utilities::DrawArea(hcal_safety_margin,4);

  cout << "------" << endl;
  cout << " Output file : " << outFile << endl;
  cout << "------" << endl << endl;

  sw->Stop();
  cout << "CPU time elapsed = " << sw->CpuTime() << " s. Real time = " << sw->RealTime() << " s. " << endl << endl;

  c1->Write();

  fout->Write();
  sw->Delete();

  return 0;
}
