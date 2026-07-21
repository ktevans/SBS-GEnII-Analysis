//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Provakar Datta
//   Modified by Sean Jeffas, sj9ry@virginia.edu
//     Last Modified August 7, 2024
//   Modified by Kate Evans, ktevans@jlab.org
//     Last Modified March 28, 2025
//   Modified by Jack Jackson, cmjackso@jlab.org
//     Last modified Aug 21, 2025
//   Modified by Kate Evans, ktevans@jlab.org
//     Last modified July 21, 2026
//   The purpose of this script is to take a configuraiton
//   file for some SBS experiment and to produced some
//   analyzed output with cuts on good elastic events.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#include <vector>
#include <iostream>
#include <cmath>
#include <cstdint>
#include <sys/resource.h>

// ROOT Headers
#include "TVirtualPad.h"
#include "TBox.h"
#include "TF1.h"
#include "TLine.h"
#include "TClassTable.h"
#include "THaRunBase.h"
#include "TCut.h"
#include "TChain.h"
#include "TVector3.h"
#include "TStopwatch.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TSystem.h"
#include "TDatime.h"
#include "TTree.h"
#include "TString.h"

#include "../include/gen-ana.h"

// Database handling from QE_ana structure
DBparse::DBInfo DBInfo;
void getDB(TString cfg){
    std::vector<DBparse::DBrequest> request = {
        {"He3 Polarization","He3 target polarization", 1},
        {"Beam Polarization","Beam Polarization values",1},
        {"Helicity Quality","Helicity readback good? (0/1 = bad/good)",1},
        {"Moller Quality","Moller measurements known? (0/1 = no/yes)",1},
        {"Field Measurement","Magnetic field direction",1}
    };
    DBInfo.cfg = cfg;
    DBInfo.var_req = request;
    DB_load(DBInfo);
}

static void EnsureSingleFileOutput() {
  rlimit rl{};
  if (getrlimit(RLIMIT_FSIZE, &rl) == 0) {
    rl.rlim_cur = RLIM_INFINITY;
    rl.rlim_max = RLIM_INFINITY;
    setrlimit(RLIMIT_FSIZE, &rl);
  }
  TTree::SetMaxTreeSize(100LL * 1024 * 1024 * 1024);
}

int QE_Parse(const std::string configfilename, std::string filebase="../outfiles/QE_Parse_pass3")
{
  if (!gClassTable->GetID("THaRunBase")) {
    gSystem->Load("libanalyzer");
  }

  std::string configdir = "../../config/";
  gErrorIgnoreLevel = kError;
  TStopwatch sw; sw.Start();

  Utilities::KinConf kin_info = Utilities::LoadKinConfig(configdir + configfilename, 1);
  getDB(kin_info.conf);

  TChain *C = LoadRawRootFiles(kin_info, 1);
  TString conf = kin_info.conf;
  SBSconfig sbsconf(conf, kin_info.sbsmag);
  sbsconf.Print();

  TCut globalcut = kin_info.globalcut;
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);

  const int maxNtr = 1000;
  const int maxhcal = 1000;
  const int kMaxHodoHits = 200;
  static const int grinchmax = 1000;

  // --- BRANCH STATUS WHIETLIST (Prevents Zeros) ---
  C->SetBranchStatus("*", 0);
  C->SetBranchStatus("g.*", 1);
  C->SetBranchStatus("bb.sh.*", 1);
  C->SetBranchStatus("bb.ps.*", 1);
  C->SetBranchStatus("bb.etot_over_p", 1);
  C->SetBranchStatus("sbs.hcal.*", 1);
  C->SetBranchStatus("bb.grinch_tdc.*", 1);
  C->SetBranchStatus("Ndata.bb.grinch_tdc.hit.pmtnum", 1);
  C->SetBranchStatus("bb.hodotdc.*", 1);
  C->SetBranchStatus("Ndata.bb.hodotdc.clus.bar.tdc.meantime", 1);
  C->SetBranchStatus("bb.tr.*", 1);
  C->SetBranchStatus("bb.tdctrig.*", 1);
  C->SetBranchStatus("scalhel.*", 1);
  C->SetBranchStatus("IGL1I00OD16_16", 1);
  C->SetBranchStatus("Lrb.BPMA.*", 1);
  C->SetBranchStatus("Lrb.Raster.*", 1);
  C->SetBranchStatus("Lrb.Raster2.*", 1);
  C->SetBranchStatus("bb.gem.track.nhits", 1);
  C->SetBranchStatus("bb.gem.track.ngoodhits", 1);
  C->SetBranchStatus("e.kine.W2", 1);

  // --- INPUT BINDINGS (Variable Structure from QE_ana) ---
  const double TI_TICK_SEC = 4e-9;
  double evtime = 0.0;
  setrootvar::setbranch(C, "g", "evtime", &evtime);
  double evnum_T;  setrootvar::setbranch(C,"g","evnum",&evnum_T);

  //double eSH, xSH, ySH, atimeSH;
  //setrootvar::setbranch(C, "bb.sh", std::vector<std::string>{"e","x","y","atimeblk"}, std::vector<void*>{&eSH, &xSH, &ySH, &atimeSH});
  // bbcal sh clus var
  double eSH,xSH,ySH,atimeSH;
  int idSH;
  std::vector<std::string> bbcalclvar = {"e","x","y","atimeblk","idblk"};
  std::vector<void*> bbcalclvar_mem = {&eSH,&xSH,&ySH,&atimeSH,&idSH};
  setrootvar::setbranch(C,"bb.sh",bbcalclvar,bbcalclvar_mem);

  //double ePS, xPS;
  //setrootvar::setbranch(C, "bb.ps", std::vector<std::string>{"e","x"}, std::vector<void*>{&ePS, &xPS});
  // bbcal ps clus var
  double ePS,xPS,yPS,atimePS;
  int idPS;
  std::vector<std::string> bbcalpsclvar = {"e","x","y","atimeblk","idblk"};
  std::vector<void*> bbcalpsclvar_mem = {&ePS,&xPS,&yPS,&atimePS,&idPS};
  setrootvar::setbranch(C,"bb.ps",bbcalpsclvar,bbcalpsclvar_mem);

  double eop;
  setrootvar::setbranch(C, "bb", std::vector<std::string>{"etot_over_p"}, std::vector<void*>{&eop});

  //double eHCAL[maxhcal], xHCAL[maxhcal], yHCAL[maxhcal], rblkHCAL[maxhcal], cblkHCAL[maxhcal], idblkHCAL[maxhcal], adctimeHCAL[maxhcal];
  //setrootvar::setbranch(C, "sbs.hcal", std::vector<std::string>{"e","x","y","rowblk","colblk","idblk","atimeblk"}, std::vector<void*>{&eHCAL, &xHCAL, &yHCAL, &rblkHCAL, &cblkHCAL, &idblkHCAL, &adctimeHCAL});

  double clusadcHCAL[maxhcal];
  setrootvar::setbranch(C, "sbs.hcal.clus", std::vector<std::string>{"adctime"}, std::vector<void*>{&clusadcHCAL});

  // hcal clus var
  double eHCAL[maxhcal], xHCAL[maxhcal], yHCAL[maxhcal], rblkHCAL[maxhcal], cblkHCAL[maxhcal], idblkHCAL[maxhcal],tdctimeHCAL[maxhcal],adctimeHCAL[maxhcal];
  std::vector<std::string> hcalclvar = {"e","x","y","rowblk","colblk","idblk","tdctimeblk","adctimeHCAL"};
  std::vector<void*> hcalclvar_mem = {&eHCAL,&xHCAL,&yHCAL,&rblkHCAL,&cblkHCAL,&idblkHCAL,&tdctimeHCAL,&atimeHCAL};
  setrootvar::setbranch(C, "sbs.hcal", hcalclvar, hcalclvar_mem);
  double hcal_blk_e[maxhcal], hcal_blk_id[maxhcal];
  setrootvar::setbranch(C, "sbs.hcal.clus_blk", "e", &hcal_blk_e);
  setrootvar::setbranch(C, "sbs.hcal.clus_blk", "id", &hcal_blk_id);

  int grinch_nhits = 0;
  double grinch_track, grinch_clus_size, grinch_tmean, grinch_amp, grinch_tmean_corr;
  std::vector<double> grinch_hit_time(grinchmax), grinch_pmtnum(grinchmax), grinch_clustindex(grinchmax), grinch_hit_time_corr(grinchmax);
  setrootvar::setbranch(C, "Ndata.bb.grinch_tdc.hit", "pmtnum", &grinch_nhits);
  setrootvar::setbranch(C, "bb.grinch_tdc.clus", {"trackindex","size","t_mean","adc","t_mean_corr"}, {&grinch_track, &grinch_clus_size, &grinch_tmean, &grinch_amp, &grinch_tmean_corr});
  setrootvar::setbranch(C, "bb.grinch_tdc.hit", {"time","pmtnum","trackindex","time_corr"}, {grinch_hit_time.data(), grinch_pmtnum.data(), grinch_clustindex.data(), grinch_hit_time_corr.data()}, grinchmax);

  double grinch_bestcluster;
  setrootvar::setbranch(C, "bb.grinch_tdc", std::vector<std::string>{"bestcluster"}, std::vector<void*>{&grinch_bestcluster});

  int nhodo_clus = 0;
  std::vector<double> HODO_bar(kMaxHodoHits), hodo_bartime(kMaxHodoHits), HODO_clus(kMaxHodoHits), hodo_clustime(kMaxHodoHits), hodo_tfinal(kMaxHodoHits);
  setrootvar::setbranch(C, "Ndata.bb.hodotdc.clus.bar.tdc", "meantime", &nhodo_clus);
  setrootvar::setbranch(C, "bb.hodotdc.clus.bar.tdc", {"id","meantime"}, {HODO_bar.data(), hodo_bartime.data()}, kMaxHodoHits);
  setrootvar::setbranch(C, "bb.hodotdc.clus", {"id","tmean","tfinal"}, {HODO_clus.data(), hodo_clustime.data(), hodo_tfinal.data()}, kMaxHodoHits);

  // hodoscope - Kate
  const int maxClus = 1000;
  double hodo_time[maxClus], hodo_id[maxClus];
  int nhodo_clus;
  setrootvar::setbranch(C,"bb.hodotdc.clus.bar.tdc","meantime",&hodo_time);
  setrootvar::setbranch(C,"bb.hodotdc.clus.bar.tdc","id",&hodo_id);
  setrootvar::setbranch(C,"bb.hodotdc.clus","tfinal",&hodo_time);
  setrootvar::setbranch(C,"Ndata.bb.hodotdc.clus.bar.tdc","meantime",&nhodo_clus);

  //double ntrack, p[maxNtr], px[maxNtr], py[maxNtr], pz[maxNtr], vx[maxNtr], vy[maxNtr], vz[maxNtr], xtgt[maxNtr], ytgt[maxNtr], thtgt[maxNtr], phtgt[maxNtr], xfp[maxNtr], yfp[maxNtr], thfp[maxNtr], phfp[maxNtr];
  //setrootvar::setbranch(C, "bb.tr", {"n","p","px","py","pz","vx","vy","vz","tg_x","tg_y","tg_th","tg_ph","r_x","r_y","r_th","r_ph"}, {&ntrack, &p, &px, &py, &pz, &vx, &vy, &vz, &xtgt, &ytgt, &thtgt, &phtgt, &xfp, &yfp, &thfp, &phfp});

  // GEMs
  double gem_nhits[maxNtr], gem_ngood[maxNtr], gem_chi2ndf[maxNtr];
  setrootvar::setbranch(C,"bb.gem.track","nhits",&gem_nhits);
  setrootvar::setbranch(C,"bb.gem.track","ngoodhits",&gem_ngood);
  setrootvar::setbranch(C,"bb.gem.track","chi2ndf",&gem_chi2ndf);

  double sbs_gem_nhits[maxNtr];
  int nsbs_gem_nhits;

  double bb_gem_nhits[maxNtr];
  double bb_gem_track_x[maxNtr];
  double bb_gem_track_y[maxNtr];
  double bb_gem_track_xp[maxNtr];
  double bb_gem_track_yp[maxNtr];

  int nbb_gem_nhits;

  //setrootvar::setbranch(C,"sbs.gem.track","nhits",&sbs_gem_nhits);
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

  //int tdcElemN; double tdcTrig[maxNtr], tdcElem[maxNtr];
  //setrootvar::setbranch(C, "bb.tdctrig", {"tdcelemID","tdcelemID","tdc"}, {&tdcElem, &tdcElemN, &tdcTrig}, 1);

  //double helicity, IHWP;
  //setrootvar::setbranch(C, "scalhel", "hel", &helicity);
  //setrootvar::setbranch(C, "IGL1I00OD16_16", "", &IHWP);

  //double BPMAx, BPMAy, rawcurx, rawcury, rawcur2x, rawcur2y;
  //setrootvar::setbranch(C, "Lrb.BPMA", {"x","y"}, {&BPMAx, &BPMAy});
  //setrootvar::setbranch(C, "Lrb.Raster.rawcur", {"x","y"}, {&rawcurx, &rawcury});
  //setrootvar::setbranch(C, "Lrb.Raster2.rawcur", {"x","y"}, {&rawcur2x, &rawcur2y});

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

  double bb_rftime[maxNtr]; C->SetBranchAddress("bb.tdctrig.rftime", bb_rftime);
  double bb_gem_goodhits[100]; C->SetBranchAddress("bb.gem.track.ngoodhits", &bb_gem_goodhits);

  // --- OUTPUT SETUP ---
  TString outDir = Form("/lustre24/expphy/volatile/halla/sbs/ktevans/Data/%s", conf.Data());
  gSystem->mkdir(outDir, kTRUE);
  EnsureSingleFileOutput();

  TString outFile = Form("%s/QE_Parse_pass3_%s_sbs%dp_nucleon_%s_model%d.root", outDir.Data(), conf.Data(), sbsconf.GetSBSmag(), kin_info.Ntype.Data(), kin_info.model);
  TFile *fout = new TFile(outFile, "RECREATE");
  TTree *Tout = new TTree("Tout", "");

  // Variables for Tout (Original List)
  int T_runnum, T_grinch_nhits, T_nhodo_clus, T_helicity, T_IHWP;
  double T_ebeam, T_trigbits, T_epoch_utc, T_nu, T_Q2, T_W2, T_dpel, T_ephi, T_etheta, T_pcentral, T_vz, T_vx, T_vy, T_xtgt, T_ytgt, T_thtgt, T_phtgt, T_thetabend, T_xfp, T_yfp, T_thfp, T_phfp, T_trP, T_trPx, T_trPy, T_trPz, T_bb_gem_ngoodhits, T_ePS, T_xPS, T_eSH, T_xSH, T_ySH, T_atimeSH, T_eop, T_eHCAL, T_xHCAL, T_yHCAL, T_xHCAL_exp, T_yHCAL_exp, T_dx, T_dy, T_theta_pq, T_adctimeHCAL, T_clusadctimeHCAL, T_grinch_track, T_grinch_clus_size, T_grinch_amp, T_grinch_bestcluster, T_grinch_tmean, T_grinch_tmean_corr, T_coin_time, T_coin_min, T_coin_max, T_hcal_time, T_bbcal_time, T_bb_rftime, T_BPMAx, T_BPMAy, T_rawcurx, T_rawcury, T_rawcur2x, T_rawcur2y, T_He3Pol;
  double T_grinch_hit_time[grinchmax], T_grinch_hit_time_corr[grinchmax], T_grinch_pmtnum[grinchmax], T_grinch_clustindex[grinchmax];
  double T_HODO_clus[kMaxHodoHits], T_HODO_bar[kMaxHodoHits], T_hodo_clustime[kMaxHodoHits], T_hodo_bartime[kMaxHodoHits];
  bool WCut, pCut, nCut, fiduCut, coinCut;

  Tout->Branch("runnum", &T_runnum, "runnum/I");
  Tout->Branch("ebeam", &T_ebeam, "ebeam/D");
  Tout->Branch("trigbits", &T_trigbits, "trigbits/D");
  Tout->Branch("epoch_utc", &T_epoch_utc, "epoch_utc/D");
  Tout->Branch("WCut", &WCut, "WCut/B");
  Tout->Branch("pCut", &pCut, "pCut/B");
  Tout->Branch("nCut", &nCut, "nCut/B");
  Tout->Branch("fiduCut", &fiduCut, "fiduCut/B");
  Tout->Branch("coinCut", &coinCut, "coinCut/B");
  Tout->Branch("nu", &T_nu, "nu/D"); Tout->Branch("Q2", &T_Q2, "Q2/D"); Tout->Branch("W2", &T_W2, "W2/D");
  Tout->Branch("dpel", &T_dpel, "dpel/D"); Tout->Branch("ephi", &T_ephi, "ephi/D"); Tout->Branch("etheta", &T_etheta, "etheta/D");
  Tout->Branch("pcentral", &T_pcentral, "pcentral/D");
  Tout->Branch("vz", &T_vz, "vz/D"); Tout->Branch("vx", &T_vx, "vx/D"); Tout->Branch("vy", &T_vy, "vy/D");
  Tout->Branch("xtgt", &T_xtgt, "xtgt/D"); Tout->Branch("ytgt", &T_ytgt, "ytgt/D"); Tout->Branch("thtgt", &T_thtgt, "thtgt/D"); Tout->Branch("phtgt", &T_phtgt, "phtgt/D");
  Tout->Branch("thetabend", &T_thetabend, "thetabend/D");
  Tout->Branch("xfp", &T_xfp, "xfp/D"); Tout->Branch("yfp", &T_yfp, "yfp/D"); Tout->Branch("thfp", &T_thfp, "thfp/D"); Tout->Branch("phfp", &T_phfp, "phfp/D");
  Tout->Branch("trP", &T_trP, "trP/D"); Tout->Branch("trPx", &T_trPx, "trPx/D"); Tout->Branch("trPy", &T_trPy, "trPy/D"); Tout->Branch("trPz", &T_trPz, "trPz/D");
  Tout->Branch("bb_gem_goodhits", &T_bb_gem_ngoodhits, "bb_gem_goodhits/D");
  Tout->Branch("ePS", &T_ePS, "ePS/D"); Tout->Branch("xPS", &T_xPS, "xPS/D");
  Tout->Branch("eSH", &T_eSH, "eSH/D"); Tout->Branch("xSH", &T_xSH, "xSH/D"); Tout->Branch("ySH", &T_ySH, "ySH/D"); Tout->Branch("atimeSH", &T_atimeSH, "atimeSH/D"); Tout->Branch("eop", &T_eop, "eop/D");
  Tout->Branch("nhodo_clus", &T_nhodo_clus, "nhodo_clus/I");
  Tout->Branch("HODO_clus", T_HODO_clus, "HODO_clus[nhodo_clus]/D");
  Tout->Branch("HODO_bar", T_HODO_bar, "HODO_bar[nhodo_clus]/D");
  Tout->Branch("hodo_clustime", T_hodo_clustime, "hodo_clustime[nhodo_clus]/D");
  Tout->Branch("hodo_bartime", T_hodo_bartime, "hodo_bartime[nhodo_clus]/D");
  Tout->Branch("eHCAL", &T_eHCAL, "eHCAL/D"); Tout->Branch("xHCAL", &T_xHCAL, "xHCAL/D"); Tout->Branch("yHCAL", &T_yHCAL, "yHCAL/D");
  Tout->Branch("xHCAL_exp", &T_xHCAL_exp, "xHCAL_exp/D"); Tout->Branch("yHCAL_exp", &T_yHCAL_exp, "yHCAL_exp/D");
  Tout->Branch("dx", &T_dx, "dx/D"); Tout->Branch("dy", &T_dy, "dy/D"); Tout->Branch("theta_pq", &T_theta_pq, "theta_pq/D");
  Tout->Branch("adctimeHCAL", &T_adctimeHCAL, "adctimeHCAL/D"); Tout->Branch("clusadctimeHCAL", &T_clusadctimeHCAL, "clusadctimeHCAL/D");
  Tout->Branch("grinch_track", &T_grinch_track, "grinch_track/D"); Tout->Branch("grinch_clus_size", &T_grinch_clus_size, "grinch_clus_size/D");
  Tout->Branch("grinch_nhits", &T_grinch_nhits, "grinch_nhits/I");
  Tout->Branch("grinch_hit_time", T_grinch_hit_time, "grinch_hit_time[grinch_nhits]/D");
  Tout->Branch("grinch_hit_time_corr", T_grinch_hit_time_corr, "grinch_hit_time_corr[grinch_nhits]/D");
  Tout->Branch("grinch_pmtnum", T_grinch_pmtnum, "grinch_pmtnum[grinch_nhits]/D");
  Tout->Branch("grinch_amp", &T_grinch_amp, "grinch_amp/D");
  Tout->Branch("grinch_clustindex", T_grinch_clustindex, "grinch_clustindex[grinch_nhits]/D");
  Tout->Branch("grinch_bestcluster", &T_grinch_bestcluster, "grinch_bestcluster/D");
  Tout->Branch("grinch_tmean", &T_grinch_tmean, "grinch_tmean/D");
  Tout->Branch("grinch_tmean_corr", &T_grinch_tmean_corr, "grinch_tmean_corr/D");
  Tout->Branch("coin_time", &T_coin_time, "coin_time/D");
  Tout->Branch("coin_min", &T_coin_min, "coin_min/D"); Tout->Branch("coin_max", &T_coin_max, "coin_max/D");
  Tout->Branch("hcal_time", &T_hcal_time, "hcal_time/D"); Tout->Branch("bbcal_time", &T_bbcal_time, "bbcal_time/D");
  Tout->Branch("bb_rftime", &T_bb_rftime, "bb_rftime/D");
  Tout->Branch("BPMAx", &T_BPMAx, "BPMAx/D"); Tout->Branch("BPMAy", &T_BPMAy, "BPMAy/D");
  Tout->Branch("Rasterx", &T_rawcurx, "Rasterx/D"); Tout->Branch("Rastery", &T_rawcury, "Rastery/D");
  Tout->Branch("Raster2x", &T_rawcur2x, "Raster2x/D"); Tout->Branch("Raster2y", &T_rawcur2y, "Raster2y/D");
  Tout->Branch("helicity", &T_helicity, "helicity/I"); Tout->Branch("IHWP", &T_IHWP, "IHWP/I");
  Tout->Branch("He3Pol", &T_He3Pol, "He3Pol/D");

  long nevent = 0, nevents = C->GetEntries();
  int treenum = 0, currenttreenum = 0, IHWP_run = -100;
  time_t run_time_unix = 0;

  // --- LOOP ---
  while (C->GetEntry(nevent++)) {
    if (nevent % 1000 == 0) { std::cout << nevent * 100.0 / nevents << "% \r"; std::cout.flush(); }

    currenttreenum = C->GetTreeNumber();
    if (nevent == 1 || currenttreenum != treenum) {
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();
      std::string s = C->GetFile()->GetName();
      T_runnum = stoi(s.substr(s.find("_stream0") - 4, 4));
      run_time_unix = C->GetFile()->Get<THaRunBase>("Run_Data")->GetDate().Convert(true);

      int start_event = nevent;
      int max_check = std::min((int)C->GetTree()->GetEntries() - 100, 8000);
      while (C->GetEntry(start_event++) && start_event < nevent + max_check) {
          if(IHWP == 1 || IHWP == -1) IHWP_run = (int)IHWP;
      }
      C->GetEntry(nevent - 1);
    }

    if (GlobalCut->EvalInstance(0) == 0) continue;

    T_epoch_utc = static_cast<double>(run_time_unix) + (evtime * TI_TICK_SEC);
    TDatime time_abs(run_time_unix + (int)(evtime * TI_TICK_SEC));

    auto it = DBInfo.He3Pol.find(time_abs);
    T_He3Pol = (it == DBInfo.He3Pol.end()) ? -1.0 : it->second;
    T_IHWP = IHWP_run;
    T_helicity = (int)helicity;

    T_ebeam = sbsconf.GetEbeam();
    T_trigbits = trigbits;
    T_ePS = ePS; T_xPS = xPS; T_eSH = eSH; T_xSH = xSH; T_ySH = ySH; T_atimeSH = atimeSH; T_eop = eop;
    T_eHCAL = eHCAL[0]; T_xHCAL = xHCAL[0]; T_yHCAL = yHCAL[0]; T_adctimeHCAL = adctimeHCAL[0];
    T_clusadctimeHCAL = clusadcHCAL[0];
    T_bb_gem_ngoodhits = bb_gem_goodhits[0];
    T_bb_rftime = bb_rftime[0];

    T_nhodo_clus = nhodo_clus;
    for(int i=0; i<nhodo_clus && i<kMaxHodoHits; ++i){
      T_HODO_bar[i] = HODO_bar[i]; T_HODO_clus[i] = HODO_clus[i];
      T_hodo_bartime[i] = hodo_bartime[i]; T_hodo_clustime[i] = hodo_clustime[i];
    }

    T_coin_time = atimeSH - adctimeHCAL[0];
    T_coin_min = kin_info.coin_time_cut[0] - kin_info.Nsigma_coin_time*kin_info.coin_time_cut[1];
    T_coin_max = kin_info.coin_time_cut[0] + kin_info.Nsigma_coin_time*kin_info.coin_time_cut[1];
    coinCut = T_coin_time > T_coin_min && T_coin_time < T_coin_max;

    T_grinch_nhits = (grinch_nhits > grinchmax) ? grinchmax : grinch_nhits;
    T_grinch_track = grinch_track; T_grinch_clus_size = grinch_clus_size;
    T_grinch_bestcluster = grinch_bestcluster; T_grinch_amp = grinch_amp;
    T_grinch_tmean = grinch_tmean - atimeSH;
    T_grinch_tmean_corr = grinch_tmean_corr - atimeSH;
    for (int i = 0; i < T_grinch_nhits; i++) {
      T_grinch_hit_time[i] = grinch_hit_time[i] - hodo_tfinal[0];
      T_grinch_hit_time_corr[i] = grinch_hit_time_corr[i];
      T_grinch_pmtnum[i] = grinch_pmtnum[i];
    }

    T_vz = vz[0]; T_vx = vx[0]; T_vy = vy[0];
    T_xtgt = xtgt[0]; T_ytgt = ytgt[0]; T_thtgt = thtgt[0]; T_phtgt = phtgt[0];
    T_xfp = xfp[0]; T_yfp = yfp[0]; T_thfp = thfp[0]; T_phfp = phfp[0];
    T_trP = p[0]; T_trPx = px[0]; T_trPy = py[0]; T_trPz = pz[0];

    TLorentzVector Pe(0, 0, T_ebeam, T_ebeam);
    TLorentzVector Peprime(px[0]*(p[0]/p[0]), py[0]*(p[0]/p[0]), pz[0]*(p[0]/p[0]), p[0]);
    TLorentzVector q = Pe - Peprime;
    T_nu = q.E(); T_Q2 = -q.M2(); T_W2 = (Pe + TLorentzVector(0,0,0,0.938) - Peprime).M2();
    T_dpel = Peprime.E()/kine::pcentral(T_ebeam, kine::etheta(Peprime), kin_info.Ntype.Data()) - 1.0;

    TVector3 pNhat = q.Vect().Unit();
    std::vector<TVector3> HCAL_axes; kine::SetHCALaxes(sbsconf.GetSBStheta_rad(), HCAL_axes);
    TVector3 HCAL_origin = sbsconf.GetHCALdist()*HCAL_axes[2] + kin_info.hcal_voffset*HCAL_axes[0];
    std::vector<double> xyHCAL_exp; kine::GetxyHCALexpect(TVector3(0,0,vz[0]), pNhat, HCAL_origin, HCAL_axes, xyHCAL_exp);
    T_dx = xHCAL[0] - xyHCAL_exp[0]; T_dy = yHCAL[0] - xyHCAL_exp[1];
    T_xHCAL_exp = xyHCAL_exp[0]; T_yHCAL_exp = xyHCAL_exp[1];

    Tout->Fill();
  }

  fout->Write(); fout->Close();
  if (GlobalCut) delete GlobalCut;
  if (C) { C->ResetBranchAddresses(); delete C; }
  return 0;
}
