//
// To execute, do the following:
//
// > root
// > .L scripts/PreShowerEnergyCut.C
// > PreShowerEnergyCut()
//
// Note: you must have a directory called "plots" created in order to save
// the canvas to a PDF properly.
//
// Kate Evans
// Last modified: July 12 2024


#include <TF1.h>

void PreShowerEnergyCut()
{

 //Add root files to the chain T, can add as many files as you want
 TChain* T = new TChain("T");
 T->Add("/volatile/halla/sbs/sbs-gen/GEN_REPLAYS/Rootfiles/pass1/GEN2/He3/rootfiles/e1209016_fullreplay_2212_stream0_2_seg11_11.root");

 int maxNtr = 200;
 
 //Set the branch address for variables you want to use
 T->SetBranchStatus("*", 0);
 // bb.ps branches
 T->SetBranchStatus("bb.ps.*", 1);
  Double_t psNclus;     T->SetBranchAddress("bb.ps.nclus", &psNclus);
  Double_t psIdblk;     T->SetBranchAddress("bb.ps.idblk", &psIdblk);
  Double_t psRowblk;    T->SetBranchAddress("bb.ps.rowblk", &psRowblk);
  Double_t psColblk;    T->SetBranchAddress("bb.ps.colblk", &psColblk);
  Double_t psNblk;      T->SetBranchAddress("bb.ps.nblk", &psNblk);
  Double_t psAtime;     T->SetBranchAddress("bb.ps.atimeblk", &psAtime);
  Double_t psE;         T->SetBranchAddress("bb.ps.e", &psE);
  Double_t psX;         T->SetBranchAddress("bb.ps.x", &psX);
  Double_t psY;         T->SetBranchAddress("bb.ps.y", &psY);
  Double_t psClBlkId[maxNtr];  T->SetBranchAddress("bb.ps.clus_blk.id", &psClBlkId);
  Double_t psClBlkE[maxNtr];   T->SetBranchAddress("bb.ps.clus_blk.e", &psClBlkE);
  Double_t psClBlkX[maxNtr];   T->SetBranchAddress("bb.ps.clus_blk.x", &psClBlkX);
  Double_t psClBlkY[maxNtr];   T->SetBranchAddress("bb.ps.clus_blk.y", &psClBlkY);
  Double_t psClBlkRow[maxNtr]; T->SetBranchAddress("bb.ps.clus_blk.row", &psClBlkRow);
  Double_t psClBlkCol[maxNtr]; T->SetBranchAddress("bb.ps.clus_blk.col", &psClBlkCol);
  Double_t psClBlkAtime[maxNtr]; T->SetBranchAddress("bb.ps.clus_blk.atime", &psClBlkAtime);
  Double_t psAgainblk;         T->SetBranchAddress("bb.ps.againblk", &psAgainblk);
 // bb.sh branches
 T->SetBranchStatus("bb.sh.*", 1);
  Double_t shNclus;            T->SetBranchAddress("bb.sh.nclus", &shNclus);
  Double_t shIdblk;            T->SetBranchAddress("bb.sh.idblk", &shIdblk);
  Double_t shRowblk;           T->SetBranchAddress("bb.sh.rowblk", &shRowblk);
  Double_t shColblk;           T->SetBranchAddress("bb.sh.colblk", &shColblk);
  Double_t shNblk;             T->SetBranchAddress("bb.sh.nblk", &shNblk);
  Double_t shAtime;            T->SetBranchAddress("bb.sh.atimeblk", &shAtime);
  Double_t shE;                T->SetBranchAddress("bb.sh.e", &shE);
  Double_t shX;                T->SetBranchAddress("bb.sh.x", &shX);
  Double_t shY;                T->SetBranchAddress("bb.sh.y", &shY);
  Double_t shClBlkId[maxNtr];  T->SetBranchAddress("bb.sh.clus_blk.id", &shClBlkId);
  Double_t shClBlkE[maxNtr];   T->SetBranchAddress("bb.sh.clus_blk.e", &shClBlkE);
  Double_t shClBlkX[maxNtr];   T->SetBranchAddress("bb.sh.clus_blk.x", &shClBlkX);
  Double_t shClBlkY[maxNtr];   T->SetBranchAddress("bb.sh.clus_blk.y", &shClBlkY);
  Double_t shClBlkRow[maxNtr]; T->SetBranchAddress("bb.sh.clus_blk.row", &shClBlkRow);
  Double_t shClBlkCol[maxNtr]; T->SetBranchAddress("bb.sh.clus_blk.col", &shClBlkCol);
  Double_t shClBlkAtime[maxNtr]; T->SetBranchAddress("bb.sh.clus_blk.atime", &shClBlkAtime);
  Double_t shAgainblk;         T->SetBranchAddress("bb.sh.againblk", &shAgainblk);
 // sbs.hcal branches
  Double_t hcalE;    T->SetBranchStatus("sbs.hcal.e",1); T->SetBranchAddress("sbs.hcal.e", &hcalE);
  Double_t hcalX;    T->SetBranchStatus("sbs.hcal.x",1); T->SetBranchAddress("sbs.hcal.x", &hcalX);
  Double_t hcalY;    T->SetBranchStatus("sbs.hcal.y",1); T->SetBranchAddress("sbs.hcal.y", &hcalY);
  Double_t hcalAtime;    T->SetBranchStatus("sbs.hcal.atimeblk",1); T->SetBranchAddress("sbs.hcal.atimeblk", &hcalAtime);
 // bb.tr branches
 T->SetBranchStatus("bb.tr.*", 1);
  Double_t trN;                T->SetBranchAddress("bb.tr.n", &trN);
  Double_t trP[maxNtr];        T->SetBranchAddress("bb.tr.p", &trP);
  Double_t trPx[maxNtr];       T->SetBranchAddress("bb.tr.px", &trPx);
  Double_t trPy[maxNtr];       T->SetBranchAddress("bb.tr.py", &trPy);
  Double_t trPz[maxNtr];       T->SetBranchAddress("bb.tr.pz", &trPz);
  Double_t trX[maxNtr];        T->SetBranchAddress("bb.tr.x", &trX);
  Double_t trY[maxNtr];        T->SetBranchAddress("bb.tr.y", &trY);
  Double_t trTh[maxNtr];       T->SetBranchAddress("bb.tr.th", &trTh);
  Double_t trPh[maxNtr];       T->SetBranchAddress("bb.tr.ph", &trPh);
  Double_t trVz[maxNtr];       T->SetBranchAddress("bb.tr.vz", &trVz);
  Double_t trVy[maxNtr];       T->SetBranchAddress("bb.tr.vy", &trVy);
  Double_t trTgth[maxNtr];     T->SetBranchAddress("bb.tr.tg_th", &trTgth);
  Double_t trTgph[maxNtr];     T->SetBranchAddress("bb.tr.tg_ph", &trTgph);
  Double_t trRx[maxNtr];       T->SetBranchAddress("bb.tr.r_x", &trRx);
  Double_t trRy[maxNtr];       T->SetBranchAddress("bb.tr.r_y", &trRy);
  Double_t trRth[maxNtr];      T->SetBranchAddress("bb.tr.r_th", &trRth);
  Double_t trRph[maxNtr];      T->SetBranchAddress("bb.tr.r_ph", &trRph);
 // bb.hodotdc branches
  Double_t thTdiff[maxNtr];    T->SetBranchStatus("bb.hodotdc.clus.tdiff",1); T->SetBranchAddress("bb.hodotdc.clus.tdiff", &thTdiff);
  Double_t thTmean[maxNtr];    T->SetBranchStatus("bb.hodotdc.clus.tmean",1); T->SetBranchAddress("bb.hodotdc.clus.tmean", &thTmean);
  Double_t thTOTmean[maxNtr];  T->SetBranchStatus("bb.hodotdc.clus.totmean",1); T->SetBranchAddress("bb.hodotdc.clus.totmean", &thTOTmean);
 // Event info
  //T->SetMakeClass(1);
  T->SetBranchStatus("fEvtHdr.*", 1);
  UInt_t rnum;           T->SetBranchAddress("fEvtHdr.fRun", &rnum);
  UInt_t trigbits;       T->SetBranchAddress("fEvtHdr.fTrigBits", &trigbits);
  ULong64_t gevnum;      T->SetBranchAddress("fEvtHdr.fEvtNum", &gevnum);
  Double_t W2;           T->SetBranchAddress("e.kine.W2", &W2);
  Double_t Q2;           T->SetBranchAddress("e.kine.Q2", &Q2);
  T->SetBranchStatus("sbs.hcal.e", 1);
  T->SetBranchStatus("bb.gem.track.nhits", 1);
  T->SetBranchStatus("bb.grinch_tdc.*", 1);
  T->SetBranchStatus("bb.gem.track.ngoodhits", 1);
  T->SetBranchStatus("bb.gem.track.chi2ndf", 1);
  T->SetBranchStatus("bb.grinch_tdc.clus.trackindex", 1);
  T->SetBranchStatus("bb.grinch_tdc.clus.size", 1);

  Double_t localMin;


 //Scan through all the entries in the TChain T
 //If the rootfiles are empty or don't exist, there will be 0 entries
 //If there are entries, then print out how many
 if(T->GetEntries()==0)
 {
   std::cerr << "\n --- No ROOT file found!! --- \n\n";
   throw;
 }
 else std::cout << "\nFound " << T->GetEntries() << " events. Starting analysis.. \n";

 //Define 1D histogram to plot ps_e
 TH1D* h_ps_e = new TH1D("h_ps_e",";ps_e",100,0.0,2.0);
 h_ps_e->GetXaxis()->SetTitle("PreShower Energy [GeV]");

 //Define 1D histogram to plot ps_e to find local min
 TH1D* h_ps_e_min = new TH1D("h_ps_e_min",";ps_e_min",70,0.1,0.8);

 //Define 1D histogram to plot ps_e with cut
 TH1D* h_ps_e_cut = new TH1D("h_ps_e_cut",";ps_e_cut",100,0.0,2.0);
 h_ps_e_cut->GetXaxis()->SetTitle("PreShower Energy [GeV] with Cut on Pion Peak");

 //Loop over all events to fill the histogram
 for (size_t iev = 0; iev < T->GetEntries(); iev++)
 {

   T->GetEntry(iev);

   h_ps_e->Fill(psE);

   if(psE >= 0.1 && psE <= 0.8)
   {
     h_ps_e_min->Fill(psE);
   }

 }//end process over events

 //localMin = h_ps_e_min->GetBinContent(h_ps_e_min->GetMinimumBin());
 localMin = h_ps_e->GetBinCenter(h_ps_e_min->GetMinimumBin());

 //Loop over all events to fill the histogram
 for (size_t iev = 0; iev < T->GetEntries(); iev++)
 {

   T->GetEntry(iev);

   if(psE >= localMin)
   {
     h_ps_e_cut->Fill(psE);
   }

 }//end process over events

 //Draw canvas, divide into two plots and draw on each one
 TCanvas *c1 = new TCanvas("c1","ps_e cut",100,100,700,700);
 c1->Divide(1,2);
 c1->cd(1);
 h_ps_e->Draw();
 c1->cd(2);
 h_ps_e_cut->Draw();

 printf("PreShower Energy Cut at %f GeV\n", localMin);

 //Save the canvas to a pdf
 c1->Print("plots/PreShowerEnergyCut.pdf");

}
