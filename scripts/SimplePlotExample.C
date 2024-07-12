// This example script will read in a chain of root files and make a 1D
// histogram of W^2 and PreShower Energy. It will then save these plots to
// a PDF. To execute, do the following:
//
// > root
// > .L scripts/SimplePlotExample.C
// > SimplePlotExample()
//
// Note: you must have a directory called "plots" created in order to save
// the canvas to a PDF properly.
//
// Kate Evans
// Last modified: July 12 2024


#include <TF1.h>

void SimplePlotExample()
{

 //Add root files to the chain T, can add as many files as you want
 TChain* T = new TChain("T");
 T->Add("/volatile/halla/sbs/sbs-gen/GEN_REPLAYS/Rootfiles/pass1/GEN2/He3/rootfiles/e1209016_fullreplay_2212_stream0_2_seg11_11.root");

 //Set the branch address for variables you want to use
 Double_t W2;           T->SetBranchAddress("e.kine.W2", &W2);
 Double_t psE;          T->SetBranchAddress("bb.ps.e", &psE);

 //Scan through all the entries in the TChain T
 //If the rootfiles are empty or don't exist, there will be 0 entries
 //If there are entries, then print out how many
 if(T->GetEntries()==0)
 {
   std::cerr << "\n --- No ROOT file found!! --- \n\n";
   throw;
 }
 else std::cout << "\nFound " << T->GetEntries() << " events. Starting analysis.. \n";

 //Define 1D histogram to plot W2
 TH1D* h_W = new TH1D("h_W",";W",150,0.0,3.0);
 h_W->GetXaxis()->SetTitle("W^2 [GeV^2]");

 //Define 1D histogram to plot ps_e
 TH1D* h_ps_e = new TH1D("h_ps_e",";ps_e",100,0.0,2.0);
 h_ps_e->GetXaxis()->SetTitle("PreShower Energy [GeV^2]");

 //Loop over all events to fill the histogram
 for (size_t iev = 0; iev < T->GetEntries(); iev++)
 {

   T->GetEntry(iev);

   h_W->Fill(W2);
   h_ps_e->Fill(psE);

 }//end process over events

 //Draw canvas, divide into two plots and draw on each one
 TCanvas *c1 = new TCanvas("c1","Example Plots",100,100,700,700);
 c1->Divide(1,2);
 c1->cd(1);
 h_W->Draw();
 c1->cd(2);
 h_ps_e->Draw();

 printf("You've completed the script!\n");

 //Save the canvas to a pdf
 c1->Print("plots/SimplePlotExample.pdf");

}
