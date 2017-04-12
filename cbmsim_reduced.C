#define cbmsim_reduced_cxx
#include "cbmsim_reduced.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void cbmsim_reduced::Loop(TString string)
{
//   In a ROOT session, you can do:
//      Root > .L cbmsim_reduced.C
//      Root > cbmsim_reduced t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
   for (Long64_t itracks=0;itracks<n_tracks_mpd;itracks++){
      if(TMath::Abs(signed_pt_mpd[itracks])<3 && TMath::Abs(eta_mpd[itracks]) < 1.5 && n_hits_mpd[itracks] > 32 && mother_ID_mc[id_from_mc_mpd[itracks]] == -1 ) {
         
         if(PDG_code_mc[id_from_mc_mpd[itracks]]==2212)//protons
         {
            h_pt_mpd_p->Fill(TMath::Abs(signed_pt_mpd[itracks]));
            h_eta_mpd_p->Fill(TMath::Abs(eta_mpd[itracks]));
            h_theta_mpd_p->Fill(TMath::Abs(theta_mpd[itracks]));
            h_phi_mpd_p->Fill(TMath::Abs(phi_mpd[itracks]));
         }
         if(PDG_code_mc[id_from_mc_mpd[itracks]]==321)//kaons
         {
            h_pt_mpd_k->Fill(TMath::Abs(signed_pt_mpd[itracks]));  
            h_eta_mpd_k->Fill(TMath::Abs(eta_mpd[itracks]));
            h_theta_mpd_k->Fill(TMath::Abs(theta_mpd[itracks]));
            h_phi_mpd_k->Fill(TMath::Abs(phi_mpd[itracks]));
         }
         if(PDG_code_mc[id_from_mc_mpd[itracks]]==211)//pions
         {
            h_pt_mpd_pi->Fill(TMath::Abs(signed_pt_mpd[itracks]));
            h_eta_mpd_pi->Fill(TMath::Abs(eta_mpd[itracks]));
            h_theta_mpd_pi->Fill(TMath::Abs(theta_mpd[itracks]));
            h_phi_mpd_pi->Fill(TMath::Abs(phi_mpd[itracks]));
         }
      

   }
   }

   for (Long64_t itracks=0;itracks<n_tracks_mc;itracks++){
      if(pt_mc[itracks] < 3 && eta_mc[itracks]<1.5 && mother_ID_mc[itracks] == -1 )  {

         if(PDG_code_mc[itracks]==2212)//protons
         {
            h_pt_mc_p->Fill(pt_mc[itracks]);
            h_eta_mc_p->Fill(eta_mc[itracks]);
            
            
         }
         if(PDG_code_mc[itracks]==130)//kaons
         {
            h_pt_mc_k->Fill(pt_mc[itracks]); 
            h_eta_mc_k->Fill(eta_mc[itracks]);
            
            
         }
         if(PDG_code_mc[itracks]==111)//pions
         {
            h_pt_mc_pi->Fill(pt_mc[itracks]);
            h_eta_mc_pi->Fill(eta_mc[itracks]);
            
            
         }
   }
   }

   }
   h_eff_pt_p->Divide(h_pt_mpd_p,h_pt_mc_p);
   h_eff_pt_k->Divide(h_pt_mpd_k,h_pt_mc_k);
   h_eff_pt_pi->Divide(h_pt_mpd_pi,h_pt_mc_pi);
   TH1F* h_pt_difference = (TH1F*)h_pt_mc->Clone("h_pt_difference");
   Float_t  diference = 0;
   // for (Int_t i =0; i<h_pt_mpd->GetEntries(); i++)
   // {
   //       diference = h_pt_mpd[i] - h_pt_mc[i];
   //       h_pt_difference->Fill(diference);
   // }
   h_pt_difference->Add(h_pt_mpd,-1);
   TFile *foutput=new TFile(string.Data(),"recreate");
   foutput->cd();
   h_pt_mpd_p->Write();
   h_pt_mc_p->Write();
   h_eff_pt_p->Write();
   h_pt_mpd_k->Write();
   h_pt_mc_k->Write();
   h_eff_pt_k->Write();
   h_pt_mpd_pi->Write();
   h_pt_mc_pi->Write();
   h_eff_pt_pi->Write();  
   foutput->Close();
   delete foutput;
   
}
