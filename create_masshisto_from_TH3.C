#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>

#include <TROOT.h>
#include <TMinuit.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TPaveText.h>
#include <TGaxis.h>
#endif

void create_masshisto_from_TH3(bool save_file = kFALSE){

  gStyle -> SetOptStat(0);
  gSystem -> CompileMacro("settings.h");
  gROOT -> ProcessLine(".x binning.C");

  //char INPUT_FILE_NAME[300] = "/home/luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/Histos_full_statistics.root"; // for ubuntu
  //char INPUT_FILE_NAME[300] = "/Users/Luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/Histos_full_statistics.root"; // for mac
  char INPUT_FILE_NAME[300] = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/Histos_full_statistics.root";
  printf("Opening %s ... \n",INPUT_FILE_NAME);
  TFile *input_file = new TFile(INPUT_FILE_NAME,"READ");
  input_file -> ls();

  return;

  char TH3_NAME[50];
  const int Npt_ranges = 12;
  TH3D *hMassCostPhiHE_2m[Npt_ranges];
  for(int i = 0;i < Npt_ranges;i++){
    sprintf(TH3_NAME,"hMassCostPhiHE_%ipt%i_2m",i,i+1);
    hMassCostPhiHE_2m[i] = (TH3D*) input_file -> Get(TH3_NAME);
  }

  sprintf(TH3_NAME,"hMassCostPhiHE_%ipt%i_2m",2,6);
  hMassCostPhiHE_2m_2pt6_prova = (TH3D*) input_file -> Get(TH3_NAME);

  //============================================================================
  printf("Defining the pT ranges ... \n");
  //============================================================================
  for(int i = 0;i < N_pt_bin;i++){printf("[%i < pT < %i GeV/c]",min_pt_bin[i],max_pt_bin[i]);}
  printf("\n");

  sprintf(TH3_NAME,"hMassCostPhiHE_%ipt%i_2m",min_pt_bin[0],max_pt_bin[0]);
  TH3D *hMassCostPhiHE_0pt2_2m = (TH3D*) hMassCostPhiHE_2m[0] -> Clone(); // 0 < pT < 1 GeV/c
  hMassCostPhiHE_0pt2_2m -> SetTitle(TH3_NAME);
  hMassCostPhiHE_0pt2_2m -> Add(hMassCostPhiHE_2m[1]); // 1 < pT < 2 GeV/c

  sprintf(TH3_NAME,"hMassCostPhiHE_%ipt%i_2m",min_pt_bin[1],max_pt_bin[1]);
  TH3D *hMassCostPhiHE_2pt6_2m = (TH3D*) hMassCostPhiHE_2m[2] -> Clone(); // 2 < pT < 3 GeV/c
  hMassCostPhiHE_2pt6_2m -> SetTitle(TH3_NAME);
  hMassCostPhiHE_2pt6_2m -> Add(hMassCostPhiHE_2m[3]); // 3 < pT < 4 GeV/c
  hMassCostPhiHE_2pt6_2m -> Add(hMassCostPhiHE_2m[4]); // 4 < pT < 5 GeV/c
  hMassCostPhiHE_2pt6_2m -> Add(hMassCostPhiHE_2m[5]); // 5 < pT < 6 GeV/c

  sprintf(TH3_NAME,"hMassCostPhiHE_%ipt%i_2m",min_pt_bin[2],max_pt_bin[2]);
  TH3D *hMassCostPhiHE_6pt12_2m = (TH3D*) hMassCostPhiHE_2m[6] -> Clone(); // 6 < pT < 7 GeV/c
  hMassCostPhiHE_6pt12_2m -> SetTitle(TH3_NAME);
  hMassCostPhiHE_6pt12_2m -> Add(hMassCostPhiHE_2m[7]); // 7 < pT < 8 GeV/c
  hMassCostPhiHE_6pt12_2m -> Add(hMassCostPhiHE_2m[8]); // 8 < pT < 9 GeV/c
  hMassCostPhiHE_6pt12_2m -> Add(hMassCostPhiHE_2m[9]); // 9 < pT < 10 GeV/c
  hMassCostPhiHE_6pt12_2m -> Add(hMassCostPhiHE_2m[10]); // 10 < pT < 11 GeV/c
  hMassCostPhiHE_6pt12_2m -> Add(hMassCostPhiHE_2m[11]); // 11 < pT < 12 GeV/c

  //============================================================================
  printf("Creating mass histos... \n");
  //============================================================================
  if(save_file){
    TFile FILE_OUT_COST_PHI_0pt2("~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/GIT_OUTPUT/mass_histos_cost_phi_0pt2.root","RECREATE");
    TFile FILE_OUT_COST_PHI_2pt6("~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/GIT_OUTPUT/mass_histos_cost_phi_2pt6.root","RECREATE");
    TFile FILE_OUT_COST_PHI_6pt12("~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/GIT_OUTPUT/mass_histos_cost_phi_6pt12.root","RECREATE");
  }

  TH1D *hist_mass_HE;
  //TH1D *hist_mass_CS;
  char hist_mass_HE_name[100];
  char hist_mass_CS_name[100];

  for(int k = 0;k < N_pt_bin;k++){
    for(int i = 0;i < N_cost_bins;i++){
      for(int j = 0;j < N_phi_bins;j++){

        sprintf(hist_mass_HE_name,"HE_%ipt%i_%icost%i_%iphi%i",min_pt_bin[k],max_pt_bin[k],min_cost_bin[i],max_cost_bin[i],min_phi_bin[j],max_phi_bin[j]);
        //sprintf(hist_mass_CS_name,"CS_%ipt%i_%icost%i_%iphi%i",min_pt_bin[k],max_pt_bin[k],min_cost_bin[i],max_cost_bin[i],min_phi_bin[j],max_phi_bin[j]);

        if(min_pt_bin[k] == 0){
          hist_mass_HE = (TH1D*) hMassCostPhiHE_0pt2_2m -> ProjectionZ(hist_mass_HE_name,min_cost_bin[i],max_cost_bin[i],min_phi_bin[j],max_phi_bin[j]);
          if(save_file){FILE_OUT_COST_PHI_0pt2.cd(); hist_mass_HE -> SetTitle(hist_mass_HE_name); hist_mass_HE -> Write();}
          //hist_mass_CS = (TH1D*) hMassCostPhiCS_0pt2_2m -> ProjectionZ(hist_mass_CS_name,min_cost_bin[i],max_cost_bin[i],min_phi_bin[j],max_phi_bin[j]);
        }
        if(min_pt_bin[k] == 2){
          hist_mass_HE = (TH1D*) hMassCostPhiHE_2pt6_2m -> ProjectionZ(hist_mass_HE_name,min_cost_bin[i],max_cost_bin[i],min_phi_bin[j],max_phi_bin[j]);
          if(save_file){FILE_OUT_COST_PHI_2pt6.cd(); hist_mass_HE -> SetTitle(hist_mass_HE_name); hist_mass_HE -> Write();}
          //hist_mass_CS = (TH1D*) hMassCostPhiCS_2pt6_2m -> ProjectionZ(hist_mass_CS_name,min_cost_bin[i],max_cost_bin[i],min_phi_bin[j],max_phi_bin[j]);
        }
        if(min_pt_bin[k] == 6){
          hist_mass_HE = (TH1D*) hMassCostPhiHE_6pt12_2m -> ProjectionZ(hist_mass_HE_name,min_cost_bin[i],max_cost_bin[i],min_phi_bin[j],max_phi_bin[j]);
          if(save_file){FILE_OUT_COST_PHI_6pt12.cd(); hist_mass_HE -> SetTitle(hist_mass_HE_name); hist_mass_HE -> Write();}
          //hist_mass_CS = (TH1D*) hMassCostPhiCS_6pt12_2m -> ProjectionZ(hist_mass_CS_name,min_cost_bin[i],max_cost_bin[i],min_phi_bin[j],max_phi_bin[j]);
        }

      }
    }
  }

  //sprintf(hist_mass_HE_name,"HE_%ipt%i",2,6);
  //hist_mass_HE = (TH1D*) hMassCostPhiHE_2m_2pt6_prova -> ProjectionZ(hist_mass_HE_name,-1,1,0,PI);
  //if(save_file){FILE_OUT_COST_PHI_2pt6.cd(); hist_mass_HE -> SetTitle(hist_mass_HE_name); hist_mass_HE -> Write();}

  FILE_OUT_COST_PHI_0pt2.Close();
  FILE_OUT_COST_PHI_2pt6.Close();
  FILE_OUT_COST_PHI_6pt12.Close();

  TCanvas *c_grid = new TCanvas("c_grid","c_grid",20,20,600,600);
  TH2D *h_grid = new TH2D("h_grid","",100,-1,1,100,0,PI);
  h_grid -> Draw();
  for(int i = 0;i < N_line_cost;i++) line_cost[i] -> Draw("same");
  for(int i = 0;i < N_line_phi;i++) line_phi[i] -> Draw("same");

}
