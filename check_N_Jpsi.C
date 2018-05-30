#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <vector>

#include <TMinuit.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TF2.h>
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

void check_N_Jpsi(){

  //gSystem -> CompileMacro("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/H_FILES/settings.h");
  //gSystem -> CompileMacro("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/H_FILES/signal_extraction.h");
  //gSystem -> CompileMacro("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/H_FILES/accXeff.h");
  gSystem -> CompileMacro("settings.h");
  gROOT -> ProcessLine(".x binning.C");

  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(1);
  TGaxis::SetMaxDigits(2);

  int matrix_N_Jpsi_HE[N_cost_bins][N_phi_bins];
  int matrix_stat_N_Jpsi_HE[N_cost_bins][N_phi_bins];

  TFile *N_Jpsi_file = new TFile("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/N_Jpsi.root","READ");
  TTree *N_Jpsi_tree = (TTree*) N_Jpsi_file -> Get("CB2_VWG");

  N_Jpsi_tree -> SetBranchAddress("N_Jpsi_HE",matrix_N_Jpsi_HE);
  N_Jpsi_tree -> SetBranchAddress("Stat_Jpsi_HE",matrix_stat_N_Jpsi_HE);

  for(int i = 0;i < N_Jpsi_tree -> GetEntries();i++){
    N_Jpsi_tree -> GetEntry(i);
  }

  for(int i = 0;i < N_cost_bins;i++){
    for(int j = 0;j < N_phi_bins;j++){
      cout << matrix_N_Jpsi_HE[i][j] << " ";
    }
    cout << endl;
  }

  return;

  //============================================================================
  // N_Jpsi HISTOGRAMS
  //============================================================================

  printf("N Jpsi matrix ... \n");

  //----------------------------------------------------------------------------
  //RANGE 0 < phi < PI ; -1 < cost < 1 [BC = Bent&Cut]
  //----------------------------------------------------------------------------

  TH2D *hist_N_Jpsi_2pt6_HE = new TH2D("hist_N_Jpsi_2pt6_HE","",N_cost_bins,value_cost,N_phi_bins,value_phi);
  hist_N_Jpsi_2pt6_HE -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  hist_N_Jpsi_2pt6_HE -> GetYaxis() -> SetTitle("#phi_{HE}");

  TH2D *hist_N_Jpsi_2pt6_HE_area_corrected = new TH2D("hist_N_Jpsi_2pt6_HE_area_corrected","",N_cost_bins,value_cost,N_phi_bins,value_phi);
  hist_N_Jpsi_2pt6_HE_area_corrected -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  hist_N_Jpsi_2pt6_HE_area_corrected -> GetYaxis() -> SetTitle("#phi_{HE}");

  for(int i = 0;i< N_cost_bins;i++){
    for(int j = 0;j < N_phi_bins;j++){
      hist_N_Jpsi_2pt6_HE -> SetBinContent(i+1,j+1,matrix_N_Jpsi_HE[i][j]);
      hist_N_Jpsi_2pt6_HE -> SetBinError(i+1,j+1,matrix_stat_N_Jpsi_HE[i][j]);
      hist_N_Jpsi_2pt6_HE_area_corrected -> SetBinContent(i+1,j+1,matrix_N_Jpsi_HE[i][j]/bin_area[i][j]);
      hist_N_Jpsi_2pt6_HE_area_corrected -> SetBinError(i+1,j+1,matrix_stat_N_Jpsi_HE[i][j]/bin_area[i][j]);
    }
  }

  TCanvas *c_N_Jpsi_2pt6_HE = new TCanvas("c_N_Jpsi_2pt6_HE","c_N_Jpsi_2pt6_HE",4,132,1024,768);
  hist_N_Jpsi_2pt6_HE -> Draw("COLZtext");

  for(int i = 0;i < N_line_cost;i++){
    if(i < N_line_phi) line_phi[i] -> Draw("same");
    line_cost[i] -> Draw("same");
  }

  TCanvas *c_N_Jpsi_2pt6_HE_area_corrected = new TCanvas("c_N_Jpsi_2pt6_HE_area_corrected","c_N_Jpsi_2pt6_HE_area_corrected",4,132,1024,768);
  hist_N_Jpsi_2pt6_HE_area_corrected -> Draw("COLZ");

  for(int i = 0;i < N_line_cost;i++){
    if(i < N_line_phi) line_phi[i] -> Draw("same");
    line_cost[i] -> Draw("same");
  }

  //============================================================================
  // RELATIVE STAT ERROR FOR EACH BIN
  //============================================================================

  printf("Relative statistic error matrix ... \n");

  TH2D *hist_rel_stat_2pt6_HE = new TH2D("hist_rel_stat_2pt6_HE","",N_cost_bins,value_cost,N_phi_bins,value_phi);
  hist_rel_stat_2pt6_HE -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  hist_rel_stat_2pt6_HE -> GetYaxis() -> SetTitle("#phi_{HE}");

  double rel_stat = 0;

  for(int i = 0;i< N_cost_bins;i++){
    for(int j = 0;j < N_phi_bins;j++){
      rel_stat = (double) matrix_stat_N_Jpsi_HE[i][j]/matrix_N_Jpsi_HE[i][j];
      printf("%3.2f ",rel_stat);
      hist_rel_stat_2pt6_HE -> SetBinContent(i+1,j+1,rel_stat);
    }
    printf("\n");
  }

  TCanvas *c_rel_stat_2pt6_HE = new TCanvas("c_rel_stat_2pt6_HE","c_rel_stat_2pt6_HE",4,132,1024,768);
  hist_rel_stat_2pt6_HE -> Draw("COLZ");

}
