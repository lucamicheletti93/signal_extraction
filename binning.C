#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>

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
#include <settings.h>
#endif

void binning(){

  gSystem -> CompileMacro("settings.h");

  //============================================================================
  //CONVERSION FACTOR cost
  //============================================================================

  double diff = 0;

  for(int i = 0;i < N_cost_bins;i++){
    width_cost[i+1] = ((max_cost - min_cost)/(double) N_TH3_cost_bins)*(bins_cost[i+1] - bins_cost[i]);
  }

  for(int i = 0;i < dim_cost;i++){
    diff += width_cost[i];
    value_cost[i] = min_cost + diff;
    cout << value_cost[i] << endl;
    if(TMath::Abs(value_cost[i]) < 0.01){value_cost[i] = 0.;}
  }

  //============================================================================
  // COMPUTE D(cost)
  //============================================================================

  for(int i = 0;i < N_cost_bins;i++){
    Dcost[i] = value_cost[i+1] - value_cost[i];
  }

  //============================================================================
  //CONVERSION FACTOR phi
  //============================================================================

  for(int i = 0;i < N_phi_bins;i++){
    width_phi[i+1] = ((max_phi - min_phi)/(double) N_TH3_phi_bins)*(bins_phi[i+1] - bins_phi[i]);
  }

  diff = 0;

  for(int i = 0;i < dim_phi;i++){
    diff += width_phi[i];
    value_phi[i] = min_phi + diff;
    cout << value_phi[i] << endl;
    if(TMath::Abs(value_phi[i]) < 0.01){value_phi[i] = 0.;}
  }

  //============================================================================
  // COMPUTE D(phi)
  //============================================================================

  double Dphi[N_phi_bins];
  for(int i = 0;i < N_phi_bins;i++){
    Dphi[i] = value_phi[i+1] - value_phi[i];
  }

  //============================================================================
  // COMPUTE THE BIN AREA
  //============================================================================

  for(int i = 0;i < N_cost_bins;i++){
    for(int j = 0;j < N_phi_bins;j++){
      bin_area[i][j] = Dcost[i]*Dphi[j];
    }
  }

  //============================================================================
  // GENERATE LINES
  //============================================================================

  for(int i = 0;i < N_line_cost;i++){
    line_cost[i] = new TLine(value_cost[i+1],0,value_cost[i+1],PI);
  }

  for(int i = 0;i < N_line_phi;i++){
    line_phi[i] = new TLine(-1,value_phi[i+1],1,value_phi[i+1]);
  }

}
