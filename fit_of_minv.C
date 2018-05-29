#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <string>
#include <vector>

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
#include <TCollection.h>
#include <TKey.h>
#include <TGaxis.h>
#include <settings.h>
#endif

Double_t Func_VWG(Double_t *, Double_t *);
Double_t Func_Jpsi_CB2(Double_t *, Double_t *);
Double_t Func_Jpsi_CB2_fix(Double_t *, Double_t *);
Double_t Func_Psi2s_CB2(Double_t *, Double_t *);
Double_t Func_Psi2s_CB2_fix(Double_t *, Double_t *);
Double_t Func_tot(Double_t *, Double_t *);

void single_histo_fit();
void loop_on_histos();
void fit_of_minv(TH1D *, int, int);

Double_t scaling_factor = 1.05154; //factor introduced to pass from the sigma of Jpsi to the sigma of Psi(2S)
double n_psi2s = 0, stat_psi2s = 0, n_jpsi = 0, stat_jpsi = 0;
char test_label[100], fit_status[10];
//double min_cost = 0, max_cost = 0, min_phi = 0, max_phi = 0;
//==============================================================================
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//==============================================================================
void single_histo_fit(){

  string const filename = "mass_histos_cost_phi_2pt6.root";
  TFile *file = new TFile(filename.c_str());

  TH1D *hist_minv = (TH1D*) file -> Get("HE_2pt6_1cost25_36phi50");
  hist_minv -> SetMarkerStyle(20);
  hist_minv -> SetMarkerSize(0.7);
  min_cost = -1, max_cost = -0.5, min_phi = 2.2, max_phi = 3.14;

  /*TH1D *hist_minv = (TH1D*) file -> Get("HE_2pt6_61cost65_26phi28");
  hist_minv -> SetMarkerStyle(20);
  hist_minv -> SetMarkerSize(0.7);
  min_cost = 0.2, max_cost = 0.3, min_phi = 1.57, max_phi = 1.76;*/

  fit_of_minv(hist_minv,0,0);
}
//==============================================================================
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//==============================================================================
void loop_on_histos(){

  gSystem -> CompileMacro("settings.h");
  gROOT -> ProcessLine(".x binning.C");

  bool save_tree = kTRUE;
  //const int N_cost_bins = 18;
  //const int N_phi_bins = 10;
  int N_Jpsi_HE[N_cost_bins][N_phi_bins];
  int Stat_Jpsi_HE[N_cost_bins][N_phi_bins];
  int counter_cost = 0; // counter to fill the matrix
  int counter_phi = 0; // counter to fill the matrix

  vector <int> N_Jpsi;
  vector <int> stat_Jpsi;

  string hist_name;

  TH1D *hist_minv_integrated = new TH1D("hist_minv_integrated","hist_minv_integrated",120,2,5);

  //string const filename = "/home/luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/mass_histos_cost_phi_2pt6.root";
  string const filename = "/Users/Luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/mass_histos_cost_phi_2pt6.root";
  TFile *file = new TFile(filename.c_str());

  TIter iter(file -> GetListOfKeys());
  TKey *key;
  while(key = (TKey*)iter()){
      TObject *obj = key -> ReadObj();
      TH1D *hist_minv = dynamic_cast<TH1D*>(obj);
      hist_name = hist_minv -> GetName();

      if(hist_name.find("HE") != std::string::npos){
        hist_minv_integrated -> Add(hist_minv);
        fit_of_minv(hist_minv,counter_cost,counter_phi);

        if(counter_phi < N_phi_bins){
          N_Jpsi_HE[counter_cost][counter_phi] = n_jpsi;
          Stat_Jpsi_HE[counter_cost][counter_phi] = stat_jpsi;
          counter_phi++;
        }
        else{
          counter_cost++;
          counter_phi = 0;
          N_Jpsi_HE[counter_cost][counter_phi] = n_jpsi;
          Stat_Jpsi_HE[counter_cost][counter_phi] = stat_jpsi;
          counter_phi++;
        }

        N_Jpsi.push_back(n_jpsi);
        stat_Jpsi.push_back(stat_jpsi);
      }
  }

  fit_of_minv(hist_minv_integrated,100,100);

  int integral = 0;
  for(int i = 0;i < N_Jpsi.size();i++){
    integral += N_Jpsi[i];
  }
  cout << "SUM OF HISTOS = " << integral << endl;

  //============================================================================
  //MATRIX OF THE SIGNAL
  //============================================================================

  printf("MATRIX OF N_Jpsi \n");
  for(int i = 0;i < N_cost_bins;i++){
    for(int j = 0;j < N_phi_bins;j++){
      cout << N_Jpsi_HE[i][j] << ",";
    }
    cout << endl;
  }


  printf("MATRIX OF STAT ERRORS \n");
  for(int i = 0;i < N_cost_bins;i++){
    for(int j = 0;j < N_phi_bins;j++){
      cout << Stat_Jpsi_HE[i][j] << ",";
    }
    cout << endl;
  }

  //============================================================================
  // SAVE RESULTS INTO A TREE
  //============================================================================

  if(save_tree){
    TTree *output_tree = new TTree("CB2_VWG","CB2_VWG");
    output_tree -> Branch("N_Jpsi_HE",N_Jpsi_HE,"N_Jpsi_HE[18][10]/I");
    output_tree -> Branch("Stat_Jpsi_HE",Stat_Jpsi_HE,"Stat_Jpsi_HE[18][10]/I");
    output_tree -> Fill();

    TFile *output_file = new TFile("/Users/Luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/N_Jpsi.root","RECREATE");
    output_tree -> Write();
    output_file -> Close();
    printf("The file is saved in /Users/Luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/ \n");
  }

}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void fit_of_minv(TH1D *hist_minv, int counter_cost, int counter_phi){

  char *sig = "CB2";
  char *bck = "VWG";
  //int min_centr_class = 2;
  //int max_centr_class = 90;
  char *fit_range = "min";
  double min_fit_range = 2;
  double max_fit_range = 5;
  char *tails = "pp8";
  char *tails_fix = "yes";

  gStyle -> SetOptStat(0);
  TGaxis::SetMaxDigits(2);

  //============================================================================
  //FIT OF THE BACKGROUND
  //============================================================================

  TH1D *histo_for_bck = (TH1D*) hist_minv -> Clone("histo_for_bck");
  double m_width = histo_for_bck -> GetBinWidth(1);
  double m_min = 2.9/m_width;
  double m_max = 3.3/m_width;

  cout << "-----------------------------------" << endl;
  cout << "Bin width : " << m_width << endl;
  cout << "bin min : " << m_min << "; bin max : " << m_max << endl;
  cout << "-----------------------------------" << endl;

  for(int i = m_min; i < m_max;i++){
    histo_for_bck -> SetBinContent(i+1,0);
    histo_for_bck -> SetBinError(i+1,0);
  }

  histo_for_bck -> Draw();

  double par_bck[4] = {500000.,0.6,0.4,0.2};
  double *par_bck_ptr = par_bck;

  TF1 *func_bck_VWG = new TF1("func_bck_VWG",Func_VWG,2.2,5.5,4);
  func_bck_VWG -> SetParameters(par_bck_ptr);
  histo_for_bck -> Fit(func_bck_VWG,"R0");

  //============================================================================
  //TAILS {alpha1,n1,alpha2,n2}
  //============================================================================
  //double tails_par[4] = {1.012,2.874,2.232,2.924}; //p-A MC COLLISIONS
  //double tails_par[4] = {1.089,3.393,2.097,8.694}; //p-p 8 TeV
  //double tails_par[4] = {0.98,6.97,1.86,14.99}; //p-p 13 TeV
  double tails_par[4] = {1.06,3.23,2.55,1.56}; //GEANT 4
  //double tails_par[4] = {0.97,3.98,2.3,3.03}; //GEANT 3
  //============================================================================
  //FIT OF THE JPSI SIGNAL
  //============================================================================

  double par_signal[8] = {50000.,3.096,7.0e-02,1.00011e+00,3.70078e+00,1.68359e+00,3.63002e+00,0.01};

  //TF1 *func_Jpsi_CB2 = new TF1("func_Jpsi_CB2",Func_Jpsi_CB2,2.9,3.3,12);
  TF1 *func_Jpsi_CB2 = new TF1("func_Jpsi_CB2",Func_Jpsi_CB2,2.9,3.3,11);
  func_Jpsi_CB2 -> SetParameter(4,par_signal[0]);
  func_Jpsi_CB2 -> SetParameter(5,par_signal[1]);
  func_Jpsi_CB2 -> FixParameter(6,par_signal[2]);
  func_Jpsi_CB2 -> FixParameter(7,tails_par[0]);
  func_Jpsi_CB2 -> FixParameter(8,tails_par[1]);
  func_Jpsi_CB2 -> FixParameter(9,tails_par[2]);
  func_Jpsi_CB2 -> FixParameter(10,tails_par[3]);
  //func_Jpsi_CB2 -> FixParameter(11,par_signal[7]);
  hist_minv -> Fit(func_Jpsi_CB2,"R0");

  //============================================================================
  //FIT OF THE TOTAL SPECTRUM
  //============================================================================
  //TF1 *func_tot = new TF1("func_tot",Func_tot,min_fit_range,max_fit_range,12);
  TF1 *func_tot = new TF1("func_tot",Func_tot,min_fit_range,max_fit_range,11);
  TFitResultPtr fit_ptr;
  for(int i = 0;i < 100;i++){
    if(i == 0){
      func_tot -> SetParameter(0,func_bck_VWG -> GetParameter(0));
      func_tot -> SetParameter(1,func_bck_VWG -> GetParameter(1));
      func_tot -> SetParameter(2,func_bck_VWG -> GetParameter(2));
      func_tot -> SetParameter(3,func_bck_VWG -> GetParameter(3));
      func_tot -> SetParameter(4,func_Jpsi_CB2 -> GetParameter(4));
      func_tot -> SetParLimits(4,0,10000000);
      func_tot -> SetParameter(5,3.096);
      func_tot -> SetParameter(6,7.0e-02);
      func_tot -> SetParLimits(6,6.0e-02,9.0e-02);

      if(strcmp(tails_fix,"yes")==0){
        func_tot -> FixParameter(7,func_Jpsi_CB2 -> GetParameter(7));
        func_tot -> FixParameter(8,func_Jpsi_CB2 -> GetParameter(8));
        func_tot -> FixParameter(9,func_Jpsi_CB2 -> GetParameter(9));
        func_tot -> FixParameter(10,func_Jpsi_CB2 -> GetParameter(10));
      }
      if(strcmp(tails_fix,"no")==0){
        func_tot -> SetParameter(7,func_Jpsi_CB2 -> GetParameter(7));
        func_tot -> SetParameter(8,func_Jpsi_CB2 -> GetParameter(8));
        func_tot -> SetParameter(9,func_Jpsi_CB2 -> GetParameter(9));
        func_tot -> SetParameter(10,func_Jpsi_CB2 -> GetParameter(10));
      }
      //func_tot -> SetParameter(11,func_Jpsi_CB2 -> GetParameter(11));
    }
    else{func_tot -> SetParameters(func_tot -> GetParameters());}
    fit_ptr = (TFitResultPtr) hist_minv -> Fit(func_tot,"RLS0");
    if(gMinuit->fCstatu.Contains("CONVERGED")) break;
  }

  double ChiSquare_NDF = func_tot -> GetChisquare()/func_tot -> GetNDF();

  //printf("\n\nfit status: %s \n\n",gMinuit.fCstatu.Data());
  if(gMinuit->fCstatu.Contains("FAILED")){
    cout << "WARNING : FIT STATUS FAILED" << endl;
    sprintf(fit_status,"FAILED");
    return;
  }
  else{sprintf(fit_status,"SUCCESS");}

  TMatrixDSym cov = fit_ptr -> GetCovarianceMatrix();

  double *fullmat;
  fullmat = cov.GetMatrixArray();

  //PSI(2S) MATRIX
  double Jpsi_mat[49];
  for(Int_t i = 0;i < 7;i++){
    for(Int_t j = 0;j < 7;j++){
        Jpsi_mat[7*i+j] = fullmat[48+j+11*i];
    }
  }

  double Jpsi_par[7];
  for(int i = 0;i < 7;i++){
    Jpsi_par[i] = func_tot -> GetParameter(4+i);
  }

  //PSI(2S) MATRIX
  //double Psi2s_mat[64];
  //for(Int_t i = 0;i < 8;i++){
    //for(Int_t j = 0;j < 8;j++){
        //Psi2s_mat[8*i+j] = fullmat[52+j+12*i];
    //}
  //}

  //double Psi2s_par[8];
  //for(int i = 0;i < 8;i++){
    //Psi2s_par[i] = func_tot -> GetParameter(4+i);
  //}

  //============================================================================
  //PLOT OF JPSI AND PSI(2S) SHAPES
  //============================================================================

  TF1 *func_bck_VWG_fix = new TF1("func_bck_VWG_fix",Func_VWG,2.,5.,4);
  func_bck_VWG_fix -> SetParameter(0,func_tot -> GetParameter(0));
  func_bck_VWG_fix -> SetParameter(1,func_tot -> GetParameter(1));
  func_bck_VWG_fix -> SetParameter(2,func_tot -> GetParameter(2));
  func_bck_VWG_fix -> SetParameter(3,func_tot -> GetParameter(3));
  func_bck_VWG_fix -> SetLineStyle(4);
  func_bck_VWG_fix -> SetLineColor(kBlue+1);
  func_bck_VWG_fix -> Draw("same");

  TF1 *func_Jpsi_CB2_fix = new TF1("func_Jpsi_CB2_fix",Func_Jpsi_CB2_fix,2.,5.,7);
  func_Jpsi_CB2_fix -> SetParameter(0,func_tot -> GetParameter(4));
  func_Jpsi_CB2_fix -> SetParameter(1,func_tot -> GetParameter(5));
  func_Jpsi_CB2_fix -> SetParameter(2,func_tot -> GetParameter(6));
  func_Jpsi_CB2_fix -> SetParameter(3,func_tot -> GetParameter(7));
  func_Jpsi_CB2_fix -> SetParameter(4,func_tot -> GetParameter(8));
  func_Jpsi_CB2_fix -> SetParameter(5,func_tot -> GetParameter(9));
  func_Jpsi_CB2_fix -> SetParameter(6,func_tot -> GetParameter(10));
  func_Jpsi_CB2_fix -> SetLineStyle(2);
  func_Jpsi_CB2_fix -> Draw("same");

  //TF1 *func_Psi2s_CB2_fix = new TF1("func_Psi2s_CB2_fix",Func_Psi2s_CB2_fix,2.,5.,8);
  //func_Psi2s_CB2_fix -> SetParameter(0,func_tot -> GetParameter(4));
  //func_Psi2s_CB2_fix -> SetParameter(1,func_tot -> GetParameter(5));
  //func_Psi2s_CB2_fix -> SetParameter(2,func_tot -> GetParameter(6));
  //func_Psi2s_CB2_fix -> SetParameter(3,func_tot -> GetParameter(7));
  //func_Psi2s_CB2_fix -> SetParameter(4,func_tot -> GetParameter(8));
  //func_Psi2s_CB2_fix -> SetParameter(5,func_tot -> GetParameter(9));
  //func_Psi2s_CB2_fix -> SetParameter(6,func_tot -> GetParameter(10));
  //func_Psi2s_CB2_fix -> SetParameter(7,func_tot -> GetParameter(11));
  //func_Psi2s_CB2_fix -> SetLineStyle(2);
  //func_Psi2s_CB2_fix -> Draw("same");

  //////////////////////////////////////////////////////////////////////////////
  //SIGNAL/BACKGROUND
  //////////////////////////////////////////////////////////////////////////////
  double sigma_min_Jpsi = func_tot -> GetParameter(5) - 3*(func_tot -> GetParameter(6));
  double sigma_max_Jpsi = func_tot -> GetParameter(5) + 3*(func_tot -> GetParameter(6));
  double N_Jpsi_3sigma = func_Jpsi_CB2_fix -> Integral(sigma_min_Jpsi,sigma_max_Jpsi)/m_width;
  double N_bck_Jpsi_3sigma = func_bck_VWG_fix -> Integral(sigma_min_Jpsi,sigma_max_Jpsi)/m_width;
  double SB_Jpsi = N_Jpsi_3sigma/N_bck_Jpsi_3sigma;
  double sigma_Jpsi = func_tot -> GetParameter(6);
  //double sigma_min_Psi2s = func_tot -> GetParameter(5) + (3.686-3.097) - 3*(func_tot -> GetParameter(6)*scaling_factor);
  //double sigma_max_Psi2s = func_tot -> GetParameter(5) + (3.686-3.097) + 3*(func_tot -> GetParameter(6)*scaling_factor);
  //double N_Psi2s_3sigma = func_Psi2s_CB2_fix -> Integral(sigma_min_Psi2s,sigma_max_Psi2s)/m_width;
  //double N_bck_Psi2s_3sigma = func_bck_VWG_fix -> Integral(sigma_min_Psi2s,sigma_max_Psi2s)/m_width;
  //double SB_Psi2s = N_Psi2s_3sigma/N_bck_Psi2s_3sigma;
  //////////////////////////////////////////////////////////////////////////////

  //============================================================================
  //NUMERICAL RESULTS
  //============================================================================

  //n_psi2s = func_Psi2s_CB2_fix -> Integral(0,5)/m_width;
  //stat_psi2s = func_Psi2s_CB2_fix -> IntegralError(0.,5.,Psi2s_par,Psi2s_mat)/m_width;
  n_jpsi = func_Jpsi_CB2_fix -> Integral(0,5)/m_width;
  stat_jpsi = func_Jpsi_CB2_fix -> IntegralError(0.,5.,Jpsi_par,Jpsi_mat)/m_width;

  //cout << "N Psi(2S) = " << n_psi2s << " +- " << stat_psi2s << endl;
  cout << "N Jpsi = " << n_jpsi << " +- " << stat_jpsi << endl;
  cout << "S/B Jpsi = " << SB_Jpsi << endl;

  //============================================================================
  //PLOT OF THE TOTAL SPECTRUM
  //============================================================================
  char title[100];

  double max_histo_value = func_bck_VWG_fix -> GetMaximum();
  max_histo_value = max_histo_value + 0.4*max_histo_value;

  TH2D *h_spectrum = new TH2D("h_spectrum","",120,2,5,100,0,max_histo_value);
  h_spectrum -> GetXaxis() -> SetTitle("#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})");
  h_spectrum -> GetXaxis() -> SetTitleSize(0.05);
  h_spectrum -> GetXaxis() -> SetTitleOffset(0.95);
  h_spectrum -> GetXaxis() -> SetLabelSize(0.05);
  h_spectrum -> GetYaxis() -> SetTitle("Counts per 25 MeV/#it{c}^{2}");
  h_spectrum -> GetYaxis() -> CenterTitle(true);
  h_spectrum -> GetYaxis() -> SetTitleSize(0.05);
  h_spectrum -> GetYaxis() -> SetTitleOffset(1.2);
  h_spectrum -> GetYaxis() -> SetLabelSize(0.05);

  TLatex *lat0 = new TLatex(0.45,0.82,"ALICE Performance");
  lat0 -> SetTextSize(0.05);
  lat0 -> SetNDC();
  lat0 -> SetTextFont(42);

  TLatex *lat1 = new TLatex(0.45,0.76,"Pb-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  lat1 -> SetTextSize(0.05);
  lat1 -> SetNDC();
  lat1 -> SetTextFont(42);

  TLatex *lat2 = new TLatex(0.45,0.70,"Inclusive J/#psi #rightarrow #mu^{#plus}#mu^{#minus}");
  lat2 -> SetTextSize(0.05);
  lat2 -> SetNDC();
  lat2 -> SetTextFont(42);

  sprintf(title,"%i < #it{p}_{T} < %i GeV/#it{c}",2,6);
  TLatex *lat3 = new TLatex(0.49,0.62,title);
  lat3 -> SetTextSize(0.05);
  lat3 -> SetNDC();
  lat3 -> SetTextFont(42);

  sprintf(title,"%2.1f < #it{y} < %i",2.5,4);
  TLatex *lat4 = new TLatex(0.49,0.55,title);
  lat4 -> SetTextSize(0.05);
  lat4 -> SetNDC();
  lat4 -> SetTextFont(42);

  sprintf(title,"%2.1f < cos#it{#theta}^{HX} < %2.1f",value_cost[counter_cost],value_cost[counter_cost+1]);
  if(counter_cost == 100) sprintf(title,"%2.1f < cos#it{#theta}^{HX} < %2.1f",-1.,1.);
  TLatex *lat5 = new TLatex(0.49,0.48,title);
  lat5 -> SetTextSize(0.05);
  lat5 -> SetNDC();
  lat5 -> SetTextFont(42);

  sprintf(title,"%3.2f < #it{#varphi}^{HX} < %3.2f rad",value_phi[counter_phi],value_phi[counter_phi+1]);
  if(counter_phi == 100) sprintf(title,"%3.2f < #it{#varphi}^{HX} < %3.2f rad",0.,PI);
  TLatex *lat6 = new TLatex(0.49,0.41,title);
  lat6 -> SetTextSize(0.05);
  lat6 -> SetNDC();
  lat6 -> SetTextFont(42);

  sprintf(title,"#chi^{2}/ndf = %3.1f",ChiSquare_NDF);
  TLatex *lat7 = new TLatex(0.55,0.34,title);
  lat7 -> SetTextSize(0.05);
  lat7 -> SetNDC();
  lat7 -> SetTextFont(42);

  TCanvas *c_spectrum = new TCanvas("c_spectrum","c_spectrum",65,73,900,806);
  c_spectrum -> Range(1.825,-6776.052,5.019444,37862.12);
  c_spectrum -> SetFillColor(0);
  c_spectrum -> SetBorderMode(0);
  c_spectrum -> SetBorderSize(0);
  c_spectrum -> SetTickx(1);
  c_spectrum -> SetTicky(1);
  c_spectrum -> SetLeftMargin(0.18);
  c_spectrum -> SetBottomMargin(0.1518219);
  c_spectrum -> SetFrameBorderMode(0);
  c_spectrum -> SetFrameBorderMode(0);

  h_spectrum -> Draw();
  hist_minv -> Draw("sameE");
  func_bck_VWG_fix -> Draw("same");
  func_tot -> Draw("same");
  func_Jpsi_CB2_fix -> Draw("same");
  lat0 -> Draw();
  lat1 -> Draw();
  lat2 -> Draw();
  lat3 -> Draw();
  lat4 -> Draw();
  lat5 -> Draw();
  lat6 -> Draw();
  lat7 -> Draw();
  //lat6 -> Draw();
  //func_Psi2s_CB2_fix -> Draw("same");
  //t_spectrum -> Draw();

  //char hist_name[30];
  //sprintf(hist_name,"FIT_PLOTS/%s.png",hist_minv -> GetName());
  //c_spectrum -> SaveAs(hist_name);
  //delete c_spectrum;
}
//==============================================================================
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//==============================================================================
Double_t Func_VWG(Double_t *x, Double_t *par){
  Double_t sigma = par[2] + par[3]*((x[0] - par[1])/par[1]);
  Double_t FitBck = par[0]*TMath::Exp(-(x[0] - par[1])*(x[0] - par[1])/(2.*sigma*sigma));
  return FitBck;
}
//==============================================================================
Double_t Func_Jpsi_CB2(Double_t *x, Double_t *par){

  Double_t t = (x[0] - par[5])/par[6];
  if (par[7] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[7]);
  Double_t absAlpha2 = fabs((Double_t)par[9]);

  if (t >= -absAlpha && t < absAlpha2) // gaussian core

  {
    return par[4]*(exp(-0.5*t*t));
  }

  if (t < -absAlpha) //left tail

  {
    Double_t a =  TMath::Power(par[8]/absAlpha,par[8])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[8]/absAlpha - absAlpha;

    return par[4]*(a/TMath::Power(b - t, par[8]));
  }

  if (t >= absAlpha2) //right tail

  {

   Double_t c =  TMath::Power(par[10]/absAlpha2,par[10])*exp(-0.5*absAlpha2*absAlpha2);
   Double_t d = par[10]/absAlpha2 - absAlpha2;

  return  par[4]*(c/TMath::Power(d + t, par[10]));
  }

  return 0. ;
}
//==============================================================================
Double_t Func_Jpsi_CB2_fix(Double_t *x, Double_t *par){
  Double_t t = (x[0]-par[1])/par[2];
  if (par[3] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);

  if (t >= -absAlpha && t < absAlpha2) // gaussian core

  {
    return par[0]*(exp(-0.5*t*t));
  }

  if (t < -absAlpha) //left tail

  {
    Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[4]/absAlpha - absAlpha;

    return par[0]*(a/TMath::Power(b - t, par[4]));
  }

  if (t >= absAlpha2) //right tail

  {

   Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
   Double_t d = par[6]/absAlpha2 - absAlpha2;

  return  par[0]*(c/TMath::Power(d + t, par[6]));
  }

  return 0. ;
}
//==============================================================================
Double_t Func_Psi2s_CB2(Double_t *x, Double_t *par){

  Double_t t = (x[0]-(par[5]+(3.686-3.097)))/(par[6]*scaling_factor);
  if (par[7] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[7]);
  Double_t absAlpha2 = fabs((Double_t)par[9]);

  if (t >= -absAlpha && t < absAlpha2) // gaussian core

  {
    return par[4]*par[11]*(exp(-0.5*t*t));
  }

  if (t < -absAlpha) //left tail

  {
    Double_t a =  TMath::Power(par[8]/absAlpha,par[8])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[8]/absAlpha - absAlpha;

    return par[4]*par[11]*(a/TMath::Power(b - t, par[8]));
  }

  if (t >= absAlpha2) //right tail

  {

   Double_t c =  TMath::Power(par[10]/absAlpha2,par[10])*exp(-0.5*absAlpha2*absAlpha2);
   Double_t d = par[10]/absAlpha2 - absAlpha2;
  return  par[4]*par[11]*(c/TMath::Power(d + t, par[10]));
  }

  return 0. ;
}
//==============================================================================
Double_t Func_Psi2s_CB2_fix(Double_t *x, Double_t *par){

  Double_t t = (x[0]-(par[1]+(3.686-3.097)))/(par[2]*scaling_factor);
  if (par[3] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);

  if (t >= -absAlpha && t < absAlpha2) // gaussian core

  {
    return par[0]*par[7]*(exp(-0.5*t*t));
  }

  if (t < -absAlpha) //left tail

  {
    Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[4]/absAlpha - absAlpha;

    return par[0]*par[7]*(a/TMath::Power(b - t, par[4]));
  }

  if (t >= absAlpha2) //right tail

  {

   Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
   Double_t d = par[6]/absAlpha2 - absAlpha2;

  return  par[0]*par[7]*(c/TMath::Power(d + t, par[6]));
  }

  return 0. ;
}
//==============================================================================
Double_t Func_tot(Double_t *x, Double_t *par){
  //return Func_VWG(x,par) + Func_Jpsi_CB2(x,par) + Func_Psi2s_CB2(x,par);
  return Func_VWG(x,par) + Func_Jpsi_CB2(x,par);
}
