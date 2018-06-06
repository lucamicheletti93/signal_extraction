#ifndef settings_h
#define settings_h

#include <TMath.h>
#include <TLine.h>

//char hist_mass_name[50];
double PI = TMath::Pi();

//==============================================================================
// VARIABLES pT
//==============================================================================

const int N_pt_bin = 3;
int min_pt_bin[N_pt_bin] = {0,2,6};
int max_pt_bin[N_pt_bin] = {2,6,12};

//==============================================================================
// VARIABLES cost
//==============================================================================

//------------------------------------------------------------------------------
// RANGE -PI < phi < PI ; -1 < cost < 1
//------------------------------------------------------------------------------
//const int N_cost_bins = 18;
const int N_cost_bins = 22;
const int N_TH3_cost_bins = 100;

//int min_cost_bin[N_cost_bins] = {1,26,31,36,41,43,45,47,49,51,53,55,57,59,61,66,71,76};
//int max_cost_bin[N_cost_bins] = {25,30,35,40,42,44,46,48,50,52,54,56,58,60,65,70,75,100};
//int min_cost_bin[N_cost_bins] = {1,19,26,31,36,41,43,45,47,49,51,53,55,57,59,61,66,71,76,81};
//int max_cost_bin[N_cost_bins] = {18,25,30,35,40,42,44,46,48,50,52,54,56,58,60,65,70,75,82,100};
int min_cost_bin[N_cost_bins] = {1,11,21,26,31,36,41,43,45,47,49,51,53,55,57,59,61,66,71,76,81,91};
int max_cost_bin[N_cost_bins] = {10,20,25,30,35,40,42,44,46,48,50,52,54,56,58,60,65,70,75,80,90,100};

const int dim_cost = N_cost_bins + 1;
//double bins_cost[dim_cost] = {0,25,30,35,40,42,44,46,48,50,52,54,56,58,60,65,70,75,100};
//double bins_cost[dim_cost] = {0,18,25,30,35,40,42,44,46,48,50,52,54,56,58,60,65,70,75,82,100};
double bins_cost[dim_cost] = {0,10,20,25,30,35,40,42,44,46,48,50,52,54,56,58,60,65,70,75,80,90,100};
double width_cost[dim_cost];

double min_cost = -1.;
double max_cost = 1.;

double value_cost[dim_cost];
double Dcost[N_cost_bins];

//==============================================================================
// VARIABLES phi
//==============================================================================

//------------------------------------------------------------------------------
// RANGE -PI < phi < PI ; -1 < cost < 1
//------------------------------------------------------------------------------
const int N_phi_bins = 10;
const int N_TH3_phi_bins = 50;
//const int N_TH3_phi_bins = 100; // double binning

int min_phi_bin[N_phi_bins] = {1,16,19,21,23,26,29,31,33,36};
int max_phi_bin[N_phi_bins] = {15,18,20,22,25,28,30,32,35,50};
//int min_phi_bin[N_phi_bins] = {1,31,37,41,45,51,57,61,65,71}; // double binning
//int max_phi_bin[N_phi_bins] = {30,36,40,44,50,56,60,64,70,100}; // double binning

const int dim_phi = N_phi_bins + 1;
double bins_phi[dim_phi] = {0,15,18,20,22,25,28,30,32,35,50};
//double bins_phi[dim_phi] = {0,30,36,40,44,50,56,60,64,70,100}; // double binning
double width_phi[dim_phi];

double min_phi = 0;
double max_phi = PI;

double value_phi[dim_phi];
double Dphi[N_phi_bins];

//==============================================================================
// COMPUTE THE BIN AREA
//==============================================================================

double bin_area[N_cost_bins][N_phi_bins];

//==============================================================================
// PLOT THE LINES
//==============================================================================

const int N_line_cost = N_cost_bins - 1;
TLine *line_cost[N_line_cost];
const int N_line_phi = N_phi_bins - 1;
TLine *line_phi[N_line_phi];

#endif
