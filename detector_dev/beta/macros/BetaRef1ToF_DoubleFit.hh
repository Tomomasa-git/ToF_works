/* * * * * * * * * * * * * 
 * BetaRef1ToF_DoubleFit.hh        *
 *                       *
 *	T. Fujiwara          *
 *	2020. 07. 24 (Fri)   *
 * * * * * * * * * * * * */
#ifndef BetaRef1ToF_DoubleFit_h
#define BetaRef1ToF_DoubleFit_h 1

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
using namespace std;

#include "TApplication.h"
#include "TApplicationImp.h"
#include "TCanvasImp.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TRint.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TObject.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TFrame.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TMath.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"

#include "./Setting.h"

#define NofCanv        14
#define NofMPPC        3
#define NofScinti      2
#define NofRef         1
#define NofCombination 3

/*-----------path------------*/
string const PATH_param_p = "/data/41a/ELS/fujiwara/ToFujiwara/beta_Aug2020/param/pede/";
string const PATH_param_q = "/data/41a/ELS/fujiwara/ToFujiwara/beta_Aug2020/param/qdc/";
string const PATH_param_t = "/data/41a/ELS/fujiwara/ToFujiwara/beta_Aug2020/param/twlk/";


TString const PATH_dat         = "/data/41a/ELS/fujiwara/ToFujiwara/beta_Aug2020/rootfiles/";
TString const PATH_fig_General = "/data/41a/ELS/fujiwara/ToFujiwara/beta_Aug2020/FIG/";
TString const PATH_fig_l       = "/data/41a/ELS/fujiwara/ToFujiwara/beta_Aug2020/FIG/leng/";
TString const PATH_fig_q       = "/data/41a/ELS/fujiwara/ToFujiwara/beta_Aug2020/FIG/qdc/";
TString const PATH_fig_p       = "/data/41a/ELS/fujiwara/ToFujiwara/beta_Aug2020/FIG/pede/";
TString const PATH_fig_g       = "/data/41a/ELS/fujiwara/ToFujiwara/beta_Aug2020/FIG/gain/";


TString const Label[NofMPPC][3] =  { {"Ref-1.QDC", "Ref-1.TDC", "Ref-1"},
                                    {"ToF-1.QDC", "ToF-1.TDC", "ToF-1"},
                                    {"ToF-2.QDC", "ToF-2.TDC", "ToF-2"} };

int const ChID[NofMPPC] = {0, 2, 3};	//Ref-1, ToF-1, ToF-2

//double const Reso_Ref_Val = 43.59;
//double const Reso_Ref_Err = 0.19;
double const Reso_Ref_Val = 34.24;
double const Reso_Ref_Err = 0.14;

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
class BetaRef1ToF_DoubleFit
{
	public:
		 BetaRef1ToF_DoubleFit();
		~BetaRef1ToF_DoubleFit();
		void SetRoot( int RNum );
		void ImportPedestal( int PNum );
		void DefineCanv();
		void DefineHistFunc();
		void FillHist();
		void GetQDCCut();
		void TwlkFit();
		void TwlkSerchBestFactorToF1vsToF2();
		void TwlkSerchBestFactorRef1vsToF1();
		void TwlkSerchBestFactorRef1vsToF2();
		void TwlkDecToF();
		void TwlkFit_Second();
		void TwlkSerchBestFactorToF1vsToF2_Second();
		void TwlkSerchBestFactorRef1vsToF1_Second();
		void TwlkSerchBestFactorRef1vsToF2_Second();
		void TwlkDecToF_Second();
		void DrawHist();
		void SaveHist(int NSave);
		void ExportTwlk();
		double Twlk_Correction( double x, double* parameter, double Coeff );
		
	private:
		Setting *Set;

		TFile *ifp;
		TTree *tree;
		TCanvas *Ca[NofCanv];
		TCanvas *Ca_dummy[6];
		TCanvas *canv_dummy_consgaus;
		TH1D *h_1D_QDCFull[NofCanv];
		TH1D *h_1D_QDCwCut[NofCanv];
		TH1D *h_1D_RawToF;
		TH1D *h_1D_DecToF;
		TH2D *h_2D_RawTQ[NofMPPC];
		TH2D *h_2D_FitTQ[NofMPPC];
		TH2D *h_2D_DecTQ[NofMPPC];
		TH2D *h_2D_FitTQ_Second[NofMPPC];
		TH2D *h_2D_DecTQ_Second[NofMPPC];
		TH2D *h_2D_CutTQ[NofMPPC];
		TH2D *h_2D_Factor[6];			//0: ToF-1 vs. ToF-2, 1: Ref-1 vs. ToF-1, 2: Ref-1 vs. ToF-2
		TH1D *h_1D_DecToF_tmp;
		TF1 *f_gaus;
		TF1 *f_Twlk[NofMPPC];
		TF1 *f_Twlk_Second[NofMPPC];
		TGraphErrors *Gr_fit[2][NofMPPC];
		TGraph *Gr_SigBest[3];
		TGraph *Gr_SigBest_Second[3];
		TF1 *f_ToF;
		TF1 *f_dTof;
		TLatex *Lat;
		TLine *Ln_Qpeak[NofMPPC];
		TLine *Ln_Qcut[NofMPPC][2];
		
		int RunNum;
		int PedNum;
		string RunName_Full;
		string RunName_str;
		TString RunName_tst;
		string Param_Ped;
		string Param_QDC;
		string Param_QCut;
		string Param_Twlk;
		TString figure;
		TString figure_open, figure_close;
		TString TextBestSig[NofMPPC];
		
		int Entry;
		int TDC[32];
		int QDC[32];

		double rawtof;
		double dectof;

		double Pede_Cnt_Val[NofMPPC];
		double Pede_Cnt_Err[NofMPPC];
		double Pede_Wid_Val[NofMPPC];
		double Pede_Wid_Err[NofMPPC];
	
		double Qpeak[NofMPPC];
		double QCut[NofMPPC][2];
		double cutfactMin[NofMPPC];
		double cutfact[NofMPPC];

		int NDiv[NofMPPC];
		double Div_min;
		double Div_max;
		double xCenter;
		double TwlkFitNDF;
		int NFitPoint;

		double TwlkFitMax[NofMPPC];
		double TwlkFitMin[NofMPPC];
		double TwlkPar_Val[NofMPPC][3];
		double TwlkPar_Err[NofMPPC][3];
		double TwlkFactor[NofMPPC];
		double Twlk_FactorBestSigma[NofMPPC];
		double Twlk_FactorBestChisNdf[NofMPPC];
		double Twlk_SigVal_tmp=0.;
		double Twlk_SigErr_tmp;
		double Twlk_SigVal_Best[3]={999., 999., 999.};
		double Twlk_SigErr_Best[3];
		double TwlkSecondPar_Val[NofMPPC][3];
		double TwlkSecondPar_Err[NofMPPC][3];
		double TwlkSecondFactor[NofMPPC];
		double TwlkSecond_FactorBestSigma[NofMPPC];
		double TwlkSecond_FactorBestChisNdf[NofMPPC];
		double TwlkSecond_SigVal_tmp=0.;
		double TwlkSecond_SigErr_tmp;
		double TwlkSecond_SigVal_Best[3]={999., 999., 999.};
		double TwlkSecond_SigErr_Best[3];
		double Twlk_Chis_tmp;
		double Twlk_Ndf_tmp;
		double Twlk_Chis_Best;
		double Twlk_Ndf_Best;
		double Twlk_ResoVal;
		double Twlk_ResoErr;
		double ToF_Adj;
		double par_tmp1=0., par_tmp2=0., par_tmp3=0.;

		double ResoToF_Val;
		double ResoToF_Err;

		int LColor = 602;
		int FColor = 422;	//kCyan-10

		double CanvHsize[NofCanv] = { 500, 500, 500, 500, 500, 500, 500, 500, 1000, 1000, 1000, 1000, 1000, 1000 };

		double HistHmax[NofMPPC];
		double HistHmin[NofMPPC];

		int NofSaveCanv;
};
#endif