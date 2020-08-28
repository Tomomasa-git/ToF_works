#ifndef TriScinti_Ref1Ref2_h
#define TriScinti_Ref1Ref2_h 1

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

#include "../Setting.h"

#define NofCanv   6
#define NofMPPC   4
#define NofScinti 3
#define NofRef    2

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
                                     {"Ref-2.QDC", "Ref-2.TDC", "Ref-2"},
                                     {"ToF-1.QDC", "ToF-1.TDC", "ToF-1"},
                                     {"ToF-2.QDC", "ToF-2.TDC", "ToF-2"} };

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
class TriScinti_Ref1Ref2
{
	public:
		 TriScinti_Ref1Ref2();
		~TriScinti_Ref1Ref2();
		void SetRoot( int RNum );
		void ImportRefPedestal( int PNum );
		void DefineCanv();
		void DefineHistFunc();
		void FillHist();
		void GetQDCCut();
		void TwlkFit();
		void TwlkSerchBestFactor();
		void TwlkDecToF();
		void DrawHist();
		void SaveHist(int NSave);
		double Twlk_Correction( double x, double* parameter, double Coeff );


	private:
		Setting *Set;

		TFile *ifp;
		TTree *tree;
		TCanvas *Ca_Ref[NofCanv];
		TCanvas *canv_dummy_consgaus;
		TH1D *h_1D_QDCFull[NofCanv];
		TH1D *h_1D_QDCCut[NofCanv];
		TH1D *h_1D_RawToF;
		TH1D *h_1D_DecToF;
		TH2D *h_2D_RawTQ[NofRef];
		TH2D *h_2D_FitTQ[NofRef];
		TH2D *h_2D_DecTQ[NofRef];
		TH2D *h_2D_CutTQ[NofRef];
		TH2D *h_2D_Factor;
		TH1D *h_1D_DecToF_tmp;
		TF1 *f_Twlk[NofRef];
		TGraph *Gr_SigBest;
		TF1 *f_ToF;
		TLatex *Lat;
		TLine *Ln_Qpeak[NofRef];
		TLine *Ln_Qcut[NofRef][2];

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
		
		int Entry;
		int TDC[32];
		int QDC[32];

		double rawtof;
		double dectof;

		double TDCcut[NofMPPC][2];	//0: Lower limit, 1: Upper limit

		double Pede_Cnt_Val[NofMPPC];
		double Pede_Cnt_Err[NofMPPC];
		double Pede_Wid_Val[NofMPPC];
		double Pede_Wid_Err[NofMPPC];
	
		double RefQDCGain[NofRef];
		double RefQDCCut[NofRef][2];

		double TwlkPar_In[3][4];
		double TwlkPar_Val[NofRef][3];
		double TwlkPar_Err[NofRef][3];
		double TwlkFitMax[NofRef];
		double TwlkFitMin[NofRef] = {5.  , 5.  };
		double TwlkFactor[NofRef];
		double Twlk_FactorBestSigma[NofRef];
		double Twlk_FactorBestChisNdf[NofRef];
		double Twlk_SigVal_tmp=0.;
		double Twlk_SigErr_tmp;
		double Twlk_SigVal_Best=999.;
		double Twlk_SigErr_Best;
		double Twlk_Chis_tmp;
		double Twlk_Ndf_tmp;
		double Twlk_Chis_Best;
		double Twlk_Ndf_Best;
		double Twlk_ResoVal;
		double Twlk_ResoErr;
		double ToF_Adj;
		double par_tmp1=0., par_tmp2=0.;

		int LColor = 602;
		int FColor = 422;	//kCyan-10

		double CanvHsize[NofCanv] = { 900., 1800., 1800., 1800., 900., 1800. };

		int ItNum;

};

#endif
