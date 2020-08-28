/* * * * * * * * * * * * * 
 * DoubleScinti_Ref2ToF.hh  *
 *                       *
 *	T. Fujiwara          *
 *	2020. 08. 25 (Wed)   *
 * * * * * * * * * * * * */
#ifndef DoubleScinti_Ref2ToF_h
#define DoubleScinti_Ref2ToF_h 1

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
#include "TPaveStats.h"
#include "TPaveText.h"

#include "../Setting.h"

#define NofCanv        9
#define NofMPPC        4
#define NofScinti      3
#define NofRef         2
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
                                    {"Ref-2.QDC", "Ref-2.TDC", "Ref-2"},
                                    {"ToF-1.QDC", "ToF-1.TDC", "ToF-1"},
                                    {"ToF-2.QDC", "ToF-2.TDC", "ToF-2"} };

int const ChID[3] = {1, 2, 3};	//Ref-2, ToF-1, ToF-2

//double const Reso_Ref_Val = 43.59;
//double const Reso_Ref_Err = 0.19;
//double const Reso_Ref_Val = 34.24;
//double const Reso_Ref_Err = 0.14;

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
class DoubleScinti_Ref2ToF
{
	public:
		 DoubleScinti_Ref2ToF();
		~DoubleScinti_Ref2ToF();
		void TwlkMan(int Type);
		void SetRoot( int RNum );
		void ImportPedestal( int PNum );
		void DefineCanv();
		void DefineHistFunc();
		void FillHist();
		void GetQDCCut();
		void TwlkFit();
		void TwlkSerchBestFactorToF1vsToF2();
		void TwlkSerchBestFactorRef2vsToF1();
		void TwlkSerchBestFactorRef2vsToF2();
		void TwlkDecToF();
		void DrawHist();
		void SaveHist(int NSave);
		void ExportTwlk();
		void ReadTwlk(int R, int I);
		void ReconstructTOF(int R, int I);
		double Twlk_Correction( double x, double* parameter, double Coeff );
		double Twlk_Correction_ToF( double x, double* parameter, double Coeff, int Mode );
		
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
		TH2D *h_2D_RawTQ[3];
		TH2D *h_2D_FitTQ[3];
		TH1D *h_2D_FitSlices[3];
		TH2D *h_2D_DecTQ[3];
		TH2D *h_2D_CutTQ[3];
		TH2D *h_2D_Factor[3];			//0: ToF-1 vs. ToF-2, 1: Ref-2 vs. ToF-1, 2: Ref-2 vs. ToF-2
		TH1D *h_1D_DecToF_tmp;
		TF1 *f_gaus;
		TF1 *f_Twlk[3];
		TGraphErrors *Gr_fit[2][3];
		TGraph *Gr_SigBest[3];
		TF1 *f_ToF;
		TF1 *f_dTof;
		TLatex *Lat;
		TLine *Ln_Qpeak[3];
		TLine *Ln_Qcut[3][2];
		TPaveText *TwlkFitStat[3];
		TPaveText *DecTOFStat;
		
		int RunNum;
		int PedNum;
		string RunName_Full;
		string RunName_str;
		TString RunName_tst;
		string Param_Ped;
		string Param_QDC;
		string Param_QCut;
		string Param_Twlk;
		string Param_Input;
		string Param_Read;
		TString figure;
		TString figure_open, figure_close;
		TString TextBestSig[3];
		
		int Entry;
		int TDC[32];
		int QDC[32];

		int QDCBin[3];
		double QDCMax[3];

		double rawtof;
		double dectof;

		double Pede_Cnt_Val[NofMPPC];
		double Pede_Cnt_Err[NofMPPC];
		double Pede_Wid_Val[NofMPPC];
		double Pede_Wid_Err[NofMPPC];
	
		double TDCcut[NofMPPC][2];

		double Qpeak[3];
		double QCut[3][2];
		double cutfactMin[NofMPPC];
		double cutfact[NofMPPC];

		int NDiv[3];
		double Div_min;
		double Div_max;
		double xCenter;
		double TwlkFitNDF;
		int NFitPoint;

		int TwlkMode;
		TString ToF_TwlkType;

		double TwlkPar_In[3][4];
		double TwlkFitMax[NofMPPC];
		double TwlkFitMin[NofMPPC];
		double TwlkPar_Val[3][4];
		double TwlkPar_Err[3][4];
		double TwlkFactor[3];
		double Twlk_FactorBestSigma[3];
		double Twlk_FactorBestChisNdf[3];
		double Twlk_SigVal_tmp=0.;
		double Twlk_SigErr_tmp;
		double Twlk_SigVal_Best[3]={999., 999., 999.};
		double Twlk_SigErr_Best[3];
		double Twlk_Chis_tmp;
		double Twlk_Ndf_tmp;
		double Twlk_Chis_Best;
		double Twlk_Ndf_Best;
		double Twlk_ResoVal;
		double Twlk_ResoErr;
		double ToF_Adj;
		double par_tmp1=0., par_tmp2=0., par_tmp3=0.;
		int FactorSearchMax[3];

		double ResoToF_Val;
		double ResoToF_Err;

		int LColor = 602;
		int FColor = 422;	//kCyan-10

		double CanvHsize[NofCanv] = { 500, 500, 500, 500, 500, 500, 1000, 1000, 1000 };

		double HistHmax[3];
		double HistHmin[3];

		int NofSaveCanv;

		int ItNum;
		bool ReadFlag=false;
};

#endif 
