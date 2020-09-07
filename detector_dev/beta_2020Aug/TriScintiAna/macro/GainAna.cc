/* * * * * * * * * * * * * 
 * GainAna.cc            *
 *                       *
 *	T. Fujiwara          *
 *	2020. 08. 25 (Wed)   *
 * * * * * * * * * * * * */
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

#define NofMPPC   4
#define NofScinti 3 

/*-----------path------------*/
string const PATH_param_p = "/data/41a/ELS/fujiwara/ToFujiwara/beta_Aug2020/param/pede/";
string const PATH_param_q = "/data/41a/ELS/fujiwara/ToFujiwara/beta_Aug2020/param/qdc/";
string const PATH_param_t = "/data/41a/ELS/fujiwara/ToFujiwara/beta_Aug2020/param/twlk/";


TString const PATH_dat         = "/data/41a/ELS/fujiwara/ToFujiwara/beta_Aug2020/rootfiles/";
TString const PATH_fig_General = "/data/41a/ELS/fujiwara/ToFujiwara/beta_Aug2020/FIG/TriScintiAna";
TString const PATH_fig_l       = "/data/41a/ELS/fujiwara/ToFujiwara/beta_Aug2020/FIG/leng/";
TString const PATH_fig_q       = "/data/41a/ELS/fujiwara/ToFujiwara/beta_Aug2020/FIG/qdc/";
TString const PATH_fig_p       = "/data/41a/ELS/fujiwara/ToFujiwara/beta_Aug2020/FIG/pede/";
TString const PATH_fig_g       = "/data/41a/ELS/fujiwara/ToFujiwara/beta_Aug2020/FIG/gain/";


TString const Label[NofMPPC][3] =  { {"Ref-1.QDC", "Ref-1.TDC", "Ref-1"},
                                     {"Ref-2.QDC", "Ref-2.TDC", "Ref-2"},
                                     {"ToF-1.QDC", "ToF-1.TDC", "ToF-1"},
                                     {"ToF-2.QDC", "ToF-2.TDC", "ToF-2"}};

int const ChID[NofMPPC] = {0, 1, 2, 3};

double LandGausSimpson( double* x, double* par );

/////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
class GainAna
{
	public:
		 GainAna();
		~GainAna();
		void SetRoot(int RNum);
		void ImportPedestal();
		void DefCanv();
		void DefObj();
		void FillHist();
		void DrawHist();
		void FitHist();
		void Export();

	private:
		Setting *Set;
		TFile *ifp;
		TTree *tree;
		TCanvas *Ca;
		TH1D *h_QDC[2];
		TF1 *f_pre[2];
		TF1 *f_QDC[2];
		TFrame *Fr[2];
		TLatex *Lat;

		int ItNum;
		int RunNum;
		string RunName_Full;
		string RunName_str;
		TString RunName_tst;
		TString figure;
		TString figure_open, figure_close;
		string param_pede;
		string param_gain;
		
		int EvN;
		int TDC[32];
		int QDC[32];
		int Entry = 0;

		int NPar;
		TString ParName[4] = {"Area", "MPV", "#sigma_{Landau}", "#sigma_{gaus}"};

		bool Trig_Clock = false;
		bool Trig_beta  = false;
	
		int NofIt;

		int NBinQDC;
		double QDCMax;
		double par_tmp;
		double par_tmp1=0., par_tmp2=0.;
		double peak[2];
		double FitMin[2];
		double FitMax[2];
		double FitPar[4];
		double ParIni[2][4];
		double PedeCenter[NofMPPC];	
		double PedeWidth[NofMPPC]; 	
		double GainVal[2];
		double GainErr[2];
		double GainAve_Val;
		double GainAve_Err;
	
		int TDCLim[2] = {1900, 2300};	//0: Lower, 1: Upper

		TString Htitle[2];
};

////////////////////////////////////////////////////////////////////////////////////////////////////
GainAna::GainAna(){
	cout<<"__________ooooooooooOOOOO000OOOOOoooooooooo__________"<<endl;
	cout<<"=====     GainAna::Constructor called     ====="<<endl;
	ItNum = 0;

	Set = new Setting();
	Set -> Setting_Gene(1);
	gStyle->SetStatFormat("6.5g");
	gStyle->SetFitFormat("6.5g");
	gStyle->SetOptFit(111);

	NBinQDC = 300;
	QDCMax  = 150.;

	cout<<"__________ooooooooooOOOOO000OOOOOoooooooooo__________"<<endl;
}

//------------------------------------------------------------------------------------------------//
GainAna::~GainAna(){
	cout<<"__________ooooooooooOOOOO000OOOOOoooooooooo__________"<<endl;
	cout<<"=====     GainAna::Destructor called     ====="<<endl;
	cout<<"__________ooooooooooOOOOO000OOOOOoooooooooo__________"<<endl;
}

//------------------------------------------------------------------------------------------------//
void GainAna::SetRoot(int RNum){
	cout<<"=====     GainAna::SetRoot     ====="<<endl;
	RunNum = RNum;

	RunName_Full = Form("run%03d.root", RunNum);
	RunName_Full = PATH_dat+RunName_Full;
	RunName_str  = Form("run%03d.root", RunNum);
	RunName_tst  = Form("run%03d"     , RunNum);

	figure = Form("Gain_run%03d_%02d.pdf", RunNum, ItNum);
	figure = PATH_fig_g+figure;
	figure_open  = figure+"[";
	figure_close = figure+"]";
	
	param_pede = Form("Pedestal_run%03d_00.dat", RunNum);
	param_pede = PATH_param_p+param_pede; 
	
	param_gain = Form("Gain_run%03d_00.dat", RunNum);
	param_gain = PATH_param_q+param_gain; 

	ifp = new TFile( RunName_Full.c_str(), "READONLY" );
	tree = (TTree*)ifp -> Get("tree");
	Entry = tree->GetEntries();
	tree -> SetBranchAddress( "tdc", &TDC );
	tree -> SetBranchAddress( "qdc", &QDC );
	cout<<"======================================"<<endl;
	cout<<"    FILE IMPORT: "<<RunName_tst<<endl;
	cout<<"    Entry: "<<Entry<<endl;
	cout<<"======================================"<<endl;
}

//------------------------------------------------------------------------------------------------//
void GainAna::ImportPedestal(){
	cout<<"=====     GainAna::ImportPedestal     ====="<<endl;

	int j=0;
	double val1, val2;

	ifstream ifs( param_pede.c_str(), ios_base::in );

	//Exception
	if(!ifs){
		cout<<"!!!WARNING!!!   There is no such file!!!"<<endl;
		return;
	}else;

	while(ifs>>val1>>val2){
		PedeCenter[j] = val1;
		PedeWidth[j]  = val2;
		j++;
	}

	cout<<"======================================"<<endl;
	cout<<"    FILE IMPORT: "<<param_pede<<endl;
	cout<<"======================================"<<endl;
}

//------------------------------------------------------------------------------------------------//
void GainAna::DefCanv(){
	cout<<"=====     GainAna::DefCanv     ====="<<endl;
	Ca = new TCanvas( "Ca", "Ca", 1200, 1200);
	Ca -> Divide(1,2);
}

//------------------------------------------------------------------------------------------------//
void GainAna::DefObj(){
	cout<<"=====     GainAna::DefObj     ====="<<endl;
	for(int i=0; i<2; i++){
		Htitle[i] = Form("run%03d: ", RunNum)+Label[i+2][2];

		h_QDC[i] = new TH1D( Form("h_QDC[%d]", i), Form("h_QDC[%d]", i), NBinQDC, 0., QDCMax );
		Set ->  Setting_Hist1D( h_QDC[i], Htitle[i], Label[i+2][0]+" [pC]", "Counts/0.5 pC", 4, 1, 62, 870, 3003 );
		
		f_pre[i] = new TF1( Form("f_pre[%d]", i), "landau(0)", 0., QDCMax );
		Set -> Setting_Func( f_pre[i], 6, 1 );
		
		f_QDC[i] = new TF1( Form("f_QDC[%d]", i), LandGausSimpson, 0., QDCMax, 4 );
		NPar = f_QDC[i] -> GetNpar();
		for(int j=0; j<NPar; j++){
			f_QDC[i]->SetParName(j, ParName[j]);
		}	
		Set -> Setting_Func( f_QDC[i], 6, 1 );
	}

	Lat = new TLatex();
	Set -> Setting_Latex( Lat, 62, 22, 2, 0.05);
}

//------------------------------------------------------------------------------------------------//
void GainAna::FillHist(){
	cout<<"=====     GainAna::FillHist     ====="<<endl;
	for(int i=0; i<Entry; i++){
		tree->GetEntry(i);
		Trig_beta = (    TDC[1]>TDCLim[0]
		              && TDC[1]<TDCLim[1]
		              && TDC[2]>TDCLim[0]
		              && TDC[2]<TDCLim[1]
		              && TDC[3]>TDCLim[0]
		              && TDC[3]<TDCLim[1] );

		if(Trig_beta){
			for(int j=0; j<2; j++){h_QDC[j]->Fill( 0.1*QDC[j+2]-0.1*PedeCenter[j+2]);}
		}else;
	}

}

//------------------------------------------------------------------------------------------------//
void GainAna::DrawHist(){
	cout<<"=====     GainAna::DrawHist     ====="<<endl;
	Ca->cd(1);
	for(int i=0; i<2; i++){
		Ca -> cd(i+1);
		gPad -> SetBottomMargin(.150);
		h_QDC[i] -> GetXaxis() -> SetNdivisions(515);
		h_QDC[i] -> Draw();
	}
}

//------------------------------------------------------------------------------------------------//
void GainAna::FitHist(){
	cout<<"=====     GainAna::FitHist     ====="<<endl;

	for(int i=0; i<2; i++){
		peak[i]   = h_QDC[i]->GetBinCenter( h_QDC[i]->GetMaximumBin() );
	//	FitMin[i] = 0.60*peak[i];
	//	FitMax[i] = 1.60*peak[i];
	}
	FitMin[0] = 0.60*peak[0];
	FitMin[1] = 0.60*peak[1];
	FitMax[0] = 1.50*peak[0];	
	FitMax[1] = 1.50*peak[1];	
//	FitMin[1] = 0.40*peak[1];
//	FitMin[1] = 0.60*peak[1];	//run099_03
//	FitMax[1] = 1.50*peak[1];	//run099_01
//	FitMax[1] = 1.75*peak[1];	//run099_02
//	FitMax[1] = 1.50*peak[1];	//run108_01
	ParIni[0][2] = 3.0;
	ParIni[0][3] = 3.0;
	ParIni[1][2] = 3.0;
	ParIni[1][3] = 3.0;


	for(int i=0; i<2; i++){
		Ca -> cd(i+1);
	//	f_pre[i] -> SetParameter(0, h_QDC[i]->GetEffectiveEntries() );
	//	f_pre[i] -> SetParameter(1, h_QDC[i]->GetBinCenter( h_QDC[i]->GetMaximumBin() ) );
	//	f_pre[i] -> SetParameter(2, 0.5);
	//	for(int j=0; j<3; j++){h_QDC[i]->Fit(f_pre[i], "", "", 0.75*peak[i], 1.25*peak[i]);}
		f_QDC[i] -> SetParameter(0, h_QDC[i]->GetEffectiveEntries() );
		f_QDC[i] -> SetParameter(1, h_QDC[i]->GetBinCenter( h_QDC[i]->GetMaximumBin() ) );
		f_QDC[i] -> SetParameter(2, ParIni[i][2]);
		f_QDC[i] -> SetParameter(3, ParIni[i][3]);
		for(int j=0; j<5; j++){h_QDC[i]->Fit( f_QDC[i], "", "", FitMin[i], FitMax[i] );}
		
		par_tmp1 = f_QDC[i]->GetParameter(1);
		par_tmp2 = f_QDC[i]->GetParError(1);

		GainVal[i] = par_tmp1;
		GainErr[i] = sqrt( par_tmp2*par_tmp2 + (0.1*PedeWidth[i+2])*(0.1*PedeWidth[i+2]) );
	}

	GainAve_Val = sqrt( GainVal[0]*GainVal[1] );
	GainAve_Err = 0.5*sqrt( (GainVal[1]/GainVal[0])*(GainErr[0]*GainErr[0]) + (GainVal[0]/GainVal[1])*(GainErr[1]*GainErr[1]) );
}

//------------------------------------------------------------------------------------------------//
void GainAna::Export(){
	cout<<"=====     GainAna::Export     ====="<<endl;

	Ca -> Print( figure_open, "pdf" );
	Ca -> Print( figure );
	Ca -> Print( figure_close, "pdf" );

	ofstream ofs( param_gain.c_str(), ios_base::trunc );
	ofs<<GainVal[0]<<" "<<GainErr[0]<<endl;
	ofs<<GainVal[1]<<" "<<GainErr[1]<<endl;
	ofs<<GainAve_Val<<" "<<GainAve_Err<<endl;

	cout<<GainVal[0]<<" "<<GainErr[0]<<endl;
	cout<<GainVal[1]<<" "<<GainErr[1]<<endl;
	cout<<GainAve_Val<<" "<<GainAve_Err<<endl;
}

//------------------------------------------------------------------------------------------------//
double LandGausSimpson( double* x, double* par ){
	//Sekibun: Simpson method
	double land;
	double ret;
	// Numeric constants
	double pi = 4.*atan(1.);
	double invsqrt2pi = 1./sqrt(2.*pi);
	double mpshift  = -0.22278298;      // Landau maximum location
	// Control constants
	double np = 1000.0;					// number of convolution steps
	double sc = 10.0;					// convolution extends to +-sc Gaussian sigmas >> +/- sc*sigma_of_gaus
	// Variables
	double xx, xxx;
	double mpc;
	double FuncBoundary[2];				//0: Upper Limit, 1: Lower Limit
	double sum_1 = 0.0;
	double sum_2 = 0.0;
	double xLowLim,xUppLim;
	double step;
	double divsize;
	double factor;
	double i;

	// MP shift correction
	mpc = par[1] - mpshift*par[2];
//	mpc = par[1];

	xLowLim = x[0] - sc*par[3];
	xUppLim = x[0] + sc*par[3];

	FuncBoundary[0] = TMath::Landau( xLowLim, mpc, par[2] )*TMath::Gaus( x[0], xLowLim, par[3] );
	FuncBoundary[1] = TMath::Landau( xUppLim, mpc, par[2] )*TMath::Gaus( x[0], xUppLim, par[3] );

	step    = (xUppLim-xLowLim)/np;
	divsize = 1./np;
	factor  = step/6.;

	for(i=1.; i<=(np-1.); i++){
		xx = xLowLim + (i*divsize)*(2.*sc*par[3]);
		land = TMath::Landau( xx, mpc, par[2] );
		sum_1 += land*TMath::Gaus( x[0], xx, par[3] );
	}

	for(i=1.; i<=np; i++){
		xxx = xLowLim +  ( ((2*i-1)*divsize)/2. )*(2.*sc*par[3]);
		land = TMath::Landau( xxx, mpc, par[2] );
		sum_2 += land*TMath::Gaus( x[0], xxx, par[3] );
	}

	ret = par[0]*factor*(invsqrt2pi/par[3])*( FuncBoundary[0] + FuncBoundary[1] + 2.*sum_1 + 4.*sum_2 );

	return ret;
}

//------------------------------------------------------------------------------------------------//
////////////////////////////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv ){
	int Flag=1;	
	int RNum=0;

	TApplication *theApp;
	GainAna *ana;

	/*=====	Set ROOT file	=====*/
	cout<<"====="<<endl;
	if(argc==2){
		RNum = atoi(argv[1]);
		cout<<"Target of Analysis >>> "<<RNum<<endl;
	}else{
		RNum = 94;
		cout<<"Target of Analysis >>> "<<RNum<<endl;
	};


	theApp = new TApplication("App", &argc, argv );
	ana = new GainAna();

	ana -> SetRoot(RNum);
	ana -> ImportPedestal();
	ana -> DefCanv();
	ana -> DefObj();
	ana -> FillHist();
	ana -> DrawHist();
	ana -> FitHist();
	ana -> Export();

	//====	END of Job	====//
	if( Flag==1 ){
		gSystem->Exit(1);
	}else;

	delete ana;

	theApp->Run();
	return 0;
}
