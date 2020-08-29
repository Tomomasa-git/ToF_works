/* * * * * * * * * * * * *
 * GetPede_TriScinti.cc  *
 * 2020. 08. 29 (Sat)    *
 * T. Fujiwara           *
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

#include "./Setting.h"


#define NofMPPC   4
#define NofScinti 3 


/*-----------path------------*/
string const PATH_param_p = "/home/fujiwara/wara_els/ToF/beta/param/pede/";
string const PATH_param_q = "/home/fujiwara/wara_els/ToF/beta/param/qdc/";
string const PATH_param_t = "/home/fujiwara/wara_els/ToF/beta/param/twlk/";


TString const PATH_dat         = "/home/fujiwara/wara_els/ToF/beta/root/";
TString const PATH_fig_General = "/home/fujiwara/wara_els/ToF/beta/FIG/";
TString const PATH_fig_l       = "/home/fujiwara/wara_els/ToF/beta/FIG/leng";
TString const PATH_fig_q       = "/home/fujiwara/wara_els/ToF/beta/FIG/qdc/";
TString const PATH_fig_p       = "/home/fujiwara/wara_els/ToF/beta/FIG/pede/";
TString const PATH_fig_g       = "/home/fujiwara/wara_els/ToF/beta/FIG/gain/";


TString const Label[NofMPPC][3] =  { {"Ref-1.QDC", "Ref-1.TDC", "Ref-1"},
                                     {"Ref-2.QDC", "Ref-2.TDC", "Ref-2"},
                                     {"ToF-1.QDC", "ToF-1.TDC", "ToF-1"},
                                     {"ToF-2.QDC", "ToF-2.TDC", "ToF-2"}};

int const ChID[NofMPPC] = {0, 1, 2, 3};

/////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
class PedeStability
{
	public:
		 PedeStability();
		~PedeStability();
		void SetRoot(int RNum);
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
		TCanvas *Ca[NofMPPC];
		TH1D *h_Pede[NofMPPC];
		TH2D *h_qdcev[NofMPPC];
		TH1D *h_sl_qdcev[NofMPPC];
		TF1 *f_Pede[NofMPPC];
		TF1 *f_Slices[NofMPPC];
		TFrame *Fr[NofMPPC];
		TFrame *Fr_Ev2D[NofMPPC];
		TLatex *Lat;

		int ItNum;
		int RunNum;
		string RunName_Full;
		string RunName_str;
		TString RunName_tst;
		TString figure;
		TString figure_open, figure_close;
		string param;
		
		int EvN;
		int TDC[32];
		int QDC[32];
		int Entry = 0;

		bool Trig_Clock;
	
		int NofIt;
		double par_tmp;
		double FitRange;
		double FitPar[3];
		double PedeCenter[NofMPPC][2];	//0: Value, 1: Error
		double PedeWidth[NofMPPC][2]; 	//0: Value, 1: Error

		int NBinEv;
		double EvMax;
		double EvMin = 0.;
		int NBinQDC = 250;
		double QDCMax = 300.;
		double QDCMin = 0.;
};
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
PedeStability::PedeStability(){
	cout<<"__________ooooooooooOOOOO000OOOOOoooooooooo__________"<<endl;
	cout<<"=====   PedeStability::Constructor is called    ====="<<endl;
	
	ItNum = 0;
	
	Set = new Setting();
	Set -> Setting_Gene(1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
PedeStability::~PedeStability(){
	cout<<"__________ooooooooooOOOOO000OOOOOoooooooooo__________"<<endl;
	cout<<"=====   PedeStability::Destructor is called    ====="<<endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void PedeStability::SetRoot(int RNum){
	cout<<"__________ooooooooooOOOOO000OOOOOoooooooooo__________"<<endl;
	cout<<"=====   PedeStability::SetRoot    ====="<<endl;
	RunNum = RNum;
	RunName_Full = Form("run%03d.root", RunNum);
	RunName_Full = PATH_dat+RunName_Full;
	RunName_str  = Form("run%03d.root", RunNum);
	RunName_tst  = Form("run%03d"     , RunNum);
	figure = Form("PedestalStability_run%03d_%02d.pdf", RunNum, ItNum);
	figure = PATH_fig_p+figure;
	figure_open  = figure+"[";
	figure_close = figure+"]";
	param = Form("Pedestal_run%03d_%02d.dat", RunNum, ItNum);
	param = PATH_param_p+param; 
	
	ifp = new TFile( RunName_Full.c_str(), "READONLY" );
	tree = (TTree*)ifp->Get("tree");
	tree -> SetBranchAddress( "evnum", &EvN );
	tree -> SetBranchAddress( "qdc"  , &QDC );
	tree -> SetBranchAddress( "tdc"  , &TDC );
	Entry = tree->GetEntries();
	cout<<"======================================"<<endl;
	cout<<"    FILE IMPORT: "<<RunName_str<<endl;
	cout<<"    Entry: "<<Entry<<endl;
	cout<<"======================================"<<endl;

	if(Entry<=25001){
		NBinEv = 125;
		EvMax  = 25000.;
	}else if(Entry<=40001){
		NBinEv = 200;
		EvMax  = 40000.;
	}else if(Entry<=50001){
		NBinEv = 250;
		EvMax  = 50000.;
	}else if(Entry<=60001){
		NBinEv = 300;
		EvMax  = 60000.;
	}else if(Entry<=75001){
		NBinEv = 375;
		EvMax  = 75000.;
	}else if(Entry<=80001){
		NBinEv = 400;
		EvMax  = 80000.;
	}else if(Entry<=100001){
		NBinEv = 500;
		EvMax  = 100000.;
	}else if(Entry<=125001){
		NBinEv = 625;
		EvMax  = 1250000.;
	}else if(Entry<=150001){
		NBinEv = 750;
		EvMax  = 150000.;
	}else if(Entry<=200001){
		NBinEv = 1000;
		EvMax  = 200000.;
	}else if(Entry<=250001){
		NBinEv = 1250;
		EvMax  = 250000.;
	}else{
		NBinEv = 400;
		EvMax  = 80000.;
	}
}

//------------------------------------------------------------------------------------------------//
void PedeStability::DefCanv(){
	cout<<"__________ooooooooooOOOOO000OOOOOoooooooooo__________"<<endl;
	cout<<"=====   PedeStability::DefCanv    ====="<<endl;
	for(int i=0; i<4; i++){
		Ca[i] = new TCanvas( Form("Ca[%d]", i), Form("Ca[%d]", i), 1600, 900);
		Ca[i] -> Divide(2,2);
	}
}

//------------------------------------------------------------------------------------------------//
void PedeStability::DefObj(){
	cout<<"__________ooooooooooOOOOO000OOOOOoooooooooo__________"<<endl;
	cout<<"=====   PedeStability::DefObj    ====="<<endl;

	for(int i=0; i<NofMPPC; i++){
		h_Pede[i] = new TH1D( Form("h_Pede[%d]",i), Form("h_Pede[%d]",i), 300, 0, 300 );
		Set->Setting_Hist1D( h_Pede[i], RunName_tst+": "+Label[i][2], Label[i][0]+" [ch/0.1 pC]", "Counts/ch", 602, 1, 62, kCyan-10, 1003 );
		
		f_Pede[i] = new TF1( Form("f_Pede[%d]",i), "gaus(0)", 0., 300. );
		Set -> Setting_Func( f_Pede[i], 6, 1 );
		
		h_qdcev[i] = new TH2D( Form("h_qdcev[%d]", i), Form("h_qdcev[%d]", i), NBinEv, EvMin, EvMax, NBinQDC, QDCMin, QDCMax );
		Set -> Setting_Hist2D( h_qdcev[i], RunName_tst+": "+Label[i][2], Label[i][0]+" [ch/0.1 pC]", "Counts/ch", "", 0. );

		f_Slices[i] = new TF1( Form("f_Slices[%d]", i), "gaus(0)", 0., 300 );
	}
	Lat = new TLatex();
	Set -> Setting_Latex( Lat, 62, 22, 6, 0.05 );
}

//------------------------------------------------------------------------------------------------//
void PedeStability::FillHist(){
	cout<<"__________ooooooooooOOOOO000OOOOOoooooooooo__________"<<endl;
	cout<<"=====   PedeStability::FillHist    ====="<<endl;
	for(int i=0; i<Entry; i++){
		tree->GetEntry(i);
		Trig_Clock = TDC[4]>0;
		if(Trig_Clock){
			for(int j=0; j<NofMPPC; j++){
				h_Pede[j]  -> Fill(QDC[j]);
				h_qdcev[j] -> Fill(EvN, QDC[j]);
			}
		}else;
	}
}

//------------------------------------------------------------------------------------------------//
void PedeStability::DrawHist(){
	cout<<"__________ooooooooooOOOOO000OOOOOoooooooooo__________"<<endl;
	cout<<"=====   PedeStability::DrawHist    ====="<<endl;
	for(int i=0; i<NofMPPC; i++){
		Ca[i]->cd(1);
		h_qdcev[i]->Draw("colz");
		Ca[i]->cd(2);
		gPad->SetLogy(1);
		h_Pede[i]->Draw();
		gPad->Update();
	}
}

//------------------------------------------------------------------------------------------------//
void PedeStability::FitHist(){
	cout<<"__________ooooooooooOOOOO000OOOOOoooooooooo__________"<<endl;
	cout<<"=====   PedeStability::FitHist    ====="<<endl;
	for(int i=0; i<NofMPPC; i++){
		Ca[i]->cd(2);
		par_tmp = h_Pede[i]->GetBinCenter( h_Pede[i]->GetMaximumBin() );
		for(int j=0; j<NofIt; j++){
			h_Pede[i] -> Fit( f_Pede[i], "", "", par_tmp-FitRange, par_tmp+FitRange );
			f_Pede[i]->GetParameters(FitPar);
			par_tmp = FitPar[1];
		}
		PedeCenter[i][0] = f_Pede[i]->GetParameter(1);
		PedeCenter[i][1] = f_Pede[i]->GetParError(1);
		PedeWidth[i][0]  = f_Pede[i]->GetParameter(2);
		PedeWidth[i][1]  = f_Pede[i]->GetParError(2);
		Lat -> DrawLatexNDC( .75, .50, Form("Pedestal: %.1lf #pm %.1lf [ch]", PedeCenter[i][0], PedeWidth[i][0] ) );
		gPad->Update();
		Fr[i] = gPad->GetFrame();
		Fr[i] -> SetFillStyle(0);
		Fr[i] -> Draw("same");
	}
	for(int i=0; i<NofMPPC; i++){
		Ca[i]->cd(1);
		h_qdcev[i]->FitSlicesY( f_Slices[i], 0, -1, 0, "QNR0" );
		h_sl_qdcev[i] = (TH1D*)gROOT->FindObject( Form("h_qdcev[%d]_1", i) );
		h_sl_qdcev[i] -> SetLineColor(6);
		h_sl_qdcev[i] -> SetMarkerStyle(29);
		h_sl_qdcev[i] -> SetMarkerColor(6);
		h_sl_qdcev[i] -> Draw("samePES");
		Ca[i]->cd(3);
		h_sl_qdcev[i] -> Draw("PE");
	}
}

//------------------------------------------------------------------------------------------------//
void PedeStability::Export(){
	cout<<"__________ooooooooooOOOOO000OOOOOoooooooooo__________"<<endl;
	cout<<"==========     PedeStability::Export     =========="<<endl;

	Ca[0]->Print( figure_open, "pdf" );
	for(int i=0; i<NofMPPC; i++){Ca[i]->Print( figure );}
	Ca[NofMPPC-1]->Print( figure_close, "pdf" );

//	ofstream ofs( param, ios_base::trunc );
//	for(int i=0; i<NofMPPC; i++){ofs<<PedeCenter[i][0]<<"   "<<PedeWidth[i][0]<<endl;}
}
//------------------------------------------------------------------------------------------------//
////////////////////////////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv){
	int Num;
	
	if(argc==2){
		Num=atoi(argv[1]);
	}else{
		cout<<"  !!!WARNING!!!   "<<endl;
		cout<<"  This program needs one argument to specify runnumber."<<endl;
		cout<<"  Please input runnumber(Only number) >>>"<<endl;
		cin>>Num;
	}
	cout<<"-+-+-+-+-+-+-+-"<<endl;
	cout<<Form("Target of Analysys >> run%03d.root", Num)<<endl;
	cout<<"-+-+-+-+-+-+-+-"<<endl;

	TApplication *theApp;
	PedeStability *ana;

	theApp = new TApplication( "App", &argc, argv );
	ana = new PedeStability();
	ana -> SetRoot(Num);
	ana -> DefCanv();
	ana -> DefObj();
	ana -> FillHist();
	ana -> DrawHist();
	ana -> FitHist();
	ana -> Export();
	delete ana;

	//====	END of Job	====//
	gSystem->Exit(1);

	theApp->Run();

	return 0;
}
