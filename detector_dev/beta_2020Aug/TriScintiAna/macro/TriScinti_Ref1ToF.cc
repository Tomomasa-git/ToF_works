/* * * * * * * * * * * * * 
 * TriScinti_Ref1ToF.cc        *
 *                       *
 *	T. Fujiwara          *
 *	2020. 07. 24 (Fri)   *
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
#include "./TriScinti_Ref1ToF.hh"

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
TriScinti_Ref1ToF::TriScinti_Ref1ToF(){
	cout<<"TriScinti_Ref1ToF::Constructer called"<<endl;

	ItNum = 11;

	cutfactMin[0] = 0.50;
	cutfactMin[1] = 0.50;
	cutfactMin[2] = 0.75;
	cutfactMin[3] = 0.75;
	cutfact[0]    = 1.75;
	cutfact[1]    = 1.75;
	cutfact[2]    = 1.50;
	cutfact[3]    = 1.50;

	TwlkFitMax[0] = 210.;
	TwlkFitMax[2] = 90.,
	TwlkFitMax[3] = 60.,
	TwlkFitMin[0] = 5.;
	TwlkFitMin[2] = 25.,
	TwlkFitMin[3] = 20.,

	NDiv[0] = 50;
	NDiv[1] = 25;
	NDiv[2] = 25;

	QDCBin[0] = 250;
	QDCBin[1] = 150;
	QDCBin[2] = 150;
	
	QDCMax[0] = 250.;
	QDCMax[1] = 150.;
	QDCMax[2] = 150.;
	

	for(int i=0; i<NofMPPC; i++){
		TDCcut[i][0] = 1900.;
		TDCcut[i][1] = 2300.;
	}

	Set = new Setting();
	Set -> Setting_Gene(1);
	gStyle -> SetNumberContours(255);
	gStyle -> SetStatX(.875);
	gStyle -> SetStatW(.175);
	gStyle -> SetPadRightMargin(.125);
	gStyle -> SetPadBottomMargin(.125);

	for(int i=0; i<3; i++){
		Ln_Qpeak[i]   = new TLine();Set->Setting_Line( Ln_Qpeak[i]  , 6  , 1, 1 );
		Ln_Qcut[i][0] = new TLine();Set->Setting_Line( Ln_Qcut[i][0], 602, 1, 7 );	//Upper-side
		Ln_Qcut[i][1] = new TLine();Set->Setting_Line( Ln_Qcut[i][1], 602, 1, 7 );	//Downer-side
	}

	Lat = new TLatex();Set -> Setting_Latex( Lat, 62, 22, 6, 0.075 );

	for(int i=0; i<3; i++){
		TwlkFitStat[i] = new TPaveText();
		TwlkFitStat[i] -> SetX1NDC(.600);
		TwlkFitStat[i] -> SetX2NDC(.825);
		TwlkFitStat[i] -> SetY1NDC(.145);
		TwlkFitStat[i] -> SetY2NDC(.425);
		TwlkFitStat[i] -> SetTextFont(42);
		TwlkFitStat[i] -> SetTextSize(0.040);
		TwlkFitStat[i] -> SetFillColorAlpha(422, 0.65);
		TwlkFitStat[i] -> SetLineColor(602);
		TwlkFitStat[i] -> SetLineWidth(1);
		TwlkFitStat[i] -> SetLineStyle(1);
	}
}

////////////////////////////////////////////////////////////////////////////////////
TriScinti_Ref1ToF::~TriScinti_Ref1ToF(){
	cout<<"TriScinti_Ref1ToF::Destructer called"<<endl;
}	

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
double TriScinti_Ref1ToF::Twlk_Correction( double x, double* parameter, double Coeff ){
	double Twlk_func=0.;
	if(RunNum==999){
		Twlk_func = parameter[0]/sqrt( x-parameter[1] ) + parameter[2]*x;
	//	Twlk_func = parameter[0]/sqrt( x-parameter[1] ) + parameter[2]*(x+parameter[3])*(x+parameter[3]);
	}else{
		Twlk_func = parameter[0]/sqrt( x-parameter[1] );
	}
	Twlk_func = Coeff*Twlk_func;
	return Twlk_func;
}
////////////////////////////////////////////////////////////////////////////////////
double TriScinti_Ref1ToF::Twlk_Correction_ToF( double x, double* parameter, double Coeff, int Mode ){
	double Twlk_func=0.;

	switch(Mode){
		case 0:
		Twlk_func = parameter[0]/sqrt( x-parameter[1] );
		break;

		case 1:
		Twlk_func = parameter[0]*exp( parameter[1]*x );
		break;

		case 2:
		Twlk_func = parameter[0]*pow( x-parameter[1], parameter[2] );
		break;

		case 3:
		Twlk_func = parameter[0]*x;
		break;

		case 4:
		Twlk_func = parameter[0]/( x-parameter[1] );
		break;

		case 5:
		Twlk_func = parameter[0]*log( x-parameter[1] );
		break;

		default:
		Twlk_func = parameter[0]/sqrt( x-parameter[1] );
		break;
	}
		
	Twlk_func = Coeff*Twlk_func;
	return Twlk_func;
}

////////////////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1ToF::TwlkMan(int Type){
	TwlkMode = Type;
//	f_Twlk[0]->SetParameters(-50. , 0., 6.E+3);
//	f_Twlk[1]->SetParameters(1.E+4, 0., 1.E+3);
//	f_Twlk[2]->SetParameters(1.E+4, 0., 1.E+3);
	
	switch(TwlkMode){
		case 0:
		ToF_TwlkType = "[0]/sqrt(x-[1])+[2]";
		TwlkPar_In[1][0] = -1.E+4;
		TwlkPar_In[1][1] = 0.;
		TwlkPar_In[1][2] = 1.E+3;
		TwlkPar_In[2][0] = -1.E+4;
		TwlkPar_In[2][1] = 0.;
		TwlkPar_In[2][2] = 1.E+3;
		break;

		case 1:
		ToF_TwlkType = "[0]*exp([1]*x)+[2]";
		TwlkPar_In[1][0] = 2.E+3;
		TwlkPar_In[1][1] = -1.E-1;
		TwlkPar_In[1][2] = 1.E+3;
		TwlkPar_In[2][0] = 2.E+3;
		TwlkPar_In[2][1] = -1.E-1;
		TwlkPar_In[2][2] = 1.E+3;
		break;

		case 2:
		ToF_TwlkType = "[0]*pow(x-[1], [2])+[3]";
		break;

		case 3:
		ToF_TwlkType = "[0]*x+[1]";
		TwlkPar_In[1][0] = -50.;
		TwlkPar_In[1][1] = 4000.;
		TwlkPar_In[2][0] = -50.;
		TwlkPar_In[2][1] = 4000.;
		break;

		case 4:
		ToF_TwlkType = "[0]/(x-[1])";
		break;

		case 5:
		ToF_TwlkType = "[0]*log(x-[1])+[2]";
		break;

		default:
		ToF_TwlkType = "[0]/sqrt(x-[1])+[2]";
		TwlkPar_In[1][0] = -1.E+4;
		TwlkPar_In[1][1] = 0.;
		TwlkPar_In[1][2] = 1.E+3;
		TwlkPar_In[2][0] = -1.E+4;
		TwlkPar_In[2][1] = 0.;
		TwlkPar_In[2][2] = 1.E+3;
		break;
	}

	if(RunNum==999){
		Ref_TwlkType = "[0]/sqrt(x-[1])+[2]*x+[3]";
		TwlkPar_In[0][0] = -50.;
		TwlkPar_In[0][1] = 0.;
		TwlkPar_In[0][2] = -10.;
		TwlkPar_In[0][3] = 2.E+3;
	//	Ref_TwlkType = "[0]/sqrt(x-[1])+[2]*(x+[3])^2+[4]";
	//	TwlkPar_In[0][0] = -1.E+3;
	//	TwlkPar_In[0][1] = 0.;
	//	TwlkPar_In[0][2] = 1.E-3;
	//	TwlkPar_In[0][3] = -60.;
	//	TwlkPar_In[0][4] = 1.E+3;
	}else{
		Ref_TwlkType = "[0]/sqrt(x-[1])+[2]";
		TwlkPar_In[0][0] = -50.;
		TwlkPar_In[0][1] = 0.;
		TwlkPar_In[0][2] = 6.E+3;
	}

	f_Twlk[0] = new TF1( "f_Twlk[0]", Ref_TwlkType, 0., 400.);
	f_Twlk[1] = new TF1( "f_Twlk[1]", ToF_TwlkType, 0., 400.);
	f_Twlk[2] = new TF1( "f_Twlk[2]", ToF_TwlkType, 0., 400.);
	//Setting for Func
	for(int i=0; i<3; i++){
		Set -> Setting_Func( f_Twlk[i], 6, 1 );
		int NPar = f_Twlk[i]->GetNpar();
		for(int j=0; j<NPar; j++){f_Twlk[i]->SetParameter(j, TwlkPar_In[i][j]);}
	}
	if(RunNum==999){
		f_Twlk[0]->SetParLimits(2,   0., -1.);
		f_Twlk[0]->SetParLimits(3, -70., -50.);
	}else;
}


////////////////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1ToF::SetRoot( int RNum ){
	RunNum = RNum;

	RunName_str  = Form("run%03d"     , RunNum);
	RunName_tst  = Form("run%03d"     , RunNum);
	RunName_Full = Form("run%03d.root", RunNum);	cout<<RunName_Full<<endl;
	RunName_Full = PATH_dat+RunName_Full;          	cout<<RunName_Full<<endl;

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

////////////////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1ToF::ImportPedestal( int PNum ){
	int i=0;
	double Val1, Val2;	

	PedNum = PNum;
	Param_Ped = Form("Pedestal_run%03d_00.dat", PedNum);
	
	ifstream ifs( PATH_param_p+Param_Ped.c_str(), ios_base::in );
	while(ifs>>Val1>>Val2){
		Pede_Cnt_Val[i]=Val1;
		Pede_Wid_Val[i]=Val2;
		cout<<Pede_Cnt_Val[i]<<"   "<<Pede_Wid_Val[i]<<endl;
		i++;
	}

	cout<<"======================================"<<endl;
	cout<<"    FILE IMPORT: "<<Param_Ped<<endl;
	cout<<"======================================"<<endl;
}

////////////////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1ToF::DefineCanv(){
	cout<<"=====   TriScinti_Ref1ToF::DefineCanv ====="<<endl;
	for(int i=0; i<NofCanv; i++){Ca[i] = new TCanvas( Form("Ca[%d]", i), Form("Ca[%d]", i), 1000, CanvHsize[i] );}
	Ca[0] -> Divide(2,2);
	Ca[1] -> Divide(2,2);
	Ca[2] -> Divide(2,2);
	Ca[3] -> Divide(2,2);
	Ca[4] -> Divide(2,2);

	for(int i=0; i<6; i++){Ca_dummy[i] = new TCanvas(Form("Ca_dummy[%d]", i), Form("Ca_dummy[%d]", i), 1200, 600);}
}

////////////////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1ToF::DefineHistFunc(){
	cout<<"=====   TriScinti_Ref1ToF::DefineHistFunc ====="<<endl;
	
	for(int i=0; i<3; i++){
		h_2D_RawTQ[i] = new TH2D( Form("h_2D_RawTQ[%d]",i), Form("h_2D_RawTQ[%d]",i), QDCBin[i], 0., QDCMax[i], 80, -4200., 7000. );
		h_2D_FitTQ[i] = new TH2D( Form("h_2D_FitTQ[%d]",i), Form("h_2D_FitTQ[%d]",i), QDCBin[i], 0., QDCMax[i], 80, -4200., 7000. );
		h_2D_DecTQ[i] = new TH2D( Form("h_2D_DecTQ[%d]",i), Form("h_2D_DecTQ[%d]",i), QDCBin[i], 0., QDCMax[i], 80, -4200., 7000. );
		h_2D_CutTQ[i] = new TH2D( Form("h_2D_CutTQ[%d]",i), Form("h_2D_CutTQ[%d]",i), QDCBin[i], 0., QDCMax[i], 80, -4200., 7000. );

		h_1D_QDCFull[i] = new TH1D( Form("h_1D_QDCFull[%d]", i), Form("h_1D_QDCFull[%d]", i), QDCBin[i], 0, QDCMax[i] );
		h_1D_QDCwCut[i] = new TH1D( Form("h_1D_QDCwCut[%d]", i), Form("h_1D_QDCwCut[%d]", i), QDCMax[i], 0, QDCMax[i] );

		//Setting for 2D hist
		Set -> Setting_Hist2D( h_2D_RawTQ[i], "RawToF vs."+Label[ChID[i]][0]                      , Label[ChID[i]][0]+" [pC]", "RawToF [psec]", "Entry/bin", 1. );
		Set -> Setting_Hist2D( h_2D_FitTQ[i], "RawToF vs."+Label[ChID[i]][0]                      , Label[ChID[i]][0]+" [pC]", "RawToF [psec]", "Entry/bin", 1. );
		Set -> Setting_Hist2D( h_2D_DecTQ[i], "DecToF vs."+Label[ChID[i]][0]+" w/ TWC"            , Label[ChID[i]][0]+" [pC]", "DecToF [psec]", "Entry/bin", 1. );
		Set -> Setting_Hist2D( h_2D_CutTQ[i], "DecToF vs."+Label[ChID[i]][0]+" w/ TWC and QDC cut", Label[ChID[i]][0]+" [pC]", "DecToF [psec]", "Entry/bin", 1. );
		//Setting for 1D hist
		Set -> Setting_Hist1D( h_1D_QDCFull[i], Label[ChID[i]][2], Label[ChID[i]][0]+" [pC]", "Counts/1.0 pC", LColor, 1, 62, FColor , 3001 );
		Set -> Setting_Hist1D( h_1D_QDCwCut[i], Label[ChID[i]][2], Label[ChID[i]][0]+" [pC]", "Counts/1.0 pC", LColor, 1, 62, kYellow, 3001 );
	}

	for(int i=0; i<NofCombination; i++){
		h_2D_Factor[i] = new TH2D( Form("h_2D_Factor[%d]", i), Form("h_2D_Factor[%d]", i), 71, -0.0125, 1.7625, 71, -0.0125, 1.7625 );
		Gr_SigBest[i] = new TGraph();
	}
	Set -> Setting_Hist2D( h_2D_Factor[0], "Factor tuning: ToF-1 vs. ToF-2", "ToF-1 Factor", "ToF-2 Factor", "#sigma_{ToF} [psec]", 1.);
	Set -> Setting_Hist2D( h_2D_Factor[1], "Factor tuning: Ref-1 vs. ToF-1", "Ref-1 Factor", "ToF-1 Factor", "#sigma_{ToF} [psec]", 1.);
	Set -> Setting_Hist2D( h_2D_Factor[2], "Factor tuning: Ref-1 vs. ToF-2", "Ref-1 Factor", "ToF-2 Factor", "#sigma_{ToF} [psec]", 1.);

	Set -> Setting_Graph( Gr_SigBest[0], "Factor tuning: ToF-1 vs. ToF-2", "ToF-1 Factor", "ToF-2 Factor", 1, 1, 62, 6, 29 );
	Set -> Setting_Graph( Gr_SigBest[1], "Factor tuning: Ref-1 vs. ToF-1", "Ref-1 Factor", "ToF-1 Factor", 1, 1, 62, 6, 29 );
	Set -> Setting_Graph( Gr_SigBest[2], "Factor tuning: Ref-1 vs. ToF-2", "Ref-1 Factor", "ToF-2 Factor", 1, 1, 62, 6, 29 );

	h_1D_DecToF = new TH1D( "h_1D_DecToF", "h_1D_DecToF", 120, -2100., 2100. );
	Set -> Setting_Hist1D( h_1D_DecToF, "dec ToF (After correction)", "ToF [psec]", "Counts/35 psec", LColor, 1, 62, 870, 3001 );
	
	f_gaus = new TF1( "f_gaus", "gaus(0)", -10000., 10000. );
	Set -> Setting_Func( f_gaus, 2, 1 );

	f_ToF = new TF1( "f_ToF", "gaus(0)", -10500., 10500. );
	Set -> Setting_Func( f_ToF, 2, 1 );
	
	f_dTof = new TF1( "f_dTof", "gaus(0)", -7000., 7000. );
	Set -> Setting_Func( f_dTof, 2, 1 );

	for(int i=0; i<2; i++){
		for(int j=0; j<3; j++){
			Gr_fit[i][j] = new TGraphErrors();
			Set -> Setting_GError( Gr_fit[i][j], "RawToF vs."+Label[j][0], Label[j][0]+" [pC]", "RawToF [psec]", 1, 1, 62, 1, 29 );
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1ToF::FillHist(){
	cout<<"=====   TriScinti_Ref1ToF::FillHist   ====="<<endl;
	for(int i=0; i<Entry; i++){
		tree->GetEntry(i);
		if(
			   TDC[0]>TDCcut[0][0]
			&& TDC[0]<TDCcut[0][1]
			&& TDC[1]>TDCcut[1][0]
			&& TDC[1]<TDCcut[1][1]
			&& TDC[2]>TDCcut[2][0]
			&& TDC[2]<TDCcut[2][1]
			&& TDC[3]>TDCcut[3][0]
			&& TDC[3]<TDCcut[3][1]
		){
			rawtof = 35.*( TDC[0] - 0.5*(TDC[2]+TDC[3]) );
			for(int j=0; j<3; j++){
				h_2D_RawTQ[j]   -> Fill(0.1*QDC[ChID[j]]-0.1*Pede_Cnt_Val[ChID[j]], rawtof);
				h_1D_QDCFull[j] -> Fill(0.1*QDC[ChID[j]]-0.1*Pede_Cnt_Val[ChID[j]]);
			}
		}else;
	}
}

////////////////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1ToF::GetQDCCut(){
	cout<<"=====   TriScinti_Ref1ToF::GetQDCCut   ====="<<endl;
	for(int i=0; i<3; i++){
		Qpeak[i] = h_1D_QDCFull[i]->GetBinCenter( h_1D_QDCFull[i]->GetMaximumBin() );
		QCut[i][0] = cutfactMin[ChID[i]]*Qpeak[i]+0.1*Pede_Cnt_Val[ChID[i]];
		QCut[i][1] = cutfact[ChID[i]]*Qpeak[i]+0.1*Pede_Cnt_Val[ChID[i]];
		cout<<Pede_Cnt_Val[i]<<endl;
		cout<<Qpeak[i]<<endl;
		cout<<QCut[i][0]<<endl;
		cout<<QCut[i][1]<<endl;

		Ln_Qpeak[i] -> SetX1(Qpeak[i]);
		Ln_Qpeak[i] -> SetX2(Qpeak[i]);
		Ln_Qcut[i][0] -> SetX1(cutfactMin[ChID[i]]*Qpeak[i]);
		Ln_Qcut[i][0] -> SetX2(cutfactMin[ChID[i]]*Qpeak[i]);
		Ln_Qcut[i][1] -> SetX1(cutfact[ChID[i]]*Qpeak[i]);
		Ln_Qcut[i][1] -> SetX2(cutfact[ChID[i]]*Qpeak[i]);
		TwlkFitMax[ChID[i]] = 1.75*Qpeak[i];
		TwlkFitMin[ChID[i]] = 0.75*Qpeak[i];
	}
//	TwlkFitMin[0] = 0.25*Qpeak[0];
//	TwlkFitMin[0] = 0.25*Qpeak[0];				//For run116_06
	TwlkFitMin[0] = 5.;				//For run116_11
	TwlkFitMax[0] = 110.;				//For run116_06
//	TwlkFitMin[0] = 12.;			//For run116_01, 02
//	if(RunNum==134){
//		TwlkFitMin[0] = 0.25*Qpeak[0];
//	}else{
//		TwlkFitMin[0] = 0.50*Qpeak[0];
//	}

	for(int i=0; i<Entry; i++){
		tree->GetEntry(i);
		for(int j=0; j<3; j++){
			if(
				   0.1*QDC[ChID[j]]>QCut[j][0]
				&& 0.1*QDC[ChID[j]]<QCut[j][1]
			){
				h_1D_QDCwCut[j] -> Fill(0.1*QDC[ChID[j]]-0.1*Pede_Cnt_Val[ChID[j]]);
			}else;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1ToF::TwlkFit(){
	cout<<"=====   TriScinti_Ref1ToF::TwlkFit   ====="<<endl;


	for(int i=0; i<Entry; i++){
		tree->GetEntry(i);
		rawtof = 35.*( TDC[ChID[0]] - 0.5*(TDC[ChID[1]]+TDC[ChID[2]]) );
		if(
			   0.1*QDC[ChID[1]]>QCut[1][0]
			&& 0.1*QDC[ChID[1]]<QCut[1][1]
			&& 0.1*QDC[ChID[2]]>QCut[2][0]
			&& 0.1*QDC[ChID[2]]<QCut[2][1]
			&& TDC[0]>TDCcut[0][0]
			&& TDC[0]<TDCcut[0][1]
			&& TDC[1]>TDCcut[1][0]
			&& TDC[1]<TDCcut[1][1]
			&& TDC[2]>TDCcut[2][0]
			&& TDC[2]<TDCcut[2][1]
			&& TDC[3]>TDCcut[3][0]
			&& TDC[3]<TDCcut[3][1]
		){
			h_2D_FitTQ[0] -> Fill(0.1*QDC[0]-0.1*Pede_Cnt_Val[0], rawtof);
		}else;
		if(
			   0.1*QDC[ChID[0]]>QCut[0][0]
			&& 0.1*QDC[ChID[0]]<QCut[0][1]
			&& TDC[0]>TDCcut[0][0]
			&& TDC[0]<TDCcut[0][1]
			&& TDC[1]>TDCcut[1][0]
			&& TDC[1]<TDCcut[1][1]
			&& TDC[2]>TDCcut[2][0]
			&& TDC[2]<TDCcut[2][1]
			&& TDC[3]>TDCcut[3][0]
			&& TDC[3]<TDCcut[3][1]
		){
			h_2D_FitTQ[1] -> Fill(0.1*QDC[2]-0.1*Pede_Cnt_Val[2], rawtof);
			h_2D_FitTQ[2] -> Fill(0.1*QDC[3]-0.1*Pede_Cnt_Val[3], rawtof);
		}else;
	}

	Ca[2]->cd(1);
	for(int i=0; i<3; i++){
		Ca[2]->cd(i+1);
		gPad->SetLogz(1);
		h_2D_FitTQ[i]->Draw("colz");
		TF1 *f_FitSlices = new TF1( "f_FitSlices", "gaus(0)", -4200., 7000. );
		h_2D_FitTQ[i]->FitSlicesY( f_FitSlices, 0, -1, 0, "QNR" );
		h_2D_FitSlices[i] = (TH1D*)gROOT->FindObject(Form("h_2D_FitTQ[%d]_1", i));
		h_2D_FitSlices[i] -> SetLineColor(1);
		Ca[2]->cd(i+1);
		h_2D_FitSlices[i] -> Draw("samePES");
		for(int j=0; j<20; j++){h_2D_FitSlices[i] -> Fit( f_Twlk[i], "", "", TwlkFitMin[ChID[i]], TwlkFitMax[ChID[i]] );}
		h_2D_FitSlices[i] -> Fit( f_Twlk[i], "", "", TwlkFitMin[ChID[i]], TwlkFitMax[ChID[i]] );
		int NPar = f_Twlk[i]->GetNpar();
		for(int j=0; j<NPar; j++){
			int par = f_Twlk[i]->GetParameter(j);
			f_Twlk[i]->SetParameter(j, par);
			cout<<j<<": "<<par<<endl;
		}
	//	cout<<"====="<<endl;
	//	g_Twlk[i] = new TGraph(f_Twlk[i]);
	//	Set -> Setting_Graph(g_Twlk[i], "RawToF vs."+Label[ChID[i]][0], Label[ChID[i]][0]+" [pC]", "RawToF [psec]", 6, 1, 62, 1, 29);
	//	g_Twlk[i] -> GetXaxis()->SetLimits(0., 250.);
	//	g_Twlk[i] -> GetYaxis()->SetRangeUser(-4200., 7000.);
	//	Ca[2]->cd(i+1);
	//	g_Twlk[i] -> Draw("sameL");
	//	Ca[2]->cd(i+1);
	//	gPad->Update();
		Ca[2]->cd(i+1);
		gPad->Update();
		TwlkFitStat[i] -> AddText( Form("#chi^{2}/NDF:   %.1lf/%d", f_Twlk[i]->GetChisquare(), f_Twlk[i]->GetNDF()) );
		for(int j=0; j<NPar; j++){
			if(j==1 && TwlkMode==1){
				TwlkFitStat[i]->AddText( Form("p_{%d} = %.4lf #pm %.4lf", j, f_Twlk[i]->GetParameter(j), f_Twlk[i]->GetParError(j)));
			}else{
				TwlkFitStat[i]->AddText( Form("p_{%d} = %.2lf #pm %.2lf", j, f_Twlk[i]->GetParameter(j), f_Twlk[i]->GetParError(j)));
			}
		}
		TwlkFitStat[i]->Draw("same");
		gPad->Update();
		for(int j=0; j<NPar; j++){TwlkPar_Val[i][j]=f_Twlk[i]->GetParameter(j);}
	}
}

////////////////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1ToF::TwlkSerchBestFactorToF1vsToF2(){
	cout<<"=====   TriScinti_Ref1ToF::TwlkSerchBestFactorToF1vsToF2   ====="<<endl;
	TwlkFactor[0] = 1.000;
	for(int i=0; i<71; i++){
		for(int j=0; j<71; j++){
			TwlkFactor[1] = 0.000+0.025*i;
			TwlkFactor[2] = 0.000+0.025*j;
		
			h_1D_DecToF_tmp = new TH1D( "h_1D_DecToF_tmp", "h_1D_DecToF_tmp", 600, -10500., 10500. );
			
			for(int k=0; k<Entry; k++){
				tree->GetEntry(k);
				dectof = 35.*( TDC[ChID[0]] - 0.5*(TDC[ChID[1]]+TDC[ChID[2]]) )
				         -Twlk_Correction( 0.1*QDC[0]-0.1*Pede_Cnt_Val[0], TwlkPar_Val[0], TwlkFactor[0] )
						 -Twlk_Correction_ToF( 0.1*QDC[2]-0.1*Pede_Cnt_Val[ChID[1]], TwlkPar_Val[1], TwlkFactor[1], TwlkMode )
				         -Twlk_Correction_ToF( 0.1*QDC[3]-0.1*Pede_Cnt_Val[ChID[2]], TwlkPar_Val[2], TwlkFactor[2], TwlkMode );
				if(
					   TDC[0]>TDCcut[0][0]
					&& TDC[0]<TDCcut[0][1]
					&& TDC[1]>TDCcut[1][0]
					&& TDC[1]<TDCcut[1][1]
					&& TDC[2]>TDCcut[2][0]
					&& TDC[2]<TDCcut[2][1]
					&& TDC[3]>TDCcut[3][0]
					&& TDC[3]<TDCcut[3][1]
				//	&& 0.1*QDC[ChID[0]]>QCut[0][0]
				//	&& 0.1*QDC[ChID[1]]>QCut[1][0]
				//	&& 0.1*QDC[ChID[2]]>QCut[2][0]
				){
					h_1D_DecToF_tmp->Fill(dectof);
				}else;
			}

			canv_dummy_consgaus = new TCanvas( "canv_dummy_consgaus", "canv_dummy_consgaus", 480, 240);
			par_tmp1 = h_1D_DecToF_tmp -> GetBinCenter( h_1D_DecToF_tmp->GetMaximumBin() );
			h_1D_DecToF_tmp  -> Fit( f_ToF, "Q", "", par_tmp1-1000., par_tmp1+1000.);
			//Resolusion
			Twlk_SigVal_tmp = f_ToF -> GetParameter(2);
			Twlk_SigErr_tmp = f_ToF -> GetParError(2);
			//Chisquare/NDF
			Twlk_Chis_tmp = f_ToF -> GetChisquare();
			Twlk_Ndf_tmp  = f_ToF -> GetNDF();
			h_1D_DecToF_tmp -> Delete();
			canv_dummy_consgaus -> Destructor();

			if( Twlk_SigVal_tmp<Twlk_SigVal_Best[0] ){
				Twlk_SigVal_Best[0] = Twlk_SigVal_tmp;
				Twlk_SigErr_Best[0] = Twlk_SigErr_tmp;
				Twlk_FactorBestSigma[1] = TwlkFactor[1];
				Twlk_FactorBestSigma[2] = TwlkFactor[2];
				cout<<"Preriminary factor ToF-1 and ToF-2 >> "<<Twlk_FactorBestSigma[1]<<",  "<<Twlk_FactorBestSigma[2]<<":  "<<"Sigma_TOF = "<<Twlk_SigVal_Best[0]<<" psec"<<endl;
			}else;
			h_2D_Factor[0] -> Fill( TwlkFactor[1], TwlkFactor[2], Twlk_SigVal_tmp );
		}
	}
	Gr_SigBest[0]  -> SetPoint( 0, Twlk_FactorBestSigma[1] , Twlk_FactorBestSigma[2]  );
	TextBestSig[0] = Form("Best Factor = (%.3lf, %.3lf)", Twlk_FactorBestSigma[1], Twlk_FactorBestSigma[2]);
	h_2D_Factor[0] -> SetMinimum( Twlk_SigVal_Best[0] );
	TwlkFactor[1] = Twlk_FactorBestSigma[1];
	TwlkFactor[2] = Twlk_FactorBestSigma[2];
	Twlk_SigVal_Best[1] = Twlk_SigVal_Best[0];

	Twlk_FactorBestSigma[0] = TwlkFactor[0];
}
 
////////////////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1ToF::TwlkSerchBestFactorRef1vsToF1(){
	cout<<"=====   TriScinti_Ref1ToF::TwlkSerchBestFactorRef1vsToF1   ====="<<endl;
	for(int i=0; i<71; i++){
		for(int j=0; j<71; j++){
			TwlkFactor[0] = 0.000+0.025*i;
			TwlkFactor[1] = 0.000+0.025*j;
		
			h_1D_DecToF_tmp = new TH1D( "h_1D_DecToF_tmp", "h_1D_DecToF_tmp", 600, -10500., 10500. );
			
			for(int k=0; k<Entry; k++){
				tree->GetEntry(k);
				dectof = 35.*( TDC[ChID[0]] - 0.5*(TDC[ChID[1]]+TDC[ChID[2]]) )
				         -Twlk_Correction( 0.1*QDC[0]-0.1*Pede_Cnt_Val[0], TwlkPar_Val[0], TwlkFactor[0] )
						 -Twlk_Correction_ToF( 0.1*QDC[2]-0.1*Pede_Cnt_Val[ChID[1]], TwlkPar_Val[1], TwlkFactor[1], TwlkMode )
				         -Twlk_Correction_ToF( 0.1*QDC[3]-0.1*Pede_Cnt_Val[ChID[2]], TwlkPar_Val[2], Twlk_FactorBestSigma[2], TwlkMode );
				if(
					   TDC[0]>TDCcut[0][0]
					&& TDC[0]<TDCcut[0][1]
					&& TDC[1]>TDCcut[1][0]
					&& TDC[1]<TDCcut[1][1]
					&& TDC[2]>TDCcut[2][0]
					&& TDC[2]<TDCcut[2][1]
					&& TDC[3]>TDCcut[3][0]
					&& TDC[3]<TDCcut[3][1]
				//	&& 0.1*QDC[ChID[0]]>QCut[0][0]
				//	&& 0.1*QDC[ChID[1]]>QCut[1][0]
				//	&& 0.1*QDC[ChID[2]]>QCut[2][0]
				){
					h_1D_DecToF_tmp->Fill(dectof);
				}else;
			}

			canv_dummy_consgaus = new TCanvas( "canv_dummy_consgaus", "canv_dummy_consgaus", 480, 240);
			par_tmp1 = h_1D_DecToF_tmp -> GetBinCenter( h_1D_DecToF_tmp->GetMaximumBin() );
			h_1D_DecToF_tmp  -> Fit( f_ToF, "Q", "", par_tmp1-800., par_tmp1+800.);
			//Resolusion
			Twlk_SigVal_tmp = f_ToF -> GetParameter(2);
			Twlk_SigErr_tmp = f_ToF -> GetParError(2);
			//Chisquare/NDF
			Twlk_Chis_tmp = f_ToF -> GetChisquare();
			Twlk_Ndf_tmp  = f_ToF -> GetNDF();
			h_1D_DecToF_tmp -> Delete();
			canv_dummy_consgaus -> Destructor();

			if( Twlk_SigVal_tmp<Twlk_SigVal_Best[1] ){
				Twlk_SigVal_Best[1] = Twlk_SigVal_tmp;
				Twlk_SigErr_Best[1] = Twlk_SigErr_tmp;
				Twlk_FactorBestSigma[0] = TwlkFactor[0];
				Twlk_FactorBestSigma[1] = TwlkFactor[1];
				cout<<"Preriminary factor Ref-1 and ToF-1 >> "<<Twlk_FactorBestSigma[0]<<",  "<<Twlk_FactorBestSigma[1]<<":  "<<"Sigma_TOF = "<<Twlk_SigVal_Best[1]<<" psec"<<endl;
			}else;
			h_2D_Factor[1] -> Fill( TwlkFactor[0], TwlkFactor[1], Twlk_SigVal_tmp );
		}
	}
	Gr_SigBest[1]  -> SetPoint( 0, Twlk_FactorBestSigma[0] , Twlk_FactorBestSigma[1]  );
	TextBestSig[1] = Form("Best Factor = (%.3lf, %.3lf)", Twlk_FactorBestSigma[0], Twlk_FactorBestSigma[1]);
	h_2D_Factor[1] -> SetMinimum( Twlk_SigVal_Best[1] );
	TwlkFactor[0] = Twlk_FactorBestSigma[0];
	TwlkFactor[1] = Twlk_FactorBestSigma[1];
	Twlk_SigVal_Best[2] = Twlk_SigVal_Best[1];
}

////////////////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1ToF::TwlkSerchBestFactorRef1vsToF2(){
	cout<<"=====   TriScinti_Ref1ToF::TwlkSerchBestFactorRef1vsToF2   ====="<<endl;
	for(int i=0; i<71; i++){
		for(int j=0; j<71; j++){
			TwlkFactor[0] = 0.000+0.025*i;
			TwlkFactor[2] = 0.000+0.025*j;
		
			h_1D_DecToF_tmp = new TH1D( "h_1D_DecToF_tmp", "h_1D_DecToF_tmp", 600, -10500., 10500. );
			
			for(int k=0; k<Entry; k++){
				tree->GetEntry(k);
				dectof = 35.*( TDC[ChID[0]] - 0.5*(TDC[ChID[1]]+TDC[ChID[2]]) )
				         -Twlk_Correction( 0.1*QDC[0]-0.1*Pede_Cnt_Val[0], TwlkPar_Val[0], TwlkFactor[0] )
						 -Twlk_Correction_ToF( 0.1*QDC[2]-0.1*Pede_Cnt_Val[ChID[1]], TwlkPar_Val[1], Twlk_FactorBestSigma[1], TwlkMode )
				         -Twlk_Correction_ToF( 0.1*QDC[3]-0.1*Pede_Cnt_Val[ChID[2]], TwlkPar_Val[2], TwlkFactor[2], TwlkMode );
				if(
					   TDC[0]>TDCcut[0][0]
					&& TDC[0]<TDCcut[0][1]
					&& TDC[1]>TDCcut[1][0]
					&& TDC[1]<TDCcut[1][1]
					&& TDC[2]>TDCcut[2][0]
					&& TDC[2]<TDCcut[2][1]
					&& TDC[3]>TDCcut[3][0]
					&& TDC[3]<TDCcut[3][1]
				//	&& 0.1*QDC[ChID[0]]>QCut[0][0]
				//	&& 0.1*QDC[ChID[1]]>QCut[1][0]
				//	&& 0.1*QDC[ChID[2]]>QCut[2][0]
				){
					h_1D_DecToF_tmp->Fill(dectof);
				}else;
			}

			canv_dummy_consgaus = new TCanvas( "canv_dummy_consgaus", "canv_dummy_consgaus", 480, 240);
			par_tmp1 = h_1D_DecToF_tmp -> GetBinCenter( h_1D_DecToF_tmp->GetMaximumBin() );
			h_1D_DecToF_tmp  -> Fit( f_ToF, "Q", "", par_tmp1-1000., par_tmp1+1000.);
			//Resolusion
			Twlk_SigVal_tmp = f_ToF -> GetParameter(2);
			Twlk_SigErr_tmp = f_ToF -> GetParError(2);
			//Chisquare/NDF
			Twlk_Chis_tmp = f_ToF -> GetChisquare();
			Twlk_Ndf_tmp  = f_ToF -> GetNDF();
			h_1D_DecToF_tmp -> Delete();
			canv_dummy_consgaus -> Destructor();

			if( Twlk_SigVal_tmp<Twlk_SigVal_Best[2] ){
				Twlk_SigVal_Best[2] = Twlk_SigVal_tmp;
				Twlk_SigErr_Best[2] = Twlk_SigErr_tmp;
				Twlk_FactorBestSigma[0] = TwlkFactor[0];
				Twlk_FactorBestSigma[2] = TwlkFactor[2];
				cout<<"Preriminary factor Ref-1 and ToF-1 >> "<<Twlk_FactorBestSigma[0]<<",  "<<Twlk_FactorBestSigma[2]<<":  "<<"Sigma_TOF = "<<Twlk_SigVal_Best[2]<<" psec"<<endl;
			}else;
			h_2D_Factor[2] -> Fill( TwlkFactor[0], TwlkFactor[2], Twlk_SigVal_tmp );
		}
	}
	Gr_SigBest[2]  -> SetPoint( 0, Twlk_FactorBestSigma[0] , Twlk_FactorBestSigma[2]  );
	TextBestSig[2] = Form("Best Factor = (%.3lf, %.3lf)", Twlk_FactorBestSigma[0], Twlk_FactorBestSigma[2]);
	h_2D_Factor[2] -> SetMinimum( Twlk_SigVal_Best[2] );
	TwlkFactor[0] = Twlk_FactorBestSigma[0];
	TwlkFactor[2] = Twlk_FactorBestSigma[2];
}

////////////////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1ToF::TwlkDecToF(){
	cout<<"=====   TriScinti_Ref1ToF::TwlkDecToF   ====="<<endl;

	h_1D_DecToF_tmp = new TH1D( "h_1D_DecToF_tmp", "h_1D_DecToF_tmp", 600, -10500., 10500. );
	canv_dummy_consgaus = new TCanvas( "canv_dummy_consgaus", "canv_dummy_consgaus", 480, 240);

	for(int i=0; i<Entry; i++){
		tree->GetEntry(i);
		dectof = 35.*( TDC[ChID[0]] - 0.5*(TDC[ChID[1]]+TDC[ChID[2]]) )
		         -Twlk_Correction( 0.1*QDC[0]-0.1*Pede_Cnt_Val[0], TwlkPar_Val[0], Twlk_FactorBestSigma[0] )
		         -Twlk_Correction_ToF( 0.1*QDC[2]-0.1*Pede_Cnt_Val[ChID[1]], TwlkPar_Val[1], Twlk_FactorBestSigma[1], TwlkMode )
		         -Twlk_Correction_ToF( 0.1*QDC[3]-0.1*Pede_Cnt_Val[ChID[2]], TwlkPar_Val[2], Twlk_FactorBestSigma[2], TwlkMode );
		if(
			   TDC[0]>TDCcut[0][0]
			&& TDC[0]<TDCcut[0][1]
			&& TDC[1]>TDCcut[1][0]
			&& TDC[1]<TDCcut[1][1]
			&& TDC[2]>TDCcut[2][0]
			&& TDC[2]<TDCcut[2][1]
			&& TDC[3]>TDCcut[3][0]
			&& TDC[3]<TDCcut[3][1]
		){
			h_1D_DecToF_tmp -> Fill(dectof);
		}else;
	}
	ToF_Adj = h_1D_DecToF_tmp->GetBinCenter( h_1D_DecToF_tmp->GetMaximumBin() );
	h_1D_DecToF_tmp -> Delete();
	canv_dummy_consgaus -> Destructor();

	for(int i=0; i<Entry; i++){
		tree->GetEntry(i);
		dectof = 35.*( TDC[ChID[0]] - 0.5*(TDC[ChID[1]]+TDC[ChID[2]]) )
		         -Twlk_Correction( 0.1*QDC[0]-0.1*Pede_Cnt_Val[0], TwlkPar_Val[0], Twlk_FactorBestSigma[0] )
		         -Twlk_Correction_ToF( 0.1*QDC[2]-0.1*Pede_Cnt_Val[ChID[1]], TwlkPar_Val[1], Twlk_FactorBestSigma[1], TwlkMode )
		         -Twlk_Correction_ToF( 0.1*QDC[3]-0.1*Pede_Cnt_Val[ChID[2]], TwlkPar_Val[2], Twlk_FactorBestSigma[2], TwlkMode )
		         -ToF_Adj;
		if(
			   TDC[0]>TDCcut[0][0]
			&& TDC[0]<TDCcut[0][1]
			&& TDC[1]>TDCcut[1][0]
			&& TDC[1]<TDCcut[1][1]
			&& TDC[2]>TDCcut[2][0]
			&& TDC[2]<TDCcut[2][1]
			&& TDC[3]>TDCcut[3][0]
		){
			for(int j=0; j<3; j++){h_2D_DecTQ[j] -> Fill(0.1*QDC[ChID[j]]-0.1*Pede_Cnt_Val[j], dectof);}
		}else;
		if(
			   TDC[0]>TDCcut[0][0]
			&& TDC[0]<TDCcut[0][1]
			&& TDC[1]>TDCcut[1][0]
			&& TDC[1]<TDCcut[1][1]
			&& TDC[2]>TDCcut[2][0]
			&& TDC[2]<TDCcut[2][1]
			&& TDC[3]>TDCcut[3][0]
			&& TDC[3]<TDCcut[3][1]
			&& 0.1*QDC[ChID[0]]>QCut[0][0]
			&& 0.1*QDC[ChID[1]]>QCut[1][0]
			&& 0.1*QDC[ChID[2]]>QCut[2][0]
		){
			for(int j=0; j<3; j++){h_2D_CutTQ[j] -> Fill(0.1*QDC[ChID[j]]-0.1*Pede_Cnt_Val[j], dectof);}
			h_1D_DecToF->Fill(dectof);
		}else;
	}

	Ca[3]->cd(1);
	for(int i=0; i<3; i++){
		Ca[3]->cd(i+1);
		gPad->SetLogz(1);
		h_2D_DecTQ[i]->Draw("colz");
	}

	Ca[4]->cd(1);
	for(int i=0; i<3; i++){
		Ca[4]->cd(i+1);
		gPad->SetLogz(1);
		h_2D_CutTQ[i]->Draw("colz");
	}

	Ca[5]->cd();
	gPad->SetRightMargin(.075);
	h_1D_DecToF->Draw("");
	par_tmp2 = h_1D_DecToF->GetBinCenter( h_1D_DecToF->GetMaximumBin() );
	h_1D_DecToF->Fit( f_dTof, "", "", par_tmp2-1000., par_tmp2+1000. );
	Twlk_ResoVal = f_dTof->GetParameter(2);
	Twlk_ResoErr = f_dTof->GetParError(2);
//	ResoToF_Val = sqrt( Twlk_ResoVal*Twlk_ResoVal - Reso_Ref_Val*Reso_Ref_Val );
//	ResoToF_Err = sqrt( pow(Twlk_ResoVal*Twlk_ResoErr,2.) + pow(Reso_Ref_Val*Reso_Ref_Err,2.) )/ResoToF_Val;

	Lat -> SetTextSize(.05);
	Lat -> DrawLatexNDC(.75, .45, Form("#sigma_{Ref1-ToF} = %.1lf #pm %.1lf psec", Twlk_ResoVal, Twlk_ResoErr) );
	
}

////////////////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1ToF::DrawHist(){
	cout<<"=====   TriScinti_Ref1ToF::DrawHist   ====="<<endl;

	Ca[0]->cd(1);
	for(int i=0; i<3; i++){
		Ca[0]->cd(i+1);
		gPad->SetLogz(1);
		h_2D_RawTQ[i]->Draw("colz");
	}

	Ca[1]->cd(1);
	for(int i=0; i<3; i++){
		Ca[1]->cd(i+1);
		gPad->SetLogy(1);
		h_1D_QDCFull[i]->Draw("");
		h_1D_QDCwCut[i] -> Draw("same");
		gPad->Update();
		HistHmin[i] = pow( 10., gPad->GetUymin()  );
		HistHmax[i] = pow( 10., gPad->GetUymax() );
		Ln_Qpeak[i] -> SetY1(HistHmin[i]);
		Ln_Qpeak[i] -> SetY2(HistHmax[i]);
		Ln_Qpeak[i] -> Draw("same");
		for(int j=0; j<2; j++){
			Ln_Qcut[i][j] -> SetY1(HistHmin[i]);
			Ln_Qcut[i][j] -> SetY2(HistHmax[i]);
			Ln_Qcut[i][j] -> Draw("same");
		}
		Lat -> DrawLatexNDC( .70, .50, Form("QDC peak: %.1lf [pC]", Qpeak[i]) );
	}

	Ca[6]->cd();
	for(int i=0; i<NofCombination; i++){
		Ca[i+6]->cd();
		gPad -> SetRightMargin(.125);
		gPad -> SetTopMargin(.090);
		gPad -> SetBottomMargin(.125);
		gPad -> SetLogz(0);
		gPad -> SetGridx(0);
		gPad -> SetGridy(0);
		h_2D_Factor[i] -> SetStats(0);
		h_2D_Factor[i] -> GetZaxis() -> SetLabelSize(.015);
		h_2D_Factor[i] -> Draw("colz");
		Gr_SigBest[i] -> SetMarkerSize(2.0);
		Gr_SigBest[i] -> Draw("same, P");
		Lat -> SetTextColor(6);
		Lat -> SetTextSize(0.025);
		Lat -> DrawLatexNDC( .60, .80, TextBestSig[i]);
	}

}

////////////////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1ToF::SaveHist( int NSave){
	cout<<"=====   TriScinti_Ref1ToF::SaveHist   ====="<<endl;
	if(NSave==0){
		NofSaveCanv = NofCanv;
	}else if( (NSave>0)&&(NSave<NofCanv) ){
		NofSaveCanv = NSave;
	}else{
		NofSaveCanv = NofCanv;
	}
	
//	figure       = Form("run%03d_TriScinti_Ref1ToF_%02d.pdf", RunNum, ItNum );
	figure       = Form("TriScintiAna/TriScinti_Ref1ToF_run%03d_%02d.pdf", RunNum, ItNum );
	figure_open  = PATH_fig_General+figure+"[";
	figure_close = PATH_fig_General+figure+"]";
	figure       = PATH_fig_General+figure;
	Ca[0] -> Print(figure_open,  "pdf");
	for(int i=0; i<NofSaveCanv; i++){Ca[i]->Print(figure);}
	Ca[NofSaveCanv-1] -> Print(figure_close,  "pdf");
}

////////////////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1ToF::ExportTwlk(){
	cout<<"=====   TriScinti_Ref1ToF::ExportTwlk   ====="<<endl;
	string ExParam = Form("TESTparamTriScintiRef1ToF_run%03d_%02d.dat", RunNum, ItNum);
	ExParam = PATH_param_t+ExParam;

	ofstream ofs( ExParam.c_str(), ios_base::trunc );
	ofs<<"BetaNum >>>"<<RunNum <<endl;
	ofs<<"PedeNum >>>"<<PedNum<<endl;
	ofs<<"====="<<endl;
	ofs<<endl;
	ofs<<Form("Ref-1 Pedestal [ch] >> %.2lf", Pede_Cnt_Val[0])<<endl;
	ofs<<Form("ToF-1 Pedestal [ch] >> %.2lf", Pede_Cnt_Val[1])<<endl;
	ofs<<Form("ToF-2 Pedestal [ch] >> %.2lf", Pede_Cnt_Val[2])<<endl;
	ofs<<"====="<<endl;
	ofs<<endl;
	ofs<<"QDC cut condition"<<endl;
	ofs<<Form("Ref-1 QDCpeak [pC] >> %.2lf", Qpeak[0]);
	ofs<<Form("Ref-1 QDC Lower Limit [pC] >> %.2lf", QCut[0][0])<<endl;
	ofs<<Form("Ref-1 QDC Upper Limit [pC] >> %.2lf", QCut[0][1])<<endl;
	ofs<<Form("ToF-1 QDCpeak [pC] >> %.2lf", Qpeak[1])<<endl;
	ofs<<Form("ToF-1 QDC Lower Limit [pC] >> %.2lf", QCut[1][0])<<endl;
	ofs<<Form("ToF-1 QDC Upper Limit [pC] >> %.2lf", QCut[1][1])<<endl;
	ofs<<Form("ToF-2 QDCpeak [pC] >> %.2lf", Qpeak[2])<<endl;
	ofs<<Form("ToF-2 QDC Lower Limit [pC] >> %.2lf", QCut[2][0])<<endl;
	ofs<<Form("ToF-2 QDC Upper Limit [pC] >> %.2lf", QCut[2][1])<<endl;
	ofs<<"====="<<endl;
	ofs<<endl;
	ofs<<"Twlk FitRange"<<endl;
	ofs<<Form("Ref-1 FitMax >> %.2lf", TwlkFitMax[0])<<endl;
	ofs<<Form("ToF-1 FitMax >> %.2lf", TwlkFitMax[2])<<endl;
	ofs<<Form("ToF-2 FitMax >> %.2lf", TwlkFitMax[3])<<endl;
	ofs<<"====="<<endl;
	ofs<<endl;
	ofs<<"Twlk Param (par[0], par[1], par[2], Factor)"<<endl;
	ofs<<"Ref-1"<<endl;
	for(int i=0; i<3; i++){ofs<<TwlkPar_Val[0][i]<<endl;}
	ofs<<Twlk_FactorBestSigma[0]<<endl;
	ofs<<endl;
	ofs<<"ToF-1"<<endl;
	for(int i=0; i<3; i++){ofs<<TwlkPar_Val[1][i]<<endl;}
	ofs<<Twlk_FactorBestSigma[1]<<endl;
	ofs<<endl;
	ofs<<"ToF-2"<<endl;
	for(int i=0; i<3; i++){ofs<<TwlkPar_Val[2][i]<<endl;}
	ofs<<Twlk_FactorBestSigma[2]<<endl;
	ofs<<"====="<<endl;
	ofs<<endl;
	ofs<<"ToF offset [psec]"<<endl;
	ofs<<ToF_Adj<<endl;
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv ){
	int Flag=1;	
	int RNum=0;
	int PNum=0;
	int FuncType=0;

	TApplication *theApp;
	TriScinti_Ref1ToF *analysis;

	theApp = new TApplication("App", &argc, argv );
	analysis = new TriScinti_Ref1ToF();

	/*=====	Set ROOT file	=====*/
	cout<<"====="<<endl;
	if(argc==3){
		RNum = atoi(argv[1]);
		PNum = atoi(argv[2]);
		FuncType = 0;
		cout<<"Target of Analysis >>> "<<RNum<<endl;
		cout<<"Pedestal Run  >>> "<<PNum<<endl;
	}else if(argc==4){
		RNum     = atoi(argv[1]);
		PNum     = atoi(argv[2]);
		FuncType = atoi(argv[3]);
		cout<<"Target of Analysis >>> "<<RNum<<endl;
		cout<<"Pedestal Run  >>> "<<PNum<<endl;
		cout<<"ToF Time walk function type  >>> "<<FuncType<<endl;
	}else{
		RNum = 220;
		PNum = 223;
		FuncType = 0;
		cout<<"Target of Analysis >>> "<<RNum<<endl;
		cout<<"Pedestal Run  >>> "<<PNum<<endl;
		cout<<"ToF Time walk function type  >>> "<<FuncType<<endl;
	};
	
	analysis->SetRoot(RNum);
	analysis->TwlkMan(FuncType);
	cout<<"-----aaaAAAAA0AAAAAaaa-----"<<endl;
	analysis->DefineCanv();
	cout<<"-----bbbBBBBB0BBBBBbbb-----"<<endl;
	analysis->DefineHistFunc();
	cout<<"-----cccCCCCC0CCCCCccc-----"<<endl;
	analysis->ImportPedestal(PNum);
	cout<<"-----dddDDDDD0DDDDDddd-----"<<endl;
	analysis->FillHist();
	analysis->GetQDCCut();
	analysis->TwlkFit();
	analysis->TwlkSerchBestFactorToF1vsToF2();
	cout<<"-----oooOOOOO0OOOOOooo-----"<<endl;
	analysis->TwlkSerchBestFactorRef1vsToF1();
	cout<<"-----oooOOOOO0OOOOOooo-----"<<endl;
	analysis->TwlkSerchBestFactorRef1vsToF2();
	cout<<"-----oooOOOOO0OOOOOooo-----"<<endl;
	analysis->TwlkDecToF();
	analysis->DrawHist();
	analysis->SaveHist(99);
//	analysis->SaveHist(3);
	analysis->ExportTwlk();

	//====	END of Job	====//
	if( Flag==1 ){
		gSystem->Exit(1);
	}else;

	theApp->Run();
	return 0;
}
/*
 */
