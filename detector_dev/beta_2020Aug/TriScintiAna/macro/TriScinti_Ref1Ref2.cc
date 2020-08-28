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

#include "./TriScinti_Ref1Ref2.hh"
#include "../Setting.h"


double TriScinti_Ref1Ref2::Twlk_Correction( double x, double* parameter, double Coeff ){
	double Twlk_func=0.;
	Twlk_func = parameter[0]/sqrt( x-parameter[1] );
	Twlk_func = Coeff*Twlk_func;
	return Twlk_func;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
TriScinti_Ref1Ref2::TriScinti_Ref1Ref2(){
	cout<<"__________ooooooooooOOOOO000OOOOOoooooooooo__________"<<endl;
	cout<<"TriScinti_Ref1Ref2::Constructer called"<<endl;
	
	ItNum = 4;

	Set = new Setting();
	Set -> Setting_Gene(1);
	gStyle -> SetStatX(.875);
//	gStyle -> SetStatW(.200);
//	gStyle -> SetStatH(.300);

	for(int i=0; i<NofRef; i++){
		Ln_Qpeak[i]   = new TLine();Set->Setting_Line( Ln_Qpeak[i]  , 6  , 1, 1 );
		Ln_Qcut[i][0] = new TLine();Set->Setting_Line( Ln_Qcut[i][0], 602, 1, 7 );	//Upper-side
		Ln_Qcut[i][1] = new TLine();Set->Setting_Line( Ln_Qcut[i][1], 602, 1, 7 );	//Downer-side
	}

	Lat = new TLatex();
	Set -> Setting_Latex( Lat, 62, 22, 6, 0.075 );

	for(int i=0; i<NofMPPC; i++){
		TDCcut[i][0] = 1900.;
		TDCcut[i][1] = 2300.;
	}
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
TriScinti_Ref1Ref2::~TriScinti_Ref1Ref2(){
	cout<<"TriScinti_Ref1Ref2::Destructer called"<<endl;
	cout<<"__________ooooooooooOOOOO000OOOOOoooooooooo__________"<<endl;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1Ref2::SetRoot( int RNum ){
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

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1Ref2::ImportRefPedestal( int PNum ){
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
/*
 */

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1Ref2::DefineCanv(){
	for(int i=0; i<NofCanv; i++){
		Ca_Ref[i] = new TCanvas( Form("Ca_Ref[%d]", i), Form("Ca_Ref[%d]", i), 1800., CanvHsize[i] );
	}
	Ca_Ref[0] -> Divide(2,2);
	Ca_Ref[1] -> Divide(1,2);
	Ca_Ref[2] -> Divide(1,2);
	Ca_Ref[3] -> Divide(1,2);
	Ca_Ref[4] -> Divide(1,1);
	Ca_Ref[5] -> Divide(1,1);
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1Ref2::DefineHistFunc(){
	for(int i=0; i<NofRef; i++){
		h_1D_QDCFull[i] = new TH1D( Form("h_1D_QDCFull[%d]", i), Form(" h_1D_QDCFull[%d]", i), 400, 0., 400.  );
		h_1D_QDCCut[i] = new TH1D( Form("h_1D_QDCCut[%d]", i), Form(" h_1D_QDCCut[%d]", i), 400, 0., 400.  );

		Set -> Setting_Hist1D( h_1D_QDCFull[i], Label[i][0], Label[i][0]+" [pC]", "Counts/0.5 pC", LColor, 1, 62, FColor, 3001 );
		Set -> Setting_Hist1D( h_1D_QDCCut[i], Label[i][0], Label[i][0]+" [pC]", "Counts/0.5 pC", LColor, 1, 62, kYellow-10, 1001 );

		h_2D_RawTQ[i] = new TH2D( Form("h_2D_RawTQ[%d]", i), Form("h_2D_RawTQ[%d]", i), 400, 0., 400., 140, -4900., 4900. );
		h_2D_FitTQ[i] = new TH2D( Form("h_2D_FitTQ[%d]", i), Form("h_2D_FitTQ[%d]", i), 400, 0., 400., 140, -4900., 4900. );
		h_2D_DecTQ[i] = new TH2D( Form("h_2D_DecTQ[%d]", i), Form("h_2D_DecTQ[%d]", i), 400, 0., 400., 140, -4900., 4900. );
		h_2D_CutTQ[i] = new TH2D( Form("h_2D_CutTQ[%d]", i), Form("h_2D_CutTQ[%d]", i), 400, 0., 400., 140, -4900., 4900. );

		Set -> Setting_Hist2D( h_2D_RawTQ[i] , Form("Ref-%d.QDC vs. RawToF", i+1)            , Label[i][0]+" [pC]", "RawToF [psec]", "", 1. );
		Set -> Setting_Hist2D( h_2D_FitTQ[i] , Form("Ref-%d.QDC vs. RawToF w/ QDC cut", i+1) , Label[i][0]+" [pC]", "RawToF [psec]", "", 1. );
		Set -> Setting_Hist2D( h_2D_DecTQ[i] , Form("Ref-%d.QDC vs. DecToF w/o QDC cut", i+1), Label[i][0]+" [pC]", "DecToF [psec] (Peak Adjusted)", "", 1. );
		Set -> Setting_Hist2D( h_2D_CutTQ[i] , Form("Ref-%d.QDC vs. RawToF w/ QDC cut", i+1) , Label[i][0]+" [pC]", "DecToF [psec] (Peak Adjusted)", "", 1. );

		f_Twlk[i] = new TF1( Form("f_Twlk[%d]", i), "[0]/sqrt(x-[1])+[2]", 0., 250.);
		Set -> Setting_Func( f_Twlk[i], 2, 1 );
	}
	h_1D_RawToF = new TH1D( "h_1D_RawToF", "h_1D_RawToF", 80 , -1400., 1400. );
	h_1D_DecToF = new TH1D( "h_1D_DecToF", "h_1D_DecToF", 80, -1400., 1400. );
	h_2D_Factor = new TH2D( "h_2D_Factor", "h_2D_Factor", 71 , -0.0125, 1.7675, 71, -0.0125, 1.7675 );
	
	Set -> Setting_Hist1D( h_1D_DecToF, "dec.ToF (After TWC)", "ToF [psec]", "Counts/35 psec", LColor, 1, 62, FColor, 3001 );
	Set -> Setting_Hist2D( h_2D_Factor, "Ref-1 vs. Ref-2", "Factor Ref-1", "Factor Ref-2", "#sigma_{Ref1-Ref2} [psec]", 1. );
	h_2D_Factor->GetXaxis()->SetNdivisions(515);
	h_2D_Factor->GetYaxis()->SetNdivisions(515);
	h_2D_Factor->GetXaxis()->SetLabelSize(.025);
	h_2D_Factor->GetYaxis()->SetLabelSize(.025);
	h_2D_Factor->GetZaxis()->SetLabelSize(.025);


	Gr_SigBest = new TGraph();
	Set -> Setting_Graph( Gr_SigBest,"Ref-1 vs. Ref-2", "Factor Ref-1", "Factor Ref-2", 1, 1, 62, 6, 29 );
	Gr_SigBest -> SetMarkerSize(4.00);

	f_ToF = new TF1( "f_ToF", "gaus(0)", -6000., 6000.);
	Set -> Setting_Func( f_ToF, 2, 1 );
	f_ToF -> SetNpx(5E+4);

	f_Twlk[0] -> SetParameters( 2.E+3, 0., -1000.);
	f_Twlk[1] -> SetParameters(-2.E+3, 0.,  1000.);
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1Ref2::FillHist(){
	for(int i=0; i<Entry; i++){
		tree->GetEntry(i);
		if(
			   TDC[0]>0
			&& TDC[1]>0
			&& TDC[2]>0
			&& TDC[3]>0
			&& TDC[0]>TDCcut[0][0]
			&& TDC[0]<TDCcut[0][1]
			&& TDC[1]>TDCcut[1][0]
			&& TDC[1]<TDCcut[1][1]
			&& TDC[2]>TDCcut[2][0]
			&& TDC[2]<TDCcut[2][1]
			&& TDC[3]>TDCcut[3][0]
			&& TDC[3]<TDCcut[3][1]
			&& QDC[0]<4095
			&& QDC[1]<4095
		){
			rawtof = 35.*((double)TDC[1]-(double)TDC[0]);
			for(int j=0; j<NofRef; j++){
				h_1D_QDCFull[j] -> Fill( 0.1*(double)QDC[j]-0.1*Pede_Cnt_Val[j] );
				h_2D_RawTQ[j]   -> Fill( 0.1*(double)QDC[j]-0.1*Pede_Cnt_Val[j], rawtof ); 
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1Ref2::GetQDCCut(){
//	if(RunNum==204){
//		RefQDCGain[0] = h_1D_QDCFull[0] -> GetBinCenter( h_1D_QDCFull[0]->GetMaximumBin() );
//		RefQDCGain[1] = 200.; 
//	}else{
//		for(int i=0; i<NofRef; i++){
//			RefQDCGain[i]   = h_1D_QDCFull[i]->GetBinCenter( h_1D_QDCFull[i]->GetMaximumBin() );
//			cout<<RefQDCGain[i]<<endl;
//		}
//	}
//	for(int i=0; i<NofRef; i++){
//		RefQDCCut[i][0] = 0.50*RefQDCGain[i]+0.1*Pede_Cnt_Val[i];
//		RefQDCCut[i][1] = 1.75*RefQDCGain[i]+0.1*Pede_Cnt_Val[i];
//	}

//	RefQDCCut[0][0] = 60.;
	RefQDCCut[0][0] = 50.;
	RefQDCCut[0][1] = 200.;
//	RefQDCCut[1][0] = 100.;
	RefQDCCut[1][0] = 50.;
//	RefQDCCut[1][1] = 300.;
	RefQDCCut[1][1] = 200.;

	for(int i=0; i<Entry; i++){
		tree->GetEntry(i);
		for(int j=0; j<NofRef; j++){
			if( (0.1*QDC[j] > RefQDCCut[j][0]) && (0.1*QDC[j] < RefQDCCut[j][1]) ){
				h_1D_QDCCut[j] -> Fill( 0.1*(double)QDC[j]-0.1*Pede_Cnt_Val[j] );
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1Ref2::TwlkFit(){
//	cout<<"___oooOOO0OOOooo___"<<endl;
	for(int i=0; i<Entry; i++){
		tree->GetEntry(i);
		for(int j=0; j<NofRef; j++){
			if( (0.1*QDC[NofRef-1-j] > RefQDCCut[NofRef-1-j][0]) && (0.1*QDC[NofRef-1-j] < RefQDCCut[NofRef-1-j][1]) ){
				rawtof = 35.*((double)TDC[1]-(double)TDC[0]);
				h_2D_FitTQ[j] -> Fill( 0.1*(double)QDC[j]-0.1*Pede_Cnt_Val[j], rawtof );
			}
		}
	}
//	double TwlkFitMax[NofRef] = {90., 100.};
	TwlkFitMax[0] = 400.;
	TwlkFitMax[1] = 400.;
	
	for(int i=0; i<NofRef; i++){
		Ca_Ref[1]->cd(i+1);
		gPad->SetLogz(1);
		gPad->SetRightMargin(.125);
		h_2D_FitTQ[i] -> Draw("colz");
		gPad -> Update();
		h_2D_FitTQ[i] -> Fit( f_Twlk[i], "", "", 10., TwlkFitMax[i] );
		gPad -> SetYstat(.80);
		gPad -> Update();
		for(int j=0; j<3; j++){
			TwlkPar_Val[i][j] = f_Twlk[i]->GetParameter(j);
			TwlkPar_Err[i][j] = f_Twlk[i]->GetParError(j);
		}
	}


}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1Ref2::TwlkSerchBestFactor(){
//	cout<<"___oooOOO0OOOooo___"<<endl;
	for(int i=0; i<71; i++){
		for(int j=0; j<71; j++){
			TwlkFactor[0] = 0.000+0.025*i;
			TwlkFactor[1] = 0.000+0.025*j;
		
			h_1D_DecToF_tmp = new TH1D( "h_1D_DecToF_tmp", "h_1D_DecToF_tmp", 300, -5250., 5250. );
			
			for(int k=0; k<Entry; k++){
				tree->GetEntry(k);
				if(
					   0.1*QDC[0] > RefQDCCut[0][0]
					&& 0.1*QDC[1] > RefQDCCut[1][0]
				){
					dectof = 35.*((double)TDC[1]-(double)TDC[0])
					         -Twlk_Correction( 0.1*QDC[0]-0.1*Pede_Cnt_Val[0], TwlkPar_Val[0], TwlkFactor[0] )
					         -Twlk_Correction( 0.1*QDC[1]-0.1*Pede_Cnt_Val[1], TwlkPar_Val[1], TwlkFactor[1] );
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

			if( Twlk_SigVal_tmp<Twlk_SigVal_Best ){
				Twlk_SigVal_Best = Twlk_SigVal_tmp;
				Twlk_SigErr_Best = Twlk_SigErr_tmp;
				Twlk_FactorBestSigma[0] = TwlkFactor[0];
				Twlk_FactorBestSigma[1] = TwlkFactor[1];
				cout<<"Preriminary factor 1 and 2 >> "<<Twlk_FactorBestSigma[0]<<",  "<<Twlk_FactorBestSigma[1]<<":  "<<"Sigma_TOF = "<<Twlk_SigVal_Best<<" psec"<<endl;
			}else;
			h_2D_Factor -> Fill( TwlkFactor[0], TwlkFactor[1], Twlk_SigVal_tmp );
		}
	}
	Gr_SigBest     -> SetPoint( 0, Twlk_FactorBestSigma[0] , Twlk_FactorBestSigma[1]  );
	h_2D_Factor -> SetMinimum( Twlk_SigVal_Best );
	for(int i=0; i<NofRef; i++){TwlkFactor[i] = Twlk_FactorBestSigma[i];}
}
/*
 
 */

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1Ref2::TwlkDecToF(){
	h_1D_DecToF_tmp = new TH1D( "h_1D_DecToF_tmp", "h_1D_DecToF_tmp", 300, -5250., 5250. );
	canv_dummy_consgaus = new TCanvas( "canv_dummy_consgaus", "canv_dummy_consgaus", 480, 240);
	canv_dummy_consgaus -> cd();
	for(int i=0; i<Entry; i++){
		tree -> GetEntry(i);
		dectof = 35.*((double)TDC[1]-(double)TDC[0])
				 -Twlk_Correction( 0.1*QDC[0]-0.1*Pede_Cnt_Val[0], TwlkPar_Val[0], TwlkFactor[0] )
				 -Twlk_Correction( 0.1*QDC[1]-0.1*Pede_Cnt_Val[1], TwlkPar_Val[1], TwlkFactor[1] );
		h_1D_DecToF_tmp -> Fill(dectof);
	}
	h_1D_DecToF_tmp -> Draw();
	ToF_Adj = h_1D_DecToF_tmp -> GetBinCenter( h_1D_DecToF_tmp->GetMaximumBin() );
	h_1D_DecToF_tmp -> Delete();
	canv_dummy_consgaus -> Destructor();

	for(int i=0; i<Entry; i++){
		tree -> GetEntry(i);
		dectof = 35.*((double)TDC[1]-(double)TDC[0])
				 -Twlk_Correction( 0.1*QDC[0]-0.1*Pede_Cnt_Val[0], TwlkPar_Val[0], TwlkFactor[0] )
				 -Twlk_Correction( 0.1*QDC[1]-0.1*Pede_Cnt_Val[1], TwlkPar_Val[1], TwlkFactor[1] )
		         -ToF_Adj;
		for(int j=0; j<NofRef; j++){h_2D_DecTQ[j]->Fill(0.1*(double)QDC[j]-0.1*Pede_Cnt_Val[j], dectof );}
		if(
			   0.1*QDC[0] > RefQDCCut[0][0]
			&& 0.1*QDC[1] > RefQDCCut[1][0]
			&& 0.1*QDC[1] < RefQDCCut[1][1]
			&& QDC[0]<4095
			&& QDC[1]<4095
			&& TDC[0]>TDCcut[0][0]
			&& TDC[0]<TDCcut[0][1]
			&& TDC[1]>TDCcut[1][0]
			&& TDC[1]<TDCcut[1][1]
			&& TDC[2]>TDCcut[2][0]
			&& TDC[2]<TDCcut[2][1]
			&& TDC[3]>TDCcut[3][0]
			&& TDC[3]<TDCcut[3][1]
		){
			for(int j=0; j<NofRef; j++){h_2D_CutTQ[j]->Fill(0.1*(double)QDC[j]-0.1*Pede_Cnt_Val[j], dectof );}
			h_1D_DecToF -> Fill(dectof);
		}else;
	}
	
	
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1Ref2::DrawHist(){
	for(int i=0; i<NofRef; i++){
		Ca_Ref[0]->cd(i+1);
		gPad->SetLogz(1);
		gPad->SetRightMargin(.125);
		h_2D_RawTQ[i]->Draw("colz");
	}
	for(int i=0; i<NofRef; i++){
		Ca_Ref[0]->cd(i+3);
		gPad->SetLogy(1);
		gPad->SetRightMargin(.125);
		h_1D_QDCFull[i]->Draw("");
		gPad -> Update();
		h_1D_QDCCut[i] ->Draw("same");
	//	Lat -> DrawLatexNDC( .70, .50, Form("QDC peak: %.1lf [pC]", RefQDCGain[i]) );
	}

	for(int i=0; i<NofRef; i++){
		Ca_Ref[2]->cd(i+1);
		gPad->SetLogz(1);
		gPad->SetRightMargin(.125);
		h_2D_DecTQ[i] -> Draw("colz");
	}

	for(int i=0; i<NofRef; i++){
		Ca_Ref[3]->cd(i+1);
		gPad->SetLogz(1);
		gPad->SetRightMargin(.125);
		h_2D_CutTQ[i] -> Draw("colz");
	}
	
	Ca_Ref[4]->cd(1);
	gPad->SetLogy(0);
	h_1D_DecToF->Draw("");
	h_1D_DecToF->Fit(f_ToF, "", "", -2100., 2100. );
	Twlk_ResoVal = f_ToF->GetParameter(2);
	Twlk_ResoErr = f_ToF->GetParError(2);
	Lat -> SetTextSize(0.045);
	Lat -> DrawLatexNDC( .75, .45, Form("#sigma_{Ref1-Ref2} = %.2lf #pm %.2lf [psec]", Twlk_ResoVal, Twlk_ResoErr ) );

	Ca_Ref[5] -> cd(1);
	gPad -> SetRightMargin(.125);
	gPad -> SetTopMargin(.090);
	gPad -> SetBottomMargin(.125);
	h_2D_Factor -> SetStats(0);
	h_2D_Factor -> Draw("colz");
	Gr_SigBest->Draw("same, P");
	Lat -> SetTextColor(6);
	Lat -> SetTextSize(0.025);
	Lat -> DrawLatex( Twlk_FactorBestSigma[0], Twlk_FactorBestSigma[1]+.050, 
		               Form("Best factor = %.3lf, %.3lf", Twlk_FactorBestSigma[0], Twlk_FactorBestSigma[1])
		            );
	cout<<"++++++++++++++++++++++++"<<endl;
	cout<<"Best Timing Resolusion: "<<Form("%.2lf", Twlk_ResoVal)<<" +/- "<<Form("%.2lf", Twlk_ResoErr)<<" psec"<<endl;
	cout<<"++++++++++++++++++++++++"<<endl;
	
}
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void TriScinti_Ref1Ref2::SaveHist( int NSave ){
	if(NSave==0){
		NSave=1;
	}else if( NSave>NofCanv ){
		NSave=NofCanv;
	}else;
	
	figure       = Form("TriScinti_Ref1Ref2_run%03d_%02d.pdf", RunNum, ItNum);
	figure_open  = PATH_fig_General+figure+"[";
	figure_close = PATH_fig_General+figure+"]";
	figure       = PATH_fig_General+figure;
	Ca_Ref[0] -> Print(figure_open,  "pdf");
//	Ca_Ref[0] -> Print(figure);
//	Ca_Ref[1] -> Print(figure);
//	Ca_Ref[0] -> Print(figure_close, "pdf");
	for(int i=0; i<NSave; i++){Ca_Ref[i] -> Print(figure);}
	Ca_Ref[NSave-1] -> Print(figure_close, "pdf");
}
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv){
	int Num;
	TApplication *theApp;
	TriScinti_Ref1Ref2 *ana;
	
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

	theApp = new TApplication("App", &argc, argv );
	ana = new TriScinti_Ref1Ref2();
	
	ana -> SetRoot(Num);
	ana -> ImportRefPedestal(Num);
	ana -> DefineCanv();
	ana -> DefineHistFunc();
	ana -> FillHist();
	ana -> GetQDCCut();
	ana -> TwlkFit();
	ana -> TwlkSerchBestFactor();
	ana -> TwlkDecToF();
	ana -> DrawHist();
	ana -> SaveHist(99);
	
	delete ana;
	//====	END of Job	====//
	gSystem->Exit(1);

	theApp->Run();
	return 0;
}

