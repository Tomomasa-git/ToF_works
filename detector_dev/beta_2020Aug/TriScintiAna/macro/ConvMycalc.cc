/* * * * * * * * * * * * *
 * ConvMyCalc.cc         *
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


double ConvLandGaus( double* x, double* par){
	//Reference: 1. Nagao-san's program in .rootmacro.cc (landaugaus)
	//           2. ROOT tutorial program $ROOTSYS/tutorials/fit/langaus.C
	//par[0] >> The area of Landau
	//par[1] >> Most probable value of Landau
	//par[2] >> sigma of Landau
	//par[3] >> sigma of Gaus: Resolusion
	double land;
	double ret;
	// Numeric constants
	double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
	double mpshift  = -0.22278298;       // Landau maximum location

	// Control constants
	double np = 1000.0;					// number of convolution steps
	double sc = 10.0;					// convolution extends to +-sc Gaussian sigmas >> +/- sc*sigma_of_gaus

	// Variables
	double xx;
	double mpc;
	double sum = 0.0;
	double xlow,xupp;
	double step;
	double i;

	// MP shift correction
	mpc = par[1] - mpshift*par[2];

	// Range of convolution integral
	xlow = x[0] - sc*par[3];
	xupp = x[0] + sc*par[3];

	step = (xupp-xlow)/np;

	// Convolution integral of Landau and Gaussian by sum
	for(i=1.0; i<=np/2.; i++){
		xx = xlow * (i-.5)*step;
		land = TMath::Landau( xx, mpc, par[2] )/par[2];
		sum += land*TMath::Gaus( x[0], xx, par[3] );

		xx = xupp * (i-.5)*step;
		land = TMath::Landau( xx, mpc, par[2] )/par[2];
		sum += land*TMath::Gaus( x[0], xx, par[3] );
	}

	ret = par[0]*step*sum*invsq2pi/par[3];

	return ret;
}

double LandGausMyInt( double* x, double* par ){
	//-sekibun
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

	xLowLim = x[0] - sc*par[3];
	xUppLim = x[0] + sc*par[3];

	FuncBoundary[0] = TMath::Landau( xLowLim, mpc, par[2] )*TMath::Gaus( x[0], xLowLim, par[3] );
	FuncBoundary[1] = TMath::Landau( xUppLim, mpc, par[2] )*TMath::Gaus( x[0], xUppLim, par[3] );

	step    = (xUppLim-xLowLim)/np;
//	divsize = 2./np;
	divsize = 1./np;
	factor  = step/6.;

	for(i=1.; i<=(np-1.); i++){
		//For Lower side
		xx = xLowLim + (i*divsize)*(2.*sc*par[3]);
		land = TMath::Landau( xx, mpc, par[2] )/par[2];
		sum_1 += land*TMath::Gaus( x[0], xx, par[3] );
	
		xx = 0.;
		
	//	//For Upper side
	//	xx = xUppLim - (i*divsize)*(2.*sc*par[3]);
	//	land = TMath::Landau( xx, mpc, par[2] )/par[2];
	//	sum_1 += land*TMath::Gaus( x[0], xx, par[3] );
	}

//	for(i=1.; i<=(np/2.); i++){
	for(i=1.; i<=np; i++){
		//For Lower side
		xxx = xLowLim +  ( ((2*i-1)*divsize)/2. )*(2.*sc*par[3]);
		land = TMath::Landau( xxx, mpc, par[2] )/par[2];
		sum_2 += land*TMath::Gaus( x[0], xxx, par[3] );
	
		xx = 0.;
		
		//For Upper side
	//	xxx = xUppLim - ( ((2*i-1)*divsize)/2. )*(2.*sc*par[3]);
	//	land = TMath::Landau( xxx, mpc, par[2] )/par[2];
	//	sum_2 += land*TMath::Gaus( x[0], xxx, par[3] );
	}

	ret = par[0]*factor*(invsqrt2pi/par[3])*( FuncBoundary[0] + FuncBoundary[1] + 2.*sum_1 + 4.*sum_2 );

	return ret;
}


void conv_test(){
	double area=10000.;
	double landmpv = 1.;
	double landsig = 1.0;
	double gaussig[5] = { 0.5, 1.0, 2.5, 5.0, 7.5 };
	int lc[6] = { 602, 2, 3, 870, 6, 800 };
	TCanvas *Ca;
	TH1D *h_raw;
	TH1D *h_con;
	TF1 *f_raw;
	TF1 *f_con[6];
	TGraph *g[6];
	TPaveText *st[6];
	double x1[6] = { .58, .72, .86, .58, .72, .86 };
	double x2[6] = { .71, .85, .99, .71, .85, .99 };
	double y1[6] = { .55, .55, .55, .05, .05, .05 };
	double y2[6] = { .95, .95, .95, .45, .45, .45 };

	Ca = new TCanvas( "Ca", "Ca", 1500, 720 );
//	Ca->Divide(2,1);
//	gStyle->SetPadRightMargin(0.5);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);

	f_raw = new TF1( "f_raw", "landau(0)", -10., 50. );
	f_raw -> SetParameters( 10000., landmpv, landsig );
	f_raw -> SetLineColor(lc[0]);
	f_raw -> SetLineWidth(1);
	f_raw -> SetNpx(1E+3);
	g[0] = new TGraph(f_raw);

	for(int i=0; i<5; i++){
		f_con[i] = new TF1( Form("f_con[%d]",i), LandGausMyInt, -10., 50., 4 );
		f_con[i] -> SetParameters( area, landmpv, landsig, gaussig[i] );
		f_con[i] -> SetLineColor(lc[i+1]);
		f_con[i] -> SetLineWidth(1);
		f_con[i] -> SetNpx(1E+3);
		g[i+1] = new TGraph(f_con[i]);
	}

	for(int i=0; i<6; i++){
		st[i] = new TPaveText();
		st[i] -> SetX1NDC(x1[i]);
		st[i] -> SetX2NDC(x2[i]);
		st[i] -> SetY1NDC(y1[i]);
		st[i] -> SetY2NDC(y2[i]);
		st[i] -> SetTextFont(42);
		st[i] -> SetTextColor(lc[i]);
		st[i] -> SetFillColor(0);
		st[i] -> SetFillStyle(0);
	//	st[i] -> SetFillColorAlpha( 19, .45);
		st[i] -> SetLineColor(lc[i]);
	}
	st[0] -> AddText( Form("p_{0} = %.1lf",f_raw->GetParameter(0)) );
	st[0] -> AddText( Form("p_{1} = %.1lf",landmpv               ) );
	st[0] -> AddText( Form("p_{2} = %.1lf",landsig               ) );
	for(int i=0; i<5; i++){
		st[i+1] -> AddText( Form("p_{0} = %.1lf", area      ) );
		st[i+1] -> AddText( Form("p_{1} = %.1lf", landmpv   ) );
		st[i+1] -> AddText( Form("p_{2} = %.1lf", landsig   ) );
		st[i+1] -> AddText( Form("p_{3} = %.1lf", gaussig[i]) );
	}

//	Ca->cd(1);
	gPad->SetTicks(1, 1);
	gPad->SetRightMargin(.45);
	gPad->SetLeftMargin(.065);
	gPad->SetLogy(0);
	gPad->Update();
	g[0]->GetXaxis()->SetNdivisions(515);
	g[0]->GetXaxis()->SetLabelFont(62);
	g[0]->GetYaxis()->SetLabelFont(62);
	g[0]->GetXaxis()->SetLimits(-10., 40. );
	g[0]->Draw("AL");
	for(int i=1; i<=5; i++){g[i]->Draw("sameL");}
	for(int i=0; i<6; i++){st[i]->Draw("same");}

	Ca->Print("../FIG/ConvMyCalc_08.pdf");
	
	gSystem->Exit(-1);
	
}

void ConvFit(){
	gStyle->SetOptFit(1111);

	TCanvas *Ca;
	TH1D *h_raw;
	TH1D *h_conv;
	TF1 *f_conv;
	int ItMax = 5E+4;
	double x;
	double landsig = 2.0;
	double landmpv = 10.;
	double gaussig = 2.5;
	
	Ca = new TCanvas( "Ca", "Ca", 1200, 720 );

	h_conv = new TH1D( "h_conv", "h_conv", 1000, -20., 80.);
	h_conv -> SetLineColor(2);
	h_raw = new TH1D( "h_raw"  , "h_raw" , 1000, -20., 80.);
	h_raw -> SetLineColor(602);

	f_conv = new TF1( "f_conv", LandGausMyInt, -10., 80., 4);
	f_conv -> SetNpx(1E+4);
	f_conv -> SetLineColor(870);
	f_conv -> SetLineWidth(2);
	f_conv -> SetParameters( (double)ItMax, landmpv, landsig, gaussig );

	for(int i=0; i<ItMax; i++){
		gRandom->SetSeed();
		x = gRandom->Landau(landmpv, landsig);
		h_raw -> Fill(x);
		x += gRandom->Gaus(0., gaussig);
		h_conv->Fill(x);
	}
	h_raw->Draw("");
	h_conv -> Draw("sameS");

	for(int i=0; i<3; i++){
		h_conv -> Fit( f_conv, "", "QN0", -10., 80.);
	}
	h_conv -> Fit( f_conv, "", "", -10., 80.);

	Ca->Print("../FIG/ConvFit_01.pdf");

	gSystem->Exit(-1);
}

int main( int argc, char** argv ){
	TApplication *theApp;
	theApp = new TApplication("App", &argc, argv );

	conv_test();
//	ConvFit();

	theApp->Run();
	return 0;
}
