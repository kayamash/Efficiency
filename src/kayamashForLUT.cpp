#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include <sstream>
#include <TF1.h>
#include <TTree.h>
#include <TStyle.h>
#include <TText.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <vector>
#include <TMath.h>
#include <TBranch.h>
#include <TObject.h>
#include <TProfile.h>
#include "../PtNewMethod/kayamashForLUT.chh"

bool kayamashForLUT::getLUTparameter(Double_t address,Double_t charge,Double_t eta,Double_t phi,Int_t (&par)[5],Double_t &tmp_phi1){
	for(Int_t i = 0; i < 4; ++i)par[i] = -1;
	if(address == 0.)par[0] = 0;//Large
	if(address == 1. && phi < -1.5)par[0] = 1;//LS sector 11
	if(address == 1. && phi > -1.5)par[0] = 2;//LS sector 15
	if(address == 2.)par[0] = 3;//Small
	if(address == 3.)par[0] = 4;//SS

	if(charge == 1.)par[1] = 0;//positive
	if(charge == -1.)par[1] = 1;//negative

	Int_t tmp_eta = static_cast<Int_t>(std::fabs(eta)*15./1.05);
    par[2] = (eta > 0) ? (tmp_eta) : (tmp_eta + 15);//Eta　divide to 30

    Double_t tmp_phi = phi;
    Int_t tmp_par = -1;
    Int_t sectorNumber = -1;
    if(address == 0 || address == 1){//Large
    	if(-0.4 < phi && 0.4 > phi)sectorNumber = 1;
    	if(0.4 < phi && 1.2 > phi){
    		tmp_phi -= TMath::Pi()/4.;
    		sectorNumber = 3;
    	}
    	if(1.2 < phi && 2.0 > phi){
    		tmp_phi -= TMath::Pi()/2.;
    		sectorNumber = 5;
    	}
    	if(2.0 < phi && 2.8 > phi){
    		tmp_phi -= 3*TMath::Pi()/4.;
    		sectorNumber = 7;
    	}
    	if(2.8 < phi){
    		tmp_phi -= TMath::Pi();
    		sectorNumber = 9;
    	}
    	if(-0.4 > phi && -1.2 < phi){
    		tmp_phi += TMath::Pi()/4.;
    		sectorNumber = 15;
    	}
    	if(-1.2 > phi && -2.0 < phi){
    		tmp_phi += TMath::Pi()/2.;
    		sectorNumber = 13;
    	}
    	if(-2.0 > phi && -2.8 < phi){
    		tmp_phi += 3*TMath::Pi()/4.;
    		sectorNumber = 11;
    	}
    	if(-2.8 > phi){
    		tmp_phi += TMath::Pi();
    		sectorNumber = 9;
    	}
    	/*
    	tmp_par = static_cast<Int_t>((tmp_phi+0.30)*28./0.50) + 1; //Phi divide to 30
    	if(tmp_phi <= -0.25)tmp_par = 0;
    	if(tmp_phi >= 0.25)tmp_par = 29;
    	*/
    	Double_t dividePhi = 0.40/30.;
    	for(Int_t loop = 0; loop < 30;++loop){
    		if(loop == 0 && tmp_phi < -0.2)tmp_par = loop;
    		if(loop != 0 && loop != 29 && tmp_phi >= -0.2 + static_cast<Double_t>(loop - 1)*dividePhi && tmp_phi < -0.2 + static_cast<Double_t>(loop)*dividePhi)tmp_par = loop;
    		if(loop == 29 && tmp_phi >= 0.2)tmp_par = loop;
    	}
    }else if(address == 2 || address == 3){//Small
    	if(0 < phi && 0.8 > phi)sectorNumber = 2;
    	if(0.8 < phi && 1.6 > phi){
    		tmp_phi -= TMath::Pi()/4.;
    		sectorNumber = 4;
    	}
    	if(1.6 < phi && 2.4 > phi){
    		tmp_phi -= TMath::Pi()/2.;
    		sectorNumber = 6;
    	}
    	if(2.4 < phi){
    		tmp_phi -= 3*TMath::Pi()/4.;
    		sectorNumber = 8;
    	}
    	if(0 > phi && -0.8 < phi){
    		tmp_phi += TMath::Pi()/4.;
    		sectorNumber = 16;
    	}
    	if(-0.8 > phi && -1.6 < phi){
    		tmp_phi += TMath::Pi()/2.;
    		sectorNumber = 14;
    	}
    	if(-1.6 > phi && -2.4 < phi){
    		tmp_phi += 3*TMath::Pi()/4.;
    		sectorNumber = 12;
    	}
    	if(-2.4 > phi){
    		tmp_phi += TMath::Pi();
    		sectorNumber = 10;
    	}
    	/*
    	tmp_par = static_cast<Int_t>((tmp_phi-0.20)*28./0.40); //Phi divide to 30
    	if(tmp_phi <= 0.20)tmp_par = 0;
    	if(tmp_phi >= 0.60)tmp_par = 29;
    	*/
    	Double_t dividePhi = 0.60/30.;
    	for(Int_t loop = 0; loop < 30;++loop){
    		if(loop == 0 && tmp_phi < 0.2)tmp_par = loop;
    		if(loop != 0 && loop != 29 && tmp_phi >= 0.2 + static_cast<Double_t>(loop - 1)*dividePhi && tmp_phi < 0.2 + static_cast<Double_t>(loop)*dividePhi)tmp_par = loop;
    		if(loop == 29 && tmp_phi >= 0.8)tmp_par = loop;
    	}
    }
    if(tmp_par >= 0 && tmp_par <= 29)par[3] = tmp_par;
    par[4] = sectorNumber;
    tmp_phi1 = tmp_phi;

    if(par[0] >= 0 && par[1] >= 0 && par[2] >= 0 && par[3] >= 0)return kTRUE;
    return kFALSE;
}

bool kayamashForLUT::WriteLUT(TProfile *prof,Int_t par1,Int_t par2,Int_t par3,Int_t par4,bool alpha,std::string filename){
	Double_t fOrder = 0;
	Double_t sOrder = 0;
	Double_t parLimits = 0;
	if(alpha){//alpha
		if(par2 == 0){//charge
			switch(par1){//Sector
				case 0:
				fOrder = -0.0184048;
				sOrder = -3.62884;
				break;
				case 1:
				fOrder = -0.201341;
				sOrder = -2.69654;
				break;
				case 2:
				fOrder = -0.118256;
				sOrder = -2.91538;
				break;
				case 3:
				fOrder = -0.373771;
				sOrder = -4.20018;
				break;
				case 4:
				fOrder = -0.276773;
				sOrder = -5.23416;
				break;
			}
		}else{
			switch(par1){//Sector
				case 0:
				fOrder = 0.147342;
				sOrder = 2.77986;
				break;
				case 1:
				fOrder = -0.0414529;
				sOrder = 4.40735;
				break;
				case 2:
				fOrder = 0.191687;
				sOrder = 2.6177;
				break;
				case 3:
				fOrder = 0.250872;
				sOrder = 3.74082;
				break;
				case 4:
				fOrder = 0.225341;
				sOrder = 3.82399;
				break;
			}
		}
	}else{//beta
		if(par2 == 0){//charge
			switch(par1){//Sector
				case 0:
				fOrder = -0.0225169;
				sOrder = 4.08316;
				break;
				case 1:
				fOrder = -0.0044969;
				sOrder = 1.75247;
				break;
				case 2:
				fOrder = 0.0835498;
				sOrder = 1.01082;
				break;
				case 3:
				fOrder = 0.207164;
				sOrder = 7.40479;
				break;
				case 4:
				fOrder = 0.275761;
				sOrder = 7.10336;
				break;
			}
		}else{
			switch(par1){//Sector
				case 0:
				fOrder = -0.156158;
				sOrder = -2.28754;
				break;
				case 1:
				fOrder = 0.319591;
				sOrder = -4.28552;
				break;
				case 2:
				fOrder = 0.149876;
				sOrder = -2.47795;
				break;
				case 3:
				fOrder = -0.489921;
				sOrder = -3.20891;
				break;
				case 4:
				fOrder = -0.575507;
				sOrder = -2.90169;
				break;
			}
		}
	}
	ofstream ofs;
	ofs.open(filename.c_str(),std::ios::app);
	prof->Draw();
	Double_t p0 = 0;
	Double_t p1 = 0;
	Double_t Ndof = 0;
	Double_t chi = 0;
	Double_t pValue = 0;
	if(prof->GetEntries() != 0){
		TF1 *fitProf = new TF1("fitProf","[0]*x + [1]*x*x",0.,0.3);
		fitProf->SetParameter(0,fOrder);
		fitProf->SetParameter(1,sOrder);
		if(alpha == kTRUE && par2 == 0)fitProf->SetParLimits(2,-10.,0.);
		if(alpha == kTRUE && par2 == 1)fitProf->SetParLimits(2,0.,10.);
		if(alpha == kFALSE && par2 == 0)fitProf->SetParLimits(2,0.,10.);
		if(alpha == kFALSE && par2 == 1)fitProf->SetParLimits(2,-10.,0.);
		prof->Fit(fitProf,"I","",0.,0.26);
		p0 = fitProf->GetParameter(0);
		p1 = fitProf->GetParameter(1);
		chi = fitProf->GetChisquare();
		Ndof = fitProf->GetNDF();
		pValue = TMath::Prob(chi, Ndof);
	}

	ofs<<par1<<"   "<<par2<<"   "<<par3<<"   "<<par4<<"   "<<p0<<"   "<<p1<<std::endl;
	ofs.close();
	return kTRUE;
}

bool kayamashForLUT::ReadLUT(Int_t (&par)[4],std::string lut,Double_t par1,Double_t par2){
	ifstream ifs(lut.c_str());
	bool check = kFALSE;
	while(!ifs.eof()){
		Int_t parameter[4] = {0,0,0,0};
		Double_t parameter1 = 0;
		Double_t parameter2 = 0;
		ifs>>parameter[0]>>parameter[1]>>parameter[2]>>parameter[3]>>parameter1>>parameter2;
		if(parameter[0] == par[0] && parameter[1] == par[1] &&parameter[2] == par[2] &&parameter[3] == par[3]){
			par1 = parameter1;
			par2 = parameter2;
			check = kTRUE;
		}
	}
	if(check)return kTRUE;
	return kFALSE;
}
