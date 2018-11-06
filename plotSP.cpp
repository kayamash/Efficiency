#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <vector>
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
#include <vector>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <TBranch.h>

const string inputfilename = "/gpfs/fs6001/kayamash/dataset/Zmumu300540_hadd.root";
const string outputfilename = "/gpfs/fs6001/kayamash/Mywork/SPplot.root";
const string mesSA = "mu26ivm";

void plotSP(){
	TChain *chain = new TChain("t_tap");
	chain->Add(inputfilename.c_str());

	vector<double>  *SPZ = 0;
	vector<double>  *SPR = 0;

	vector<string>  *mes_name = 0;
	vector<double>  *probe_mesSA_superPointR_BI = 0;
	vector<double>  *probe_mesSA_superPointR_BM = 0;
	vector<double>  *probe_mesSA_superPointR_BO = 0;
	vector<double>  *probe_mesSA_superPointR_EI = 0;
	vector<double>  *probe_mesSA_superPointR_EM = 0;
	vector<double>  *probe_mesSA_superPointR_EO = 0;
	vector<double>  *probe_mesSA_superPointR_EE = 0;
	vector<double>  *probe_mesSA_superPointR_CSC = 0;
	vector<double>  *probe_mesSA_superPointR_BEE = 0;
	vector<double>  *probe_mesSA_superPointR_BME = 0;
	vector<double>  *probe_mesSA_superPointZ_BI = 0;
	vector<double>  *probe_mesSA_superPointZ_BM = 0;
	vector<double>  *probe_mesSA_superPointZ_BO = 0;
	vector<double>  *probe_mesSA_superPointZ_EI = 0;
	vector<double>  *probe_mesSA_superPointZ_EM = 0;
	vector<double>  *probe_mesSA_superPointZ_EO = 0;
	vector<double>  *probe_mesSA_superPointZ_EE = 0;
	vector<double>  *probe_mesSA_superPointZ_CSC = 0;
	vector<double>  *probe_mesSA_superPointZ_BEE = 0;
	vector<double>  *probe_mesSA_superPointZ_BME = 0;
	vector<double>  *probe_mesSA_superPointSlope_BI = 0;
	vector<double>  *probe_mesSA_superPointSlope_BM = 0;
	vector<double>  *probe_mesSA_superPointSlope_BO = 0;
	vector<double>  *probe_mesSA_superPointSlope_EI = 0;
	vector<double>  *probe_mesSA_superPointSlope_EM = 0;
	vector<double>  *probe_mesSA_superPointSlope_EO = 0;
	vector<double>  *probe_mesSA_superPointSlope_EE = 0;
	vector<double>  *probe_mesSA_superPointSlope_CSC = 0;
	vector<double>  *probe_mesSA_superPointSlope_BEE = 0;
	vector<double>  *probe_mesSA_superPointSlope_BME = 0;
	vector<double>  *probe_mesSA_superPointIntercept_BI = 0;
	vector<double>  *probe_mesSA_superPointIntercept_BM = 0;
	vector<double>  *probe_mesSA_superPointIntercept_BO = 0;
	vector<double>  *probe_mesSA_superPointIntercept_EI = 0;
	vector<double>  *probe_mesSA_superPointIntercept_EM = 0;
	vector<double>  *probe_mesSA_superPointIntercept_EO = 0;
	vector<double>  *probe_mesSA_superPointIntercept_EE = 0;
	vector<double>  *probe_mesSA_superPointIntercept_CSC = 0;
	vector<double>  *probe_mesSA_superPointIntercept_BEE = 0;
	vector<double>  *probe_mesSA_superPointIntercept_BME = 0;
	vector<double>  *probe_mesSA_superPointChi2_BI = 0;
	vector<double>  *probe_mesSA_superPointChi2_BM = 0;
	vector<double>  *probe_mesSA_superPointChi2_BO = 0;
	vector<double>  *probe_mesSA_superPointChi2_EI = 0;
	vector<double>  *probe_mesSA_superPointChi2_EM = 0;
	vector<double>  *probe_mesSA_superPointChi2_EO = 0;
	vector<double>  *probe_mesSA_superPointChi2_EE = 0;
	vector<double>  *probe_mesSA_superPointChi2_CSC = 0;
	vector<double>  *probe_mesSA_superPointChi2_BEE = 0;
	vector<double>  *probe_mesSA_superPointChi2_BME = 0;

	TBranch *b_mes_name;
	TBranch *b_probe_mesSA_superPointR_BI;
	TBranch *b_probe_mesSA_superPointR_BM;
	TBranch *b_probe_mesSA_superPointR_BO;
	TBranch *b_probe_mesSA_superPointR_EI;
	TBranch *b_probe_mesSA_superPointR_EM;
	TBranch *b_probe_mesSA_superPointR_EO;
	TBranch *b_probe_mesSA_superPointR_EE;
	TBranch *b_probe_mesSA_superPointR_CSC;
	TBranch *b_probe_mesSA_superPointR_BEE;
	TBranch *b_probe_mesSA_superPointR_BME;
	TBranch *b_probe_mesSA_superPointZ_BI;
	TBranch *b_probe_mesSA_superPointZ_BM;
	TBranch *b_probe_mesSA_superPointZ_BO;
	TBranch *b_probe_mesSA_superPointZ_EI;
	TBranch *b_probe_mesSA_superPointZ_EM;
	TBranch *b_probe_mesSA_superPointZ_EO;
	TBranch *b_probe_mesSA_superPointZ_EE;
	TBranch *b_probe_mesSA_superPointZ_CSC;
	TBranch *b_probe_mesSA_superPointZ_BEE;
	TBranch *b_probe_mesSA_superPointZ_BME;
	TBranch *b_probe_mesSA_superPointSlope_BI;
	TBranch *b_probe_mesSA_superPointSlope_BM;
	TBranch *b_probe_mesSA_superPointSlope_BO;
	TBranch *b_probe_mesSA_superPointSlope_EI;
	TBranch *b_probe_mesSA_superPointSlope_EM;
	TBranch *b_probe_mesSA_superPointSlope_EO;
	TBranch *b_probe_mesSA_superPointSlope_EE;
	TBranch *b_probe_mesSA_superPointSlope_CSC;
	TBranch *b_probe_mesSA_superPointSlope_BEE;
	TBranch *b_probe_mesSA_superPointSlope_BME;
	TBranch *b_probe_mesSA_superPointIntercept_BI;
	TBranch *b_probe_mesSA_superPointIntercept_BM;
	TBranch *b_probe_mesSA_superPointIntercept_BO;
	TBranch *b_probe_mesSA_superPointIntercept_EI;
	TBranch *b_probe_mesSA_superPointIntercept_EM;
	TBranch *b_probe_mesSA_superPointIntercept_EO;
	TBranch *b_probe_mesSA_superPointIntercept_EE;
	TBranch *b_probe_mesSA_superPointIntercept_CSC;
	TBranch *b_probe_mesSA_superPointIntercept_BEE;
	TBranch *b_probe_mesSA_superPointIntercept_BME;
	TBranch *b_probe_mesSA_superPointChi2_BI;
	TBranch *b_probe_mesSA_superPointChi2_BM;
	TBranch *b_probe_mesSA_superPointChi2_BO;
	TBranch *b_probe_mesSA_superPointChi2_EI;
	TBranch *b_probe_mesSA_superPointChi2_EM;
	TBranch *b_probe_mesSA_superPointChi2_EO;
	TBranch *b_probe_mesSA_superPointChi2_EE;
	TBranch *b_probe_mesSA_superPointChi2_CSC;
	TBranch *b_probe_mesSA_superPointChi2_BEE;
	TBranch *b_probe_mesSA_superPointChi2_BME;

	chain->SetBranchAddress("mes_name", &mes_name, &b_mes_name);
	chain->SetBranchAddress("probe_mesSA_superPointR_BI", &probe_mesSA_superPointR_BI, &b_probe_mesSA_superPointR_BI);
	chain->SetBranchAddress("probe_mesSA_superPointR_BM", &probe_mesSA_superPointR_BM, &b_probe_mesSA_superPointR_BM);
	chain->SetBranchAddress("probe_mesSA_superPointR_BO", &probe_mesSA_superPointR_BO, &b_probe_mesSA_superPointR_BO);
	chain->SetBranchAddress("probe_mesSA_superPointR_EI", &probe_mesSA_superPointR_EI, &b_probe_mesSA_superPointR_EI);
	chain->SetBranchAddress("probe_mesSA_superPointR_EM", &probe_mesSA_superPointR_EM, &b_probe_mesSA_superPointR_EM);
	chain->SetBranchAddress("probe_mesSA_superPointR_EO", &probe_mesSA_superPointR_EO, &b_probe_mesSA_superPointR_EO);
	chain->SetBranchAddress("probe_mesSA_superPointR_EE", &probe_mesSA_superPointR_EE, &b_probe_mesSA_superPointR_EE);
	chain->SetBranchAddress("probe_mesSA_superPointR_CSC", &probe_mesSA_superPointR_CSC, &b_probe_mesSA_superPointR_CSC);
	chain->SetBranchAddress("probe_mesSA_superPointR_BEE", &probe_mesSA_superPointR_BEE, &b_probe_mesSA_superPointR_BEE);
	chain->SetBranchAddress("probe_mesSA_superPointR_BME", &probe_mesSA_superPointR_BME, &b_probe_mesSA_superPointR_BME);
	chain->SetBranchAddress("probe_mesSA_superPointZ_BI", &probe_mesSA_superPointZ_BI, &b_probe_mesSA_superPointZ_BI);
	chain->SetBranchAddress("probe_mesSA_superPointZ_BM", &probe_mesSA_superPointZ_BM, &b_probe_mesSA_superPointZ_BM);
	chain->SetBranchAddress("probe_mesSA_superPointZ_BO", &probe_mesSA_superPointZ_BO, &b_probe_mesSA_superPointZ_BO);
	chain->SetBranchAddress("probe_mesSA_superPointZ_EI", &probe_mesSA_superPointZ_EI, &b_probe_mesSA_superPointZ_EI);
	chain->SetBranchAddress("probe_mesSA_superPointZ_EM", &probe_mesSA_superPointZ_EM, &b_probe_mesSA_superPointZ_EM);
	chain->SetBranchAddress("probe_mesSA_superPointZ_EO", &probe_mesSA_superPointZ_EO, &b_probe_mesSA_superPointZ_EO);
	chain->SetBranchAddress("probe_mesSA_superPointZ_EE", &probe_mesSA_superPointZ_EE, &b_probe_mesSA_superPointZ_EE);
	chain->SetBranchAddress("probe_mesSA_superPointZ_CSC", &probe_mesSA_superPointZ_CSC, &b_probe_mesSA_superPointZ_CSC);
	chain->SetBranchAddress("probe_mesSA_superPointZ_BEE", &probe_mesSA_superPointZ_BEE, &b_probe_mesSA_superPointZ_BEE);
	chain->SetBranchAddress("probe_mesSA_superPointZ_BME", &probe_mesSA_superPointZ_BME, &b_probe_mesSA_superPointZ_BME);
	chain->SetBranchAddress("probe_mesSA_superPointSlope_BI", &probe_mesSA_superPointSlope_BI, &b_probe_mesSA_superPointSlope_BI);
	chain->SetBranchAddress("probe_mesSA_superPointSlope_BM", &probe_mesSA_superPointSlope_BM, &b_probe_mesSA_superPointSlope_BM);
	chain->SetBranchAddress("probe_mesSA_superPointSlope_BO", &probe_mesSA_superPointSlope_BO, &b_probe_mesSA_superPointSlope_BO);
	chain->SetBranchAddress("probe_mesSA_superPointSlope_EI", &probe_mesSA_superPointSlope_EI, &b_probe_mesSA_superPointSlope_EI);
	chain->SetBranchAddress("probe_mesSA_superPointSlope_EM", &probe_mesSA_superPointSlope_EM, &b_probe_mesSA_superPointSlope_EM);
	chain->SetBranchAddress("probe_mesSA_superPointSlope_EO", &probe_mesSA_superPointSlope_EO, &b_probe_mesSA_superPointSlope_EO);
	chain->SetBranchAddress("probe_mesSA_superPointSlope_EE", &probe_mesSA_superPointSlope_EE, &b_probe_mesSA_superPointSlope_EE);
	chain->SetBranchAddress("probe_mesSA_superPointSlope_CSC", &probe_mesSA_superPointSlope_CSC, &b_probe_mesSA_superPointSlope_CSC);
	chain->SetBranchAddress("probe_mesSA_superPointSlope_BEE", &probe_mesSA_superPointSlope_BEE, &b_probe_mesSA_superPointSlope_BEE);
	chain->SetBranchAddress("probe_mesSA_superPointSlope_BME", &probe_mesSA_superPointSlope_BME, &b_probe_mesSA_superPointSlope_BME);
	chain->SetBranchAddress("probe_mesSA_superPointIntercept_BI", &probe_mesSA_superPointIntercept_BI, &b_probe_mesSA_superPointIntercept_BI);
	chain->SetBranchAddress("probe_mesSA_superPointIntercept_BM", &probe_mesSA_superPointIntercept_BM, &b_probe_mesSA_superPointIntercept_BM);
	chain->SetBranchAddress("probe_mesSA_superPointIntercept_BO", &probe_mesSA_superPointIntercept_BO, &b_probe_mesSA_superPointIntercept_BO);
	chain->SetBranchAddress("probe_mesSA_superPointIntercept_EI", &probe_mesSA_superPointIntercept_EI, &b_probe_mesSA_superPointIntercept_EI);
	chain->SetBranchAddress("probe_mesSA_superPointIntercept_EM", &probe_mesSA_superPointIntercept_EM, &b_probe_mesSA_superPointIntercept_EM);
	chain->SetBranchAddress("probe_mesSA_superPointIntercept_EO", &probe_mesSA_superPointIntercept_EO, &b_probe_mesSA_superPointIntercept_EO);
	chain->SetBranchAddress("probe_mesSA_superPointIntercept_EE", &probe_mesSA_superPointIntercept_EE, &b_probe_mesSA_superPointIntercept_EE);
	chain->SetBranchAddress("probe_mesSA_superPointIntercept_CSC", &probe_mesSA_superPointIntercept_CSC, &b_probe_mesSA_superPointIntercept_CSC);
	chain->SetBranchAddress("probe_mesSA_superPointIntercept_BEE", &probe_mesSA_superPointIntercept_BEE, &b_probe_mesSA_superPointIntercept_BEE);
	chain->SetBranchAddress("probe_mesSA_superPointIntercept_BME", &probe_mesSA_superPointIntercept_BME, &b_probe_mesSA_superPointIntercept_BME);
	chain->SetBranchAddress("probe_mesSA_superPointChi2_BI", &probe_mesSA_superPointChi2_BI, &b_probe_mesSA_superPointChi2_BI);
	chain->SetBranchAddress("probe_mesSA_superPointChi2_BM", &probe_mesSA_superPointChi2_BM, &b_probe_mesSA_superPointChi2_BM);
	chain->SetBranchAddress("probe_mesSA_superPointChi2_BO", &probe_mesSA_superPointChi2_BO, &b_probe_mesSA_superPointChi2_BO);
	chain->SetBranchAddress("probe_mesSA_superPointChi2_EI", &probe_mesSA_superPointChi2_EI, &b_probe_mesSA_superPointChi2_EI);
	chain->SetBranchAddress("probe_mesSA_superPointChi2_EM", &probe_mesSA_superPointChi2_EM, &b_probe_mesSA_superPointChi2_EM);
	chain->SetBranchAddress("probe_mesSA_superPointChi2_EO", &probe_mesSA_superPointChi2_EO, &b_probe_mesSA_superPointChi2_EO);
	chain->SetBranchAddress("probe_mesSA_superPointChi2_EE", &probe_mesSA_superPointChi2_EE, &b_probe_mesSA_superPointChi2_EE);
	chain->SetBranchAddress("probe_mesSA_superPointChi2_CSC", &probe_mesSA_superPointChi2_CSC, &b_probe_mesSA_superPointChi2_CSC);
	chain->SetBranchAddress("probe_mesSA_superPointChi2_BEE", &probe_mesSA_superPointChi2_BEE, &b_probe_mesSA_superPointChi2_BEE);
	chain->SetBranchAddress("probe_mesSA_superPointChi2_BME", &probe_mesSA_superPointChi2_BME, &b_probe_mesSA_superPointChi2_BME);

	for(Int_t event = 0;event < chain->GetEntries();event++){
		chain->GetEntry(event);
		Double_t buf_probe_mesSA_superPointR_BI = -99999;
		Double_t buf_probe_mesSA_superPointR_BM = -99999;
		Double_t buf_probe_mesSA_superPointR_BO = -99999;
		Double_t buf_probe_mesSA_superPointR_EI = -99999;
		Double_t buf_probe_mesSA_superPointR_EM = -99999;
		Double_t buf_probe_mesSA_superPointR_EO = -99999;
		Double_t buf_probe_mesSA_superPointR_EE = -99999;
		Double_t buf_probe_mesSA_superPointR_CSC = -99999;
		Double_t buf_probe_mesSA_superPointR_BEE = -99999;
		Double_t buf_probe_mesSA_superPointR_BME = -99999;
		Double_t buf_probe_mesSA_superPointZ_BI = -99999;
		Double_t buf_probe_mesSA_superPointZ_BM = -99999;
		Double_t buf_probe_mesSA_superPointZ_BO = -99999;
		Double_t buf_probe_mesSA_superPointZ_EI = -99999;
		Double_t buf_probe_mesSA_superPointZ_EM = -99999;
		Double_t buf_probe_mesSA_superPointZ_EO = -99999;
		Double_t buf_probe_mesSA_superPointZ_EE = -99999;
		Double_t buf_probe_mesSA_superPointZ_CSC = -99999;
		Double_t buf_probe_mesSA_superPointZ_BEE = -99999;
		Double_t buf_probe_mesSA_superPointZ_BME = -99999;
		for(Int_t mes = 0;mes < probe_mesSA_superPointR_BI->size();mes++){
			if(mes_name->at(mes) == mesSA){
				buf_probe_mesSA_superPointR_BI = probe_mesSA_superPointR_BI->at(mes);
				buf_probe_mesSA_superPointR_BM = probe_mesSA_superPointR_BM->at(mes);
				buf_probe_mesSA_superPointR_BO = probe_mesSA_superPointR_BO->at(mes);
				buf_probe_mesSA_superPointR_EI = probe_mesSA_superPointR_EI->at(mes);
				buf_probe_mesSA_superPointR_EM = probe_mesSA_superPointR_EM->at(mes);
				buf_probe_mesSA_superPointR_EO = probe_mesSA_superPointR_EO->at(mes);
				buf_probe_mesSA_superPointR_EE = probe_mesSA_superPointR_EE->at(mes);
				buf_probe_mesSA_superPointR_CSC = probe_mesSA_superPointR_CSC->at(mes);
				buf_probe_mesSA_superPointR_BEE = probe_mesSA_superPointR_BEE->at(mes);
				buf_probe_mesSA_superPointR_BME = probe_mesSA_superPointR_BME->at(mes);
				buf_probe_mesSA_superPointZ_BI = probe_mesSA_superPointZ_BI->at(mes);
				buf_probe_mesSA_superPointZ_BM = probe_mesSA_superPointZ_BM->at(mes);
				buf_probe_mesSA_superPointZ_BO = probe_mesSA_superPointZ_BO->at(mes);
				buf_probe_mesSA_superPointZ_EI = probe_mesSA_superPointZ_EI->at(mes);
				buf_probe_mesSA_superPointZ_EM = probe_mesSA_superPointZ_EM->at(mes);
				buf_probe_mesSA_superPointZ_EO = probe_mesSA_superPointZ_EO->at(mes);
				buf_probe_mesSA_superPointZ_EE = probe_mesSA_superPointZ_EE->at(mes);
				buf_probe_mesSA_superPointZ_CSC = probe_mesSA_superPointZ_CSC->at(mes);
				buf_probe_mesSA_superPointZ_BEE = probe_mesSA_superPointZ_BEE->at(mes);
				buf_probe_mesSA_superPointZ_BME = probe_mesSA_superPointZ_BME->at(mes);
			}
		}
		if(buf_probe_mesSA_superPointR_BI != -99999 && buf_probe_mesSA_superPointZ_BI != -99999){
			SPR.push_back(buf_probe_mesSA_superPointR_BI);
			SPZ.push_back(buf_probe_mesSA_superPointZ_BI);
		}
		if(buf_probe_mesSA_superPointR_BM != -99999 && buf_probe_mesSA_superPointZ_BM != -99999){
			SPR.push_back(buf_probe_mesSA_superPointR_BM);
			SPZ.push_back(buf_probe_mesSA_superPointZ_BM);
		}
		if(buf_probe_mesSA_superPointR_BO != -99999 && buf_probe_mesSA_superPointZ_BO != -99999){
			SPR.push_back(buf_probe_mesSA_superPointR_BO);
			SPZ.push_back(buf_probe_mesSA_superPointZ_BO);
		}
		if(buf_probe_mesSA_superPointR_EI != -99999 && buf_probe_mesSA_superPointZ_EI != -99999){
			SPR.push_back(buf_probe_mesSA_superPointR_EI);
			SPZ.push_back(buf_probe_mesSA_superPointZ_EI);
		}
		if(buf_probe_mesSA_superPointR_EM != -99999 && buf_probe_mesSA_superPointZ_EM != -99999){
			SPR.push_back(buf_probe_mesSA_superPointR_EM);
			SPZ.push_back(buf_probe_mesSA_superPointZ_EM);
		}
		if(buf_probe_mesSA_superPointR_EO != -99999 && buf_probe_mesSA_superPointZ_EO != -99999){
			SPR.push_back(buf_probe_mesSA_superPointR_EO);
			SPZ.push_back(buf_probe_mesSA_superPointZ_EO);
		}
		if(buf_probe_mesSA_superPointR_EE != -99999 && buf_probe_mesSA_superPointZ_EE != -99999){
			SPR.push_back(buf_probe_mesSA_superPointR_EE);
			SPZ.push_back(buf_probe_mesSA_superPointZ_EE);
		}
		if(buf_probe_mesSA_superPointR_CSC != -99999 && buf_probe_mesSA_superPointZ_CSC != -99999){
			SPR.push_back(buf_probe_mesSA_superPointR_CSC);
			SPZ.push_back(buf_probe_mesSA_superPointZ_CSC);
		}
		if(buf_probe_mesSA_superPointR_BEE != -99999 && buf_probe_mesSA_superPointZ_BEE != -99999){
			SPR.push_back(buf_probe_mesSA_superPointR_BEE);
			SPZ.push_back(buf_probe_mesSA_superPointZ_BEE);
		}
		if(buf_probe_mesSA_superPointR_BME != -99999 && buf_probe_mesSA_superPointZ_BME != -99999){
			SPR.push_back(buf_probe_mesSA_superPointR_BME);
			SPZ.push_back(buf_probe_mesSA_superPointZ_BME);
		}
	}

	TCanvas *c1 = new TCanvas("c1","c1",1600,900);
	TGraph *gr = new TGraph(SPR.size(),&(SPR.at(0)),&(SPZ.at(0)));
	gr->Draw("P");
	gr->SetName("R_Zplot");
	TFile *output_file = new TFile(outputfilename.c_str(),"RECREATE");
	gr->Write();
	delete gr;
	delete c1;


}
