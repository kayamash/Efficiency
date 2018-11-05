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

	for(Int_t event = 0;event < chain->GetEntries();event++){
		chain->GetEntry(event);
		for(Int_t mes = 0;mes < probe_mesSAsuperPointR_BI;mes++){
			
		}
	}


}