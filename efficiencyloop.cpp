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
#include <Efficiency.cpp>

const string trigger = "mu26ivm";
//const string trigger = "mu4";
//const string inputfilename = "~/dataset/efficiencysample2.root";
//const string inputfilename = "/gpfs/fs6001/kayamash/data18_physics_Main_Ztap_hadd.root";
//const string inputfilename = "/gpfs/fs6001/kayamash/efficiency_output/mc16_13TeVZmumu070.root";
const string inputfilename = "/gpfs/fs6001/kayamash/dataset/Zmumu300540_hadd.root";
const string outputfilename = "/gpfs/fs6001/kayamash/output/" + trigger + "/plottest.root";
const Int_t efficiency_maxenergy = 101;
const Double_t efficiency_x_err = 0.25;
const Int_t nhist = 13;
const Int_t thpitch = 2;

//main function
void efficiencyloop(){
	std::cout<<"start!"<<std::endl;
	Efficiency eff;
	TFile *output_file = new TFile(outputfilename.c_str(),"RECREATE");
	TFile *tf1 = TFile::Open(inputfilename.c_str(),"read");
	TTree *tr1 = dynamic_cast<TTree*>(tf1->Get("t_tap"));
	eff.Init(tr1,trigger,24,20,3.0,2.5,0.08,efficiency_maxenergy,efficiency_x_err,nhist,thpitch);
	for(Int_t event = 0;event < tr1->GetEntries(); event++){
		eff.Execute(event);
	}
	cout<<"finalize"<<endl;
	eff.Finalize(output_file);
}
