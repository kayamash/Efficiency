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
//const string inputfilename = "/gpfs/fs6001/kayamash/data18_physics_Main_Ztap_hadd.root";
//const string inputfilename = "/gpfs/fs6001/kayamash/efficiency_output/mc16_13TeVZmumu070.root";
const string outputfilename = "/gpfs/fs6001/kayamash/dataset/Zmumu300540_hadd.root";
const string outputfilename = "/gpfs/fs6001/kayamash/output/" + trigger + "/plottest.root";
const Int_t efficiency_maxenergy = 101;
const Int_t efficiency_nbin = 300;
const Double_t efficiency_x_err = 0.25;

//main function
void efficiencyloop(){
	std::cout<<"start!"<<std::endl;
  //gROOT->LoadMacro("Efficiency.cpp");
	Efficiency eff;
  eff.SetNbin(efficiency_nbin);
	TFile *output_file = new TFile(outputfilename.c_str(),"RECREATE");
	TFile *tf1 = TFile::Open(inputfilename.c_str(),"read");
	TTree *tr1 = dynamic_cast<TTree*>(tf1->Get("t_tap"));
	for(Int_t i = 0;i < 82;i += 2){
		eff.Init(tr1,trigger,i,24,20,3.0,2.5,0.08,efficiency_maxenergy,efficiency_x_err);
		cout<<"start init"<<endl;
		for(Int_t event = 0;event < tr1->GetEntries(); event++){
			eff.Execute(event);
		}
		cout<<"start final"<<endl;
		eff.Final(output_file);
	}
}
