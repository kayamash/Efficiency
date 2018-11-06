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
#include "Efficiency.cpp"

const string trigger = "mu26ivm";
//const string trigger = "mu4";
//const string inputfilename = "~/dataset/efficiencysample2.root";
//const string inputfilename = "/gpfs/fs6001/kayamash/data18_physics_Main_Ztap_hadd.root";
//const string inputfilename = "/gpfs/fs6001/kayamash/efficiency_output/mc16_13TeVZmumu070.root";
//const string inputfilename = "/gpfs/fs6001/kayamash/dataset/Zmumu300540_hadd.root";
const string inputfilelist = "/home/kayamash/list/data16_grid.list";
const string outputfilename = "/gpfs/fs6001/kayamash/Mywork/efficiencyloopoutput/plotadd4.root";
const Int_t efficiency_maxenergy = 101;
const Double_t efficiency_x_err = 0.25;
const Int_t nhist = 8;
const Int_t thpitch = 4;

//main function
void efficiencyloop(){
	Efficiency eff;
	
	std::ofstream ofs("LargeSpecialEvent.dat");
	ofs.close();

	TFile *output_file = new TFile(outputfilename.c_str(),"RECREATE");
    	/*
	TFile *tf1 = TFile::Open(inputfilename.c_str(),"read");
	TTree *tr1 = dynamic_cast<TTree*>(tf1->Get("t_tap"));
	*/
        
	TChain *tr1 = new TChain("t_tap");
	std::ifstream ifs(inputfilelist.c_str());
	std::string str;
	while(getline(ifs,str)){
		tr1->Add(str.c_str());
	}
	
	if(!tr1)cout<<"read failed!"<<endl; 
	cout<<"Initialize"<<endl;
	eff.Init(tr1,trigger,48,80,3.0,2.5,0.08,efficiency_maxenergy,efficiency_x_err,nhist,thpitch);
	cout<<tr1->GetEntries()<<endl;
	cout<<"Execute"<<endl;
	for(Int_t event = 0;event < tr1->GetEntries(); event++){
	//for(Int_t event = 0;event < 100000; event++){
		eff.Execute(event);
	}
	cout<<"Finalize"<<endl;
	eff.Finalize(output_file);

	delete output_file;
	delete tr1;
}
