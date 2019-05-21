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
#include "../src/Efficiency.cpp"

//const string trigger = "mu26ivm";
//const string trigger = "data18mu26ivm";
const string trigger = "mu4";
//Jpsitap == 1,Ztap == 3
Int_t proc = 1;

//const string inputfilelist = "/home/kayamash/efflist/Zmumu364160.list";
const string inputfilelist = "/home/kayamash/efflist/data18_physics_Main_Ztap.list";
//const string inputfilelist = "/home/kayamash/efflist/Jpsi_noMdtCsm1k.list";
//const string inputfilelist = "/home/kayamash/efflist/newmc16345099.list"; 
//const string outputfilename = "/gpfs/fs6001/kayamash/Mywork/efficiencyloopoutput/Jpsi_noMdtCsm1k.root";
const string outputfilename = "/gpfs/fs6001/kayamash/Mywork/efficiencyloopoutput/20190510/newdata18MU4JPZtest.root";
//const string outputfilename = "/gpfs/fs6001/kayamash/Mywork/efficiencyloopoutput/Zmumu364160.root";
const Int_t efficiency_maxenergy = 61;
const Double_t efficiency_x_err = 0.25;
const bool EventFullScan = kTRUE;//true = full scan,false = sample scan(1000000 event)

//main function
void efficiencyloop(){
	cout<<"start!"<<endl;
	TChain *tr1 = new TChain("t_tap");
	
	std::ifstream ifs(inputfilelist.c_str());
	std::string str;
	while(getline(ifs,str)){
		tr1->Add(str.c_str());
	}
	if(!tr1)cout<<"tree failed"<<endl;

	Efficiency *eff = new Efficiency(tr1);
	TFile *output_file = new TFile(outputfilename.c_str(),"RECREATE");
	cout<<"Initialize"<<endl;
	eff->Init(trigger,48,80,3.0,2.5,0.08,efficiency_maxenergy,efficiency_x_err,proc);
	cout<<tr1->GetEntries()<<endl;
	cout<<"Execute"<<endl;
	Int_t nEvent = (EventFullScan)?(tr1->GetEntries()):(1000000);
	for(Int_t event = 0;event < nEvent; ++event){
		eff->Execute(event);
		if(event % 1000000 == 0)cout<<"The event is "<<event<<endl;
	}
	cout<<"Finalize"<<endl;
	eff->Finalize(output_file);

	delete output_file;
	delete tr1;
}
