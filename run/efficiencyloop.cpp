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
//const string trigger = "mu6";
//Jpsitap == 1,Ztap == 3
Int_t proc = 1;

const string inputfilelist = "/home/kayamash/efflist/20190416data18_physics_Main_Ztap.list";
//const string inputfilelist = "/home/kayamash/efflist/mc16Jpsi_tsakai.list";
//const string inputfilelist = "/home/kayamash/efflist/mc16_410472.list";
const string outputfilename = "/gpfs/fs6001/kayamash/Mywork/efficiencyloopoutput/20190510/data18_physics_Main_ZtapMU4JPZ.root";
//const string outputfilename = "/gpfs/fs6001/kayamash/Mywork/efficiencyloopoutput/20190318/mc16Jpsi_tsakaiMU20.root";
//const string outputfilename = "/gpfs/fs6001/kayamash/Mywork/efficiencyloopoutput/20190328/mc16_410472.root";
const Int_t efficiency_maxenergy = 61;
const Double_t efficiency_x_err = 0.25;
const Int_t eventmode = 0;//eventmode = 0,full scan eventmode = 1,sample scan

//main function
void efficiencyloop(){
	cout<<"start!"<<endl;
	TChain *tr1 = new TChain("t_tap");
	
	std::ifstream ifs(inputfilelist.c_str());
	std::string str;
	while(getline(ifs,str)){
		tr1->Add(str.c_str());
	}
	
	//tr1->Add(inputfilename.c_str());

	if(!tr1)cout<<"tree failed"<<endl;
	Efficiency *eff = new Efficiency(tr1);
	
	TFile *output_file = new TFile(outputfilename.c_str(),"RECREATE");

	cout<<"Initialize"<<endl;
	eff->Init(trigger,48,80,3.0,2.5,0.08,efficiency_maxenergy,efficiency_x_err,proc);
	cout<<tr1->GetEntries()<<endl;
	cout<<"Execute"<<endl;
	Int_t nevent = 0;
	if(eventmode == 0){
		nevent = tr1->GetEntries();
	}else if(eventmode == 1){
		nevent = 100000;
	}
	for(Int_t event = 0;event < nevent; event++){
		//if(event%100000 == 0)printf("%d\r",event);
		if(event%1000000 == 0)cout<<"The number of entry is "<<event<<" now!"<<endl;
		eff->Execute(event);
	}
	cout<<"Finalize"<<endl;
	eff->Finalize(output_file);

	delete output_file;
	delete tr1;
}
