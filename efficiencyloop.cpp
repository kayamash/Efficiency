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
//const string trigger = "data18mu26ivm";
//const string trigger = "mu4";
//Jpsitap == 1,Ztap == 3
Int_t proc = 3;

//const string inputfilelist = "/home/kayamash/efflist/data18_physics_Main_Ztap.list";
//const string inputfilelist = "/home/kayamash/efflist/mc16Jpsi_tsakai.list";
const string inputfilelist = "/home/kayamash/efflist/data15284484.list";
//const string outputfilename = "/gpfs/fs6001/kayamash/Mywork/efficiencyloopoutput/data18Ztap_kayamashMU20.root";
//const string outputfilename = "/gpfs/fs6001/kayamash/Mywork/efficiencyloopoutput/mc16Jpsi_tsakaiMU20.root";
const string outputfilename = "/gpfs/fs6001/kayamash/Mywork/efficiencyloopoutput/mdata15_284484.root";
const Int_t efficiency_maxenergy = 101;
const Double_t efficiency_x_err = 0.25;
const Int_t thmin = 0;
const Int_t nhist = 1;
const Int_t thpitch = 4;
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
	Efficiency *eff = new Efficiency(nhist,thpitch,thmin,tr1);
	
	TFile *output_file = new TFile(outputfilename.c_str(),"RECREATE");

	cout<<"Initialize"<<endl;
	eff->Init(trigger,48,80,3.0,2.5,0.08,efficiency_maxenergy,efficiency_x_err,nhist,proc);
	cout<<tr1->GetEntries()<<endl;
	cout<<"Execute"<<endl;
	Int_t nevent = 0;
	if(eventmode == 0){
		nevent = tr1->GetEntries();
	}else if(eventmode == 1){
		nevent = 100000;
	}
	for(Int_t event = 0;event < nevent; event++){
		if(event%100000 == 0)printf("%d\r",event);
		eff->Execute(event);
	}
	cout<<"Finalize"<<endl;
	eff->Finalize(output_file);

	delete output_file;
	delete tr1;
}
