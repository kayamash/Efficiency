#include "CalcEfficiency.chh"
#include <TObject.h>
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
#include <TROOT.h>
#include <TChain.h>

void CalcEfficiency::DrawEfficiency(TH1D *h1,TH1D *h2,Double_t max,Int_t nbin,Double_t err){
  TCanvas *c1 = new TCanvas("c1","c1",1600,900);
  gStyle->SetTitleYOffset(m_yoffset);
  TH1F *frame = gPad->DrawFrame(-2.0,0,max,1.1);
  std::vector<Double_t> eff_x;
  std::vector<Double_t> eff_y;
  std::vector<Double_t> eff_x_err;
  std::vector<Double_t> eff_y_err;
  for(Int_t i = 0;i <= nbin;i++){
    Double_t buf_m = static_cast<Double_t>(h1->GetBinContent(i+1));
    Double_t buf_n = static_cast<Double_t>(h2->GetBinContent(i+1));
    if(buf_m != 0){
    eff_y.push_back(buf_n/buf_m);
    eff_y_err.push_back(sqrt(buf_n*(1 - buf_n/buf_m))/buf_m);
    eff_x.push_back(h1->GetBinCenter(i+1));
   eff_x_err.push_back(err);
   }
  }

  TGraphErrors *tg1 = new TGraphErrors(eff_x.size(),&(eff_x.at(0)),&(eff_y.at(0)),&(eff_x_err.at(0)),&(eff_y_err.at(0)));
  tg1->SetMarkerStyle(20);
  tg1->SetMarkerSize(0.4);
  frame->SetTitle(m_title.c_str());
  tg1->Draw("P");
  tg1->SetName(m_name.c_str());
  tg1->Write();
  eff_x.clear();
  eff_x_err.clear();
  eff_y.clear();
  eff_y_err.clear();
  delete tg1;
  delete c1; 
}

void CalcEfficiency::DrawEfficiencyeta(TH1D *h1,TH1D *h2){
  TCanvas *c1 = new TCanvas("c1","c1",1600,900);
  gStyle->SetTitleYOffset(m_yoffset);
  TH1F *frame = gPad->DrawFrame(-2.5,0,2.5,1);
  std::vector<Double_t> eff_x;
  std::vector<Double_t> eff_y;
  std::vector<Double_t> eff_x_err;
  std::vector<Double_t> eff_y_err;
  for(Int_t i = 0;i <= 50;i++){
    Double_t buf_m = static_cast<Double_t>(h1->GetBinContent(i+1));
    Double_t buf_n = static_cast<Double_t>(h2->GetBinContent(i+1));
    if(buf_m != 0){
    eff_y.push_back(buf_n/buf_m);
    eff_y_err.push_back(sqrt(buf_n*(1 - buf_n/buf_m))/buf_m);
    eff_x.push_back(static_cast<Double_t>(h1->GetBinCenter(i+1)));
    eff_x_err.push_back(0.05);
   }
  }

  TGraphErrors *tg1 = new TGraphErrors(eff_x.size(),&(eff_x.at(0)),&(eff_y.at(0)),&(eff_x_err.at(0)),&(eff_y_err.at(0)));
  tg1->SetMarkerStyle(20);
  tg1->SetMarkerSize(1.0);
  frame->SetTitle(m_title.c_str());
  tg1->Draw("P");
  tg1->SetName(m_name.c_str());
  tg1->Write();
  eff_x.clear();
  eff_x_err.clear();
  eff_y.clear();
  eff_y_err.clear();
  delete tg1;
  delete c1; 
}

void CalcEfficiency::DrawEfficiencyphi(TH1D *h1,TH1D *h2){
  TCanvas *c1 = new TCanvas("c1","c1",1600,900);
  gStyle->SetTitleYOffset(m_yoffset);
  TH1F *frame = gPad->DrawFrame(-2.5,0,2.5,1);
  std::vector<Double_t> eff_x;
  std::vector<Double_t> eff_y;
  std::vector<Double_t> eff_x_err;
  std::vector<Double_t> eff_y_err;
  for(Int_t i = 0;i <= 48;i++){
    Double_t buf_m = static_cast<Double_t>(h1->GetBinContent(i+1));
    Double_t buf_n = static_cast<Double_t>(h2->GetBinContent(i+1));
    if(buf_m != 0){
    eff_y.push_back(buf_n/buf_m);
    eff_y_err.push_back(sqrt(buf_n*(1 - buf_n/buf_m))/buf_m);
    eff_x.push_back(static_cast<Double_t>(h1->GetBinCenter(i+1)));
    eff_x_err.push_back(1.0/16.0);
   }
  }

  TGraphErrors *tg1 = new TGraphErrors(eff_x.size(),&(eff_x.at(0)),&(eff_y.at(0)),&(eff_x_err.at(0)),&(eff_y_err.at(0)));
  tg1->SetMarkerStyle(20);
  tg1->SetMarkerSize(1.0);
  frame->SetTitle(m_title.c_str());
  tg1->Draw("P");
  tg1->SetName(m_name.c_str());
  tg1->Write();
  eff_x.clear();
  eff_x_err.clear();
  eff_y.clear();
  eff_y_err.clear();
  delete tg1;
  delete c1; 
}

void CalcEfficiency::DrawEfficiencypileup(TH1D *h1,TH1D *h2){
  TCanvas *c1 = new TCanvas("c1","c1",1600,900);
  gStyle->SetTitleYOffset(m_yoffset);
  TH1F *frame = gPad->DrawFrame(-2.5,0,2.5,1);
  std::vector<Double_t> eff_x;
  std::vector<Double_t> eff_y;
  std::vector<Double_t> eff_x_err;
  std::vector<Double_t> eff_y_err;
  for(Int_t i = 0;i <= 50;i++){
    Double_t buf_m = static_cast<Double_t>(h1->GetBinContent(i+1));
    Double_t buf_n = static_cast<Double_t>(h2->GetBinContent(i+1));
    if(buf_m != 0){
    eff_y.push_back(buf_n/buf_m);
    eff_y_err.push_back(sqrt(buf_n*(1 - buf_n/buf_m))/buf_m);
    eff_x.push_back(static_cast<Double_t>(h1->GetBinCenter(i+1)));
    eff_x_err.push_back(0.5);
   }
  }

  TGraphErrors *tg1 = new TGraphErrors(eff_x.size(),&(eff_x.at(0)),&(eff_y.at(0)),&(eff_x_err.at(0)),&(eff_y_err.at(0)));
  tg1->SetMarkerStyle(20);
  tg1->SetMarkerSize(1.0);
  frame->SetTitle(m_title.c_str());
  tg1->Draw("P");
  tg1->SetName(m_name.c_str());
  tg1->Write();
  eff_x.clear();
  eff_x_err.clear();
  eff_y.clear();
  eff_y_err.clear();
  delete tg1;
  delete c1; 
}

void CalcEfficiency::DrawEfficiency2D(TH2F *h1,TH2F *h2){
	TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
	gStyle->SetTitleOffset(m_yoffset);
	TH2F *h3 = new TH2F(m_name.c_str(),m_title.c_str(),m_nbineta,-1*m_etamax,m_etamax,m_nbinphi,-1*m_phimax,m_phimax);
  gStyle->SetPalette(1);
	for(Int_t i = 0;i < m_nbineta;i++){
		for(Int_t j = 0;j < m_nbinphi;j++){
			Double_t tmp_eff1 = h1->GetBinContent(i+1,j+1);
			Double_t tmp_eff2 = h2->GetBinContent(i+1,j+1);
      if(tmp_eff1 == 0){
      h3->SetBinContent(i+1,j+1,0);
      }else{
			h3->SetBinContent(i+1,j+1,tmp_eff2/tmp_eff1);
      }
		}
	}
	h3->GetYaxis()->SetTitleOffset(m_yoffset);
	h3->SetStats(0);
	h3->Draw("colz");
	h3->Write();
	delete h3;
	delete c1;
}

void CalcEfficiency::DrawResidualplot(TH1D *h1,TH1D *h2,TH1D *h3,TH1D *h4,TH1D *h5,TH1D *h6,TH1D *h7,TH1D *h8,TH1D *h9,string type,Int_t num,Int_t pitch,Int_t thmin){
  TCanvas *c1 = new TCanvas("c1","c1",1600,900);
  gStyle->SetTitleOffset(m_yoffset);
  std::vector<Double_t> res_x;
  std::vector<Double_t> res_y;
  std::vector<Double_t> res_x_err;
  std::vector<Double_t> res_y_err;
  TH1D hist[9] = {*h1,*h2,*h3,*h4,*h5,*h6,*h7,*h8,*h9};
  string title = type + Form("_%d",num*pitch + thmin);
  for(Int_t i = 0;i < 9;i++){
    res_x.push_back((i + 1)*5.0);
    res_x_err.push_back(0);
    hist[i].Draw();
    TF1 *f1 = new TF1("f1","gaus",-1.0,1.0);
    hist[i].Fit("f1","","",-0.1,0.1);
    res_y.push_back(f1->GetParameter(1));
    res_y_err.push_back(f1->GetParameter(2));
    c1->SaveAs(("~/plot/" + title + Form("_%d.png",i)).c_str());
    c1->Clear();
  }
  TGraphErrors *tg1 = new TGraphErrors(res_x.size(),&(res_x.at(0)),&(res_y.at(0)),&(res_x_err.at(0)),&(res_y_err.at(0)));
  tg1->SetMarkerStyle(20);
  tg1->SetMarkerSize(1.0);
  tg1->Draw("ap");
  tg1->SetName(title.c_str());
  tg1->Write();
  res_x.clear();
  res_x_err.clear();
  res_y.clear();
  res_y_err.clear();
  delete tg1;
  delete c1;
}

void CalcEfficiency::SetCondition(string title,Double_t offset,Double_t tmargin,Double_t bmargin,Double_t lmargin,Double_t rmargin){
	m_title = title;
	m_yoffset = offset;
	m_topmargin = tmargin;
	m_bottommargin = bmargin;
	m_leftmargin = lmargin;
	m_rightmargin = rmargin;
}

void CalcEfficiency::SetConditionlabel(string label1,string label2,string label3,string label4,string label5){
	m_label[0] = label1;
	m_label[1] = label2;
	m_label[2] = label3;
	m_label[3] = label4;
	m_label[4] = label5;
}

void CalcEfficiency::SetConditionbin(Double_t etabin,Double_t phibin,Double_t etamax,Double_t phimax){
	m_nbineta = etabin;
	m_nbinphi = phibin;
	m_etamax = etamax;
	m_phimax = phimax;
}

void CalcEfficiency::SetConditionName(string name){
	m_name = name;
}

