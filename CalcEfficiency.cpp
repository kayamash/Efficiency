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

void CalcEfficiency::DrawHist1D(TH1D *h1){
	TCanvas *c1 = new TCanvas("c1","c1",1600,900);
	gStyle->SetTitleYOffset(m_yoffset);
	if(m_title != "test")h1->SetTitle(m_title.c_str());
	h1->SetFillColor(kRed);
	h1->SetLineColor(kRed);
	h1->Draw();
	h1->Write();
	delete c1;
}

void CalcEfficiency::DrawHist2D(TH2F *h1){
	TCanvas *c1 = new TCanvas("c1","c1",1200,1200);
	gStyle->SetPadTopMargin(m_topmargin);
	gStyle->SetPadBottomMargin(m_bottommargin);
	gStyle->SetPadLeftMargin(m_leftmargin);
	gStyle->SetPadRightMargin(m_rightmargin);
	h1->GetYaxis()->SetTitleOffset(m_yoffset);
	h1->SetStats(0);
	h1->Draw("colz");
	h1->Write();
	delete c1;
}

void CalcEfficiency::DrawHistAll(TH1D *h1,TH1D *h2,TH1D *h3,TH1D *h4,TH1D *h5){
	TCanvas *c1 = new TCanvas("c1","c1",1600,900);
	gStyle->SetTitleYOffset(m_yoffset);
	h1->SetStats(0);
	h2->SetStats(0);
	h3->SetStats(0);
	h4->SetStats(0);
	h5->SetStats(0);
	h1->SetFillColor(kRed);
	h1->SetLineColor(kRed);
	h2->SetFillColor(kGreen);
	h2->SetLineColor(kGreen);
	h3->SetFillColor(kBlue);
	h3->SetLineColor(kBlue);
	h4->SetFillColor(kYellow);
	h4->SetLineColor(kYellow);
	h5->SetFillColor(kMagenta);
	h5->SetLineColor(kMagenta);
	h1->Draw();
	h2->Draw("same");
	h3->Draw("same");
	h4->Draw("same");
	h5->Draw("same");
	c1->RedrawAxis();
	TLegend *legend = new TLegend(0.7,0.1,0.9,0.4,"");
	legend->AddEntry(h1,m_label[0].c_str(),"f");
	legend->AddEntry(h2,m_label[1].c_str(),"f");
	legend->AddEntry(h3,m_label[2].c_str(),"f");
	legend->AddEntry(h4,m_label[3].c_str(),"f");
	legend->AddEntry(h5,m_label[4].c_str(),"f");
	legend->Draw();
	delete c1;
}

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
  tg1->SetMarkerSize(0.3);
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

void CalcEfficiency::DrawEfficiencyAll(TH1D *h1,TH1D *h2,TH1D *h3,TH1D *h4,TH1D *h5,Double_t max,Int_t nbin,Double_t err){
  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  gStyle->SetTitleYOffset(m_yoffset);
  std::vector<Double_t> eff_x1;
  std::vector<Double_t> eff_y1;
  std::vector<Double_t> eff_x1_err;
  std::vector<Double_t> eff_y1_err;
  std::vector<Double_t> eff_x2;
  std::vector<Double_t> eff_y2;
  std::vector<Double_t> eff_x2_err;
  std::vector<Double_t> eff_y2_err;
  std::vector<Double_t> eff_x3;
  std::vector<Double_t> eff_y3;
  std::vector<Double_t> eff_x3_err;
  std::vector<Double_t> eff_y3_err;
  std::vector<Double_t> eff_x4;
  std::vector<Double_t> eff_y4;
  std::vector<Double_t> eff_x4_err;
  std::vector<Double_t> eff_y4_err;
  for(Int_t i = 0;i <= nbin;i++){
    Double_t buf_1 = static_cast<Double_t>(h1->GetBinContent(i+1));
    Double_t buf_2 = static_cast<Double_t>(h2->GetBinContent(i+1));
    Double_t buf_3 = static_cast<Double_t>(h3->GetBinContent(i+1));
    Double_t buf_4 = static_cast<Double_t>(h4->GetBinContent(i+1));
    Double_t buf_5 = static_cast<Double_t>(h5->GetBinContent(i+1));

    if(buf_1 != 0){
      eff_y1.push_back(buf_2/buf_1);
      eff_y1_err.push_back(sqrt(buf_2*(1 - buf_2/buf_1))/buf_1);
    }else{
      eff_y1.push_back(0);
      eff_y1_err.push_back(0);
    }
    if(buf_2 != 0){
      eff_y2.push_back(buf_3/buf_2);
      eff_y2_err.push_back(sqrt(buf_3*(1 - buf_3/buf_2))/buf_2);
    }else{
      eff_y2.push_back(0);
      eff_y2_err.push_back(0);
    }
    if(buf_3 != 0){
      eff_y3.push_back(buf_4/buf_3);
      eff_y3_err.push_back(sqrt(buf_4*(1 - buf_4/buf_3))/buf_3);
    }else{
      eff_y3.push_back(0);
      eff_y3_err.push_back(0);
    }
    if(buf_4 != 0){
      eff_y4.push_back(buf_5/buf_4);
      eff_y4_err.push_back(sqrt(buf_5*(1 - buf_5/buf_4))/buf_4);
   }else{
      eff_y4.push_back(0);
      eff_y4_err.push_back(0);
    }

    eff_x1.push_back(h1->GetBinCenter(i+1));
    eff_x1_err.push_back(err);
    eff_x2.push_back(h1->GetBinCenter(i+1));
    eff_x2_err.push_back(err);
    eff_x3.push_back(h1->GetBinCenter(i+1));
    eff_x3_err.push_back(err);
    eff_x4.push_back(h1->GetBinCenter(i+1));
    eff_x4_err.push_back(err);
  }

  TH1F *frame = gPad->DrawFrame(-2.0,0,max,1.1);
  frame->SetTitle(m_title.c_str());

  TGraphErrors *tg1 = new TGraphErrors(eff_x1.size(),&(eff_x1.at(0)),&(eff_y1.at(0)),&(eff_x1_err.at(0)),&(eff_y1_err.at(0)));
  tg1->SetMarkerStyle(20);
  tg1->SetMarkerSize(0.3);
  tg1->SetMarkerColor(1);
  tg1->SetLineColor(1);
  tg1->Draw("P");
  TGraphErrors *tg2 = new TGraphErrors(eff_x2.size(),&(eff_x2.at(0)),&(eff_y2.at(0)),&(eff_x2_err.at(0)),&(eff_y2_err.at(0)));
  tg2->SetMarkerStyle(20);
  tg2->SetMarkerSize(0.3);
  tg2->SetMarkerColor(2);
  tg2->SetLineColor(2);
  tg2->Draw("P");
  TGraphErrors *tg3 = new TGraphErrors(eff_x3.size(),&(eff_x3.at(0)),&(eff_y3.at(0)),&(eff_x3_err.at(0)),&(eff_y3_err.at(0)));
  tg3->SetMarkerStyle(20);
  tg3->SetMarkerSize(0.3);
  tg3->SetMarkerColor(3);
  tg3->SetLineColor(3);
  tg3->Draw("P");
  TGraphErrors *tg4 = new TGraphErrors(eff_x4.size(),&(eff_x4.at(0)),&(eff_y4.at(0)),&(eff_x4_err.at(0)),&(eff_y4_err.at(0)));
  tg4->SetMarkerStyle(20);
  tg4->SetMarkerSize(0.3);
  tg4->SetMarkerColor(6);
  tg4->SetLineColor(6);
  tg4->Draw("P");

  TLegend *legend = new TLegend(0.7,0.1,0.9,0.4,"");
  legend->AddEntry(tg1,m_label[0].c_str(),"lep");
  legend->AddEntry(tg2,m_label[1].c_str(),"lep");
  legend->AddEntry(tg3,m_label[2].c_str(),"lep");
  legend->AddEntry(tg4,m_label[3].c_str(),"lep");
  legend->Draw();

  eff_x1.clear();
  eff_x1_err.clear();
  eff_y1.clear();
  eff_y1_err.clear();
  eff_x2.clear();
  eff_x2_err.clear();
  eff_y2.clear();
  eff_y2_err.clear();
  eff_x3.clear();
  eff_x3_err.clear();
  eff_y3.clear();
  eff_y3_err.clear();
  eff_x4.clear();
  eff_x4_err.clear();
  eff_y4.clear();
  eff_y4_err.clear();
  delete tg1;
  delete tg2;
  delete tg3;
  delete tg4;
  delete c1; 
}

void CalcEfficiency::DrawEfficiencyeta(TH1D *h1,TH1D *h2,Double_t max){
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
    eff_x.push_back(static_cast<Double_t>(-2.5 + 0.1*i));
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

void CalcEfficiency::SetCondition(string title,string png,Double_t offset,Double_t tmargin,Double_t bmargin,Double_t lmargin,Double_t rmargin){
	m_title = title;
	m_png = png;
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

void CalcEfficiency::SetDirectry(string dir){
	m_dir = dir;
}

