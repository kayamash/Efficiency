#include "Efficiency.chh"
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
#include <CalcEfficiency.cpp>

  void Efficiency::Init(TTree *tree,std::string name,const Int_t np,const Int_t ne,const Double_t mp,const Double_t me,Double_t req,Int_t max,Double_t err,const Int_t nh,const Int_t th){
   if (tree){
     tChain = tree;
     m_nbin_phi = np;
     m_nbin_eta = ne;
     m_phi_max = mp;
     m_eta_max = me;
     m_method_name = name;
     m_binmax = max;
     m_efficiency_xerr = err;
     m_nhist = nh;
     m_thpitch = th;
    
     //initialize
     m_toff_pt = 0;
     m_toff_eta = 0;
     m_toff_exteta = 0;
     m_toff_extinneta = 0;
     m_toff_phi = 0;
     m_toff_extphi = 0;
     m_toff_extinnphi = 0;
     m_tag_charge = 0;
     m_tag_d0 = 0;
     m_tag_z0 = 0;
     m_poff_pt = 0;
     m_poff_eta = 0;
     m_poff_exteta = 0;
     m_poff_extinneta = 0;
     m_poff_phi = 0;
     m_poff_extphi = 0;
     m_poff_extinnphi = 0;
     m_probe_charge = 0;
     m_probe_d0 = 0;
     m_probe_z0 = 0;
     m_tL1_pt = 0;
     m_tL1_eta = 0;
     m_tL1_phi = 0;
     m_tSA_pt = 0;
     m_tSA_eta = 0;
     m_tSA_phi = 0;
     m_tCB_pt = 0;
     m_tCB_eta = 0;
     m_tCB_phi = 0;
     m_tEF_pt = 0;
     m_tEF_eta = 0;
     m_tEF_phi = 0;
     m_mes_name = 0;
     m_pL1_pt = 0;
     m_pL1_eta = 0;
     m_pL1_phi = 0;
     m_pL1_dR = 0;
     m_pL1_pass = 0;
     m_pSA_pt = 0;
     m_pSA_eta = 0;
     m_pSA_phi = 0;
     m_pSA_dR = 0;
     m_pSA_pass = 0;
     m_pCB_pt = 0;
     m_pCB_eta = 0;
     m_pCB_phi = 0;
     m_pCB_dR = 0;
     m_pCB_pass = 0;
     m_pEF_pt = 0;
     m_pEF_eta = 0;
     m_pEF_phi = 0;
     m_pEF_dR = 0;
     m_pEF_pass = 0;
     m_reqL1dR = req;

    //active only need branch 
     tChain->SetBranchStatus("*",0);
     tChain->SetBranchStatus("mes_name",1);
     tChain->SetBranchStatus("tag_pt",1);
     tChain->SetBranchStatus("tag_eta",1);
     tChain->SetBranchStatus("tag_exteta",1);
     tChain->SetBranchStatus("tag_extinneta",1);
     tChain->SetBranchStatus("tag_phi",1);
     tChain->SetBranchStatus("tag_extphi",1);
     tChain->SetBranchStatus("tag_extinnphi",1);
     tChain->SetBranchStatus("tag_charge",1);
     tChain->SetBranchStatus("tag_d0",1);
     tChain->SetBranchStatus("tag_z0",1);
     tChain->SetBranchStatus("probe_pt",1);
     tChain->SetBranchStatus("probe_eta",1);
     tChain->SetBranchStatus("probe_exteta",1);
     tChain->SetBranchStatus("probe_extinneta",1);
     tChain->SetBranchStatus("probe_phi",1);
     tChain->SetBranchStatus("probe_extphi",1);
     tChain->SetBranchStatus("probe_extinnphi",1);
     tChain->SetBranchStatus("probe_charge",1);
     tChain->SetBranchStatus("probe_d0",1);
     tChain->SetBranchStatus("probe_z0",1);
     tChain->SetBranchStatus("tag_L1_pt",1);
     tChain->SetBranchStatus("tag_L1_eta",1);
     tChain->SetBranchStatus("tag_L1_phi",1);
     tChain->SetBranchStatus("tag_SA_pt",1);
     tChain->SetBranchStatus("tag_SA_eta",1);
     tChain->SetBranchStatus("tag_SA_phi",1);
     tChain->SetBranchStatus("tag_CB_pt",1);
     tChain->SetBranchStatus("tag_CB_eta",1);
     tChain->SetBranchStatus("tag_CB_phi",1);
     tChain->SetBranchStatus("tag_EF_pt",1);
     tChain->SetBranchStatus("tag_EF_eta",1);
     tChain->SetBranchStatus("tag_EF_phi",1);
     tChain->SetBranchStatus("probe_mesL1_pt",1);
     tChain->SetBranchStatus("probe_mesL1_eta",1);
     tChain->SetBranchStatus("probe_mesL1_phi",1);
     tChain->SetBranchStatus("probe_mesL1_pass",1);
     tChain->SetBranchStatus("probe_mesL1_dR",1);
     tChain->SetBranchStatus("probe_mesSA_pt",1);
     tChain->SetBranchStatus("probe_mesSA_eta",1);
     tChain->SetBranchStatus("probe_mesSA_phi",1);
     tChain->SetBranchStatus("probe_mesSA_pass",1);
     tChain->SetBranchStatus("probe_mesSA_dR",1);
     tChain->SetBranchStatus("probe_mesCB_pt",1);
     tChain->SetBranchStatus("probe_mesCB_eta",1);
     tChain->SetBranchStatus("probe_mesCB_phi",1);  
     tChain->SetBranchStatus("probe_mesCB_pass",1);
     tChain->SetBranchStatus("probe_mesCB_dR",1);
     tChain->SetBranchStatus("probe_mesEF_pt",1);
     tChain->SetBranchStatus("probe_mesEF_eta",1);
     tChain->SetBranchStatus("probe_mesEF_phi",1);
     tChain->SetBranchStatus("probe_mesEF_pass",1);
     tChain->SetBranchStatus("probe_mesEF_dR",1);
     //setting each branch address
     tChain->SetBranchAddress("mes_name",&m_mes_name,&b_mes_name);
     tChain->SetBranchAddress("tag_pt",&m_toff_pt,&b_tag_pt);
     tChain->SetBranchAddress("tag_eta",&m_toff_eta,&b_tag_eta);
     tChain->SetBranchAddress("tag_exteta",&m_toff_exteta,&b_tag_exteta);
     tChain->SetBranchAddress("tag_extinneta",&m_toff_extinneta,&b_tag_extinneta);
     tChain->SetBranchAddress("tag_phi",&m_toff_phi,&b_tag_phi);
     tChain->SetBranchAddress("tag_extphi",&m_toff_extphi,&b_tag_extphi);
     tChain->SetBranchAddress("tag_extinnphi",&m_toff_extinnphi,&b_tag_extinnphi);
     tChain->SetBranchAddress("tag_charge",&m_tag_charge,&b_tag_charge);
     tChain->SetBranchAddress("tag_d0",&m_tag_d0,&b_tag_d0);
     tChain->SetBranchAddress("tag_z0",&m_tag_z0,&b_tag_z0);
     tChain->SetBranchAddress("probe_pt",&m_poff_pt,&b_probe_pt);
     tChain->SetBranchAddress("probe_eta",&m_poff_eta,&b_probe_eta);
     tChain->SetBranchAddress("probe_exteta",&m_poff_exteta,&b_probe_exteta);
     tChain->SetBranchAddress("probe_extinneta",&m_poff_extinneta,&b_probe_extinneta);
     tChain->SetBranchAddress("probe_phi",&m_poff_phi,&b_probe_phi);
     tChain->SetBranchAddress("probe_extphi",&m_poff_extphi,&b_probe_extphi);
     tChain->SetBranchAddress("probe_extinnphi",&m_poff_extinnphi,&b_probe_extinnphi);
     tChain->SetBranchAddress("probe_charge",&m_probe_charge,&b_probe_charge);
     tChain->SetBranchAddress("probe_d0",&m_probe_d0,&b_probe_d0);
     tChain->SetBranchAddress("probe_z0",&m_probe_z0,&b_probe_z0);
     tChain->SetBranchAddress("tag_L1_pt",&m_tL1_pt,&b_tL1_pt);
     tChain->SetBranchAddress("tag_L1_eta",&m_tL1_eta,&b_tL1_eta);
     tChain->SetBranchAddress("tag_L1_phi",&m_tL1_phi,&b_tL1_phi);
     tChain->SetBranchAddress("tag_SA_pt",&m_tSA_pt,&b_tSA_pt);
     tChain->SetBranchAddress("tag_SA_eta",&m_tSA_eta,&b_tSA_eta);
     tChain->SetBranchAddress("tag_SA_phi",&m_tSA_phi,&b_tSA_phi);
     tChain->SetBranchAddress("tag_CB_pt",&m_tCB_pt,&b_tCB_pt);
     tChain->SetBranchAddress("tag_CB_eta",&m_tCB_eta,&b_tCB_eta);
     tChain->SetBranchAddress("tag_CB_phi",&m_tCB_phi,&b_tCB_phi);
     tChain->SetBranchAddress("tag_EF_pt",&m_tEF_pt,&b_tEF_pt);
     tChain->SetBranchAddress("tag_EF_eta",&m_tEF_eta,&b_tEF_eta);
     tChain->SetBranchAddress("tag_EF_phi",&m_tEF_phi,&b_tEF_phi);
     tChain->SetBranchAddress("probe_mesL1_pt",&m_pL1_pt,&b_pL1_pt);
     tChain->SetBranchAddress("probe_mesL1_eta",&m_pL1_eta,&b_pL1_eta);
     tChain->SetBranchAddress("probe_mesL1_phi",&m_pL1_phi,&b_pL1_phi);
     tChain->SetBranchAddress("probe_mesL1_pass",&m_pL1_pass,&b_pL1_pass);
     tChain->SetBranchAddress("probe_mesL1_dR",&m_pL1_dR,&b_pL1_dR);
     tChain->SetBranchAddress("probe_mesSA_pt",&m_pSA_pt,&b_pSA_pt);
     tChain->SetBranchAddress("probe_mesSA_eta",&m_pSA_eta,&b_pSA_eta);
     tChain->SetBranchAddress("probe_mesSA_phi",&m_pSA_phi,&b_pSA_phi);
     tChain->SetBranchAddress("probe_mesSA_pass",&m_pSA_pass,&b_pSA_pass);
     tChain->SetBranchAddress("probe_mesSA_dR",&m_pSA_dR,&b_pSA_dR);
     tChain->SetBranchAddress("probe_mesCB_pt",&m_pCB_pt,&b_pCB_pt);
     tChain->SetBranchAddress("probe_mesCB_eta",&m_pCB_eta,&b_pCB_eta);
     tChain->SetBranchAddress("probe_mesCB_phi",&m_pCB_phi,&b_pCB_phi);
     tChain->SetBranchAddress("probe_mesCB_pass",&m_pCB_pass,&b_pCB_pass);
     tChain->SetBranchAddress("probe_mesCB_dR",&m_pCB_dR,&b_pCB_dR);
     tChain->SetBranchAddress("probe_mesEF_pt",&m_pEF_pt,&b_pEF_pt);
     tChain->SetBranchAddress("probe_mesEF_eta",&m_pEF_eta,&b_pEF_eta);
     tChain->SetBranchAddress("probe_mesEF_phi",&m_pEF_phi,&b_pEF_phi);
     tChain->SetBranchAddress("probe_mesEF_pass",&m_pEF_pass,&b_pEF_pass);
     tChain->SetBranchAddress("probe_mesEF_dR",&m_pEF_dR,&b_pEF_dR);
     //define each histgram
     for(Int_t i = 0;i < m_nhist;i++){
          m_h_poff_pt.push_back(new TH1D(Form("h_poff_pt_%dGeV",i*thpitch),"probe offline pt;offline pt[GeV];Entries",150,0,150));
          m_h_pL1_pt.push_back(new TH1D(Form("h_pL1_pt_%dGeV",i*thpitch),"probe L1 pt;L1 pt[GeV];Entries",150,0,150));
          m_h_pSA_pt.push_back(new TH1D(Form("h_pSA_pt_%dGeV",i*thpitch),"probe L2MuonSA pt;L2MuonSA pt[GeV];Entries",150,0,150));
          m_h_pCB_pt.push_back(new TH1D(Form("h_pCB_pt_%dGeV",i*thpitch),"probe muComb pt;muComb pt[GeV];Entries",150,0,150));
          m_h_pEF_pt.push_back(new TH1D(Form("h_pEF_pt_%dGeV",i*thpitch),"probe EventFilter pt;EventFilter pt[GeV];Entries",150,0,150));
          m_h_pL1_dR.push_back(new TH1D(Form("h_L1 dR_%dGeV",i*thpitch),"L1 dR;dR;Entries",1000,0,0.1));
          m_h_pSA_dR.push_back(new TH1D(Form("h_SA dR_%dGeV",i*thpitch),"L2MuonSA dR;dR;Entries",500,0,0.05));
          m_h_pCB_dR.push_back(new TH1D(Form("h_CB dR_%dGeV",i*thpitch),"muComb dR;dR;Entries",200,0,0.002));
          m_h_pEF_dR.push_back(new TH1D(Form("h_EF dR_%dGeV",i*thpitch),"EventFilter dR;dR;Entries",100,0,0.001));
          m_h_textL1_dR.push_back(new TH1D(Form("h_textL1_dR_%dGeV",i*thpitch),"tag extL1 dR;dR;Entries",1000,0,0.1));
          m_h_textSA_dR.push_back(new TH1D(Form("h_textSA_dR_%dGeV",i*thpitch),"tag extL2MuonSA dR;dR;Entries",1000,0,0.1));
          m_h_textCB_dR.push_back(new TH1D(Form("h_textCB_dR_%dGeV",i*thpitch),"tag extmuComb dR;dR;Entries",1000,0,0.01));
          m_h_textEF_dR.push_back(new TH1D(Form("h_textEF_dR_%dGeV",i*thpitch),"tag extEventFilter dR;dR;Entries",1000,0,0.01));
          m_h_pextL1_dR.push_back(new TH1D(Form("h_pextL1_dR_%dGeV",i*thpitch),"probe extL1 dR;dR;Entries",1000,0,0.1));
          m_h_pextSA_dR.push_back(new TH1D(Form("h_pextSA_dR_%dGeV",i*thpitch),"probe extL2MuonSA dR;dR;Entries",1000,0,0.1));
          m_h_pextCB_dR.push_back(new TH1D(Form("h_pextCB_dR_%dGeV",i*thpitch),"probe extmuComb dR;dR;Entries",1000,0,0.01));
          m_h_pextEF_dR.push_back(new TH1D(Form("h_pextEF_dR_%dGeV",i*thpitch),"probe extEventFilter dR;dR;Entries",1000,0,0.01));
          m_h_pSA_respt.push_back(new TH1D(Form("h_pSA_respt_%dGeV",i*thpitch),"probe L2MuonSA residual pt;residual pt;Entries",300,-1,1));
          m_h_pCB_respt.push_back(new TH1D(Form("h_pCB_respt_%dGeV",i*thpitch),"probe muComb residual pt;residual pt;Entries",640,-0.8,0.8));
          m_h_pEF_respt.push_back(new TH1D(Form("h_pEF_respt_%dGeV",i*thpitch),"probe EventFilter residual pt;residual pt;Entries",1800,-0.3,0.3));
          m_h_eoff_pt.push_back(new TH1D(Form("h_eoff_pt_%dGeV",i*thpitch),"mesoff_pt;offline pt[GeV];Entries",300,-0.25,149.75));
          m_h_eL1_pt.push_back(new TH1D(Form("h_eL1_pt_%dGeV",i*thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
          m_h_eSA_pt.push_back(new TH1D(Form("h_eSA_pt_%dGeV",i*thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
          m_h_eCB_pt.push_back(new TH1D(Form("h_eCB_pt_%dGeV",i*thpitch),"mesCB_pt;CB pt[GeV];Entries",300,-0.25,149.75));
          m_h_eEF_pt.push_back(new TH1D(Form("h_eEF_pt_%dGeV",i*thpitch),"mesEF_pt;EF pt[GeV];Entries",300,-0.25,149.75));
          m_h_eoff_eta.push_back(new TH1D(Form("h_eoff_eta_%dGeV",i*thpitch),"off eta;off eta;Entries",50,-2.5,2.5));
          m_h_eL1_eta.push_back(new TH1D(Form("h_eL1_eta_%dGeV",i*thpitch),"L1 eta;L1 eta;Entries",50,-2.5,2.5));
          m_h_eSA_eta.push_back(new TH1D(Form("h_eSA_eta_%dGeV",i*thpitch),"SA eta;SA eta;Entries",50,-2.5,2.5));
          m_h_eCB_eta.push_back(new TH1D(Form("h_eCB_eta_%dGeV",i*thpitch),"CB eta;CB eta;Entries",50,-2.5,2.5));
          m_h_eEF_eta.push_back(new TH1D(Form("h_eEF_eta_%dGeV",i*thpitch),"EF eta;EF eta;Entries",50,-2.5,2.5));
          m_h_eoff_pt_barrel.push_back(new TH1D(Form("h_eoff_pt_barrel_%dGeV",i*thpitch),"off_pt;off pt[GeV];Entries",300,-0.25,149.75));
          m_h_eL1_pt_barrel.push_back(new TH1D(Form("h_eL1_pt_barrel_%dGeV",i*thpitch),"L1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
          m_h_eSA_pt_barrel.push_back(new TH1D(Form("h_eSA_pt_barrel_%dGeV",i*thpitch),"SA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
          m_h_eCB_pt_barrel.push_back(new TH1D(Form("h_eCB_pt_barrel_%dGeV",i*thpitch),"CB_pt;CB pt[GeV];Entries",300,-0.25,149.75));
          m_h_eEF_pt_barrel.push_back(new TH1D(Form("h_eEF_pt_barrel_%dGeV",i*thpitch),"EF_pt;EF pt[GeV];Entries",300,-0.25,149.75));
          m_h_eoff_pt_end.push_back(new TH1D(Form("h_eoff_pt_end_%dGeV",i*thpitch),"off_pt;off pt[GeV];Entries",300,-0.25,149.75));
          m_h_eL1_pt_end.push_back(new TH1D(Form("h_eL1_pt_end_%dGeV",i*thpitch),"L1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
          m_h_eSA_pt_end.push_back(new TH1D(Form("h_eSA_pt_end_%dGeV",i*thpitch),"SA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
          m_h_eCB_pt_end.push_back(new TH1D(Form("h_eCB_pt_end_%dGeV",i*thpitch),"CB_pt;CB pt[GeV];Entries",300,-0.25,149.75));
          m_h_eEF_pt_end.push_back(new TH1D(Form("h_eEF_pt_end_%dGeV",i*thpitch),"EF_pt;EF pt[GeV];Entries",300,-0.25,149.75));
          m_h_eff_poff_etaphi.push_back(new TH2F(Form("h_eff_poff_etaphi_%dGeV",i*thpitch),"offlineeta vs offlinephi;offline eta;offline phi",m_nbin_eta,-m_eta_max,m_eta_max,m_nbin_phi,-m_phi_max,m_phi_max));
          m_h_eff_pL1_etaphi.push_back(new TH2F(Form("h_eff_pL1_etaphi_%dGeV",i*thpitch),"L1eta vs L1phi;L1 eta;L1 phi",m_nbin_eta,-m_eta_max,m_eta_max,m_nbin_phi,-m_phi_max,m_phi_max));
          m_h_eff_pSA_etaphi.push_back(new TH2F(Form("h_eff_pSA_etaphi_%dGeV",i*thpitch),"L2MuonSAeta vs L2MuonSAphi;L2MuonSA eta;L2MuonSA phi",m_nbin_eta,-m_eta_max,m_eta_max,m_nbin_phi,-m_phi_max,m_phi_max));
          m_h_poffvsSA_pt.push_back(new TH2F(Form("h_poffvsSA_pt_%dGeV",i*thpitch),"probe offline pt vs probe L2MuonSA pt@mu26ivm;probe offline pt[GeV];probe L2MuonSA pt[GeV]",150,0,150,150,0,150));
     }
   }
 }

 bool Efficiency::Dicision_barrel(Double_t eta){
     if(std::fabs(eta) <= 1.05){
          return kTRUE;
     }else{
          return kFALSE;
     }
}

bool Efficiency::Cut_L1(Int_t pass){
  if(pass > -1){
          return kTRUE;
     }else{
          return kFALSE;
     }
}

bool Efficiency::Cut_SA(Int_t pass,Double_t pt,Double_t th){
     if(th == 0){
          if(pass == 1){
               return kTRUE;
          }else{
               return kFALSE;
          }
     }else{
          if(pass == 1 && std::fabs(pt) > th){
               return kTRUE;
          }else{
               return kFALSE;
          }
     }
}

bool Efficiency::Cut_CB(Int_t pass){
     if(pass == 1){
          return kTRUE;
     }else{
          return kFALSE;
     }
}

bool Efficiency::Cut_EF(Int_t pass){
     if(pass == 1){
          return kTRUE;
     }else{
          return kFALSE;
     }
}

void Efficiency::Execute(Int_t ev){
     tChain->GetEntry(ev);
     for(Int_t i = 0;m_nhist;i++){
          Double_t pextL1_dR = 1; 
          Double_t pextSA_dR = 1; 
          Double_t pextCB_dR = 1; 
          Double_t pextEF_dR = 1; 
          Double_t pL1_pt = -99999;
          Double_t pL1_eta = 0;
          Double_t pL1_phi = 0;
          Double_t pL1_dR = 1;
          Int_t pL1_pass = 0;
          Double_t pSA_pt = -99999;
          Double_t pSA_eta = 0;
          Double_t pSA_phi = 0;
          Double_t pSA_dR = 1;
          Int_t pSA_pass = 0;
          Double_t pCB_pt = -99999;
          Double_t pCB_eta = 0;
          Double_t pCB_phi = 0;
          Double_t pCB_dR = 1;
          Int_t pCB_pass = 0;
          Double_t pEF_pt = -99999;
          Double_t pEF_eta = 0;
          Double_t pEF_phi = 0;
          Double_t pEF_dR = 1;
          Int_t pEF_pass = 0;
          for(Int_t method = 0;method < 25;method++){
               if(m_mes_name->at(method) == m_method_name){
                    pL1_pt = m_pL1_pt->at(method);
                    pSA_pt = m_pSA_pt->at(method);
                    pCB_pt = m_pCB_pt->at(method);
                    pEF_pt = m_pEF_pt->at(method);
                    pL1_eta = m_pL1_eta->at(method);
                    pSA_eta = m_pSA_eta->at(method);
                    pCB_eta = m_pCB_eta->at(method);
                    pEF_eta = m_pEF_eta->at(method);
                    pL1_phi = m_pL1_phi->at(method);
                    pSA_phi = m_pSA_phi->at(method);
                    pCB_phi = m_pCB_phi->at(method);
                    pEF_phi = m_pEF_phi->at(method);
                    pL1_pass = m_pL1_pass->at(method);
                    pSA_pass = m_pSA_pass->at(method);
                    pCB_pass = m_pCB_pass->at(method);
                    pEF_pass = m_pEF_pass->at(method);
                    pL1_dR = m_pL1_dR->at(method);
                    pSA_dR = m_pSA_dR->at(method);
                    pCB_dR = m_pCB_dR->at(method);
                    pEF_dR = m_pEF_dR->at(method);
               }
          }
          pL1_dR = TMath::Sqrt(pow(pL1_eta - m_poff_eta,2) + pow(pL1_phi - m_poff_phi,2) );
          if(std::fabs(m_poff_pt)*0.001 < 10.0)m_reqL1dR = -0.00001*std::fabs(m_poff_pt) + 0.18;
          //offline
          h_poff_pt->at(i)->Fill(m_poff_pt*0.001);
          h_eoff_pt->at(i)->Fill(std::fabs(m_poff_pt*0.001));
          if(std::fabs(m_poff_pt*0.001) > 40)h_eoff_eta->Fill(m_poff_eta);
          if(Dicision_barrel(m_poff_eta)){
               h_eoff_pt_barrel->at(i)->Fill(std::fabs(m_poff_pt*0.001));
          }else{
               h_eoff_pt_end->at(i)->Fill(std::fabs(m_poff_pt*0.001));
          }
          if(std::fabs(m_poff_pt*0.001) > 40)h_eff_poff_etaphi->at(i)->Fill(m_poff_eta,m_poff_phi);

          //L1
          if(Cut_L1(pL1_pass)){
               Double_t textL1_dR = TMath::Sqrt(pow(m_tL1_eta - m_toff_exteta,2) + pow(m_tL1_phi - m_toff_extphi,2));
               pextL1_dR = TMath::Sqrt(pow(pL1_eta - m_poff_exteta,2) + pow(pL1_phi - m_poff_extphi,2));

               h_pL1_pt->at(i)->Fill(std::fabs(pL1_pt*0.001));
               h_pL1_dR->at(i)->Fill(pL1_dR);
               h_textL1_dR->at(i)->Fill(textL1_dR);
               h_pextL1_dR->at(i)->Fill(pextL1_dR);
               h_eL1_pt->at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(Dicision_barrel(m_poff_eta)){
                    h_eL1_pt_barrel->at(i)->Fill(std::fabs(m_poff_pt*0.001));
               }else{
                    h_eL1_pt_end->at(i)->Fill(std::fabs(m_poff_pt*0.001));
               }
               if(std::fabs(m_poff_pt*0.001) > 40){
                    h_eff_pL1_etaphi->at(i)->Fill(m_poff_eta,m_poff_phi);
                    h_eL1_eta->at(i)->Fill(m_poff_eta);
               }

               //SA
               if(Cut_SA(pSA_pass,pSA_pt,i*thpitch)){
                    Double_t textSA_dR = TMath::Sqrt(pow(m_tSA_eta - m_toff_exteta,2) + pow(m_tSA_phi - m_toff_extphi,2));
                    pextSA_dR = TMath::Sqrt(pow(pSA_eta - m_poff_exteta,2) + pow(pSA_phi - m_poff_extphi,2));
                    Double_t resSA_pt = std::fabs(m_poff_pt*0.001)/std::fabs(pSA_pt) - 1.0;
                    Double_t buf_pSA_dR = TMath::Sqrt(pow(pSA_eta - m_poff_eta,2) + pow(pSA_phi - m_poff_phi,2));

                    h_pSA_pt->at(i)->Fill(std::fabs(pSA_pt));
                    h_pSA_dR->at(i)->Fill(buf_pSA_dR);
                    h_textSA_dR->at(i)->Fill(textSA_dR);
                    h_pextSA_dR->at(i)->Fill(pextSA_dR);
                    h_eSA_pt->at(i)->Fill(std::fabs(m_poff_pt*0.001));
                    h_pSA_respt->at(i)->Fill(resSA_pt);
                    if(Dicision_barrel(m_poff_eta)){
                         h_eSA_pt_barrel->at(i)->Fill(std::fabs(m_poff_pt*0.001));
                    }else{
                         h_eSA_pt_end->at(i)->Fill(std::fabs(m_poff_pt*0.001));
                    }
                    if(std::fabs(m_poff_pt*0.001) > 40){
                         h_eff_pSA_etaphi->at(i)->Fill(m_poff_eta,m_poff_phi);
                         h_eSA_eta->at(i)->Fill(m_poff_eta);
                    }
                    h_poffvsSA_pt->at(i)->Fill(std::fabs(m_poff_pt*0.001),std::fabs(pSA_pt));

                    //CB
                    if(Cut_CB(pCB_pass)){
                         Double_t textCB_dR = TMath::Sqrt(pow(m_tCB_eta - m_toff_exteta,2) + pow(m_tCB_phi - m_toff_extphi,2));
                         pextCB_dR = TMath::Sqrt(pow(pCB_eta - m_poff_exteta,2) + pow(pCB_phi - m_poff_extphi,2));
                         Double_t resCB_pt = std::fabs(m_poff_pt)/std::fabs(pCB_pt) - 1.0;
                         h_pCB_pt->at(i)->Fill(std::fabs(pCB_pt*0.001));
                         h_pCB_dR->at(i)->Fill(pCB_dR);
                         h_textCB_dR->at(i)->Fill(textCB_dR);
                         h_pextCB_dR->at(i)->Fill(pextCB_dR);
                         h_eCB_pt->at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         if(std::fabs(m_poff_pt*0.001) > 40)h_eCB_eta->Fill(m_poff_eta);
                         h_pCB_respt->at(i)->Fill(resCB_pt);
                         if(Dicision_barrel(m_poff_eta)){
                              h_eCB_pt_barrel->at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         }else{
                              h_eCB_pt_end->at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         }

                         //EF
                         if(Cut_EF(pEF_pass)){
                              Double_t textEF_dR = TMath::Sqrt(pow(m_tEF_eta - m_toff_exteta,2) + pow(m_tEF_phi - m_toff_extphi,2));
                              pextEF_dR = TMath::Sqrt(pow(pEF_eta - m_poff_exteta,2) + pow(pEF_phi - m_poff_extphi,2));
                              Double_t resEF_pt = std::fabs(m_poff_pt)/std::fabs(pEF_pt) - 1.0;
                              h_pEF_pt->at(i)->Fill(std::fabs(pEF_pt*0.001));
                              h_pEF_dR->at(i)->Fill(pEF_dR);
                              h_textEF_dR->at(i)->Fill(textEF_dR);
                              h_pextEF_dR->at(i)->Fill(pextEF_dR);
                              h_eEF_pt->at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              if(std::fabs(m_poff_pt*0.001) > 40)h_eEF_eta->Fill(m_poff_eta);
                              h_pEF_respt->at(i)->Fill(resEF_pt);
                              if(Dicision_barrel(m_poff_eta)){
                                   h_eEF_pt_barrel->at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              }else{
                                   h_eEF_pt_end->at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              }
                         }
                    }
               }
          }
     }
}

void Efficiency::Finalize(TFile *tf1){
     CalcEfficiency ceff;
     tf1->cd();
     //SetCondition
     //title,file title,yoffset,top margin,bottom margin,left margin,right margin
     for(Int_t i = 0;i < m_nhist;i++){
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_poff_pt->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_pL1_pt->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_pSA_pt->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_pCB_pt->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_pEF_pt->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_pL1_dR->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_pSA_dR->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_pCB_dR->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_pEF_dR->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_textL1_dR->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_textSA_dR->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_textCB_dR->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_textEF_dR->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_pextL1_dR->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_pextSA_dR->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_pextCB_dR->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_pextEF_dR->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_pSA_respt->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_pCB_respt->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_pEF_respt->at(i));
          ceff.SetCondition("trigger;offline eta;count",1.5,0,0,0,0);
          ceff.DrawHist1D(h_eL1_eta->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_eSA_eta->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_eCB_eta->at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(h_eEF_eta->at(i));

          ceff.SetCondition("test",1.5,0.1,0.1,0.105,0.165);
          ceff.DrawHist2D(h_eff_poff_etaphi->at(i));
          ceff.SetCondition("test",1.5,0.1,0.1,0.105,0.165);
          ceff.DrawHist2D(h_eff_pL1_etaphi->at(i));
          ceff.SetCondition("test",1.5,0.1,0.1,0.105,0.165);
          ceff.DrawHist2D(h_eff_pSA_etaphi->at(i));
          ceff.SetCondition("test",1.5,0.1,0.1,0.105,0.165);
          ceff.DrawHist2D(h_poffvsSA_pt->at(i));

          //base,target
          ceff.SetConditionName(Form("L1Efficiency_%dGeV",i*thpitch));
          ceff.SetCondition("L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(h_eoff_pt->at(i),h_eL1_pt->at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiency_%dGeV",i*thpitch));
          ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(h_eL1_pt->at(i),h_eSA_pt->at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("CBEfficiency_%dGeV",i*thpitch));
          ceff.SetCondition("muComb Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(h_eSA_pt->at(i),h_eCB_pt->at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("EFEfficiency_%dGeV",i*thpitch));
          ceff.SetCondition("EventFilter Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(h_eCB_pt->at(i),h_eEF_pt->at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1Efficiency_eta_%dGeV",i*thpitch));
          ceff.SetCondition("L1 Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyeta(h_eoff_eta->at(i),h_eL1_eta->at(i));
          ceff.SetConditionName(Form("SAEfficiency_eta_%dGeV",i*thpitch));
          ceff.SetCondition("L2MuonSA Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyeta(h_eL1_eta->at(i),h_eSA_eta->at(i));
          ceff.SetConditionName(Form("CBEfficiency_eta_%dGeV",i*thpitch));
          ceff.SetCondition("muComb Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyeta(h_eSA_eta->at(i),h_eCB_eta->at(i));
          ceff.SetConditionName(Form("EFEfficiency_eta_%dGeV",i*thpitch));
          ceff.SetCondition("EventFilter Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyeta(h_eCB_eta->at(i),h_eEF_eta->at(i));
          ceff.SetConditionName(Form("L1Efficiency_barrel_%dGeV",i*thpitch));
          ceff.SetCondition("L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(h_eoff_pt_barrel->at(i),h_eL1_pt_barrel->at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiency_barrel_%dGeV",i*thpitch));
          ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(h_eL1_pt_barrel->at(i),h_eSA_pt_barrel->at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("CBEfficiency_barrel_%dGeV",i*thpitch));
          ceff.SetCondition("muComb Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(h_eSA_pt_barrel->at(i),h_eCB_pt_barrel->at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("EFEfficiency_barrel_%dGeV",i*thpitch));
          ceff.SetCondition("EventFilter Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(h_eCB_pt_barrel->at(i),h_eEF_pt_barrel->at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1Efficiency_end_%dGeV",i*thpitch));
          ceff.SetCondition("L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(h_eoff_pt_end->at(i),h_eL1_pt_end->at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiency_end_%dGeV",i*thpitch));
          ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(h_eL1_pt_end->at(i),h_eSA_pt_end->at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("CBEfficiency_end_%dGeV",i*thpitch));
          ceff.SetCondition("muComb Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(h_eSA_pt_end->at(i),h_eCB_pt_end->at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("EFEfficiency_end_%dGeV",i*thpitch));
          ceff.SetCondition("EventFilter Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(h_eCB_pt_end->at(i),h_eEF_pt_end->at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SA2DEfficiency_%dGeV",i*thpitch));
          ceff.SetCondition("L1vsL2MuonSA Efficiency;offline eta;offline phi",1.5,0.1,0.1,0.105,0.165);
          ceff.SetConditionbin(m_nbin_eta,m_nbin_phi,m_eta_max,m_phi_max);
          ceff.DrawEfficiency2D(h_eff_pL1_etaphi,h_eff_pSA_etaphi);
          ceff.SetConditionName(Form("L12DEfficiency_%dGeV",i*thpitch));
          ceff.SetCondition("L1 Efficiency;offline eta;offline phi",1.5,0.1,0.1,0.105,0.165);
          ceff.SetConditionbin(m_nbin_eta,m_nbin_phi,m_eta_max,m_phi_max);
          ceff.DrawEfficiency2D(h_eff_poff_etaphi->at(i),h_eff_pL1_etaphi->at(i));
     }

     h_poff_pt->clear();
     h_pL1_pt->clear();
     h_pSA_pt->clear();
     h_pCB_pt->clear();
     h_pEF_pt->clear();
     h_pL1_dR->clear();
     h_pSA_dR->clear();
     h_pCB_dR->clear();
     h_pEF_dR->clear();
     h_textL1_dR->clear();
     h_textSA_dR->clear();
     h_textCB_dR->clear();
     h_textEF_dR->clear();
     h_pextL1_dR->clear();
     h_pextSA_dR->clear();
     h_pextCB_dR->clear();
     h_pextEF_dR->clear();
     h_pSA_respt->clear();
     h_pCB_respt->clear();
     h_pEF_respt->clear();
     h_eoff_pt->clear();
     h_eL1_pt->clear();
     h_eSA_pt->clear();
     h_eCB_pt->clear();
     h_eEF_pt->clear();
     h_eoff_eta->clear();
     h_eL1_eta->clear();
     h_eSA_eta->clear();
     h_eCB_eta->clear();
     h_eEF_eta->clear();
     h_eoff_pt_barrel->clear();
     h_eL1_pt_barrel->clear();
     h_eSA_pt_barrel->clear();
     h_eCB_pt_barrel->clear();
     h_eEF_pt_barrel->clear();
     h_eoff_pt_end->clear();
     h_eL1_pt_end->clear();
     h_eSA_pt_end->clear();
     h_eCB_pt_end->clear();
     h_eEF_pt_end->clear();
     h_eff_poff_etaphi->clear();
     h_eff_pL1_etaphi->clear();
     h_eff_pSA_etaphi->clear();
     h_poffvsSA_pt->clear();
     m_mes_name->clear();
     m_pL1_pt->clear();
     m_pL1_eta->clear();
     m_pL1_phi->clear();
     m_pL1_pass->clear();
     m_pL1_dR->clear();
     m_pSA_pt->clear();
     m_pSA_eta->clear();
     m_pSA_phi->clear();
     m_pSA_pass->clear();
     m_pSA_dR->clear();
     m_pCB_pt->clear();
     m_pCB_eta->clear();
     m_pCB_phi->clear();
     m_pCB_pass->clear();
     m_pCB_dR->clear();
     m_pEF_pt->clear();
     m_pEF_eta->clear();
     m_pEF_phi->clear();
     m_pEF_pass->clear();
     m_pEF_dR->clear();
}
