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

 void Efficiency::SetNbin(Int_t nbin){
     m_efficiency_nbin = nbin;
 }

 void Efficiency::Init(TTree *tree,std::string name,const Int_t threshold_SA,const Int_t np,const Int_t ne,const Double_t mp,const Double_t me,Double_t req,Int_t max,Double_t err){
   if (tree){
     tChain = tree;
     m_nbin_phi = np;
     m_nbin_eta = ne;
     m_phi_max = mp;
     m_eta_max = me;
     m_method_name = name;
     m_threshold_SA = threshold_SA;
     m_binmax = max;
     m_efficiency_xerr = err;
    
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
     m_countSA = 0;
     m_countEF = 0;
     m_countSA_barrel = 0;
     m_countEF_barrel = 0;
     m_countSA_endcap = 0;
     m_countEF_endcap = 0;
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
     h_poff_pt = new TH1D(Form("h_poff_pt_%dGeV",m_threshold_SA),"probe offline pt;offline pt[GeV];Entries",150,0,150);
     h_pL1_pt = new TH1D(Form("h_pL1_pt_%dGeV",m_threshold_SA),"probe L1 pt;L1 pt[GeV];Entries",150,0,150);
     h_pSA_pt = new TH1D(Form("h_pSA_pt_%dGeV",m_threshold_SA),"probe L2MuonSA pt;L2MuonSA pt[GeV];Entries",150,0,150);
     h_pCB_pt = new TH1D(Form("h_pCB_pt_%dGeV",m_threshold_SA),"probe muComb pt;muComb pt[GeV];Entries",150,0,150);
     h_pEF_pt = new TH1D(Form("h_pEF_pt_%dGeV",m_threshold_SA),"probe EventFilter pt;EventFilter pt[GeV];Entries",150,0,150);
     h_pL1_dR = new TH1D(Form("h_L1 dR_%dGeV",m_threshold_SA),"L1 dR;dR;Entries",1000,0,0.1);
     h_pSA_dR = new TH1D(Form("h_SA dR_%dGeV",m_threshold_SA),"L2MuonSA dR;dR;Entries",500,0,0.05);
     h_pCB_dR = new TH1D(Form("h_CB dR_%dGeV",m_threshold_SA),"muComb dR;dR;Entries",200,0,0.002);
     h_pEF_dR = new TH1D(Form("h_EF dR_%dGeV",m_threshold_SA),"EventFilter dR;dR;Entries",100,0,0.001);
     h_textL1_dR = new TH1D(Form("h_textL1_dR_%dGeV",m_threshold_SA),"tag extL1 dR;dR;Entries",1000,0,0.1);
     h_textSA_dR = new TH1D(Form("h_textSA_dR_%dGeV",m_threshold_SA),"tag extL2MuonSA dR;dR;Entries",1000,0,0.1);
     h_textCB_dR = new TH1D(Form("h_textCB_dR_%dGeV",m_threshold_SA),"tag extmuComb dR;dR;Entries",1000,0,0.01);
     h_textEF_dR = new TH1D(Form("h_textEF_dR_%dGeV",m_threshold_SA),"tag extEventFilter dR;dR;Entries",1000,0,0.01);
     h_pextL1_dR = new TH1D(Form("h_pextL1_dR_%dGeV",m_threshold_SA),"probe extL1 dR;dR;Entries",1000,0,0.1);
     h_pextSA_dR = new TH1D(Form("h_pextSA_dR_%dGeV",m_threshold_SA),"probe extL2MuonSA dR;dR;Entries",1000,0,0.1);
     h_pextCB_dR = new TH1D(Form("h_pextCB_dR_%dGeV",m_threshold_SA),"probe extmuComb dR;dR;Entries",1000,0,0.01);
     h_pextEF_dR = new TH1D(Form("h_pextEF_dR_%dGeV",m_threshold_SA),"probe extEventFilter dR;dR;Entries",1000,0,0.01);
     h_td0 = new TH1D(Form("h_td0_%dGeV",m_threshold_SA),"tag d0;tag d0;Entries",1000,0,0.1);
     h_tz0 = new TH1D(Form("h_tz0_%dGeV",m_threshold_SA),"tag z0;tag z0;Entries",10000,0,1000);
     h_pd0 = new TH1D(Form("h_pd0_%dGeV",m_threshold_SA),"probe d0;probe d0;Entries",500,0,0.05);
     h_pz0 = new TH1D(Form("h_pz0_%dGeV",m_threshold_SA),"probe z0;probe z0;Entries",10000,0,1000);
     h_pSA_respt = new TH1D(Form("h_pSA_respt_%dGeV",m_threshold_SA),"probe L2MuonSA residual pt;residual pt;Entries",300,-1,1);
     h_pCB_respt = new TH1D(Form("h_pCB_respt_%dGeV",m_threshold_SA),"probe muComb residual pt;residual pt;Entries",640,-0.8,0.8);
     h_pEF_respt = new TH1D(Form("h_pEF_respt_%dGeV",m_threshold_SA),"probe EventFilter residual pt;residual pt;Entries",1800,-0.3,0.3);
     h_invmass_off = new TH1D(Form("h_invmass_off_%dGeV",m_threshold_SA),"offline invariant mass;invariant mass[GeV];Entries",800,0,200);
     h_invmass_L1 = new TH1D(Form("h_invmass_L1_%dGeV",m_threshold_SA),"L1 invariant mass;invariant mass[GeV];Entries",800,0,200);
     h_invmass_SA = new TH1D(Form("h_invmass_SA_%dGeV",m_threshold_SA),"L2MuonSA invariant mass;invariant mass[GeV];Entries",800,0,200);
     h_invmass_CB = new TH1D(Form("h_invmass_CB_%dGeV",m_threshold_SA),"muComb invariant mass;invariant mass[GeV];Entries",800,0,200);
     h_invmass_EF = new TH1D(Form("h_invmass_EF_%dGeV",m_threshold_SA),"EventFilter invariant mass;invariant mass[GeV];Entries",800,0,200);
     h_eoff_pt = new TH1D(Form("h_eoff_pt_%dGeV",m_threshold_SA),"mesoff_pt;offline pt[GeV];Entries",m_efficiency_nbin,-0.25,150-0.25);
     h_eL1_pt = new TH1D(Form("h_eL1_pt_%dGeV",m_threshold_SA),"mesL1_pt;L1 pt[GeV];Entries",m_efficiency_nbin,-0.25,150-0.25);
     h_eSA_pt = new TH1D(Form("h_eSA_pt_%dGeV",m_threshold_SA),"mesSA_pt;SA pt[GeV];Entries",m_efficiency_nbin,-0.25,150-0.25);
     h_eCB_pt = new TH1D(Form("h_eCB_pt_%dGeV",m_threshold_SA),"mesCB_pt;CB pt[GeV];Entries",m_efficiency_nbin,-0.25,150-0.25);
     h_eEF_pt = new TH1D(Form("h_eEF_pt_%dGeV",m_threshold_SA),"mesEF_pt;EF pt[GeV];Entries",m_efficiency_nbin,-0.25,150-0.25);
     h_eoff_eta = new TH1D(Form("h_eoff_eta_%dGeV",m_threshold_SA),"off eta;off eta;Entries",50,-2.5,2.5);
     h_eL1_eta = new TH1D(Form("h_eL1_eta_%dGeV",m_threshold_SA),"L1 eta;L1 eta;Entries",50,-2.5,2.5);
     h_eSA_eta = new TH1D(Form("h_eSA_eta_%dGeV",m_threshold_SA),"SA eta;SA eta;Entries",50,-2.5,2.5);
     h_eCB_eta = new TH1D(Form("h_eCB_eta_%dGeV",m_threshold_SA),"CB eta;CB eta;Entries",50,-2.5,2.5);
     h_eEF_eta = new TH1D(Form("h_eEF_eta_%dGeV",m_threshold_SA),"EF eta;EF eta;Entries",50,-2.5,2.5);
     h_eoff_pt_barrel = new TH1D(Form("h_eoff_pt_barrel_%dGeV",m_threshold_SA),"off_pt;off pt[GeV];Entries",m_efficiency_nbin,-0.25,150-0.25);
     h_eL1_pt_barrel = new TH1D(Form("h_eL1_pt_barrel_%dGeV",m_threshold_SA),"L1_pt;L1 pt[GeV];Entries",m_efficiency_nbin,-0.25,150-0.25);
     h_eSA_pt_barrel = new TH1D(Form("h_eSA_pt_barrel_%dGeV",m_threshold_SA),"SA_pt;SA pt[GeV];Entries",m_efficiency_nbin,-0.25,150-0.25);
     h_eCB_pt_barrel = new TH1D(Form("h_eCB_pt_barrel_%dGeV",m_threshold_SA),"CB_pt;CB pt[GeV];Entries",m_efficiency_nbin,-0.25,150-0.25);
     h_eEF_pt_barrel = new TH1D(Form("h_eEF_pt_barrel_%dGeV",m_threshold_SA),"EF_pt;EF pt[GeV];Entries",m_efficiency_nbin,-0.25,150-0.25);
     h_eoff_pt_end = new TH1D(Form("h_eoff_pt_end_%dGeV",m_threshold_SA),"off_pt;off pt[GeV];Entries",m_efficiency_nbin,-0.25,150-0.25);
     h_eL1_pt_end = new TH1D(Form("h_eL1_pt_end_%dGeV",m_threshold_SA),"L1_pt;L1 pt[GeV];Entries",m_efficiency_nbin,-0.25,150-0.25);
     h_eSA_pt_end = new TH1D(Form("h_eSA_pt_end_%dGeV",m_threshold_SA),"SA_pt;SA pt[GeV];Entries",m_efficiency_nbin,-0.25,150-0.25);
     h_eCB_pt_end = new TH1D(Form("h_eCB_pt_end_%dGeV",m_threshold_SA),"CB_pt;CB pt[GeV];Entries",m_efficiency_nbin,-0.25,150-0.25);
     h_eEF_pt_end = new TH1D(Form("h_eEF_pt_end_%dGeV",m_threshold_SA),"EF_pt;EF pt[GeV];Entries",m_efficiency_nbin,-0.25,150-0.25);
     h_eff_poff_etaphi = new TH2F(Form("h_eff_poff_etaphi_%dGeV",m_threshold_SA),"offlineeta vs offlinephi;offline eta;offline phi",m_nbin_eta,-m_eta_max,m_eta_max,m_nbin_phi,-m_phi_max,m_phi_max);
     h_eff_pL1_etaphi = new TH2F(Form("h_eff_pL1_etaphi_%dGeV",m_threshold_SA),"L1eta vs L1phi;L1 eta;L1 phi",m_nbin_eta,-m_eta_max,m_eta_max,m_nbin_phi,-m_phi_max,m_phi_max);
     h_eff_pSA_etaphi = new TH2F(Form("h_eff_pSA_etaphi_%dGeV",m_threshold_SA),"L2MuonSAeta vs L2MuonSAphi;L2MuonSA eta;L2MuonSA phi",m_nbin_eta,-m_eta_max,m_eta_max,m_nbin_phi,-m_phi_max,m_phi_max);
     h_poffvsSA_pt = new TH2F(Form("h_poffvsSA_pt_%dGeV",m_threshold_SA),"probe offline pt vs probe L2MuonSA pt@mu26ivm;probe offline pt[GeV];probe L2MuonSA pt[GeV]",150,0,150,150,0,150);
   }
 }

 bool Efficiency::Dicision_barrel(Double_t eta){
     if(std::fabs(eta) <= 1.05){
          return kTRUE;
     }else{
          return kFALSE;
     }
}

bool Efficiency::Cut_L1(Int_t pass,Double_t dr){
  if(pass > -1){
          return kTRUE;
     }else{
          return kFALSE;
     }
}

bool Efficiency::Cut_SA(Int_t pass,Double_t pt,Double_t th){
     if(th == 0){
          if(pass == 1){
               m_countSA++;
               return kTRUE;
          }else{
               return kFALSE;
          }
     }else{
          if(pass == 1 && std::fabs(pt) > th){
               m_countSA++;
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
          m_countEF++;
          return kTRUE;
     }else{
          return kFALSE;
     }
}

void Efficiency::Execute(Int_t ev){
     tChain->GetEntry(ev);
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
       h_poff_pt->Fill(m_poff_pt*0.001);
       h_invmass_off->Fill(0.001*TMath::Sqrt(2.0*m_toff_pt*m_poff_pt*(TMath::CosH(m_toff_eta - m_poff_eta) - TMath::Cos(m_toff_phi - m_poff_phi))));
       h_td0->Fill(m_tag_d0);
       h_pd0->Fill(m_probe_d0);
       h_tz0->Fill(m_tag_z0);
       h_pz0->Fill(m_probe_z0);
       h_eoff_pt->Fill(std::fabs(m_poff_pt*0.001));
        if(std::fabs(m_poff_pt*0.001) > 40)h_eoff_eta->Fill(m_poff_eta);
       if(Dicision_barrel(m_poff_eta)){
               h_eoff_pt_barrel->Fill(std::fabs(m_poff_pt*0.001));
       }else{
            h_eoff_pt_end->Fill(std::fabs(m_poff_pt*0.001));
       }
    if(std::fabs(m_poff_pt*0.001) > 40)h_eff_poff_etaphi->Fill(m_poff_eta,m_poff_phi);

     //L1
       if(Cut_L1(pL1_pass,pL1_dR)){
     Double_t textL1_dR = TMath::Sqrt(pow(m_tL1_eta - m_toff_exteta,2) + pow(m_tL1_phi - m_toff_extphi,2));
       pextL1_dR = TMath::Sqrt(pow(pL1_eta - m_poff_exteta,2) + pow(pL1_phi - m_poff_extphi,2));
       Double_t invmass_L1 = 0.001*TMath::Sqrt(2.0*m_tL1_pt*pL1_pt*(TMath::CosH(m_tL1_eta - pL1_eta) - TMath::Cos(m_tL1_phi - pL1_phi)));

       h_pL1_pt->Fill(std::fabs(pL1_pt*0.001));
       h_pL1_dR->Fill(pL1_dR);
       h_textL1_dR->Fill(textL1_dR);
       h_pextL1_dR->Fill(pextL1_dR);
       h_invmass_L1->Fill(invmass_L1);
       h_eL1_pt->Fill(std::fabs(m_poff_pt*0.001));
       if(Dicision_barrel(m_poff_eta)){
               h_eL1_pt_barrel->Fill(std::fabs(m_poff_pt*0.001));
       }else{
            h_eL1_pt_end->Fill(std::fabs(m_poff_pt*0.001));
       }
       if(std::fabs(m_poff_pt*0.001) > 40){h_eff_pL1_etaphi->Fill(m_poff_eta,m_poff_phi);
       h_eL1_eta->Fill(m_poff_eta);
     }

       //SA
       if(Cut_SA(pSA_pass,pSA_pt,m_threshold_SA)){
       Double_t textSA_dR = TMath::Sqrt(pow(m_tSA_eta - m_toff_exteta,2) + pow(m_tSA_phi - m_toff_extphi,2));
       pextSA_dR = TMath::Sqrt(pow(pSA_eta - m_poff_exteta,2) + pow(pSA_phi - m_poff_extphi,2));
       Double_t invmass_SA = TMath::Sqrt(2.0*m_tSA_pt*pSA_pt*(TMath::CosH(m_tSA_eta - pSA_eta) - TMath::Cos(m_tSA_phi - pSA_phi)));
       Double_t resSA_pt = std::fabs(m_poff_pt*0.001)/std::fabs(pSA_pt) - 1.0;
    Double_t buf_pSA_dR = TMath::Sqrt(pow(pSA_eta - m_poff_eta,2) + pow(pSA_phi - m_poff_phi,2));

       h_pSA_pt->Fill(std::fabs(pSA_pt));
       h_pSA_dR->Fill(buf_pSA_dR);
       h_textSA_dR->Fill(textSA_dR);
       h_pextSA_dR->Fill(pextSA_dR);
       h_invmass_SA->Fill(invmass_SA);
       h_eSA_pt->Fill(std::fabs(m_poff_pt*0.001));
       h_pSA_respt->Fill(resSA_pt);
       if(Dicision_barrel(m_poff_eta)){
           h_eSA_pt_barrel->Fill(std::fabs(m_poff_pt*0.001));
           m_countSA_barrel++;
       }else{
           h_eSA_pt_end->Fill(std::fabs(m_poff_pt*0.001));
           m_countSA_endcap++;
       }
       if(std::fabs(m_poff_pt*0.001) > 40){h_eff_pSA_etaphi->Fill(m_poff_eta,m_poff_phi);
       h_eSA_eta->Fill(m_poff_eta);
}
       h_poffvsSA_pt->Fill(std::fabs(m_poff_pt*0.001),std::fabs(pSA_pt));

       //CB
       if(Cut_CB(pCB_pass)){
       Double_t textCB_dR = TMath::Sqrt(pow(m_tCB_eta - m_toff_exteta,2) + pow(m_tCB_phi - m_toff_extphi,2));
       pextCB_dR = TMath::Sqrt(pow(pCB_eta - m_poff_exteta,2) + pow(pCB_phi - m_poff_extphi,2));
       Double_t invmass_CB = TMath::Sqrt(2.0*m_tSA_pt*pSA_pt*(TMath::CosH(m_tCB_eta - pCB_eta) - TMath::Cos(m_tCB_phi - pCB_phi)));
       Double_t resCB_pt = std::fabs(m_poff_pt)/std::fabs(pCB_pt) - 1.0;

       h_pCB_pt->Fill(std::fabs(pCB_pt*0.001));
       h_pCB_dR->Fill(pCB_dR);
       h_textCB_dR->Fill(textCB_dR);
       h_pextCB_dR->Fill(pextCB_dR);
       h_invmass_CB->Fill(invmass_CB);
       h_eCB_pt->Fill(std::fabs(m_poff_pt*0.001));
       if(std::fabs(m_poff_pt*0.001) > 40)h_eCB_eta->Fill(m_poff_eta);
       h_pCB_respt->Fill(resCB_pt);
       if(Dicision_barrel(m_poff_eta)){
               h_eCB_pt_barrel->Fill(std::fabs(m_poff_pt*0.001));
       }else{
           h_eCB_pt_end->Fill(std::fabs(m_poff_pt*0.001));
       }

       //EF
       if(Cut_EF(pEF_pass)){
       Double_t textEF_dR = TMath::Sqrt(pow(m_tEF_eta - m_toff_exteta,2) + pow(m_tEF_phi - m_toff_extphi,2));
       pextEF_dR = TMath::Sqrt(pow(pEF_eta - m_poff_exteta,2) + pow(pEF_phi - m_poff_extphi,2));
       Double_t invmass_EF = TMath::Sqrt(2.0*m_tSA_pt*pSA_pt*(TMath::CosH(m_tEF_eta - pEF_eta) - TMath::Cos(m_tEF_phi - pEF_phi)));
       Double_t resEF_pt = std::fabs(m_poff_pt)/std::fabs(pEF_pt) - 1.0;

       h_pEF_pt->Fill(std::fabs(pEF_pt*0.001));
       h_pEF_dR->Fill(pEF_dR);
       h_textEF_dR->Fill(textEF_dR);
       h_pextEF_dR->Fill(pextEF_dR);
       h_invmass_EF->Fill(invmass_EF);
       h_eEF_pt->Fill(std::fabs(m_poff_pt*0.001));
       if(std::fabs(m_poff_pt*0.001) > 40)h_eEF_eta->Fill(m_poff_eta);
       h_pEF_respt->Fill(resEF_pt);
       if(Dicision_barrel(m_poff_eta)){
           h_eEF_pt_barrel->Fill(std::fabs(m_poff_pt*0.001));
           m_countEF_barrel++;
       }else{
           h_eEF_pt_end->Fill(std::fabs(m_poff_pt*0.001));
           m_countEF_endcap++;
       }
    }
   }
  }
 }
}

void Efficiency::Final(TFile *tf1){
     gROOT->LoadMacro("CalcEfficiency.cpp");
     CalcEfficiency ceff;
     tf1->cd();
     //SetCondition
     //title,file title,yoffset,top margin,bottom margin,left margin,right margin

     ceff.SetDirectry("/gpfs/fs6001/kayamash/output/" + m_method_name + Form("/%dGeV/",m_threshold_SA) );
     ceff.GetDirectry();
     ceff.SetCondition("test","probe_offline_pt.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_poff_pt);
     ceff.SetCondition("test","probeL1_pt.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_pL1_pt);
     ceff.SetCondition("test","probeSA_pt.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_pSA_pt);
     ceff.SetCondition("test","probeCB_pt.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_pCB_pt);
     ceff.SetCondition("test","probeEF_pt.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_pEF_pt);
     ceff.SetCondition("test","probeL1_dR.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_pL1_dR);
     ceff.SetCondition("test","probeSA_dR.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_pSA_dR);
     ceff.SetCondition("test","probeCB_dR.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_pCB_dR);
     ceff.SetCondition("test","probeEF_dR.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_pEF_dR);
     ceff.SetCondition("test","textL1_dR.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_textL1_dR);
     ceff.SetCondition("test","textSA_dR.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_textSA_dR);
     ceff.SetCondition("test","textCB_dR.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_textCB_dR);
     ceff.SetCondition("test","textEF_dR.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_textEF_dR);
     ceff.SetCondition("test","pextL1_dR.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_pextL1_dR);
     ceff.SetCondition("test","pextSA_dR.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_pextSA_dR);
     ceff.SetCondition("test","pextCB_dR.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_pextCB_dR);
     ceff.SetCondition("test","pextEF_dR.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_pextEF_dR);
     ceff.SetCondition("test","tag_d0.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_td0);
     ceff.SetCondition("test","tag_z0.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_tz0);
     ceff.SetCondition("test","probe_d0.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_pd0);
     ceff.SetCondition("test","probe_d0.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_pz0);
     ceff.SetCondition("test","probeSA_respt.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_pSA_respt);
     ceff.SetCondition("test","probeCB_respt.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_pCB_respt);
     ceff.SetCondition("test","probeEF_respt.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_pEF_respt);
     ceff.SetCondition("test","invmass_offline.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_invmass_off);
     ceff.SetCondition("test","invmass_L1.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_invmass_L1);
     ceff.SetCondition("test","invmass_SA.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_invmass_SA);
     ceff.SetCondition("test","invmass_CB.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_invmass_CB);
     ceff.SetCondition("test","invmass_EF.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_invmass_EF);
     ceff.SetCondition("trigger;offline eta;count","eL1_eta.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_eL1_eta);
     ceff.SetCondition("test","eSA_eta.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_eSA_eta);
     ceff.SetCondition("test","eCB_eta.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_eCB_eta);
     ceff.SetCondition("test","eEF_eta.png",1.5,0,0,0,0);
     ceff.DrawHist1D(h_eEF_eta);


     ceff.SetCondition("test","offline_eta_phi.png",1.5,0.1,0.1,0.105,0.165);
     ceff.DrawHist2D(h_eff_poff_etaphi);
     ceff.SetCondition("test","L1_eta_phi.png",1.5,0.1,0.1,0.105,0.165);
     ceff.DrawHist2D(h_eff_pL1_etaphi);
     ceff.SetCondition("test","SA_eta_phi.png",1.5,0.1,0.1,0.105,0.165);
     ceff.DrawHist2D(h_eff_pSA_etaphi);
     ceff.SetCondition("test","probe_mesSA_pt.png",1.5,0.1,0.1,0.105,0.165);
     ceff.DrawHist2D(h_poffvsSA_pt);

     ceff.SetConditionlabel("offline","L1","L2MuonSA","muComb","EventFilter");
     ceff.SetCondition("z invariant mass;Mmumu[GeV];Entries","z_invariant_mass_All.png",1.5,0.1,0.1,0.105,0.165);
     ceff.DrawHistAll(h_invmass_off,h_invmass_L1,h_invmass_SA,h_invmass_CB,h_invmass_EF);
     ceff.SetConditionlabel("offline pt","L1 pt","L2MuonSA pt","muComb pt","EventFilter pt");
     ceff.SetCondition("pt distribution;pt[GeV];Entries","pt_distribution_All.png",1.5,0.1,0.1,0.105,0.165);
     ceff.DrawHistAll(h_poff_pt,h_pL1_pt,h_pSA_pt,h_pCB_pt,h_pEF_pt);

     //base,target
     ceff.SetConditionName(Form("L1Efficiency_%dGeV",m_threshold_SA));
     ceff.SetCondition("L1 Efficiency;offline pt[GeV];Efficiency","L1efficiency_pt.png",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(h_eoff_pt,h_eL1_pt,m_binmax,m_efficiency_nbin);
     ceff.SetConditionName(Form("SAEfficiency_%dGeV",m_threshold_SA));
     ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency","SAefficiency_pt.png",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(h_eL1_pt,h_eSA_pt,m_binmax,m_efficiency_nbin);
     ceff.SetConditionName(Form("CBEfficiency_%dGeV",m_threshold_SA));
     ceff.SetCondition("muComb Efficiency;offline pt[GeV];Efficiency","CBefficiency_pt.png",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(h_eSA_pt,h_eCB_pt,m_binmax,m_efficiency_nbin);
     ceff.SetConditionName(Form("EFEfficiency_%dGeV",m_threshold_SA));
     ceff.SetCondition("EventFilter Efficiency;offline pt[GeV];Efficiency","EFefficiency_pt.png",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(h_eCB_pt,h_eEF_pt,m_binmax,m_efficiency_nbin);
     ceff.SetConditionName(Form("L1Efficiency_eta_%dGeV",m_threshold_SA));
     ceff.SetCondition("L1 Efficiency;offline eta;Efficiency","L1efficiency_eta.png",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyeta(h_eoff_eta,h_eL1_eta,m_binmax);
     ceff.SetConditionName(Form("SAEfficiency_eta_%dGeV",m_threshold_SA));
     ceff.SetCondition("L2MuonSA Efficiency;offline eta;Efficiency","SAefficiency_eta.png",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyeta(h_eL1_eta,h_eSA_eta,m_binmax);
     ceff.SetConditionName(Form("CBEfficiency_eta_%dGeV",m_threshold_SA));
     ceff.SetCondition("muComb Efficiency;offline eta;Efficiency","CBefficiency_eta.png",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyeta(h_eSA_eta,h_eCB_eta,m_binmax);
     ceff.SetConditionName(Form("EFEfficiency_eta_%dGeV",m_threshold_SA));
     ceff.SetCondition("EventFilter Efficiency;offline eta;Efficiency","EFefficiency_eta.png",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyeta(h_eCB_eta,h_eEF_eta,m_binmax);
     ceff.SetConditionName(Form("L1Efficiency_barrel_%dGeV",m_threshold_SA));
     ceff.SetCondition("L1 Efficiency;offline pt[GeV];Efficiency","L1efficiency_pt_barrel.png",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(h_eoff_pt_barrel,h_eL1_pt_barrel,m_binmax,m_efficiency_nbin);
     ceff.SetConditionName(Form("SAEfficiency_barrel_%dGeV",m_threshold_SA));
     ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency","SAefficiency_pt_barrel.png",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(h_eL1_pt_barrel,h_eSA_pt_barrel,m_binmax,m_efficiency_nbin);
     ceff.SetConditionName(Form("CBEfficiency_barrel_%dGeV",m_threshold_SA));
     ceff.SetCondition("muComb Efficiency;offline pt[GeV];Efficiency","CBefficiency_pt_barrel.png",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(h_eSA_pt_barrel,h_eCB_pt_barrel,m_binmax,m_efficiency_nbin);
     ceff.SetConditionName(Form("EFEfficiency_barrel_%dGeV",m_threshold_SA));
     ceff.SetCondition("EventFilter Efficiency;offline pt[GeV];Efficiency","EFefficiency_pt_barrel.png",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(h_eCB_pt_barrel,h_eEF_pt_barrel,m_binmax,m_efficiency_nbin);
     ceff.SetConditionName(Form("L1Efficiency_end_%dGeV",m_threshold_SA));
     ceff.SetCondition("L1 Efficiency;offline pt[GeV];Efficiency","L1efficiency_pt_end.png",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(h_eoff_pt_end,h_eL1_pt_end,m_binmax,m_efficiency_nbin);
     ceff.SetConditionName(Form("SAEfficiency_end_%dGeV",m_threshold_SA));
     ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency","SAefficiency_pt_end.png",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(h_eL1_pt_end,h_eSA_pt_end,m_binmax,m_efficiency_nbin);
     ceff.SetConditionName(Form("CBEfficiency_end_%dGeV",m_threshold_SA));
     ceff.SetCondition("muComb Efficiency;offline pt[GeV];Efficiency","CBefficiency_pt_end.png",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(h_eSA_pt_end,h_eCB_pt_end,m_binmax,m_efficiency_nbin);
     ceff.SetConditionName(Form("EFEfficiency_end_%dGeV",m_threshold_SA));
     ceff.SetCondition("EventFilter Efficiency;offline pt[GeV];Efficiency","EFefficiency_pt_end.png",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(h_eCB_pt_end,h_eEF_pt_end,m_binmax,m_efficiency_nbin);
     ceff.SetConditionName(Form("SA2DEfficiency_%dGeV",m_threshold_SA));
     ceff.SetCondition("L1vsL2MuonSA Efficiency;offline eta;offline phi","2DSAefficiency.png",1.5,0.1,0.1,0.105,0.165);
     ceff.SetConditionbin(m_nbin_eta,m_nbin_phi,m_eta_max,m_phi_max);
     ceff.DrawEfficiency2D(h_eff_pL1_etaphi,h_eff_pSA_etaphi);
     ceff.SetConditionName(Form("L12DEfficiency_%dGeV",m_threshold_SA));
     ceff.SetCondition("L1 Efficiency;offline eta;offline phi","2DL1efficiency.png",1.5,0.1,0.1,0.105,0.165);
     ceff.SetConditionbin(m_nbin_eta,m_nbin_phi,m_eta_max,m_phi_max);
     ceff.DrawEfficiency2D(h_eff_poff_etaphi,h_eff_pL1_etaphi);
     ceff.SetConditionlabel("Off/L1 efficiency","L1/SA efficiency","SA/CB efficiency","CB/EF efficiency","error");
     ceff.SetCondition("Efficiency;offline pt[GeV];Efficiency","efficiencyAll.png",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyAll(h_eoff_pt,h_eL1_pt,h_eSA_pt,h_eCB_pt,h_eEF_pt,m_binmax,m_efficiency_nbin);
     ceff.SetCondition("Efficiency;offline pt[GeV];Efficiency","efficiencyAll_barrel.png",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyAll(h_eoff_pt_barrel,h_eL1_pt_barrel,h_eSA_pt_barrel,h_eCB_pt_barrel,h_eEF_pt_barrel,m_binmax,m_efficiency_nbin);
     ceff.SetCondition("Efficiency;offline pt[GeV];Efficiency","efficiencyAll_end.png",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyAll(h_eoff_pt_end,h_eL1_pt_end,h_eSA_pt_end,h_eCB_pt_end,h_eEF_pt_end,m_binmax,m_efficiency_nbin);

     delete h_poff_pt;
     delete h_pL1_pt;
     delete h_pSA_pt;
     delete h_pCB_pt;
     delete h_pEF_pt;
     delete h_pL1_dR;
     delete h_pSA_dR;
     delete h_pCB_dR;
     delete h_pEF_dR;
     delete h_textL1_dR;
     delete h_textSA_dR;
     delete h_textCB_dR;
     delete h_textEF_dR;
     delete h_pextL1_dR;
     delete h_pextSA_dR;
     delete h_pextCB_dR;
     delete h_pextEF_dR;
     delete h_td0;
     delete h_tz0;
     delete h_pd0;
     delete h_pz0;
     delete h_pSA_respt;
     delete h_pCB_respt;
     delete h_pEF_respt;
     delete h_invmass_off;
     delete h_invmass_L1;
     delete h_invmass_SA;
     delete h_invmass_CB;
     delete h_invmass_EF;
     delete h_eoff_pt;
     delete h_eL1_pt;
     delete h_eSA_pt;
     delete h_eCB_pt;
     delete h_eEF_pt;
     delete h_eoff_eta;
     delete h_eL1_eta;
     delete h_eSA_eta;
     delete h_eCB_eta;
     delete h_eEF_eta;
     delete h_eoff_pt_barrel;
     delete h_eL1_pt_barrel;
     delete h_eSA_pt_barrel;
     delete h_eCB_pt_barrel;
     delete h_eEF_pt_barrel;
     delete h_eoff_pt_end;
     delete h_eL1_pt_end;
     delete h_eSA_pt_end;
     delete h_eCB_pt_end;
     delete h_eEF_pt_end;
     delete h_eff_poff_etaphi;
     delete h_eff_pL1_etaphi;
     delete h_eff_pSA_etaphi;
     delete h_poffvsSA_pt;
     delete h_poffvsSA_pt_bad;
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
