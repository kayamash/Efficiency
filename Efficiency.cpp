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
#include "CalcEfficiency.cpp"

void Efficiency::Init(TTree *tree,std::string name,const Int_t np,const Int_t ne,const Double_t mp,const Double_t me,Double_t req,Int_t max,Double_t err,const Int_t nh,const Int_t th){
     if(tree){
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
     	m_pSA_sAddress = 0;
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
     	m_pEFTAG_pass = 0;
     	m_reqL1dR = req;
     	m_tp_extdR = 0;
     	m_sumReqdRL1 = 0;
     	m_sumReqdREF = 0;
     	m_tp_dR = 0;
     	m_pSA_phims = 0;
     	m_pSA_roiphi = 0;
     	m_countall = 0;
     	m_countoff = 0;
     	m_countL1 = 0;
     	m_countSA = 0;
     	m_countCB = 0;
     	m_countEF = 0;
     	m_count = 0;

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
     	tChain->SetBranchStatus("probe_mesSA_sAddress",1);
     	tChain->SetBranchStatus("probe_mesSA_phims",1);
     	tChain->SetBranchStatus("probe_mesSA_roiPhi",1);
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
     	tChain->SetBranchStatus("probe_mesEFTAG_pass",1);
     	tChain->SetBranchStatus("tp_extdR",1);
     	tChain->SetBranchStatus("sumReqdRL1",1);
     	tChain->SetBranchStatus("sumReqdREF",1);
     	tChain->SetBranchStatus("tp_dR",1);
     	//setting each branch address
    	tChain->SetBranchAddress("mes_name",&m_mes_name,&b_mes_name);
     	tChain->SetBranchAddress("sumReqdRL1",&m_sumReqdRL1,&b_sumReqdRL1);
     	tChain->SetBranchAddress("sumReqdREF",&m_sumReqdREF,&b_sumReqdREF);
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
     	tChain->SetBranchAddress("tp_dR",&m_tp_dR,&b_tp_dR);
     	tChain->SetBranchAddress("tp_extdR",&m_tp_extdR,&b_tp_extdR);
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
     	tChain->SetBranchAddress("probe_mesSA_sAddress",&m_pSA_sAddress,&b_pSA_sAddress);
     	tChain->SetBranchAddress("probe_mesSA_phims",&m_pSA_phims,&b_pSA_phims);
     	tChain->SetBranchAddress("probe_mesSA_roiPhi",&m_pSA_roiphi,&b_pSA_roiphi);
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
     	tChain->SetBranchAddress("probe_mesEFTAG_pass",&m_pEFTAG_pass,&b_pEFTAG_pass);

     	//define each histgram
     	m_h_offphi_LargeSpecial = new TH1D("h_offphi_LargeSpecial_0GeV","offline phi;offline phi;Entries",600,-3.0,3.0);
     	m_h_saphims_LargeSpecial = new TH1D("h_saphims_LargeSpecial_0GeV","L2MuonSA phi;L2MuonSA phims;Entries",600,-3.0,3.0);
     	m_h_saroiphi_LargeSpecial = new TH1D("h_saroiphi_LargeSpecial_0GeV","RoI phi;RoI phi;Entries",600,-3.0,3.0);
          m_h_saroiphi_SmallSpecial = new TH1D("h_saroiphi_SmallSpecial_0GeV","RoI phi;RoI phi;Entries",600,-3.0,3.0);
     	for(Int_t i = 0;i <= m_nhist;i++){
     		m_h_poff_pt.push_back(new TH1D(Form("h_poff_pt_%dGeV",i*m_thpitch),"probe offline pt;offline pt[GeV];Entries",150,0,150));
        	m_h_pL1_pt.push_back(new TH1D(Form("h_pL1_pt_%dGeV",i*m_thpitch),"probe L1 pt;L1 pt[GeV];Entries",150,0,150));
        	m_h_pSA_pt.push_back(new TH1D(Form("h_pSA_pt_%dGeV",i*m_thpitch),"probe L2MuonSA pt;L2MuonSA pt[GeV];Entries",150,0,150)); 
        	m_h_pCB_pt.push_back(new TH1D(Form("h_pCB_pt_%dGeV",i*m_thpitch),"probe muComb pt;muComb pt[GeV];Entries",150,0,150));
        	m_h_pEF_pt.push_back(new TH1D(Form("h_pEF_pt_%dGeV",i*m_thpitch),"probe EventFilter pt;EventFilter pt[GeV];Entries",150,0,150));
        	m_h_pL1_dR.push_back(new TH1D(Form("h_L1 dR_%dGeV",i*m_thpitch),"L1 dR;dR;Entries",1000,0,0.1));
        	m_h_pSA_dR.push_back(new TH1D(Form("h_SA dR_%dGeV",i*m_thpitch),"L2MuonSA dR;dR;Entries",500,0,0.05));
          	m_h_pCB_dR.push_back(new TH1D(Form("h_CB dR_%dGeV",i*m_thpitch),"muComb dR;dR;Entries",200,0,0.002));
          	m_h_pEF_dR.push_back(new TH1D(Form("h_EF dR_%dGeV",i*m_thpitch),"EventFilter dR;dR;Entries",100,0,0.001));
          	m_h_textL1_dR.push_back(new TH1D(Form("h_textL1_dR_%dGeV",i*m_thpitch),"tag extL1 dR;dR;Entries",1000,0,0.1));
          	m_h_textSA_dR.push_back(new TH1D(Form("h_textSA_dR_%dGeV",i*m_thpitch),"tag extL2MuonSA dR;dR;Entries",1000,0,0.1));
          	m_h_textCB_dR.push_back(new TH1D(Form("h_textCB_dR_%dGeV",i*m_thpitch),"tag extmuComb dR;dR;Entries",1000,0,0.01));
          	m_h_textEF_dR.push_back(new TH1D(Form("h_textEF_dR_%dGeV",i*m_thpitch),"tag extEventFilter dR;dR;Entries",1000,0,0.01));
          	m_h_pextL1_dR.push_back(new TH1D(Form("h_pextL1_dR_%dGeV",i*m_thpitch),"probe extL1 dR;dR;Entries",1000,0,0.1));
          	m_h_pextSA_dR.push_back(new TH1D(Form("h_pextSA_dR_%dGeV",i*m_thpitch),"probe extL2MuonSA dR;dR;Entries",1000,0,0.1));
          	m_h_pextCB_dR.push_back(new TH1D(Form("h_pextCB_dR_%dGeV",i*m_thpitch),"probe extmuComb dR;dR;Entries",1000,0,0.01));
          	m_h_pextEF_dR.push_back(new TH1D(Form("h_pextEF_dR_%dGeV",i*m_thpitch),"probe extEventFilter dR;dR;Entries",1000,0,0.01));
          	m_h_pSA_respt.push_back(new TH1D(Form("h_pSA_respt_%dGeV",i*m_thpitch),"probe L2MuonSA residual pt;residual pt;Entries",300,-1,1));
          	m_h_pCB_respt.push_back(new TH1D(Form("h_pCB_respt_%dGeV",i*m_thpitch),"probe muComb residual pt;residual pt;Entries",640,-0.8,0.8));
          	m_h_pEF_respt.push_back(new TH1D(Form("h_pEF_respt_%dGeV",i*m_thpitch),"probe EventFilter residual pt;residual pt;Entries",1800,-0.3,0.3));
          	m_h_eoff_pt.push_back(new TH1D(Form("h_eoff_pt_%dGeV",i*m_thpitch),"mesoff_pt;offline pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eL1_pt.push_back(new TH1D(Form("h_eL1_pt_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eSA_pt.push_back(new TH1D(Form("h_eSA_pt_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eCB_pt.push_back(new TH1D(Form("h_eCB_pt_%dGeV",i*m_thpitch),"mesCB_pt;CB pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eEF_pt.push_back(new TH1D(Form("h_eEF_pt_%dGeV",i*m_thpitch),"mesEF_pt;EF pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eoff_eta.push_back(new TH1D(Form("h_eoff_eta_%dGeV",i*m_thpitch),"off eta;off eta;Entries",50,-2.5,2.5));
          	m_h_eL1_eta.push_back(new TH1D(Form("h_eL1_eta_%dGeV",i*m_thpitch),"L1 eta;L1 eta;Entries",50,-2.5,2.5));
          	m_h_eSA_eta.push_back(new TH1D(Form("h_eSA_eta_%dGeV",i*m_thpitch),"SA eta;SA eta;Entries",50,-2.5,2.5));
          	m_h_eCB_eta.push_back(new TH1D(Form("h_eCB_eta_%dGeV",i*m_thpitch),"CB eta;CB eta;Entries",50,-2.5,2.5));
          	m_h_eEF_eta.push_back(new TH1D(Form("h_eEF_eta_%dGeV",i*m_thpitch),"EF eta;EF eta;Entries",50,-2.5,2.5));
          	m_h_eoff_pt_barrel.push_back(new TH1D(Form("h_eoff_pt_barrel_%dGeV",i*m_thpitch),"off_pt;off pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eL1_pt_barrel.push_back(new TH1D(Form("h_eL1_pt_barrel_%dGeV",i*m_thpitch),"L1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eSA_pt_barrel.push_back(new TH1D(Form("h_eSA_pt_barrel_%dGeV",i*m_thpitch),"SA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eCB_pt_barrel.push_back(new TH1D(Form("h_eCB_pt_barrel_%dGeV",i*m_thpitch),"CB_pt;CB pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eEF_pt_barrel.push_back(new TH1D(Form("h_eEF_pt_barrel_%dGeV",i*m_thpitch),"EF_pt;EF pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eoff_pt_end.push_back(new TH1D(Form("h_eoff_pt_end_%dGeV",i*m_thpitch),"off_pt;off pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eL1_pt_end.push_back(new TH1D(Form("h_eL1_pt_end_%dGeV",i*m_thpitch),"L1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eSA_pt_end.push_back(new TH1D(Form("h_eSA_pt_end_%dGeV",i*m_thpitch),"SA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eCB_pt_end.push_back(new TH1D(Form("h_eCB_pt_end_%dGeV",i*m_thpitch),"CB_pt;CB pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eEF_pt_end.push_back(new TH1D(Form("h_eEF_pt_end_%dGeV",i*m_thpitch),"EF_pt;EF pt[GeV];Entries",300,-0.25,149.75));
          	m_h_SA_respt0.push_back(new TH1D(Form("h_SA_respt0_%dGeV",i*m_thpitch),"SAresidualpt_Large;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
          	m_h_SA_respt1.push_back(new TH1D(Form("h_SA_respt1_%dGeV",i*m_thpitch),"SAresidualpt_LargeSpecial;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
          	m_h_SA_respt2.push_back(new TH1D(Form("h_SA_respt2_%dGeV",i*m_thpitch),"SAresidualpt_Small;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
          	m_h_SA_respt3.push_back(new TH1D(Form("h_SA_respt3_%dGeV",i*m_thpitch),"SAresidualpt_SmallSpecial;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
          	m_h_eff_poff_etaphi.push_back(new TH2F(Form("h_eff_poff_etaphi_%dGeV",i*m_thpitch),"offlineeta vs offlinephi;offline eta;offline phi",m_nbin_eta,-m_eta_max,m_eta_max,m_nbin_phi,-m_phi_max,m_phi_max));
          	m_h_eff_pL1_etaphi.push_back(new TH2F(Form("h_eff_pL1_etaphi_%dGeV",i*m_thpitch),"L1eta vs L1phi;L1 eta;L1 phi",m_nbin_eta,-m_eta_max,m_eta_max,m_nbin_phi,-m_phi_max,m_phi_max));
          	m_h_eff_pSA_etaphi.push_back(new TH2F(Form("h_eff_pSA_etaphi_%dGeV",i*m_thpitch),"L2MuonSAeta vs L2MuonSAphi;L2MuonSA eta;L2MuonSA phi",m_nbin_eta,-m_eta_max,m_eta_max,m_nbin_phi,-m_phi_max,m_phi_max));
          	m_h_poffvsSA_pt.push_back(new TH2F(Form("h_poffvsSA_pt_%dGeV",i*m_thpitch),"probe offline pt vs probe L2MuonSA pt@mu26ivm;probe offline pt[GeV];probe L2MuonSA pt[GeV]",150,0,150,150,0,150));
          	m_h_off_ptvsSA_resptplus0.push_back(new TH2F(Form("h_off_ptvsSA_resptplus0_%dGeV",i*m_thpitch),"Large Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
          	m_h_off_ptvsSA_resptplus1.push_back(new TH2F(Form("h_off_ptvsSA_resptplus1_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
          	m_h_off_ptvsSA_resptplus2.push_back(new TH2F(Form("h_off_ptvsSA_resptplus2_%dGeV",i*m_thpitch),"Small Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
          	m_h_off_ptvsSA_resptplus3.push_back(new TH2F(Form("h_off_ptvsSA_resptplus3_%dGeV",i*m_thpitch),"SmallSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
          	m_h_off_ptvsSA_resptminus0.push_back(new TH2F(Form("h_off_ptvsSA_resptminus0_%dGeV",i*m_thpitch),"Large Qeta/|eta|=-1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
          	m_h_off_ptvsSA_resptminus1.push_back(new TH2F(Form("h_off_ptvsSA_resptminus1_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=-1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
          	m_h_off_ptvsSA_resptminus2.push_back(new TH2F(Form("h_off_ptvsSA_resptminus2_%dGeV",i*m_thpitch),"Small Qeta/|eta|=-1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
          	m_h_off_ptvsSA_resptminus3.push_back(new TH2F(Form("h_off_ptvsSA_resptminus3_%dGeV",i*m_thpitch),"SmallSpecial Qeta/|eta|=-1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
          	m_h_offphivsSA_sAddress.push_back(new TH2F(Form("h_offphivsSA_sAddress_%dGeV",i*m_thpitch),"offline phi vs sAddress;offline phi;sAddress",140,-3.5,3.5,4,0.0,4.0));
          	m_h_offphivsSA_respt0.push_back(new TH2F(Form("h_offphivsSA_respt0_%dGeV",i*m_thpitch),"Large;offline phi;SApt residual",140,-3.5,3.5,300,-2.0,1.0));
          	m_h_offphivsSA_respt1.push_back(new TH2F(Form("h_offphivsSA_respt1_%dGeV",i*m_thpitch),"LargeSpecial;offline phi;SApt residual",140,-3.5,3.5,300,-2.0,1.0));
          	m_h_offphivsSA_respt2.push_back(new TH2F(Form("h_offphivsSA_respt2_%dGeV",i*m_thpitch),"Small;offline phi;SApt residual",140,-3.5,3.5,300,-2.0,1.0));
          	m_h_offphivsSA_respt3.push_back(new TH2F(Form("h_offphivsSA_respt3_%dGeV",i*m_thpitch),"SmallSpecial;offline phi;SApt residual",140,-3.5,3.5,300,-2.0,1.0));
          	m_h_eoff_pt_Large.push_back(new TH1D(Form("h_eoff_ptLarge_%dGeV",i*m_thpitch),"mesoff_pt;offline pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eL1_pt_Large.push_back(new TH1D(Form("h_eL1_ptLarge_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eSA_pt_Large.push_back(new TH1D(Form("h_eSA_ptLarge_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eoff_pt_LargeSpecial.push_back(new TH1D(Form("h_eoff_ptLargeSpecial_%dGeV",i*m_thpitch),"mesoff_pt;offline pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eL1_pt_LargeSpecial.push_back(new TH1D(Form("h_eL1_ptLargeSpecial_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eSA_pt_LargeSpecial.push_back(new TH1D(Form("h_eSA_ptLargeSpecial_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eL1_pt_LargeSpecial11.push_back(new TH1D(Form("h_eL1_ptLargeSpecial11_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eSA_pt_LargeSpecial11.push_back(new TH1D(Form("h_eSA_ptLargeSpecial11_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eL1_pt_LargeSpecial15.push_back(new TH1D(Form("h_eL1_ptLargeSpecial15_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eSA_pt_LargeSpecial15.push_back(new TH1D(Form("h_eSA_ptLargeSpecial15_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eL1_pt_LargeSpecialplus11.push_back(new TH1D(Form("h_eL1_ptLargeSpecialplus11_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eSA_pt_LargeSpecialplus11.push_back(new TH1D(Form("h_eSA_ptLargeSpecialplus11_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eL1_pt_LargeSpecialminus11.push_back(new TH1D(Form("h_eL1_ptLargeSpecialminus11_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eSA_pt_LargeSpecialminus11.push_back(new TH1D(Form("h_eSA_ptLargeSpecialminus11_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eL1_pt_LargeSpecialplus15.push_back(new TH1D(Form("h_eL1_ptLargeSpecialplus15_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eSA_pt_LargeSpecialplus15.push_back(new TH1D(Form("h_eSA_ptLargeSpecialplus15_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eL1_pt_LargeSpecialminus15.push_back(new TH1D(Form("h_eL1_ptLargeSpecialminus15_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eSA_pt_LargeSpecialminus15.push_back(new TH1D(Form("h_eSA_ptLargeSpecialminus15_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eoff_pt_Small.push_back(new TH1D(Form("h_eoff_ptSmall_%dGeV",i*m_thpitch),"mesoff_pt;offline pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eL1_pt_Small.push_back(new TH1D(Form("h_eL1_ptSmall_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eSA_pt_Small.push_back(new TH1D(Form("h_eSA_ptSmall_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eoff_pt_SmallSpecial.push_back(new TH1D(Form("h_eoff_ptSmallSpecial_%dGeV",i*m_thpitch),"mesoff_pt;offline pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eL1_pt_SmallSpecial.push_back(new TH1D(Form("h_eL1_ptSmallSpecial_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
          	m_h_eSA_pt_SmallSpecial.push_back(new TH1D(Form("h_eSA_ptSmallSpecial_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
          	m_h_offphivsSAphims.push_back(new TH2F(Form("h_offphivsSAphims_%dGeV",i*m_thpitch),"offphi vs phims;offline phi;L2MuonSA phims",140,-3.5,3.5,140,-3.5,3.5));
          	m_h_off_ptvsSA_resptplusLS11.push_back(new TH2F(Form("h_off_ptvsSA_resptplusLS11_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
          	m_h_off_ptvsSA_resptplusLS15.push_back(new TH2F(Form("h_off_ptvsSA_resptplusLS15_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
          	m_h_off_ptvsSA_resptplusLSplus11.push_back(new TH2F(Form("h_off_ptvsSA_resptplusLSplus11_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
          	m_h_off_ptvsSA_resptplusLSplus15.push_back(new TH2F(Form("h_off_ptvsSA_resptplusLSplus15_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
          	m_h_off_ptvsSA_resptplusLSminus11.push_back(new TH2F(Form("h_off_ptvsSA_resptplusLSminus11_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
          	m_h_off_ptvsSA_resptplusLSminus15.push_back(new TH2F(Form("h_off_ptvsSA_resptplusLSminus15_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
          	m_h_off_ptvsSA_resptminusLS11.push_back(new TH2F(Form("h_off_ptvsSA_resptminusLS11_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
          	m_h_off_ptvsSA_resptminusLS15.push_back(new TH2F(Form("h_off_ptvsSA_resptminusLS15_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
          	m_h_off_ptvsSA_resptminusLSplus11.push_back(new TH2F(Form("h_off_ptvsSA_resptminusLSplus11_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
          	m_h_off_ptvsSA_resptminusLSplus15.push_back(new TH2F(Form("h_off_ptvsSA_resptminusLSplus15_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
          	m_h_off_ptvsSA_resptminusLSminus11.push_back(new TH2F(Form("h_off_ptvsSA_resptminusLSminus11_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
          	m_h_off_ptvsSA_resptminusLSminus15.push_back(new TH2F(Form("h_off_ptvsSA_resptminusLSminus15_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));	
          	m_h_SA_resptplus11.push_back(new TH1D(Form("h_SA_resptLS11+_%dGeV",i*m_thpitch),"SAresidualpt_LS11+;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_SA_resptminus11.push_back(new TH1D(Form("h_SA_resptLS11-_%dGeV",i*m_thpitch),"SAresidualpt_LS11-;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_SA_resptplus15.push_back(new TH1D(Form("h_SA_resptLS15+_%dGeV",i*m_thpitch),"SAresidualpt_LS15+;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_SA_resptminus15.push_back(new TH1D(Form("h_SA_resptLS15-_%dGeV",i*m_thpitch),"SAresidualpt_LS15-;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_offetavsSA_respt0.push_back(new TH2F(Form("h_offetavsSA_respt0_%dGeV",i*m_thpitch),"Large;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_respt1.push_back(new TH2F(Form("h_offetavsSA_respt1_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_respt2.push_back(new TH2F(Form("h_offetavsSA_respt2_%dGeV",i*m_thpitch),"Small;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_respt3.push_back(new TH2F(Form("h_offetavsSA_respt3_%dGeV",i*m_thpitch),"SmallSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_resptLSplus11.push_back(new TH2F(Form("h_offetavsSA_resptLS11+_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_resptLSminus11.push_back(new TH2F(Form("h_offetavsSA_resptLS11-_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_resptLSplus15.push_back(new TH2F(Form("h_offetavsSA_resptLS15+_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_resptLSminus15.push_back(new TH2F(Form("h_offetavsSA_resptLS15-_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_highoffetavsSA_resptLSplus11.push_back(new TH2F(Form("h_highoffetavsSA_resptLS11+_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_highoffetavsSA_resptLSminus11.push_back(new TH2F(Form("h_highoffetavsSA_resptLS11-_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_highoffetavsSA_resptLSplus15.push_back(new TH2F(Form("h_highoffetavsSA_resptLS15+_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_highoffetavsSA_resptLSminus15.push_back(new TH2F(Form("h_highoffetavsSA_resptLS15-_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_highoffphivsSA_resptLSplus11.push_back(new TH2F(Form("h_highoffphivsSA_resptLS11+_%dGeV",i*m_thpitch),"LargeSpecial;offline phi;SApt residual",140,-3.5,3.5,300,-2.0,1.0));
               m_h_highoffphivsSA_resptLSminus11.push_back(new TH2F(Form("h_highoffphivsSA_resptLS11-_%dGeV",i*m_thpitch),"LargeSpecial;offline phi;SApt residual",140,-3.5,3.5,300,-2.0,1.0));
               m_h_highoffphivsSA_resptLSplus15.push_back(new TH2F(Form("h_highoffphivsSA_resptLS15+_%dGeV",i*m_thpitch),"LargeSpecial;offline phi;SApt residual",140,-3.5,3.5,300,-2.0,1.0));
               m_h_highoffphivsSA_resptLSminus15.push_back(new TH2F(Form("h_highoffphivsSA_resptLS15-_%dGeV",i*m_thpitch),"LargeSpecial;offline phi;SApt residual",140,-3.5,3.5,300,-2.0,1.0));

               m_countLarge.push_back(0);
          	m_countLargeSpecial.push_back(0);
          	m_countSmall.push_back(0);
          	m_countSmallSpecial.push_back(0);
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

bool Efficiency::Cut_tagprobe(Int_t pass){
     if(m_sumReqdRL1 < m_tp_extdR && 0.2 < m_tp_extdR && m_sumReqdREF < m_tp_dR && pass > -1){
          return kTRUE;
     }else{
          return kFALSE;
     }
}

bool Efficiency::Cut_L1(Int_t pass){
  if(pass > -1){
          m_countL1++;
	  return kTRUE;
     }else{
          return kFALSE;
     }
}

bool Efficiency::Cut_SA(Int_t pass,Double_t pt,Double_t th){
     if(pass == 1 && std::fabs(pt) > th){
 	  m_countSA++;
          return kTRUE;
     }else{
          return kFALSE;
     }
}

bool Efficiency::Cut_CB(Int_t pass){
     if(pass == 1){
	  m_countCB++;
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
     for(Int_t i = 0;i <= m_nhist;i++){     
          m_count++;
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
          Double_t pSA_sAddress = -1;
          Double_t pSA_phims = -99999;
          float pSA_roiphi = -99999;
          Double_t pCB_pt = -99999;
          Double_t pCB_eta = 0;
          Double_t pCB_phi = 0;
          Double_t pCB_dR = 1;
          Int_t pCB_pass = 0;
          Double_t pEF_pt = -99999;
          Double_t pEF_eta = 0;
          Double_t pEF_phi = 0;
          Double_t pEF_dR = 1;
          Double_t tL1_dR = 0;
          Double_t tEF_dR = 0;
          Int_t pEF_pass = 0;
          Int_t pEFTAG_pass = -1;
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
                    pEFTAG_pass = m_pEFTAG_pass->at(method);
                    pSA_sAddress = m_pSA_sAddress->at(method);
                    pSA_phims = m_pSA_phims->at(method);
                    pSA_roiphi = m_pSA_roiphi->at(method);
               }
          }
          m_countall++;
          tL1_dR = TMath::Sqrt(pow(m_tL1_eta - m_toff_eta,2) + pow(m_tL1_phi - m_toff_phi,2) );
          tEF_dR = TMath::Sqrt(pow(m_tEF_eta - m_toff_eta,2) + pow(m_tEF_phi - m_toff_phi,2) );
          if(std::fabs(m_toff_pt)*0.001 < 10.0)m_reqL1dR = -0.00001*std::fabs(m_toff_pt) + 0.18;
          if(Cut_tagprobe(pEFTAG_pass)){
               m_countoff++;
	          //offline
               if(i == 0 && static_cast<Int_t>(pSA_sAddress) == 1)m_h_offphi_LargeSpecial->Fill(m_poff_phi);
               m_countoff++;
               m_h_poff_pt.at(i)->Fill(m_poff_pt*0.001);
               m_h_eoff_pt.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(std::fabs(m_poff_pt*0.001) > 40)m_h_eoff_eta.at(i)->Fill(m_poff_eta);
               if(Dicision_barrel(m_poff_eta)){
                    m_h_eoff_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               }else{
                    m_h_eoff_pt_end.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               }
               if(std::fabs(m_poff_pt*0.001) > 40)m_h_eff_poff_etaphi.at(i)->Fill(m_poff_eta,m_poff_phi);
               switch(static_cast<Int_t>(pSA_sAddress)){
                    case 0:
                         m_h_eoff_pt_Large.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                    case 1:
                         m_h_eoff_pt_LargeSpecial.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                    case 2:
                         m_h_eoff_pt_Small.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                    case 3:
                         m_h_eoff_pt_SmallSpecial.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                    default:
                         break;
               }

               //L1
               if(Cut_L1(pL1_pass)){
                    m_countL1++;
                    Double_t textL1_dR = TMath::Sqrt(pow(m_tL1_eta - m_toff_exteta,2) + pow(m_tL1_phi - m_toff_extphi,2));
                    pextL1_dR = TMath::Sqrt(pow(pL1_eta - m_poff_exteta,2) + pow(pL1_phi - m_poff_extphi,2));

                    m_h_pL1_pt.at(i)->Fill(std::fabs(pL1_pt*0.001));
                    m_h_pL1_dR.at(i)->Fill(pL1_dR);
                    m_h_textL1_dR.at(i)->Fill(textL1_dR);
                    m_h_pextL1_dR.at(i)->Fill(pextL1_dR);
                    m_h_eL1_pt.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                    if(Dicision_barrel(m_poff_eta)){
                         m_h_eL1_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                    }else{
                         m_h_eL1_pt_end.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                    }
                    if(std::fabs(m_poff_pt*0.001) > 40){
                         m_h_eff_pL1_etaphi.at(i)->Fill(m_poff_eta,m_poff_phi);
                         m_h_eL1_eta.at(i)->Fill(m_poff_eta);
                    }
                    switch(static_cast<Int_t>(pSA_sAddress)){
                         case 0:
                              m_h_eL1_pt_Large.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                         case 1:
                              m_h_eL1_pt_LargeSpecial.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              if(pSA_roiphi < -0.6 && pSA_roiphi > -1.0){//15
                                   m_h_eL1_pt_LargeSpecial15.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                   if(pSA_roiphi > -0.8){
                                        m_h_eL1_pt_LargeSpecialplus15.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                   }else{
                                        m_h_eL1_pt_LargeSpecialminus15.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                   }
                              }
                              if(pSA_roiphi < -2.0 && pSA_roiphi > -2.6){//11
                                   m_h_eL1_pt_LargeSpecial11.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                   if(pSA_roiphi > -2.4){
                                        m_h_eL1_pt_LargeSpecialplus11.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                   }else{
                                        m_h_eL1_pt_LargeSpecialminus11.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                   }
                              }
                              break;
                         case 2:
                              m_h_eL1_pt_Small.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                         case 3:
                              m_h_eL1_pt_SmallSpecial.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                         default:
                              break;
                    }

                    //SA
                    if(Cut_SA(pSA_pass,pSA_pt,i*m_thpitch)){
                         m_countSA++;
                         Double_t textSA_dR = TMath::Sqrt(pow(m_tSA_eta - m_toff_exteta,2) + pow(m_tSA_phi - m_toff_extphi,2));
                         pextSA_dR = TMath::Sqrt(pow(pSA_eta - m_poff_exteta,2) + pow(pSA_phi - m_poff_extphi,2));
                         Double_t resSA_pt = std::fabs(m_poff_pt*0.001)/std::fabs(pSA_pt) - 1.0;
                         Double_t buf_pSA_dR = TMath::Sqrt(pow(pSA_eta - m_poff_eta,2) + pow(pSA_phi - m_poff_phi,2));

                         m_h_pSA_pt.at(i)->Fill(std::fabs(pSA_pt));
                         m_h_pSA_dR.at(i)->Fill(buf_pSA_dR);
                         m_h_textSA_dR.at(i)->Fill(textSA_dR);
                         m_h_pextSA_dR.at(i)->Fill(pextSA_dR);
                         m_h_eSA_pt.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_pSA_respt.at(i)->Fill(resSA_pt);
                         if(Dicision_barrel(m_poff_eta)){
                              m_h_eSA_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         }else{
                              m_h_eSA_pt_end.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         }
                         if(std::fabs(m_poff_pt*0.001) > 40){
                              m_h_eff_pSA_etaphi.at(i)->Fill(m_poff_eta,m_poff_phi);
                              m_h_eSA_eta.at(i)->Fill(m_poff_eta);
                         }
                         m_h_poffvsSA_pt.at(i)->Fill(std::fabs(m_poff_pt*0.001),std::fabs(pSA_pt));
                         switch(static_cast<Int_t>(pSA_sAddress)){
                              case 0:
                                   if(m_probe_charge*m_poff_eta/std::fabs(m_poff_eta)==1)m_h_off_ptvsSA_resptplus0.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                   if(m_probe_charge*m_poff_eta/std::fabs(m_poff_eta)==-1)m_h_off_ptvsSA_resptminus0.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                   m_h_SA_respt0.at(i)->Fill(resSA_pt);
                                   m_h_eSA_pt_Large.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                   m_h_offphivsSA_respt0.at(i)->Fill(m_poff_phi,resSA_pt);
                                   m_h_offetavsSA_respt0.at(i)->Fill(m_poff_eta,resSA_pt);
                                   m_countLarge.at(i)++;
                                   break;
                              case 1:
                                   if(m_probe_charge*m_poff_eta/std::fabs(m_poff_eta)==1){
                                        m_h_off_ptvsSA_resptplus1.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                        if(pSA_roiphi < -0.6 && pSA_roiphi > -1.0){//15
                                             m_h_off_ptvsSA_resptplusLS15.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                             if(pSA_roiphi > -0.8){
                                                  m_h_off_ptvsSA_resptplusLSplus15.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                                  m_h_SA_resptplus15.at(i)->Fill(resSA_pt);
                                                  m_h_offetavsSA_resptLSplus15.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  if(std::fabs(m_poff_pt*0.001) > 30.0 && std::fabs(m_poff_pt*0.001) < 50.0){
                                                       m_h_highoffetavsSA_resptLSplus15.at(i)->Fill(m_poff_eta,resSA_pt);
                                                       m_h_highoffphivsSA_resptLSplus15.at(i)->Fill(m_poff_phi,resSA_pt);
                                                  }
                                             }else{
                                                  m_h_off_ptvsSA_resptplusLSminus15.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                                  m_h_offetavsSA_resptLSminus15.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  m_h_SA_resptminus15.at(i)->Fill(resSA_pt);
                                                  if(std::fabs(m_poff_pt*0.001) > 30.0 && std::fabs(m_poff_pt*0.001) < 50.0){
                                                       m_h_highoffetavsSA_resptLSminus15.at(i)->Fill(m_poff_eta,resSA_pt);
                                                       m_h_highoffphivsSA_resptLSminus15.at(i)->Fill(m_poff_phi,resSA_pt);
                                                  }
                                             }
                                        }
                                        if(pSA_roiphi < -2.0 && pSA_roiphi > -2.6){//11
                                             m_h_off_ptvsSA_resptplusLS11.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                             if(pSA_roiphi > -2.4){
                                                  m_h_off_ptvsSA_resptplusLSplus11.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                                  m_h_offetavsSA_resptLSplus11.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  m_h_SA_resptplus11.at(i)->Fill(resSA_pt);
                                                  if(std::fabs(m_poff_pt*0.001) > 30.0 && std::fabs(m_poff_pt*0.001) < 50.0){
                                                       m_h_highoffetavsSA_resptLSplus11.at(i)->Fill(m_poff_eta,resSA_pt);
                                                       m_h_highoffphivsSA_resptLSplus11.at(i)->Fill(m_poff_phi,resSA_pt);
                                                  }
                                             }else{
                                                  m_h_off_ptvsSA_resptplusLSminus11.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                                  m_h_offetavsSA_resptLSminus11.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  m_h_SA_resptminus11.at(i)->Fill(resSA_pt);
                                                  if(std::fabs(m_poff_pt*0.001) > 30.0 && std::fabs(m_poff_pt*0.001) < 50.0){
                                                       m_h_highoffetavsSA_resptLSminus11.at(i)->Fill(m_poff_eta,resSA_pt);
                                                       m_h_highoffphivsSA_resptLSminus11.at(i)->Fill(m_poff_phi,resSA_pt);
                                                  }
                                             }
                                        }
                                   }
                                   if(m_probe_charge*m_poff_eta/std::fabs(m_poff_eta)==-1){
                                        m_h_off_ptvsSA_resptminus1.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                        if(pSA_roiphi < -0.6 && pSA_roiphi > -1.0){//15
                                             m_h_off_ptvsSA_resptminusLS15.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                             if(pSA_roiphi > -0.8){
                                                  m_h_off_ptvsSA_resptminusLSplus15.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                             }else{
                                                  m_h_off_ptvsSA_resptminusLSminus15.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                             }
                                        }
                                        if(pSA_roiphi < -2.0 && pSA_roiphi > -2.6){//11
                                             m_h_off_ptvsSA_resptminusLS11.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                             if(pSA_roiphi > -2.4){
                                                  m_h_off_ptvsSA_resptminusLSplus11.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                             }else{
                                                  m_h_off_ptvsSA_resptminusLSminus11.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                             }
                                        }
                                   }
                                   m_h_SA_respt1.at(i)->Fill(resSA_pt);
                                   m_h_eSA_pt_LargeSpecial.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                   m_h_offphivsSA_respt1.at(i)->Fill(m_poff_phi,resSA_pt);
                                   m_h_offetavsSA_respt1.at(i)->Fill(m_poff_eta,resSA_pt);
                                   if(i == 0){
                                        m_h_saphims_LargeSpecial->Fill(pSA_phims);
                                        m_h_saroiphi_LargeSpecial->Fill(pSA_roiphi);
                                   }

                                   if(pSA_roiphi < -0.6 && pSA_roiphi > -1.0){//15
                                        m_h_eSA_pt_LargeSpecial15.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        if(pSA_roiphi > -0.8){
                                             m_h_eSA_pt_LargeSpecialplus15.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        }else{
                                             m_h_eSA_pt_LargeSpecialminus15.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        }
                                   }
                                   if(pSA_roiphi < -2.0 && pSA_roiphi > -2.6){//11
                                        m_h_eSA_pt_LargeSpecial11.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        if(pSA_roiphi > -2.4){
                                             m_h_eSA_pt_LargeSpecialplus11.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        }else{
                                             m_h_eSA_pt_LargeSpecialminus11.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        }
                                   }
                                   m_countLargeSpecial.at(i)++;
                                   break;
                              case 2:
                                   if(m_probe_charge*m_poff_eta/std::fabs(m_poff_eta)==1)m_h_off_ptvsSA_resptplus2.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                   if(m_probe_charge*m_poff_eta/std::fabs(m_poff_eta)==-1)m_h_off_ptvsSA_resptminus2.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                   m_h_SA_respt2.at(i)->Fill(resSA_pt);
                                   m_h_eSA_pt_Small.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                   m_h_offphivsSA_respt2.at(i)->Fill(m_poff_phi,resSA_pt);
                                   m_h_offetavsSA_respt2.at(i)->Fill(m_poff_eta,resSA_pt);
                                   m_countSmall.at(i)++;
                                   break;
                              case 3:
                                   if(m_probe_charge*m_poff_eta/std::fabs(m_poff_eta)==1)m_h_off_ptvsSA_resptplus3.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                   if(m_probe_charge*m_poff_eta/std::fabs(m_poff_eta)==-1)m_h_off_ptvsSA_resptminus3.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                   m_h_SA_respt3.at(i)->Fill(resSA_pt);
                                   m_h_eSA_pt_SmallSpecial.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                   m_h_offphivsSA_respt3.at(i)->Fill(m_poff_phi,resSA_pt);
                                   m_h_offetavsSA_respt3.at(i)->Fill(m_poff_eta,resSA_pt);
                                   if(i == 0)m_h_saroiphi_SmallSpecial->Fill(pSA_roiphi);
                                   m_countSmallSpecial.at(i)++;
                                   break;
                              default:
                                   break;
                         }
                         if(static_cast<Int_t>(pSA_sAddress) == 0 || static_cast<Int_t>(pSA_sAddress) == 1 || static_cast<Int_t>(pSA_sAddress) == 2 || static_cast<Int_t>(pSA_sAddress) == 3){
                              m_h_offphivsSA_sAddress.at(i)->Fill(pSA_phims,pSA_sAddress);
                              m_h_offphivsSA_respt.at(i)->Fill(m_poff_phi,resSA_pt);
                         }
                         m_h_offphivsSAphims.at(i)->Fill(m_poff_phi,pSA_phims);

                         //CB
                         if(Cut_CB(pCB_pass)){
                              m_countCB++;
                              Double_t textCB_dR = TMath::Sqrt(pow(m_tCB_eta - m_toff_exteta,2) + pow(m_tCB_phi - m_toff_extphi,2));
                              pextCB_dR = TMath::Sqrt(pow(pCB_eta - m_poff_exteta,2) + pow(pCB_phi - m_poff_extphi,2));
                              Double_t resCB_pt = std::fabs(m_poff_pt)/std::fabs(pCB_pt) - 1.0;
                              m_h_pCB_pt.at(i)->Fill(std::fabs(pCB_pt*0.001));
                              m_h_pCB_dR.at(i)->Fill(pCB_dR);
                              m_h_textCB_dR.at(i)->Fill(textCB_dR);
                              m_h_pextCB_dR.at(i)->Fill(pextCB_dR);
                              m_h_eCB_pt.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              if(std::fabs(m_poff_pt*0.001) > 40)m_h_eCB_eta.at(i)->Fill(m_poff_eta);
                              m_h_pCB_respt.at(i)->Fill(resCB_pt);
                              if(Dicision_barrel(m_poff_eta)){
                                   m_h_eCB_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              }else{
                                   m_h_eCB_pt_end.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              }

                              //EF
                              if(Cut_EF(pEF_pass)){
                                   m_countEF++;
                                   Double_t textEF_dR = TMath::Sqrt(pow(m_tEF_eta - m_toff_exteta,2) + pow(m_tEF_phi - m_toff_extphi,2));
                                   pextEF_dR = TMath::Sqrt(pow(pEF_eta - m_poff_exteta,2) + pow(pEF_phi - m_poff_extphi,2));
                                   Double_t resEF_pt = std::fabs(m_poff_pt)/std::fabs(pEF_pt) - 1.0;
                                   m_h_pEF_pt.at(i)->Fill(std::fabs(pEF_pt*0.001));
                                   m_h_pEF_dR.at(i)->Fill(pEF_dR);
                                   m_h_textEF_dR.at(i)->Fill(textEF_dR);
                                   m_h_pextEF_dR.at(i)->Fill(pextEF_dR);
                                   m_h_eEF_pt.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                   if(std::fabs(m_poff_pt*0.001) > 40)m_h_eEF_eta.at(i)->Fill(m_poff_eta);
                                   m_h_pEF_respt.at(i)->Fill(resEF_pt);
                                   if(Dicision_barrel(m_poff_eta)){
                                        m_h_eEF_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                   }else{
                                        m_h_eEF_pt_end.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                   }
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
     cout<<"ptSAth   nLarge   nLargeS   nSmall   nSmallS"<<endl;
     //SetCondition
     //title,file title,yoffset,top margin,bottom margin,left margin,right margin
     ceff.SetCondition("test",1.5,0,0,0,0);
     ceff.DrawHist1D(m_h_offphi_LargeSpecial);
     ceff.SetCondition("test",1.5,0,0,0,0);
     ceff.DrawHist1D(m_h_saphims_LargeSpecial);
     ceff.SetCondition("test",1.5,0,0,0,0);
     ceff.DrawHist1D(m_h_saroiphi_LargeSpecial);
     m_h_saroiphi_SmallSpecial->Write();
     for(Int_t i = 0;i <= m_nhist;i++){
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_poff_pt.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_pL1_pt.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_pSA_pt.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_pCB_pt.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_pEF_pt.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_pL1_dR.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_pSA_dR.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_pCB_dR.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_pEF_dR.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_textL1_dR.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_textSA_dR.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_textCB_dR.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_textEF_dR.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_pextL1_dR.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_pextSA_dR.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_pextCB_dR.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_pextEF_dR.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_pSA_respt.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_pCB_respt.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_pEF_respt.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_SA_respt0.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_SA_respt1.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_SA_respt2.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_SA_respt3.at(i));
          ceff.SetCondition("trigger;offline eta;count",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_eL1_eta.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_eSA_eta.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_eCB_eta.at(i));
          ceff.SetCondition("test",1.5,0,0,0,0);
          ceff.DrawHist1D(m_h_eEF_eta.at(i));

          ceff.SetCondition("test",1.5,0.1,0.1,0.105,0.165);
          ceff.DrawHist2D(m_h_eff_poff_etaphi.at(i));
          ceff.SetCondition("test",1.5,0.1,0.1,0.105,0.165);
          ceff.DrawHist2D(m_h_eff_pL1_etaphi.at(i));
          ceff.SetCondition("test",1.5,0.1,0.1,0.105,0.165);
          ceff.DrawHist2D(m_h_eff_pSA_etaphi.at(i));
          ceff.SetCondition("test",1.5,0.1,0.1,0.105,0.165);
          ceff.DrawHist2D(m_h_poffvsSA_pt.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptplus0.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptplus1.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptplus2.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptplus3.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptminus0.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptminus1.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptminus2.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptminus3.at(i));
          ceff.DrawHist2D(m_h_offphivsSA_sAddress.at(i));
          ceff.DrawHist2D(m_h_offphivsSA_respt0.at(i));
          ceff.DrawHist2D(m_h_offphivsSA_respt1.at(i));
          ceff.DrawHist2D(m_h_offphivsSA_respt2.at(i));
          ceff.DrawHist2D(m_h_offphivsSA_respt3.at(i));
          ceff.DrawHist2D(m_h_offphivsSAphims.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptplusLS11.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptplusLS15.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptplusLSplus11.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptplusLSplus15.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptplusLSminus11.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptplusLSminus15.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptminusLS11.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptminusLS15.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptminusLSplus11.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptminusLSplus15.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptminusLSminus11.at(i));
          ceff.DrawHist2D(m_h_off_ptvsSA_resptminusLSminus15.at(i));

          //base,target
          ceff.SetConditionName(Form("L1Efficiency_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt.at(i),m_h_eL1_pt.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiency_%dGeV",i*m_thpitch));
          ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt.at(i),m_h_eSA_pt.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("CBEfficiency_%dGeV",i*m_thpitch));
          ceff.SetCondition("muComb Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eSA_pt.at(i),m_h_eCB_pt.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("EFEfficiency_%dGeV",i*m_thpitch));
          ceff.SetCondition("EventFilter Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eCB_pt.at(i),m_h_eEF_pt.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1Efficiency_eta_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyeta(m_h_eoff_eta.at(i),m_h_eL1_eta.at(i));
          ceff.SetConditionName(Form("SAEfficiency_eta_%dGeV",i*m_thpitch));
          ceff.SetCondition("L2MuonSA Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyeta(m_h_eL1_eta.at(i),m_h_eSA_eta.at(i));
          ceff.SetConditionName(Form("CBEfficiency_eta_%dGeV",i*m_thpitch));
          ceff.SetCondition("muComb Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyeta(m_h_eSA_eta.at(i),m_h_eCB_eta.at(i));
          ceff.SetConditionName(Form("EFEfficiency_eta_%dGeV",i*m_thpitch));
          ceff.SetCondition("EventFilter Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyeta(m_h_eCB_eta.at(i),m_h_eEF_eta.at(i));
          ceff.SetConditionName(Form("L1Efficiency_barrel_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_barrel.at(i),m_h_eL1_pt_barrel.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiency_barrel_%dGeV",i*m_thpitch));
          ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_barrel.at(i),m_h_eSA_pt_barrel.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("CBEfficiency_barrel_%dGeV",i*m_thpitch));
          ceff.SetCondition("muComb Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eSA_pt_barrel.at(i),m_h_eCB_pt_barrel.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("EFEfficiency_barrel_%dGeV",i*m_thpitch));
          ceff.SetCondition("EventFilter Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eCB_pt_barrel.at(i),m_h_eEF_pt_barrel.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1Efficiency_end_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_end.at(i),m_h_eL1_pt_end.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiency_end_%dGeV",i*m_thpitch));
          ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_end.at(i),m_h_eSA_pt_end.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("CBEfficiency_end_%dGeV",i*m_thpitch));
          ceff.SetCondition("muComb Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eSA_pt_end.at(i),m_h_eCB_pt_end.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("EFEfficiency_end_%dGeV",i*m_thpitch));
          ceff.SetCondition("EventFilter Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eCB_pt_end.at(i),m_h_eEF_pt_end.at(i),m_binmax,300,m_efficiency_xerr);

          ceff.SetConditionName(Form("L1EfficiencyLarge_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 Large Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_Large.at(i),m_h_eL1_pt_Large.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLarge_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA Large Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_Large.at(i),m_h_eSA_pt_Large.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1EfficiencyLargeSpecial_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_LargeSpecial.at(i),m_h_eL1_pt_LargeSpecial.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecial_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecial.at(i),m_h_eSA_pt_LargeSpecial.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecial11_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecial11.at(i),m_h_eSA_pt_LargeSpecial11.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecial15_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecial15.at(i),m_h_eSA_pt_LargeSpecial15.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialplus11_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus11.at(i),m_h_eSA_pt_LargeSpecialplus11.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialplus15_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus15.at(i),m_h_eSA_pt_LargeSpecialplus15.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialminus11_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus11.at(i),m_h_eSA_pt_LargeSpecialminus11.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialminus15_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus15.at(i),m_h_eSA_pt_LargeSpecialminus15.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1EfficiencySmall_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 Small Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_Small.at(i),m_h_eL1_pt_Small.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencySmall_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA Small Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_Small.at(i),m_h_eSA_pt_Small.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1EfficiencySmallSpecial_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 SmallSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_SmallSpecial.at(i),m_h_eL1_pt_SmallSpecial.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencySmallSpecial_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA SmallSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_SmallSpecial.at(i),m_h_eSA_pt_SmallSpecial.at(i),m_binmax,300,m_efficiency_xerr);

          ceff.SetConditionName(Form("SA2DEfficiency_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1vsL2MuonSA Efficiency;offline eta;offline phi",1.5,0.1,0.1,0.105,0.165);
          ceff.SetConditionbin(m_nbin_eta,m_nbin_phi,m_eta_max,m_phi_max);
          ceff.DrawEfficiency2D(m_h_eff_pL1_etaphi.at(i),m_h_eff_pSA_etaphi.at(i));
          ceff.SetConditionName(Form("L12DEfficiency_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 Efficiency;offline eta;offline phi",1.5,0.1,0.1,0.105,0.165);
          ceff.SetConditionbin(m_nbin_eta,m_nbin_phi,m_eta_max,m_phi_max);
          ceff.DrawEfficiency2D(m_h_eff_poff_etaphi.at(i),m_h_eff_pL1_etaphi.at(i));


	  m_h_eoff_pt.at(i)->Write();
	  m_h_eL1_pt.at(i)->Write();
	  m_h_eSA_pt.at(i)->Write();
	  m_h_eCB_pt.at(i)->Write();
	  m_h_eEF_pt.at(i)->Write();
	  m_h_eoff_eta.at(i)->Write();
	  m_h_eoff_pt_barrel.at(i)->Write();
	  m_h_eL1_pt_barrel.at(i)->Write();
	  m_h_eSA_pt_barrel.at(i)->Write();
	  m_h_eCB_pt_barrel.at(i)->Write();
	  m_h_eEF_pt_barrel.at(i)->Write();
	  m_h_eoff_pt_end.at(i)->Write();
	  m_h_eL1_pt_end.at(i)->Write();
	  m_h_eSA_pt_end.at(i)->Write();
	  m_h_eCB_pt_end.at(i)->Write();
	  m_h_eEF_pt_end.at(i)->Write();
	  m_h_eoff_pt_LargeSpecial.at(i)->Write();
	  m_h_eL1_pt_LargeSpecial.at(i)->Write();
	  m_h_eSA_pt_LargeSpecial.at(i)->Write();
	  m_h_eL1_pt_LargeSpecial11.at(i)->Write();
	  m_h_eSA_pt_LargeSpecial11.at(i)->Write();
	  m_h_eL1_pt_LargeSpecial15.at(i)->Write();
	  m_h_eSA_pt_LargeSpecial15.at(i)->Write();
	  m_h_eL1_pt_LargeSpecialplus11.at(i)->Write();
	  m_h_eSA_pt_LargeSpecialplus11.at(i)->Write();
	  m_h_eL1_pt_LargeSpecialplus15.at(i)->Write();
	  m_h_eSA_pt_LargeSpecialplus15.at(i)->Write();
	  m_h_eL1_pt_LargeSpecialminus11.at(i)->Write();
	  m_h_eSA_pt_LargeSpecialminus11.at(i)->Write();
	  m_h_eL1_pt_LargeSpecialminus15.at(i)->Write();
	  m_h_eSA_pt_LargeSpecialminus15.at(i)->Write();
	  m_h_eoff_pt_Small.at(i)->Write();
	  m_h_eL1_pt_Small.at(i)->Write();
	  m_h_eSA_pt_Small.at(i)->Write();
	  m_h_eoff_pt_SmallSpecial.at(i)->Write();
	  m_h_eL1_pt_SmallSpecial.at(i)->Write();
	  m_h_eSA_pt_SmallSpecial.at(i)->Write();
       m_h_SA_resptplus11.at(i)->Write();
       m_h_SA_resptminus11.at(i)->Write();
       m_h_SA_resptplus15.at(i)->Write();
       m_h_SA_resptminus15.at(i)->Write();
       m_h_offphivsSAphims.at(i)->Write();
       m_h_offphivsSA_respt.at(i)->Write();
       m_h_offetavsSA_respt0.at(i)->Write();
       m_h_offetavsSA_respt1.at(i)->Write();
       m_h_offetavsSA_respt2.at(i)->Write();
       m_h_offetavsSA_respt3.at(i)->Write();
       m_h_offetavsSA_resptLSplus11.at(i)->Write();
       m_h_offetavsSA_resptLSminus11.at(i)->Write();
       m_h_offetavsSA_resptLSplus15.at(i)->Write();
       m_h_offetavsSA_resptLSminus15.at(i)->Write();
       m_h_highoffetavsSA_resptLSplus11.at(i)->Write();
       m_h_highoffetavsSA_resptLSminus11.at(i)->Write();
       m_h_highoffetavsSA_resptLSplus15.at(i)->Write();
       m_h_highoffetavsSA_resptLSminus15.at(i)->Write();
       m_h_highoffphivsSA_resptLSplus11.at(i)->Write();
       m_h_highoffphivsSA_resptLSminus11.at(i)->Write();
       m_h_highoffphivsSA_resptLSplus15.at(i)->Write();
       m_h_highoffphivsSA_resptLSminus15.at(i)->Write();


         cout<<i*m_thpitch<<"      "<<m_countLarge.at(i)<<"      "<<m_countLargeSpecial.at(i)<<"      "<<m_countSmall.at(i)<<"      "<<m_countSmallSpecial.at(i)<<endl;
     }
     
     cout<<m_count<<"   "<<m_countall<<"   "<<m_countoff<<"   "<<m_countL1<<"   "<<m_countSA<<"   "<<m_countCB<<"   "<<m_countEF<<endl;
    
     delete m_h_offphi_LargeSpecial;
     delete m_h_saphims_LargeSpecial;
     delete m_h_saroiphi_LargeSpecial;

     m_h_poff_pt.clear();
     m_h_pL1_pt.clear();
     m_h_pSA_pt.clear();
     m_h_pCB_pt.clear();
     m_h_pEF_pt.clear();
     m_h_pL1_dR.clear();
     m_h_pSA_dR.clear();
     m_h_pCB_dR.clear();
     m_h_pEF_dR.clear();
     m_h_textL1_dR.clear();
     m_h_textSA_dR.clear();
     m_h_textCB_dR.clear();
     m_h_textEF_dR.clear();
     m_h_pextL1_dR.clear();
     m_h_pextSA_dR.clear();
     m_h_pextCB_dR.clear();
     m_h_pextEF_dR.clear();
     m_h_pSA_respt.clear();
     m_h_pCB_respt.clear();
     m_h_pEF_respt.clear();
     m_h_eoff_pt.clear();
     m_h_eL1_pt.clear();
     m_h_eSA_pt.clear();
     m_h_eCB_pt.clear();
     m_h_eEF_pt.clear();
     m_h_eoff_eta.clear();
     m_h_eL1_eta.clear();
     m_h_eSA_eta.clear();
     m_h_eCB_eta.clear();
     m_h_eEF_eta.clear();
     m_h_eoff_pt_barrel.clear();
     m_h_eL1_pt_barrel.clear();
     m_h_eSA_pt_barrel.clear();
     m_h_eCB_pt_barrel.clear();
     m_h_eEF_pt_barrel.clear();
     m_h_eoff_pt_end.clear();
     m_h_eL1_pt_end.clear();
     m_h_eSA_pt_end.clear();
     m_h_eCB_pt_end.clear();
     m_h_eEF_pt_end.clear();
     m_h_SA_respt0.clear();
     m_h_SA_respt1.clear();
     m_h_SA_respt2.clear();
     m_h_SA_respt3.clear();
     m_h_eff_poff_etaphi.clear();
     m_h_eff_pL1_etaphi.clear();
     m_h_eff_pSA_etaphi.clear();
     m_h_poffvsSA_pt.clear();
     m_h_off_ptvsSA_resptplus0.clear();
     m_h_off_ptvsSA_resptplus1.clear();
     m_h_off_ptvsSA_resptplus2.clear();
     m_h_off_ptvsSA_resptplus3.clear();
     m_h_off_ptvsSA_resptminus0.clear();
     m_h_off_ptvsSA_resptminus1.clear();
     m_h_off_ptvsSA_resptminus2.clear();
     m_h_off_ptvsSA_resptminus3.clear();
     m_h_offphivsSA_sAddress.clear();
     m_h_offphivsSA_respt0.clear();
     m_h_offphivsSA_respt1.clear();
     m_h_offphivsSA_respt2.clear();
     m_h_offphivsSA_respt3.clear();
     m_h_eoff_pt_Large.clear();
     m_h_eL1_pt_Large.clear();
     m_h_eSA_pt_Large.clear();
     m_h_eoff_pt_LargeSpecial.clear();
     m_h_eL1_pt_LargeSpecial.clear();
     m_h_eSA_pt_LargeSpecial.clear();
     m_h_eoff_pt_Small.clear();
     m_h_eL1_pt_Small.clear();
     m_h_eSA_pt_Small.clear();
     m_h_eoff_pt_SmallSpecial.clear();
     m_h_eL1_pt_SmallSpecial.clear();
     m_h_eSA_pt_SmallSpecial.clear();
     m_h_offphivsSAphims.clear();
     m_h_off_ptvsSA_resptplusLS11.clear();
     m_h_off_ptvsSA_resptplusLS15.clear();
     m_h_off_ptvsSA_resptplusLSplus11.clear();
     m_h_off_ptvsSA_resptplusLSplus15.clear();
     m_h_off_ptvsSA_resptplusLSminus11.clear();
     m_h_off_ptvsSA_resptplusLSminus15.clear();
     m_h_off_ptvsSA_resptminusLS11.clear();
     m_h_off_ptvsSA_resptminusLS15.clear();
     m_h_off_ptvsSA_resptminusLSplus11.clear();
     m_h_off_ptvsSA_resptminusLSplus15.clear();
     m_h_off_ptvsSA_resptminusLSminus11.clear();
     m_h_off_ptvsSA_resptminusLSminus15.clear();
     m_h_eL1_pt_LargeSpecial11.clear();
     m_h_eSA_pt_LargeSpecial11.clear();
     m_h_eL1_pt_LargeSpecial15.clear();
     m_h_eSA_pt_LargeSpecial15.clear();
     m_h_eL1_pt_LargeSpecialplus11.clear();
     m_h_eSA_pt_LargeSpecialplus11.clear();
     m_h_eL1_pt_LargeSpecialplus15.clear();
     m_h_eSA_pt_LargeSpecialplus15.clear();
     m_h_eL1_pt_LargeSpecialminus11.clear();
     m_h_eSA_pt_LargeSpecialminus11.clear();
     m_h_eL1_pt_LargeSpecialminus15.clear();
     m_h_eSA_pt_LargeSpecialminus15.clear();
     m_h_offphivsSA_respt.clear();

     m_mes_name->clear();
     m_pEFTAG_pass->clear();
     m_pL1_pt->clear();
     m_pL1_eta->clear();
     m_pL1_phi->clear();
     m_pL1_pass->clear();
     m_pL1_dR->clear();
     m_pSA_pt->clear();
     m_pSA_eta->clear();
     m_pSA_phi->clear();
     m_pSA_phims->clear();
     m_pSA_roiphi->clear();
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
     m_countLarge.clear();
     m_countLargeSpecial.clear();
     m_countSmall.clear();
     m_countSmallSpecial.clear();
}
