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

void Efficiency::Init(TTree *tree,std::string name,const Int_t np,const Int_t ne,const Double_t mp,const Double_t me,Double_t req,Int_t max,Double_t err,const Int_t nh,const Int_t th,Int_t proc){
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
          m_proc = proc;
    
     	//initialize
          m_tag_proc = 0;
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
     	m_poff_charge = 0;
     	m_poff_d0 = 0;
     	m_poff_z0 = 0;
          m_tp_dR = 0;
          m_tp_extdR = 0;
     	m_tL1_pt = 0;
     	m_tL1_eta = 0;
     	m_tL1_phi = 0;
          m_reqL1dR = req;
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
          m_sumReqdRL1 = 0;
     	m_pSA_pt = 0;
     	m_pSA_eta = 0;
     	m_pSA_phi = 0;
     	m_pSA_dR = 0;
     	m_pSA_pass = 0;
     	m_pSA_sAddress = 0;
          m_pSA_rpcX = 0;
          m_pSA_rpcY = 0;
          m_pSA_mdtZ = 0;
          m_pSA_mdtR = 0;
          m_pSA_mdtPhi = 0;
          m_pSA_phims = 0;
          m_pSA_roiphi = 0;
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
     	m_sumReqdREF = 0;

     	//active only need branch 
     	tChain->SetBranchStatus("*",0);
     	tChain->SetBranchStatus("mes_name",1);
          tChain->SetBranchStatus("tag_proc",1);
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
          tChain->SetBranchStatus("probe_segment_etaIndex",1);
          tChain->SetBranchStatus("tp_dR",1);
          tChain->SetBranchStatus("tp_extdR",1);
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
          tChain->SetBranchStatus("sumReqdRL1",1);
     	tChain->SetBranchStatus("probe_mesSA_pt",1);
     	tChain->SetBranchStatus("probe_mesSA_eta",1);
     	tChain->SetBranchStatus("probe_mesSA_phi",1);
     	tChain->SetBranchStatus("probe_mesSA_pass",1);
     	tChain->SetBranchStatus("probe_mesSA_dR",1);
     	tChain->SetBranchStatus("probe_mesSA_sAddress",1);
     	tChain->SetBranchStatus("probe_mesSA_phims",1);
     	tChain->SetBranchStatus("probe_mesSA_roiPhi",1);
          tChain->SetBranchStatus("probe_mesSA_rpcHitX",1);
          tChain->SetBranchStatus("probe_mesSA_rpcHitY",1);
          tChain->SetBranchStatus("probe_mesSA_mdtHitZ",1);
          tChain->SetBranchStatus("probe_mesSA_mdtHitR",1);
          tChain->SetBranchStatus("probe_mesSA_mdtHitPhi",1);
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
          tChain->SetBranchStatus("sumReqdREF",1);
     	tChain->SetBranchStatus("probe_mesEFTAG_pass",1);

     	//setting each branch address
    	     tChain->SetBranchAddress("mes_name",&m_mes_name,&b_mes_name);
          tChain->SetBranchAddress("tag_proc",&m_tag_proc,&b_tag_proc);
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
     	tChain->SetBranchAddress("probe_charge",&m_poff_charge,&b_probe_charge);
     	tChain->SetBranchAddress("probe_d0",&m_poff_d0,&b_probe_d0);
     	tChain->SetBranchAddress("probe_z0",&m_poff_z0,&b_probe_z0);
          tChain->SetBranchAddress("probe_segment_etaIndex",m_probe_segment_etaIndex,&b_probe_segment_etaIndex);
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
          tChain->SetBranchAddress("sumReqdRL1",&m_sumReqdRL1,&b_sumReqdRL1);
     	tChain->SetBranchAddress("probe_mesSA_pt",&m_pSA_pt,&b_pSA_pt);
     	tChain->SetBranchAddress("probe_mesSA_eta",&m_pSA_eta,&b_pSA_eta);
     	tChain->SetBranchAddress("probe_mesSA_phi",&m_pSA_phi,&b_pSA_phi);
     	tChain->SetBranchAddress("probe_mesSA_pass",&m_pSA_pass,&b_pSA_pass);
     	tChain->SetBranchAddress("probe_mesSA_dR",&m_pSA_dR,&b_pSA_dR);
     	tChain->SetBranchAddress("probe_mesSA_sAddress",&m_pSA_sAddress,&b_pSA_sAddress);
          tChain->SetBranchAddress("probe_mesSA_rpcHitX",&m_pSA_rpcX,&b_pSA_rpcX);
          tChain->SetBranchAddress("probe_mesSA_rpcHitY",&m_pSA_rpcY,&b_pSA_rpcY);
          tChain->SetBranchAddress("probe_mesSA_mdtHitZ",&m_pSA_mdtZ,&b_pSA_mdtZ);
          tChain->SetBranchAddress("probe_mesSA_mdtHitR",&m_pSA_mdtR,&b_pSA_mdtR);
          tChain->SetBranchAddress("probe_mesSA_mdtHitPhi",&m_pSA_mdtPhi,&b_pSA_mdtPhi);
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
          tChain->SetBranchAddress("sumReqdREF",&m_sumReqdREF,&b_sumReqdREF);
     	tChain->SetBranchAddress("probe_mesEFTAG_pass",&m_pEFTAG_pass,&b_pEFTAG_pass);

     	//define each histgram
     	m_h_offphi_LargeSpecial = new TH1D("h_offphi_LargeSpecial_0GeV","offline phi;offline phi;Entries",600,-3.0,3.0);
     	m_h_saphims_LargeSpecial = new TH1D("h_saphims_LargeSpecial_0GeV","L2MuonSA phi;L2MuonSA phims;Entries",600,-3.0,3.0);
     	m_h_saroiphi_LargeSpecial = new TH1D("h_saroiphi_LargeSpecial_0GeV","RoI phi;RoI phi;Entries",600,-3.0,3.0);
          m_h_saroiphi_SmallSpecial = new TH1D("h_saroiphi_SmallSpecial_0GeV","RoI phi;RoI phi;Entries",600,-3.0,3.0);
     	for(Int_t i = 0;i <= m_nhist;i++){
               //Standard
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
               m_h_poffvsSA_pt.push_back(new TH2F(Form("h_poffvsSA_pt_%dGeV",i*m_thpitch),"probe offline pt vs probe L2MuonSA pt@mu26ivm;probe offline pt[GeV];probe L2MuonSA pt[GeV]",150,0,150,150,0,150));
               m_h_offphivsSA_sAddress.push_back(new TH2F(Form("h_offphivsSA_sAddress_%dGeV",i*m_thpitch),"offline phi vs sAddress;offline phi;sAddress",140,-3.5,3.5,4,0.0,4.0));
               m_h_offphivsSAphims.push_back(new TH2F(Form("h_offphivsSAphims_%dGeV",i*m_thpitch),"offphi vs phims;offline phi;L2MuonSA phims",140,-3.5,3.5,140,-3.5,3.5));
               m_h_rpchitXY.push_back(new TH2F(Form("h_rpchitXvsrpchitY_%dGeV",i*m_thpitch),"RPC hit distribution;RPC hit X[cm];RPC hit Y[cm]",2000,-10000.0,10000.0,2000,-10000.0,10000.0));
               m_h_rpchitXYLargeSpecial.push_back(new TH2F(Form("h_rpchitXvsrpchitYLS_%dGeV",i*m_thpitch),"RPCLS hit distribution;RPC hit X[cm];RPC hit Y[cm]",2000,-10000.0,10000.0,2000,-10000.0,10000.0));
               m_h_rpchitXYLargeSpecialplus11out.push_back(new TH2F(Form("h_rpchitXvsrpchitYLargeSpecialplus11out_%dGeV",i*m_thpitch),"RPCLS hit distribution;RPC hit X[cm];RPC hit Y[cm]",2000,-10000.0,10000.0,2000,-10000.0,10000.0));
               m_h_rpchitXYLargeSpecialplus11in.push_back(new TH2F(Form("h_rpchitXvsrpchitYLargeSpecialplus11in_%dGeV",i*m_thpitch),"RPCLS hit distribution;RPC hit X[cm];RPC hit Y[cm]",2000,-10000.0,10000.0,2000,-10000.0,10000.0));
               m_h_rpchitXYLargeSpecialplus15out.push_back(new TH2F(Form("h_rpchitXvsrpchitYLargeSpecialplus15out_%dGeV",i*m_thpitch),"RPCLS hit distribution;RPC hit X[cm];RPC hit Y[cm]",2000,-10000.0,10000.0,2000,-10000.0,10000.0));
               m_h_rpchitXYLargeSpecialplus15in.push_back(new TH2F(Form("h_rpchitXvsrpchitYLargeSpecialplus15in_%dGeV",i*m_thpitch),"RPCLS hit distribution;RPC hit X[cm];RPC hit Y[cm]",2000,-10000.0,10000.0,2000,-10000.0,10000.0));
               m_h_rpchitXYLargeSpecialminus11out.push_back(new TH2F(Form("h_rpchitXvsrpchitYLargeSpecialminus11out_%dGeV",i*m_thpitch),"RPCLS hit distribution;RPC hit X[cm];RPC hit Y[cm]",2000,-10000.0,10000.0,2000,-10000.0,10000.0));
               m_h_rpchitXYLargeSpecialminus11in.push_back(new TH2F(Form("h_rpchitXvsrpchitYLargeSpecialminus11in_%dGeV",i*m_thpitch),"RPCLS hit distribution;RPC hit X[cm];RPC hit Y[cm]",2000,-10000.0,10000.0,2000,-10000.0,10000.0));
               m_h_rpchitXYLargeSpecialminus15out.push_back(new TH2F(Form("h_rpchitXvsrpchitYLargeSpecialminus15out_%dGeV",i*m_thpitch),"RPCLS hit distribution;RPC hit X[cm];RPC hit Y[cm]",2000,-10000.0,10000.0,2000,-10000.0,10000.0));
               m_h_rpchitXYLargeSpecialminus15in.push_back(new TH2F(Form("h_rpchitXvsrpchitYLargeSpecialminus15in_%dGeV",i*m_thpitch),"RPCLS hit distribution;RPC hit X[cm];RPC hit Y[cm]",2000,-10000.0,10000.0,2000,-10000.0,10000.0));
               m_h_mdthitXY.push_back(new TH2F(Form("h_mdthitXvsmdthitY_%dGeV",i*m_thpitch),"MDT hit distribution;MDT hit X[cm];MDT hit Y[cm]",5000,-25000.0,25000.0,5000,-25000.0,25000.0));
               m_h_mdthitXYLargeSpecial.push_back(new TH2F(Form("h_mdthitXvsmdthitYLargeSpecial_%dGeV",i*m_thpitch),"MDT hit distribution;MDT hit X[cm];MDT hit Y[cm]",5000,-25000.0,25000.0,5000,-25000.0,25000.0));
               m_h_mdthitXYLargeSpecialplus11out.push_back(new TH2F(Form("h_mdthitXvsmdthitYLargeSpecialplus11out_%dGeV",i*m_thpitch),"MDT hit distribution;MDT hit X[cm];MDT hit Y[cm]",5000,-25000.0,25000.0,5000,-25000.0,25000.0));
               m_h_mdthitXYLargeSpecialplus11in.push_back(new TH2F(Form("h_mdthitXvsmdthitYLargeSpecialplus11in_%dGeV",i*m_thpitch),"MDT hit distribution;MDT hit X[cm];MDT hit Y[cm]",5000,-25000.0,10000.0,5000,-25000.0,25000.0));
               m_h_mdthitXYLargeSpecialplus15out.push_back(new TH2F(Form("h_mdthitXvsmdthitYLargeSpecialplus15out_%dGeV",i*m_thpitch),"MDT hit distribution;MDT hit X[cm];MDT hit Y[cm]",5000,-25000.0,10000.0,5000,-25000.0,25000.0));
               m_h_mdthitXYLargeSpecialplus15in.push_back(new TH2F(Form("h_mdthitXvsmdthitYLargeSpecialplus15in_%dGeV",i*m_thpitch),"MDT hit distribution;MDT hit X[cm];MDT hit Y[cm]",5000,-25000.0,10000.0,5000,-25000.0,25000.0));
               m_h_mdthitXYLargeSpecialminus11out.push_back(new TH2F(Form("h_mdthitXvsmdthitYLargeSpecialminus11out_%dGeV",i*m_thpitch),"MDT hit distribution;MDT hit X[cm];MDT hit Y[cm]",5000,-25000.0,25000.0,5000,-25000.0,25000.0));
               m_h_mdthitXYLargeSpecialminus11in.push_back(new TH2F(Form("h_mdthitXvsmdthitYLargeSpecialminus11in_%dGeV",i*m_thpitch),"MDT hit distribution;MDT hit X[cm];MDT hit Y[cm]",5000,-25000.0,10000.0,5000,-25000.0,25000.0));
               m_h_mdthitXYLargeSpecialminus15out.push_back(new TH2F(Form("h_mdthitXvsmdthitYLargeSpecialminus15out_%dGeV",i*m_thpitch),"MDT hit distribution;MDT hit X[cm];MDT hit Y[cm]",5000,-25000.0,10000.0,5000,-25000.0,25000.0));
               m_h_mdthitXYLargeSpecialminus15in.push_back(new TH2F(Form("h_mdthitXvsmdthitYLargeSpecialminus15in_%dGeV",i*m_thpitch),"MDT hit distribution;MDT hit X[cm];MDT hit Y[cm]",5000,-25000.0,10000.0,5000,-25000.0,25000.0));

               //Efficiency
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
               m_h_eoff_pt_Largeplus.push_back(new TH1D(Form("h_eoff_ptLargeplus_%dGeV",i*m_thpitch),"mesoff_pt;offline pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_Largeplus.push_back(new TH1D(Form("h_eL1_ptLargeplus_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_Largeplus.push_back(new TH1D(Form("h_eSA_ptLargeplus_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
               m_h_eoff_pt_Largeminus.push_back(new TH1D(Form("h_eoff_ptLargeminus_%dGeV",i*m_thpitch),"mesoff_pt;offline pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_Largeminus.push_back(new TH1D(Form("h_eL1_ptLargeminus_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_Largeminus.push_back(new TH1D(Form("h_eSA_ptLargeminus_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
               m_h_eoff_pt_LargeSpecialplus.push_back(new TH1D(Form("h_eoff_ptLargeSpecialplus_%dGeV",i*m_thpitch),"mesoff_pt;offline pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_LargeSpecialplus.push_back(new TH1D(Form("h_eL1_ptLargeSpecialplus_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_LargeSpecialplus.push_back(new TH1D(Form("h_eSA_ptLargeSpecialplus_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_LargeSpecialplus11.push_back(new TH1D(Form("h_eL1_ptLargeSpecial11_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_LargeSpecialplus11.push_back(new TH1D(Form("h_eSA_ptLargeSpecial11_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_LargeSpecialplus15.push_back(new TH1D(Form("h_eL1_ptLargeSpecial15_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_LargeSpecialplus15.push_back(new TH1D(Form("h_eSA_ptLargeSpecial15_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_LargeSpecialplus11in.push_back(new TH1D(Form("h_eL1_ptLargeSpecialplus11in_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_LargeSpecialplus11out.push_back(new TH1D(Form("h_eL1_ptLargeSpecialplus11out_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_LargeSpecialplus11in.push_back(new TH1D(Form("h_eSA_ptLargeSpecialplus11in_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_LargeSpecialplus11out.push_back(new TH1D(Form("h_eSA_ptLargeSpecialplus11out_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_LargeSpecialplus15in.push_back(new TH1D(Form("h_eL1_ptLargeSpecialplus15in_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_LargeSpecialplus15out.push_back(new TH1D(Form("h_eL1_ptLargeSpecialplus15out_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_LargeSpecialplus15in.push_back(new TH1D(Form("h_eSA_ptLargeSpecialplus15in_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_LargeSpecialplus15out.push_back(new TH1D(Form("h_eSA_ptLargeSpecialplus15out_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
               m_h_eoff_pt_LargeSpecialminus.push_back(new TH1D(Form("h_eoff_ptLargeSpecialminus_%dGeV",i*m_thpitch),"mesoff_pt;offline pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_LargeSpecialminus.push_back(new TH1D(Form("h_eL1_ptLargeSpecialminus_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_LargeSpecialminus.push_back(new TH1D(Form("h_eSA_ptLargeSpecialminus_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_LargeSpecialminus11.push_back(new TH1D(Form("h_eL1_ptLargeSpecialminus11_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_LargeSpecialminus11.push_back(new TH1D(Form("h_eSA_ptLargeSpecialminus11_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_LargeSpecialminus15.push_back(new TH1D(Form("h_eL1_ptLargeSpecialminus15_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_LargeSpecialminus15.push_back(new TH1D(Form("h_eSA_ptLargeSpecialminus15_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_LargeSpecialminus11in.push_back(new TH1D(Form("h_eL1_ptLargeSpecialminus11in_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_LargeSpecialminus11in.push_back(new TH1D(Form("h_eSA_ptLargeSpecialminus11in_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_LargeSpecialminus15in.push_back(new TH1D(Form("h_eL1_ptLargeSpecialminus15in_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_LargeSpecialminus15in.push_back(new TH1D(Form("h_eSA_ptLargeSpecialminus15in_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_LargeSpecialminus11out.push_back(new TH1D(Form("h_eL1_ptLargeSpecialminus11out_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_LargeSpecialminus11out.push_back(new TH1D(Form("h_eSA_ptLargeSpecialminus11out_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_LargeSpecialminus15out.push_back(new TH1D(Form("h_eL1_ptLargeSpecialminus15out_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_LargeSpecialminus15out.push_back(new TH1D(Form("h_eSA_ptLargeSpecialminus15out_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));

               m_h_eoff_pt_Smallplus.push_back(new TH1D(Form("h_eoff_ptSmallplus_%dGeV",i*m_thpitch),"mesoff_pt;offline pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_Smallplus.push_back(new TH1D(Form("h_eL1_ptSmallplus_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_Smallplus.push_back(new TH1D(Form("h_eSA_ptSmallplus_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
               m_h_eoff_pt_Smallminus.push_back(new TH1D(Form("h_eoff_ptSmallminus_%dGeV",i*m_thpitch),"mesoff_pt;offline pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_Smallminus.push_back(new TH1D(Form("h_eL1_ptSmallminus_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_Smallminus.push_back(new TH1D(Form("h_eSA_ptSmallminus_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
               m_h_eoff_pt_SmallSpecialplus.push_back(new TH1D(Form("h_eoff_ptSmallSpecialplus_%dGeV",i*m_thpitch),"mesoff_pt;offline pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_SmallSpecialplus.push_back(new TH1D(Form("h_eL1_ptSmallSpecialplus_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_SmallSpecialplus.push_back(new TH1D(Form("h_eSA_ptSmallSpecialplus_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
               m_h_eoff_pt_SmallSpecialminus.push_back(new TH1D(Form("h_eoff_ptSmallSpecialminus_%dGeV",i*m_thpitch),"mesoff_pt;offline pt[GeV];Entries",300,-0.25,149.75));
               m_h_eL1_pt_SmallSpecialminus.push_back(new TH1D(Form("h_eL1_ptSmallSpecialminus_%dGeV",i*m_thpitch),"mesL1_pt;L1 pt[GeV];Entries",300,-0.25,149.75));
               m_h_eSA_pt_SmallSpecialminus.push_back(new TH1D(Form("h_eSA_ptSmallSpecialminus_%dGeV",i*m_thpitch),"mesSA_pt;SA pt[GeV];Entries",300,-0.25,149.75));
               m_h_eff_poff_etaphi.push_back(new TH2F(Form("h_eff_poff_etaphi_%dGeV",i*m_thpitch),"offlineeta vs offlinephi;offline eta;offline phi",m_nbin_eta,-m_eta_max,m_eta_max,m_nbin_phi,-m_phi_max,m_phi_max));
               m_h_eff_pL1_etaphi.push_back(new TH2F(Form("h_eff_pL1_etaphi_%dGeV",i*m_thpitch),"L1eta vs L1phi;L1 eta;L1 phi",m_nbin_eta,-m_eta_max,m_eta_max,m_nbin_phi,-m_phi_max,m_phi_max));
               m_h_eff_pSA_etaphi.push_back(new TH2F(Form("h_eff_pSA_etaphi_%dGeV",i*m_thpitch),"L2MuonSAeta vs L2MuonSAphi;L2MuonSA eta;L2MuonSA phi",m_nbin_eta,-m_eta_max,m_eta_max,m_nbin_phi,-m_phi_max,m_phi_max));

               //Residual
               m_h_pSA_respt.push_back(new TH1D(Form("h_pSA_respt_%dGeV",i*m_thpitch),"probe L2MuonSA residual pt;residual pt;Entries",300,-1,1));
               m_h_pCB_respt.push_back(new TH1D(Form("h_pCB_respt_%dGeV",i*m_thpitch),"probe muComb residual pt;residual pt;Entries",640,-0.8,0.8));
               m_h_pEF_respt.push_back(new TH1D(Form("h_pEF_respt_%dGeV",i*m_thpitch),"probe EventFilter residual pt;residual pt;Entries",1800,-0.3,0.3));
               m_h_SA_resptLargeplus.push_back(new TH1D(Form("h_SA_resptLargeplus_%dGeV",i*m_thpitch),"SAresidualpt_Large;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_SA_resptLargeSpecialplus.push_back(new TH1D(Form("h_SA_resptLargeSpecialplus_%dGeV",i*m_thpitch),"SAresidualpt_LargeSpecial;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_SA_resptSmallplus.push_back(new TH1D(Form("h_SA_resptSmallplus_%dGeV",i*m_thpitch),"SAresidualpt_Small;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_SA_resptSmallSpecialplus.push_back(new TH1D(Form("h_SA_resptSmallSpecialplus_%dGeV",i*m_thpitch),"SAresidualpt_SmallSpecial;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_SA_resptLargeminus.push_back(new TH1D(Form("h_SA_resptLargeminus_%dGeV",i*m_thpitch),"SAresidualpt_Large;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_SA_resptLargeSpecialminus.push_back(new TH1D(Form("h_SA_resptLargeSpecialminus_%dGeV",i*m_thpitch),"SAresidualpt_LargeSpecial;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_SA_resptSmallminus.push_back(new TH1D(Form("h_SA_resptSmallminus_%dGeV",i*m_thpitch),"SAresidualpt_Small;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_SA_resptSmallSpecialminus.push_back(new TH1D(Form("h_SA_resptSmallSpecialminus_%dGeV",i*m_thpitch),"SAresidualpt_SmallSpecial;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_SA_resptLargeSpecialplus11out.push_back(new TH1D(Form("h_SA_resptLargeSpecialplus11out_%dGeV",i*m_thpitch),"SAresidualpt_LS11+;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_SA_resptLargeSpecialplus11in.push_back(new TH1D(Form("h_SA_resptLargeSpecialplus11in_%dGeV",i*m_thpitch),"SAresidualpt_LS11+;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_SA_resptLargeSpecialplus15out.push_back(new TH1D(Form("h_SA_resptLargeSpecialplus15out_%dGeV",i*m_thpitch),"SAresidualpt_LS15+;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_SA_resptLargeSpecialplus15in.push_back(new TH1D(Form("h_SA_resptLargeSpecialplus15in_%dGeV",i*m_thpitch),"SAresidualpt_LS15+;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_SA_resptLargeSpecialminus11out.push_back(new TH1D(Form("h_SA_resptLargeSpecialsminus11out_%dGeV",i*m_thpitch),"SAresidualpt_LS11-;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_SA_resptLargeSpecialminus11in.push_back(new TH1D(Form("h_SA_resptLargeSpecialsminus11in_%dGeV",i*m_thpitch),"SAresidualpt_LS11-;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_SA_resptLargeSpecialminus15out.push_back(new TH1D(Form("h_SA_resptLargeSpecialminus15out_%dGeV",i*m_thpitch),"SAresidualpt_LS15-;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_SA_resptLargeSpecialminus15in.push_back(new TH1D(Form("h_SA_resptLargeSpecialminus15in_%dGeV",i*m_thpitch),"SAresidualpt_LS15-;L2MuonSA residual pt;Entries",1000,-5.0,5.0));
               m_h_off_ptvsSA_resptLargeplus.push_back(new TH2F(Form("h_off_ptvsSA_resptLargeplus_%dGeV",i*m_thpitch),"Large Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
               m_h_off_ptvsSA_resptLargeSpecialplus.push_back(new TH2F(Form("h_off_ptvsSA_resptLargeSpecialplus_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
               m_h_off_ptvsSA_resptSmallplus.push_back(new TH2F(Form("h_off_ptvsSA_resptSmallplus_%dGeV",i*m_thpitch),"Small Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
               m_h_off_ptvsSA_resptSmallSpecialplus.push_back(new TH2F(Form("h_off_ptvsSA_resptSmallSpecialplus_%dGeV",i*m_thpitch),"SmallSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
               m_h_off_ptvsSA_resptLargeminus.push_back(new TH2F(Form("h_off_ptvsSA_resptLargeminus_%dGeV",i*m_thpitch),"Large Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
               m_h_off_ptvsSA_resptLargeSpecialminus.push_back(new TH2F(Form("h_off_ptvsSA_resptLargeSpecialminus_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
               m_h_off_ptvsSA_resptSmallminus.push_back(new TH2F(Form("h_off_ptvsSA_resptSmallminus_%dGeV",i*m_thpitch),"Small Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
               m_h_off_ptvsSA_resptSmallSpecialminus.push_back(new TH2F(Form("h_off_ptvsSA_resptSmallSpecialminus_%dGeV",i*m_thpitch),"SmallSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
               m_h_offphivsSA_resptLargeplus.push_back(new TH2F(Form("h_offphivsSA_resptLargeplus_%dGeV",i*m_thpitch),"Large;offline phi;SApt residual",140,-3.5,3.5,300,-2.0,1.0));
               m_h_offphivsSA_resptLargeSpecialplus.push_back(new TH2F(Form("h_offphivsSA_resptLargeSpecialplus_%dGeV",i*m_thpitch),"LargeSpecial;offline phi;SApt residual",140,-3.5,3.5,300,-2.0,1.0));
               m_h_offphivsSA_resptSmallplus.push_back(new TH2F(Form("h_offphivsSA_resptSmallplus_%dGeV",i*m_thpitch),"Small;offline phi;SApt residual",140,-3.5,3.5,300,-2.0,1.0));
               m_h_offphivsSA_resptSmallSpecialplus.push_back(new TH2F(Form("h_offphivsSA_resptSmallSpecialplus_%dGeV",i*m_thpitch),"SmallSpecial;offline phi;SApt residual",140,-3.5,3.5,300,-2.0,1.0));
               m_h_offphivsSA_resptLargeminus.push_back(new TH2F(Form("h_offphivsSA_resptLargeminus_%dGeV",i*m_thpitch),"Large;offline phi;SApt residual",140,-3.5,3.5,300,-2.0,1.0));
               m_h_offphivsSA_resptLargeSpecialminus.push_back(new TH2F(Form("h_offphivsSA_resptLargeSpecialminus_%dGeV",i*m_thpitch),"LargeSpecial;offline phi;SApt residual",140,-3.5,3.5,300,-2.0,1.0));
               m_h_offphivsSA_resptSmallminus.push_back(new TH2F(Form("h_offphivsSA_resptSmallminus_%dGeV",i*m_thpitch),"Small;offline phi;SApt residual",140,-3.5,3.5,300,-2.0,1.0));
               m_h_offphivsSA_resptSmallSpecialminus.push_back(new TH2F(Form("h_offphivsSA_resptSmallSpecialminus_%dGeV",i*m_thpitch),"SmallSpecial;offline phi;SApt residual",140,-3.5,3.5,300,-2.0,1.0));
               m_h_off_ptvsSA_resptLargeSpecialplus11out.push_back(new TH2F(Form("h_off_ptvsSA_resptLargeSpcialplus11out_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
               m_h_off_ptvsSA_resptLargeSpecialplus15out.push_back(new TH2F(Form("h_off_ptvsSA_resptLargeSpecialplus15out_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
               m_h_off_ptvsSA_resptLargeSpecialplus11in.push_back(new TH2F(Form("h_off_ptvsSA_resptplusLargeSpecialplus11in_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
               m_h_off_ptvsSA_resptLargeSpecialplus15in.push_back(new TH2F(Form("h_off_ptvsSA_resptplusLargeSpecialplus15in_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
               m_h_off_ptvsSA_resptLargeSpecialminus11out.push_back(new TH2F(Form("h_off_ptvsSA_resptLargeSpcialminus11out_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
               m_h_off_ptvsSA_resptLargeSpecialminus15out.push_back(new TH2F(Form("h_off_ptvsSA_resptLargeSpecialminus15out_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
               m_h_off_ptvsSA_resptLargeSpecialminus11in.push_back(new TH2F(Form("h_off_ptvsSA_resptplusLargeSpecialminus11in_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
               m_h_off_ptvsSA_resptLargeSpecialminus15in.push_back(new TH2F(Form("h_off_ptvsSA_resptplusLargeSpecialminus15in_%dGeV",i*m_thpitch),"LargeSpecial Qeta/|eta|=+1;probe offline pt[GeV];SApt residual",35,0,70,40,-2.0,2.0));
               m_h_offetavsSA_resptLargeplus.push_back(new TH2F(Form("h_offetavsSA_resptLargeplus_%dGeV",i*m_thpitch),"Large;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_resptLargeSpecialplus.push_back(new TH2F(Form("h_offetavsSA_resptLargeSpecialplus_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_resptSmallplus.push_back(new TH2F(Form("h_offetavsSA_resptSmallplus_%dGeV",i*m_thpitch),"Small;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_resptSmallSpecialplus.push_back(new TH2F(Form("h_offetavsSA_resptSmallSpecialplus_%dGeV",i*m_thpitch),"SmallSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_resptLargeminus.push_back(new TH2F(Form("h_offetavsSA_resptLargeminus_%dGeV",i*m_thpitch),"Large;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_resptLargeSpecialminus.push_back(new TH2F(Form("h_offetavsSA_resptLargeSpecialminus_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_resptSmallminus.push_back(new TH2F(Form("h_offetavsSA_resptSmallminus_%dGeV",i*m_thpitch),"Small;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_resptSmallSpecialminus.push_back(new TH2F(Form("h_offetavsSA_resptSmallSpecialminus_%dGeV",i*m_thpitch),"SmallSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_resptLargeSpecialplus11out.push_back(new TH2F(Form("h_offetavsSA_resptLargeSpecialplus11out_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_resptLargeSpecialplus15out.push_back(new TH2F(Form("h_offetavsSA_resptLargeSpecialplus15out_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_resptLargeSpecialplus11in.push_back(new TH2F(Form("h_offetavsSA_resptLargeSpecialplus11in_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_resptLargeSpecialplus15in.push_back(new TH2F(Form("h_offetavsSA_resptLargeSpecialplus15in_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_resptLargeSpecialminus11out.push_back(new TH2F(Form("h_offetavsSA_resptLargeSpecialminus11out_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_resptLargeSpecialminus15out.push_back(new TH2F(Form("h_offetavsSA_resptLargeSpecialminus15out_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_resptLargeSpecialminus11in.push_back(new TH2F(Form("h_offetavsSA_resptLargeSpecialminus11in_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_offetavsSA_resptLargeSpecialminus15in.push_back(new TH2F(Form("h_offetavsSA_resptLargeSpecialminus15in_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));

               m_h_mdtetavsSA_resptLargeplus.push_back(new TH2F(Form("h_mdtetavsSA_resptLargeplus_%dGeV",i*m_thpitch),"Large;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_mdtetavsSA_resptLargeSpecialplus.push_back(new TH2F(Form("h_mdtetavsSA_resptLargeSpecialplus_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_mdtetavsSA_resptSmallplus.push_back(new TH2F(Form("h_mdtetavsSA_resptSmallplus_%dGeV",i*m_thpitch),"Small;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_mdtetavsSA_resptSmallSpecialplus.push_back(new TH2F(Form("h_mdtetavsSA_resptSmallSpecialplus_%dGeV",i*m_thpitch),"SmallSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_mdtetavsSA_resptLargeminus.push_back(new TH2F(Form("h_mdtetavsSA_resptLargeminus_%dGeV",i*m_thpitch),"Large;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_mdtetavsSA_resptLargeSpecialminus.push_back(new TH2F(Form("h_mdtetavsSA_resptLargeSpecialminus_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_mdtetavsSA_resptSmallminus.push_back(new TH2F(Form("h_mdtetavsSA_resptSmallminus_%dGeV",i*m_thpitch),"Small;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_mdtetavsSA_resptSmallSpecialminus.push_back(new TH2F(Form("h_mdtetavsSA_resptSmallSpecialminus_%dGeV",i*m_thpitch),"SmallSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_mdtetavsSA_resptLargeSpecialplus11out.push_back(new TH2F(Form("h_mdtetavsSA_resptLargeSpecialplus11out_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_mdtetavsSA_resptLargeSpecialplus15out.push_back(new TH2F(Form("h_mdtetavsSA_resptLargeSpecialplus15out_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_mdtetavsSA_resptLargeSpecialplus11in.push_back(new TH2F(Form("h_mdtetavsSA_resptLargeSpecialplus11in_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_mdtetavsSA_resptLargeSpecialplus15in.push_back(new TH2F(Form("h_mdtetavsSA_resptLargeSpecialplus15in_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_mdtetavsSA_resptLargeSpecialminus11out.push_back(new TH2F(Form("h_mdtetavsSA_resptLargeSpecialminus11out_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_mdtetavsSA_resptLargeSpecialminus15out.push_back(new TH2F(Form("h_mdtetavsSA_resptLargeSpecialminus15out_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_mdtetavsSA_resptLargeSpecialminus11in.push_back(new TH2F(Form("h_mdtetavsSA_resptLargeSpecialminus11in_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_mdtetavsSA_resptLargeSpecialminus15in.push_back(new TH2F(Form("h_mdtetavsSA_resptLargeSpecialminus15in_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));

               m_h_highoffetavsSA_resptLargeSpecialplus11out.push_back(new TH2F(Form("h_highoffetavsSA_resptLargeSpecialplus11out_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_highoffetavsSA_resptLargeSpecialplus15out.push_back(new TH2F(Form("h_highoffetavsSA_resptLargeSpecialplus15out_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_highoffetavsSA_resptLargeSpecialplus11in.push_back(new TH2F(Form("h_highoffetavsSA_resptLargeSpecialplus11in_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_highoffetavsSA_resptLargeSpecialplus15in.push_back(new TH2F(Form("h_highoffetavsSA_resptLargeSpecialplus15in_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_highoffetavsSA_resptLargeSpecialminus11out.push_back(new TH2F(Form("h_highoffetavsSA_resptLargeSpecialminus11out_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_highoffetavsSA_resptLargeSpecialminus15out.push_back(new TH2F(Form("h_highoffetavsSA_resptLargeSpecialminus15out_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_highoffetavsSA_resptLargeSpecialminus11in.push_back(new TH2F(Form("h_highoffetavsSA_resptLargeSpecial,minus11in_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_highoffetavsSA_resptLargeSpecialminus15in.push_back(new TH2F(Form("h_highoffetavsSA_resptLargeSpecialminus15in_%dGeV",i*m_thpitch),"LargeSpecial;offline eta;SApt residual",100,-2.5,2.5,300,-2.0,1.0));
               m_h_etaIndexvsSA_respt.push_back(new TH2F(Form("h_etaIndexvsSA_respt_%dGeV",i*m_thpitch),"etaIndex vs respt;offline eta;SApt residual",17,-8.5,8.5,300,-2.0,1.0));
               m_h_etaIndexvsSA_resptLargeplus.push_back(new TH2F(Form("h_etaIndexvsSA_resptLargeplus_%dGeV",i*m_thpitch),"etaIndex vs respt;offline eta;SApt residual",17,-8.5,8.5,300,-2.0,1.0));
               m_h_etaIndexvsSA_resptLargeSpecialplus.push_back(new TH2F(Form("h_etaIndexvsSA_resptLargeSpecialplus_%dGeV",i*m_thpitch),"etaIndex vs respt;offline eta;SApt residual",17,-8.5,8.5,300,-2.0,1.0));
               m_h_etaIndexvsSA_resptSmallplus.push_back(new TH2F(Form("h_etaIndexvsSA_resptSmallplus_%dGeV",i*m_thpitch),"etaIndex vs respt;offline eta;SApt residual",17,-8.5,8.5,300,-2.0,1.0));
               m_h_etaIndexvsSA_resptSmallSpecialplus.push_back(new TH2F(Form("h_etaIndexvsSA_resptSmallSpecialplus_%dGeV",i*m_thpitch),"etaIndex vs respt;offline eta;SApt residual",17,-8.5,8.5,300,-2.0,1.0));
               m_h_etaIndexvsSA_resptLargeminus.push_back(new TH2F(Form("h_etaIndexvsSA_resptLargeminus_%dGeV",i*m_thpitch),"etaIndex vs respt;offline eta;SApt residual",17,-8.5,8.5,300,-2.0,1.0));
               m_h_etaIndexvsSA_resptLargeSpecialminus.push_back(new TH2F(Form("h_etaIndexvsSA_resptLargeSpecialminus_%dGeV",i*m_thpitch),"etaIndex vs respt;offline eta;SApt residual",17,-8.5,8.5,300,-2.0,1.0));
               m_h_etaIndexvsSA_resptSmallminus.push_back(new TH2F(Form("h_etaIndexvsSA_resptSmallminus_%dGeV",i*m_thpitch),"etaIndex vs respt;offline eta;SApt residual",17,-8.5,8.5,300,-2.0,1.0));
               m_h_etaIndexvsSA_resptSmallSpecialminus.push_back(new TH2F(Form("h_etaIndexvsSA_resptSmallSpecialminus_%dGeV",i*m_thpitch),"etaIndex vs respt;offline eta;SApt residual",17,-8.5,8.5,300,-2.0,1.0));
               m_h_etaIndexvsSA_resptLargeSpecialplus11out.push_back(new TH2F(Form("h_etaIndexvsSA_resptLargeSpecialplus11out_%dGeV",i*m_thpitch),"etaIndex vs respt;offline eta;SApt residual",17,-8.5,8.5,300,-2.0,1.0));
               m_h_etaIndexvsSA_resptLargeSpecialplus11in.push_back(new TH2F(Form("h_etaIndexvsSA_resptLargeSpecialplus11in_%dGeV",i*m_thpitch),"etaIndex vs respt;offline eta;SApt residual",17,-8.5,8.5,300,-2.0,1.0));
               m_h_etaIndexvsSA_resptLargeSpecialplus15out.push_back(new TH2F(Form("h_etaIndexvsSA_resptLargeSpecialplus15out_%dGeV",i*m_thpitch),"etaIndex vs respt;offline eta;SApt residual",17,-8.5,8.5,300,-2.0,1.0));
               m_h_etaIndexvsSA_resptLargeSpecialplus15in.push_back(new TH2F(Form("h_etaIndexvsSA_resptLargeSpecialplus15in_%dGeV",i*m_thpitch),"etaIndex vs respt;offline eta;SApt residual",17,-8.5,8.5,300,-2.0,1.0));
               m_h_etaIndexvsSA_resptLargeSpecialminus11out.push_back(new TH2F(Form("h_etaIndexvsSA_resptLargeSpecialminus11out_%dGeV",i*m_thpitch),"etaIndex vs respt;offline eta;SApt residual",17,-8.5,8.5,300,-2.0,1.0));
               m_h_etaIndexvsSA_resptLargeSpecialminus11in.push_back(new TH2F(Form("h_etaIndexvsSA_resptLargeSpecialminus11in_%dGeV",i*m_thpitch),"etaIndex vs respt;offline eta;SApt residual",17,-8.5,8.5,300,-2.0,1.0));
               m_h_etaIndexvsSA_resptLargeSpecialminus15out.push_back(new TH2F(Form("h_etaIndexvsSA_resptLargeSpecialminus15out_%dGeV",i*m_thpitch),"etaIndex vs respt;offline eta;SApt residual",17,-8.5,8.5,300,-2.0,1.0));
               m_h_etaIndexvsSA_resptLargeSpecialminus15in.push_back(new TH2F(Form("h_etaIndexvsSA_resptLargeSpecialminus15in_%dGeV",i*m_thpitch),"etaIndex vs respt;offline eta;SApt residual",17,-8.5,8.5,300,-2.0,1.0));

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
     if(m_sumReqdRL1 < m_tp_extdR && 0.2 < m_tp_extdR && m_sumReqdREF < m_tp_dR && pass > -1 && m_tag_proc == m_proc){
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
     if(pass == 1 && std::fabs(pt) > th){
          return kTRUE;
     }else{
          return kFALSE;
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
     for(Int_t i = 0;i <= m_nhist;i++){
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
          vector<float> *pSA_rpcX = 0;
          vector<float> *pSA_rpcY = 0;
          vector<float> *pSA_mdtZ = 0;
          vector<float> *pSA_mdtR = 0;
          vector<float> *pSA_mdtPhi = 0;
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
                    pSA_rpcX = &(m_pSA_rpcX->at(method));
                    pSA_rpcY = &(m_pSA_rpcY->at(method));
                    /*vector<float> buf_pSA_mdtZ;
                    if((signed int)m_pSA_mdtZ->at(method).size() != 0)cout<<m_pSA_mdtZ->at(method).size()<<endl;
                    for(Int_t j = 0;j > (signed int)m_pSA_mdtZ->at(method).size();j++){
                         buf_pSA_mdtZ.push_back(m_pSA_mdtZ->at(method).at(j));
                         cout<<m_pSA_mdtZ->at(method).at(j)<<endl;
                    }
                    pSA_mdtZ = &(buf_pSA_mdtZ);*/
                    pSA_mdtZ = &(m_pSA_mdtZ->at(method));
                    pSA_mdtR = &(m_pSA_mdtR->at(method));
                    pSA_mdtPhi = &(m_pSA_mdtPhi->at(method));
               }
          }

          tL1_dR = TMath::Sqrt(pow(m_tL1_eta - m_toff_eta,2) + pow(m_tL1_phi - m_toff_phi,2) );
          tEF_dR = TMath::Sqrt(pow(m_tEF_eta - m_toff_eta,2) + pow(m_tEF_phi - m_toff_phi,2) );
          if(std::fabs(m_toff_pt)*0.001 < 10.0)m_reqL1dR = -0.00001*std::fabs(m_toff_pt) + 0.18;

          //offline
          if(Cut_tagprobe(pEFTAG_pass)){
               if(i == 0 && static_cast<Int_t>(pSA_sAddress) == 1)m_h_offphi_LargeSpecial->Fill(m_poff_phi);
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
                         if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1)m_h_eoff_pt_Largeplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1)m_h_eoff_pt_Largeminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                    case 1:
                         if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1)m_h_eoff_pt_LargeSpecialplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1)m_h_eoff_pt_LargeSpecialminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                    case 2:
                         if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1)m_h_eoff_pt_Smallplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1)m_h_eoff_pt_Smallminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                    case 3:
                         if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1)m_h_eoff_pt_SmallSpecialplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1)m_h_eoff_pt_SmallSpecialminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                    default:
                         break;
               }

               //L1
               if(Cut_L1(pL1_pass)){
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
                              if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1)m_h_eL1_pt_Largeplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1)m_h_eL1_pt_Largeminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                         case 1:
                              if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1){
                                   m_h_eL1_pt_LargeSpecialplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                   if(pSA_roiphi < -0.6 && pSA_roiphi > -1.0){//15
                                        m_h_eL1_pt_LargeSpecialplus15.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        if(pSA_roiphi > -0.8){
                                             m_h_eL1_pt_LargeSpecialplus15out.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        }else{
                                             m_h_eL1_pt_LargeSpecialplus15in.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        }
                                   }
                                   if(pSA_roiphi < -2.0 && pSA_roiphi > -2.6){//11
                                        m_h_eL1_pt_LargeSpecialplus11.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        if(pSA_roiphi > -2.4){
                                             m_h_eL1_pt_LargeSpecialplus11in.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        }else{
                                             m_h_eL1_pt_LargeSpecialplus11out.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        }
                                   }
                              }
                              if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1){
                                   m_h_eL1_pt_LargeSpecialminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                   if(pSA_roiphi < -0.6 && pSA_roiphi > -1.0){//15
                                        m_h_eL1_pt_LargeSpecialminus15.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        if(pSA_roiphi > -0.8){
                                             m_h_eL1_pt_LargeSpecialminus15out.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        }else{
                                             m_h_eL1_pt_LargeSpecialminus15in.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        }
                                   }
                                   if(pSA_roiphi < -2.0 && pSA_roiphi > -2.6){//11
                                        m_h_eL1_pt_LargeSpecialminus11.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        if(pSA_roiphi > -2.4){
                                             m_h_eL1_pt_LargeSpecialminus11in.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        }else{
                                             m_h_eL1_pt_LargeSpecialminus11out.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        }
                                   }
                              }

                              break;
                         case 2:
                              if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1)m_h_eL1_pt_Smallplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1)m_h_eL1_pt_Smallminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                         case 3:
                              if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1)m_h_eL1_pt_SmallSpecialplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1)m_h_eL1_pt_SmallSpecialminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                         default:
                              break;
                    }

                    //SA
                    if(Cut_SA(pSA_pass,pSA_pt,i*m_thpitch)){
                         Double_t textSA_dR = TMath::Sqrt(pow(m_tSA_eta - m_toff_exteta,2) + pow(m_tSA_phi - m_toff_extphi,2));
                         pextSA_dR = TMath::Sqrt(pow(pSA_eta - m_poff_exteta,2) + pow(pSA_phi - m_poff_extphi,2));
                         Double_t resSA_pt = std::fabs(m_poff_pt*0.001)/std::fabs(pSA_pt) - 1.0;
                         Double_t buf_pSA_dR = TMath::Sqrt(pow(pSA_eta - m_poff_eta,2) + pow(pSA_phi - m_poff_phi,2));
                         Double_t buf_eta = 0;
                         for(Int_t size = 0;size < (signed int)pSA_mdtZ->size();size++){
                              buf_eta += -TMath::Log((sqrt(pow(pSA_mdtZ->at(size),2) + pow(pSA_mdtR->at(size),2)) - pSA_mdtZ->at(size))/(sqrt(pow(pSA_mdtZ->at(size),2) + pow(pSA_mdtR->at(size),2)) + pow(pSA_mdtZ->at(size),2)))/2.0;
                         }
                         Double_t ave_mdteta = buf_eta/static_cast<Double_t>(pSA_mdtZ->size());

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
                                   if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1){
                                        m_h_off_ptvsSA_resptLargeplus.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                        m_h_SA_resptLargeplus.at(i)->Fill(resSA_pt);
                                        m_h_eSA_pt_Largeplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        m_h_offphivsSA_resptLargeplus.at(i)->Fill(m_poff_phi,resSA_pt);
                                        m_h_offetavsSA_resptLargeplus.at(i)->Fill(m_poff_eta,resSA_pt);
                                        m_h_mdtetavsSA_resptLargeplus.at(i)->Fill(ave_mdteta,resSA_pt);
                                        for(Int_t index = 0;index < 10;index++){
                                             if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeplus.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                        }
                                   }
                                   if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1){
                                        m_h_off_ptvsSA_resptLargeminus.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                        m_h_SA_resptLargeminus.at(i)->Fill(resSA_pt);
                                        m_h_eSA_pt_Largeminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        m_h_offphivsSA_resptLargeminus.at(i)->Fill(m_poff_phi,resSA_pt);
                                        m_h_offetavsSA_resptLargeminus.at(i)->Fill(m_poff_eta,resSA_pt);
                                        m_h_mdtetavsSA_resptLargeminus.at(i)->Fill(ave_mdteta,resSA_pt);
                                        for(Int_t index = 0;index < 10;index++){
                                             if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeminus.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                        }
                                   }
                                   m_countLarge.at(i)++;
                                   break;

                              case 1:
                                   //rpc,mdthitXY
                                   for(Int_t rpc = 0;rpc < (signed int)pSA_rpcX->size();rpc++){
                                        m_h_rpchitXYLargeSpecial.at(i)->Fill(pSA_rpcX->at(rpc),pSA_rpcY->at(rpc));
                                   }
                                   for(Int_t mdt = 0;mdt < (signed int)pSA_mdtZ->size();mdt++){
                                        m_h_mdthitXYLargeSpecial.at(i)->Fill(pSA_mdtR->at(mdt)*cos(pSA_mdtPhi->at(mdt)),pSA_mdtR->at(mdt)*sin(pSA_mdtPhi->at(mdt)));
                                   }

                                   if(i == 0){
                                        m_h_saphims_LargeSpecial->Fill(pSA_phims);
                                        m_h_saroiphi_LargeSpecial->Fill(pSA_roiphi);
                                   }

                                   if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1){
                                        m_h_eSA_pt_LargeSpecialplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        m_h_SA_resptLargeSpecialplus.at(i)->Fill(resSA_pt);
                                        m_h_off_ptvsSA_resptLargeSpecialplus.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                        m_h_offphivsSA_resptLargeSpecialplus.at(i)->Fill(m_poff_phi,resSA_pt);
                                        m_h_offetavsSA_resptLargeSpecialplus.at(i)->Fill(m_poff_eta,resSA_pt);
                                        m_h_mdtetavsSA_resptLargeSpecialplus.at(i)->Fill(ave_mdteta,resSA_pt);
                                        for(Int_t index = 0;index < 10;index++){
                                             if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialplus.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                        }

                                        if(pSA_roiphi < -0.6 && pSA_roiphi > -1.0){//15
                                             if(pSA_roiphi > -0.8){
                                                  m_h_off_ptvsSA_resptLargeSpecialplus15out.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                                  m_h_SA_resptLargeSpecialplus15out.at(i)->Fill(resSA_pt);
                                                  m_h_offetavsSA_resptLargeSpecialplus15out.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  m_h_mdtetavsSA_resptLargeSpecialplus15out.at(i)->Fill(ave_mdteta,resSA_pt);
                                                  m_h_eSA_pt_LargeSpecialplus15out.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                                  if(std::fabs(m_poff_pt*0.001) > 30.0 && std::fabs(m_poff_pt*0.001) < 50.0){
                                                       m_h_highoffetavsSA_resptLargeSpecialplus15out.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  }
                                                  for(Int_t index = 0;index < 10;index++){
                                                      if(m_probe_segment_etaIndex[index] >= -8 && m_probe_segment_etaIndex[index] <= 8)m_h_etaIndexvsSA_resptLargeSpecialplus15out.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                                  }
                                                  //rpc,mdthitXY
                                                  for(Int_t rpc = 0;rpc < (signed int)pSA_rpcX->size();rpc++){
                                                       m_h_rpchitXYLargeSpecialplus15out.at(i)->Fill(pSA_rpcX->at(rpc),pSA_rpcY->at(rpc));
                                                  }
                                                  for(Int_t mdt = 0;mdt < (signed int)pSA_mdtZ->size();mdt++){
                                                       m_h_mdthitXYLargeSpecialplus15out.at(i)->Fill(pSA_mdtR->at(mdt)*cos(pSA_mdtPhi->at(mdt)),pSA_mdtR->at(mdt)*sin(pSA_mdtPhi->at(mdt)));
                                                  }

                                             }else{
                                                  m_h_off_ptvsSA_resptLargeSpecialplus15in.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                                  m_h_SA_resptLargeSpecialplus15in.at(i)->Fill(resSA_pt);
                                                  m_h_offetavsSA_resptLargeSpecialplus15in.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  m_h_mdtetavsSA_resptLargeSpecialplus15in.at(i)->Fill(ave_mdteta,resSA_pt);
                                                  m_h_eSA_pt_LargeSpecialplus15in.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                                  for(Int_t index = 0;index < 10;index++){
                                                      if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialplus15in.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                                  }
                                                  if(std::fabs(m_poff_pt*0.001) > 30.0 && std::fabs(m_poff_pt*0.001) < 50.0){
                                                       m_h_highoffetavsSA_resptLargeSpecialplus15in.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  }
                                             //rpc,mdthitXY
                                                  for(Int_t rpc = 0;rpc < (signed int)pSA_rpcX->size();rpc++){
                                                       m_h_rpchitXYLargeSpecialplus15in.at(i)->Fill(pSA_rpcX->at(rpc),pSA_rpcY->at(rpc));
                                                  }
                                                  for(Int_t mdt = 0;mdt < (signed int)pSA_mdtZ->size();mdt++){
                                                       m_h_mdthitXYLargeSpecialplus15in.at(i)->Fill(pSA_mdtR->at(mdt)*cos(pSA_mdtPhi->at(mdt)),pSA_mdtR->at(mdt)*sin(pSA_mdtPhi->at(mdt)));
                                                  }
                                             }
                                        }
                                        if(pSA_roiphi < -2.0 && pSA_roiphi > -2.6){//11
                                             if(pSA_roiphi > -2.4){
                                                  m_h_off_ptvsSA_resptLargeSpecialplus11in.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                                  m_h_SA_resptLargeSpecialplus11in.at(i)->Fill(resSA_pt);
                                                  m_h_offetavsSA_resptLargeSpecialplus11in.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  m_h_mdtetavsSA_resptLargeSpecialplus11in.at(i)->Fill(ave_mdteta,resSA_pt);
                                                  m_h_eSA_pt_LargeSpecialplus11in.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                                  for(Int_t index = 0;index < 10;index++){
                                                      if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialplus11in.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                                  }
                                                  if(std::fabs(m_poff_pt*0.001) > 30.0 && std::fabs(m_poff_pt*0.001) < 50.0){
                                                       m_h_highoffetavsSA_resptLargeSpecialplus11in.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  }
                                                  //rpc,mdthitXY
                                                  for(Int_t rpc = 0;rpc < (signed int)pSA_rpcX->size();rpc++){
                                                       m_h_rpchitXYLargeSpecialplus11in.at(i)->Fill(pSA_rpcX->at(rpc),pSA_rpcY->at(rpc));
                                                  }
                                                  for(Int_t mdt = 0;mdt < (signed int)pSA_mdtZ->size();mdt++){
                                                       m_h_mdthitXYLargeSpecialplus11in.at(i)->Fill(pSA_mdtR->at(mdt)*cos(pSA_mdtPhi->at(mdt)),pSA_mdtR->at(mdt)*sin(pSA_mdtPhi->at(mdt)));
                                                  }

                                             }else{
                                                  m_h_off_ptvsSA_resptLargeSpecialplus11out.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                                  m_h_SA_resptLargeSpecialplus11out.at(i)->Fill(resSA_pt);
                                                  m_h_offetavsSA_resptLargeSpecialplus11out.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  m_h_mdtetavsSA_resptLargeSpecialplus11out.at(i)->Fill(ave_mdteta,resSA_pt);
                                                  for(Int_t index = 0;index < 10;index++){
                                                      if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialplus11out.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                                  }
                                                  if(std::fabs(m_poff_pt*0.001) > 30.0 && std::fabs(m_poff_pt*0.001) < 50.0){
                                                       m_h_highoffetavsSA_resptLargeSpecialplus11out.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  }
                                             //rpc,mdthitXY
                                                  for(Int_t rpc = 0;rpc < (signed int)pSA_rpcX->size();rpc++){
                                                       m_h_rpchitXYLargeSpecialplus11out.at(i)->Fill(pSA_rpcX->at(rpc),pSA_rpcY->at(rpc));
                                                  }
                                                  for(Int_t mdt = 0;mdt < (signed int)pSA_mdtZ->size();mdt++){
                                                       m_h_mdthitXYLargeSpecialplus11out.at(i)->Fill(pSA_mdtR->at(mdt)*cos(pSA_mdtPhi->at(mdt)),pSA_mdtR->at(mdt)*sin(pSA_mdtPhi->at(mdt)));
                                                  }
                                             }
                                        }
                                   }

                                   if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1){
                                        m_h_eSA_pt_LargeSpecialminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        m_h_SA_resptLargeSpecialminus.at(i)->Fill(resSA_pt);
                                        m_h_off_ptvsSA_resptLargeSpecialminus.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                        m_h_offphivsSA_resptLargeSpecialminus.at(i)->Fill(m_poff_phi,resSA_pt);
                                        m_h_offetavsSA_resptLargeSpecialminus.at(i)->Fill(m_poff_eta,resSA_pt);
                                        m_h_mdtetavsSA_resptLargeSpecialminus.at(i)->Fill(ave_mdteta,resSA_pt);
                                        for(Int_t index = 0;index < 10;index++){
                                             if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialminus.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                        }

                                        if(pSA_roiphi < -0.6 && pSA_roiphi > -1.0){//15
                                             if(pSA_roiphi > -0.8){
                                                  m_h_off_ptvsSA_resptLargeSpecialminus15out.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                                  m_h_SA_resptLargeSpecialminus15out.at(i)->Fill(resSA_pt);
                                                  m_h_offetavsSA_resptLargeSpecialminus15out.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  m_h_mdtetavsSA_resptLargeSpecialminus15out.at(i)->Fill(ave_mdteta,resSA_pt);
                                                  m_h_eSA_pt_LargeSpecialminus15out.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                                  for(Int_t index = 0;index < 10;index++){
                                                      if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialminus15out.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                                  }
                                                  if(std::fabs(m_poff_pt*0.001) > 30.0 && std::fabs(m_poff_pt*0.001) < 50.0){
                                                       m_h_highoffetavsSA_resptLargeSpecialminus15out.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  }
                                                  //rpc,mdthitXY
                                                  for(Int_t rpc = 0;rpc < (signed int)pSA_rpcX->size();rpc++){
                                                       m_h_rpchitXYLargeSpecialminus15out.at(i)->Fill(pSA_rpcX->at(rpc),pSA_rpcY->at(rpc));
                                                  }
                                                  for(Int_t mdt = 0;mdt < (signed int)pSA_mdtZ->size();mdt++){
                                                       m_h_mdthitXYLargeSpecialminus15out.at(i)->Fill(pSA_mdtR->at(mdt)*cos(pSA_mdtPhi->at(mdt)),pSA_mdtR->at(mdt)*sin(pSA_mdtPhi->at(mdt)));
                                                  }

                                             }else{
                                                  m_h_off_ptvsSA_resptLargeSpecialminus15in.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                                  m_h_SA_resptLargeSpecialminus15in.at(i)->Fill(resSA_pt);
                                                  m_h_offetavsSA_resptLargeSpecialminus15in.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  m_h_mdtetavsSA_resptLargeSpecialminus15in.at(i)->Fill(ave_mdteta,resSA_pt);
                                                  m_h_eSA_pt_LargeSpecialminus15in.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                                  for(Int_t index = 0;index < 10;index++){
                                                      if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialminus15in.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                                  }
                                                  if(std::fabs(m_poff_pt*0.001) > 30.0 && std::fabs(m_poff_pt*0.001) < 50.0){
                                                       m_h_highoffetavsSA_resptLargeSpecialminus15in.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  }
                                             //rpc,mdthitXY
                                                  for(Int_t rpc = 0;rpc < (signed int)pSA_rpcX->size();rpc++){
                                                       m_h_rpchitXYLargeSpecialminus15in.at(i)->Fill(pSA_rpcX->at(rpc),pSA_rpcY->at(rpc));
                                                  }
                                                  for(Int_t mdt = 0;mdt < (signed int)pSA_mdtZ->size();mdt++){
                                                       m_h_mdthitXYLargeSpecialminus15in.at(i)->Fill(pSA_mdtR->at(mdt)*cos(pSA_mdtPhi->at(mdt)),pSA_mdtR->at(mdt)*sin(pSA_mdtPhi->at(mdt)));
                                                  }
                                             }
                                        }

                                        if(pSA_roiphi < -2.0 && pSA_roiphi > -2.6){//11
                                             if(pSA_roiphi > -2.4){
                                                  m_h_off_ptvsSA_resptLargeSpecialminus11in.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                                  m_h_SA_resptLargeSpecialminus11in.at(i)->Fill(resSA_pt);
                                                  m_h_offetavsSA_resptLargeSpecialminus11in.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  m_h_mdtetavsSA_resptLargeSpecialminus11in.at(i)->Fill(ave_mdteta,resSA_pt);
                                                  m_h_eSA_pt_LargeSpecialminus11in.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                                  for(Int_t index = 0;index < 10;index++){
                                                      if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialminus11in.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                                  }
                                                  if(std::fabs(m_poff_pt*0.001) > 30.0 && std::fabs(m_poff_pt*0.001) < 50.0){
                                                       m_h_highoffetavsSA_resptLargeSpecialminus11in.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  }
                                                  //rpc,mdthitXY
                                                  for(Int_t rpc = 0;rpc < (signed int)pSA_rpcX->size();rpc++){
                                                       m_h_rpchitXYLargeSpecialminus11in.at(i)->Fill(pSA_rpcX->at(rpc),pSA_rpcY->at(rpc));
                                                  }
                                                  for(Int_t mdt = 0;mdt < (signed int)pSA_mdtZ->size();mdt++){
                                                       m_h_mdthitXYLargeSpecialminus11in.at(i)->Fill(pSA_mdtR->at(mdt)*cos(pSA_mdtPhi->at(mdt)),pSA_mdtR->at(mdt)*sin(pSA_mdtPhi->at(mdt)));
                                                  }
                                             }else{
                                                  m_h_off_ptvsSA_resptLargeSpecialminus11out.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                                  m_h_SA_resptLargeSpecialminus11out.at(i)->Fill(resSA_pt);
                                                  m_h_offetavsSA_resptLargeSpecialminus11out.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  m_h_mdtetavsSA_resptLargeSpecialminus11out.at(i)->Fill(ave_mdteta,resSA_pt);
                                                  for(Int_t index = 0;index < 10;index++){
                                                      if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialminus11out.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                                  }
                                                  if(std::fabs(m_poff_pt*0.001) > 30.0 && std::fabs(m_poff_pt*0.001) < 50.0){
                                                       m_h_highoffetavsSA_resptLargeSpecialminus11out.at(i)->Fill(m_poff_eta,resSA_pt);
                                                  }
                                             //rpc,mdthitXY
                                                  for(Int_t rpc = 0;rpc < (signed int)pSA_rpcX->size();rpc++){
                                                       m_h_rpchitXYLargeSpecialminus11out.at(i)->Fill(pSA_rpcX->at(rpc),pSA_rpcY->at(rpc));
                                                  }
                                                  for(Int_t mdt = 0;mdt < (signed int)pSA_mdtZ->size();mdt++){
                                                       m_h_mdthitXYLargeSpecialminus11out.at(i)->Fill(pSA_mdtR->at(mdt)*cos(pSA_mdtPhi->at(mdt)),pSA_mdtR->at(mdt)*sin(pSA_mdtPhi->at(mdt)));
                                                  }
                                             }
                                        }
                                   }
                                   m_countLargeSpecial.at(i)++;
                                   break;

                              case 2:
                                   if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1){
                                        m_h_off_ptvsSA_resptSmallplus.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                        m_h_SA_resptSmallplus.at(i)->Fill(resSA_pt);
                                        m_h_eSA_pt_Smallplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        m_h_offphivsSA_resptSmallplus.at(i)->Fill(m_poff_phi,resSA_pt);
                                        m_h_offetavsSA_resptSmallplus.at(i)->Fill(m_poff_eta,resSA_pt);
                                        m_h_mdtetavsSA_resptSmallplus.at(i)->Fill(ave_mdteta,resSA_pt);
                                        for(Int_t index = 0;index < 10;index++){
                                             if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptSmallplus.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                        }
                                   }
                                   if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1){
                                        m_h_off_ptvsSA_resptSmallminus.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                        m_h_SA_resptSmallminus.at(i)->Fill(resSA_pt);
                                        m_h_eSA_pt_Smallminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        m_h_offphivsSA_resptSmallminus.at(i)->Fill(m_poff_phi,resSA_pt);
                                        m_h_offetavsSA_resptSmallminus.at(i)->Fill(m_poff_eta,resSA_pt);
                                        m_h_mdtetavsSA_resptSmallminus.at(i)->Fill(ave_mdteta,resSA_pt);
                                        for(Int_t index = 0;index < 10;index++){
                                             if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptSmallminus.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                        }
                                   }
                                   m_countSmall.at(i)++;
                                   break;

                              case 3:
                                   if(i == 0)m_h_saroiphi_SmallSpecial->Fill(pSA_roiphi);
                                   if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1){
                                        m_h_off_ptvsSA_resptSmallSpecialplus.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                        m_h_SA_resptSmallSpecialplus.at(i)->Fill(resSA_pt);
                                        m_h_eSA_pt_SmallSpecialplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        m_h_offphivsSA_resptSmallSpecialplus.at(i)->Fill(m_poff_phi,resSA_pt);
                                        m_h_offetavsSA_resptSmallSpecialplus.at(i)->Fill(m_poff_eta,resSA_pt);
                                        m_h_mdtetavsSA_resptSmallSpecialplus.at(i)->Fill(ave_mdteta,resSA_pt);
                                        for(Int_t index = 0;index < 10;index++){
                                             if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptSmallSpecialplus.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                        }
                                   }
                                   if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1){
                                        m_h_off_ptvsSA_resptSmallSpecialminus.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                        m_h_SA_resptSmallSpecialminus.at(i)->Fill(resSA_pt);
                                        m_h_eSA_pt_SmallSpecialminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                        m_h_offphivsSA_resptSmallSpecialminus.at(i)->Fill(m_poff_phi,resSA_pt);
                                        m_h_offetavsSA_resptSmallSpecialminus.at(i)->Fill(m_poff_eta,resSA_pt);
                                        for(Int_t index = 0;index < 10;index++){
                                             if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptSmallSpecialminus.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                        }
                                   }
                                   m_countSmallSpecial.at(i)++;
                                   break;
                              default:
                                   break;
                         }

                         if(static_cast<Int_t>(pSA_sAddress) == 0 || static_cast<Int_t>(pSA_sAddress) == 1 || static_cast<Int_t>(pSA_sAddress) == 2 || static_cast<Int_t>(pSA_sAddress) == 3){
                              m_h_offphivsSA_sAddress.at(i)->Fill(pSA_phims,pSA_sAddress);
                              for(Int_t size = 0;size < (signed int)pSA_rpcX->size();size++){
                                   m_h_rpchitXY.at(i)->Fill(pSA_rpcX->at(size),pSA_rpcY->at(size));
                              }
                              for(Int_t size = 0;size < (signed int)pSA_mdtZ->size();size++){
                                   m_h_mdthitXY.at(i)->Fill(pSA_mdtR->at(size)*cos(pSA_mdtR->at(size)),pSA_mdtR->at(size)*sin(pSA_mdtR->at(size)));
                              }
                              for(Int_t index = 0;index < 10;index++){
                                   if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_respt.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              }
                         }
                         m_h_offphivsSAphims.at(i)->Fill(m_poff_phi,pSA_phims);

                         //CB
                         if(Cut_CB(pCB_pass)){
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
                              }//EF

                         }//CB

                    }//SA

               }//L1

          }//tag

     }//for

}//Exe

void Efficiency::Finalize(TFile *tf1){
     CalcEfficiency ceff;
     tf1->cd();
     cout<<"ptSAth   nLarge   nLargeS   nSmall   nSmallS"<<endl;
     //SetCondition
     //title,file title,yoffset,top margin,bottom margin,left margin,right margin

     m_h_saroiphi_LargeSpecial->Write();
     m_h_saroiphi_SmallSpecial->Write();
     m_h_saphims_LargeSpecial->Write();
     m_h_offphi_LargeSpecial->Write();
     for(Int_t i = 0;i <= m_nhist;i++){
          /*ceff.SetCondition("test",1.5,0,0,0,0);
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
          ceff.DrawHist2D(m_h_off_ptvsSA_resptminusLSminus15.at(i));*/

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

          ceff.SetConditionName(Form("L1EfficiencyLargeplus_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 Large Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_Largeplus.at(i),m_h_eL1_pt_Largeplus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeplus_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA Large Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_Largeplus.at(i),m_h_eSA_pt_Largeplus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1EfficiencyLargeSpecialplus_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_LargeSpecialplus.at(i),m_h_eL1_pt_LargeSpecialplus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialplus_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus.at(i),m_h_eSA_pt_LargeSpecialplus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1EfficiencyLargeminus_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 Large Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_Largeminus.at(i),m_h_eL1_pt_Largeminus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeminus_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA Large Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_Largeminus.at(i),m_h_eSA_pt_Largeminus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1EfficiencyLargeSpecialminus_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_LargeSpecialminus.at(i),m_h_eL1_pt_LargeSpecialminus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialminus_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus.at(i),m_h_eSA_pt_LargeSpecialminus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialplus11_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus11.at(i),m_h_eSA_pt_LargeSpecialplus11.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialminus11_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus11.at(i),m_h_eSA_pt_LargeSpecialminus11.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialplus15_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus15.at(i),m_h_eSA_pt_LargeSpecialplus15.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialminus15_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus15.at(i),m_h_eSA_pt_LargeSpecialminus15.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialplus11out_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus11out.at(i),m_h_eSA_pt_LargeSpecialplus11out.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialminus11out_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus11out.at(i),m_h_eSA_pt_LargeSpecialminus11out.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialplus15out_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus15out.at(i),m_h_eSA_pt_LargeSpecialplus15out.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialminus15out_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus15out.at(i),m_h_eSA_pt_LargeSpecialminus15out.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialplus11in_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus11in.at(i),m_h_eSA_pt_LargeSpecialplus11in.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialminus11in_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus11in.at(i),m_h_eSA_pt_LargeSpecialminus11in.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialplus15in_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus15in.at(i),m_h_eSA_pt_LargeSpecialplus15in.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialminus15in_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus15in.at(i),m_h_eSA_pt_LargeSpecialminus15in.at(i),m_binmax,300,m_efficiency_xerr);

          ceff.SetConditionName(Form("L1EfficiencySmallplus_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 Small Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_Smallplus.at(i),m_h_eL1_pt_Smallplus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencySmallplus_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA Small Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_Smallplus.at(i),m_h_eSA_pt_Smallplus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1EfficiencySmallminus_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 Small Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_Smallminus.at(i),m_h_eL1_pt_Smallminus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencySmallminus_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA Small Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_Smallminus.at(i),m_h_eSA_pt_Smallminus.at(i),m_binmax,300,m_efficiency_xerr);

          ceff.SetConditionName(Form("L1EfficiencySmallSpecialplus_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 SmallSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_SmallSpecialplus.at(i),m_h_eL1_pt_SmallSpecialplus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencySmallSpecialplus_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA SmallSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_SmallSpecialplus.at(i),m_h_eSA_pt_SmallSpecialplus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1EfficiencySmallSpecialminus_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 SmallSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_SmallSpecialminus.at(i),m_h_eL1_pt_SmallSpecialminus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencySmallSpecialminus_%dGeV",i*m_thpitch));
          ceff.SetCondition("SA SmallSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_SmallSpecialminus.at(i),m_h_eSA_pt_SmallSpecialminus.at(i),m_binmax,300,m_efficiency_xerr);

          ceff.SetConditionName(Form("SA2DEfficiency_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1vsL2MuonSA Efficiency;offline eta;offline phi",1.5,0.1,0.1,0.105,0.165);
          ceff.SetConditionbin(m_nbin_eta,m_nbin_phi,m_eta_max,m_phi_max);
          ceff.DrawEfficiency2D(m_h_eff_pL1_etaphi.at(i),m_h_eff_pSA_etaphi.at(i));
          ceff.SetConditionName(Form("L12DEfficiency_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 Efficiency;offline eta;offline phi",1.5,0.1,0.1,0.105,0.165);
          ceff.SetConditionbin(m_nbin_eta,m_nbin_phi,m_eta_max,m_phi_max);
          ceff.DrawEfficiency2D(m_h_eff_poff_etaphi.at(i),m_h_eff_pL1_etaphi.at(i));

          m_h_poff_pt.at(i)->Write();
          m_h_pL1_pt.at(i)->Write();
          m_h_pSA_pt.at(i)->Write();
          m_h_pCB_pt.at(i)->Write();
          m_h_pEF_pt.at(i)->Write();
          m_h_pL1_dR.at(i)->Write();
          m_h_pSA_dR.at(i)->Write();
          m_h_pCB_dR.at(i)->Write();
          m_h_pEF_dR.at(i)->Write();
          m_h_textL1_dR.at(i)->Write();
          m_h_textSA_dR.at(i)->Write();
          m_h_textCB_dR.at(i)->Write();
          m_h_textEF_dR.at(i)->Write();
          m_h_pextL1_dR.at(i)->Write();
          m_h_pextSA_dR.at(i)->Write();
          m_h_pextCB_dR.at(i)->Write();
          m_h_pextEF_dR.at(i)->Write();
          m_h_poffvsSA_pt.at(i)->Write();
          m_h_offphivsSA_sAddress.at(i)->Write();
          m_h_offphivsSAphims.at(i)->Write();
          m_h_rpchitXY.at(i)->Write();
          m_h_rpchitXYLargeSpecial.at(i)->Write();
          m_h_rpchitXYLargeSpecialplus11out.at(i)->Write();
          m_h_rpchitXYLargeSpecialplus11in.at(i)->Write();
          m_h_rpchitXYLargeSpecialplus15out.at(i)->Write();
          m_h_rpchitXYLargeSpecialplus15in.at(i)->Write();
          m_h_rpchitXYLargeSpecialminus11out.at(i)->Write();
          m_h_rpchitXYLargeSpecialminus11in.at(i)->Write();
          m_h_rpchitXYLargeSpecialminus15out.at(i)->Write();
          m_h_rpchitXYLargeSpecialminus15in.at(i)->Write();
          m_h_mdthitXY.at(i)->Write();
          m_h_mdthitXYLargeSpecial.at(i)->Write();
          m_h_mdthitXYLargeSpecialplus11out.at(i)->Write();
          m_h_mdthitXYLargeSpecialplus11in.at(i)->Write();
          m_h_mdthitXYLargeSpecialplus15out.at(i)->Write();
          m_h_mdthitXYLargeSpecialplus15in.at(i)->Write();
          m_h_mdthitXYLargeSpecialminus11out.at(i)->Write();
          m_h_mdthitXYLargeSpecialminus11in.at(i)->Write();
          m_h_mdthitXYLargeSpecialminus15out.at(i)->Write();
          m_h_mdthitXYLargeSpecialminus15in.at(i)->Write();
          m_h_eoff_pt.at(i)->Write();
          m_h_eL1_pt.at(i)->Write();
          m_h_eSA_pt.at(i)->Write();
          m_h_eCB_pt.at(i)->Write();
          m_h_eEF_pt.at(i)->Write();
          m_h_eoff_eta.at(i)->Write();
          m_h_eL1_eta.at(i)->Write();
          m_h_eSA_eta.at(i)->Write();
          m_h_eCB_eta.at(i)->Write();
          m_h_eEF_eta.at(i)->Write();
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
          m_h_eoff_pt_Largeplus.at(i)->Write();
          m_h_eL1_pt_Largeplus.at(i)->Write();
          m_h_eSA_pt_Largeplus.at(i)->Write();
          m_h_eoff_pt_Largeminus.at(i)->Write();
          m_h_eL1_pt_Largeminus.at(i)->Write();
          m_h_eSA_pt_Largeminus.at(i)->Write();
          m_h_eoff_pt_LargeSpecialplus.at(i)->Write();
          m_h_eL1_pt_LargeSpecialplus.at(i)->Write();
          m_h_eSA_pt_LargeSpecialplus.at(i)->Write();
          m_h_eL1_pt_LargeSpecialplus11.at(i)->Write();
          m_h_eSA_pt_LargeSpecialplus11.at(i)->Write();
          m_h_eL1_pt_LargeSpecialplus15.at(i)->Write();
          m_h_eSA_pt_LargeSpecialplus15.at(i)->Write();
          m_h_eL1_pt_LargeSpecialplus11in.at(i)->Write();
          m_h_eSA_pt_LargeSpecialplus11in.at(i)->Write();
          m_h_eL1_pt_LargeSpecialplus11out.at(i)->Write();
          m_h_eSA_pt_LargeSpecialplus11out.at(i)->Write();
          m_h_eL1_pt_LargeSpecialplus15out.at(i)->Write();
          m_h_eSA_pt_LargeSpecialplus15out.at(i)->Write();
          m_h_eL1_pt_LargeSpecialplus15in.at(i)->Write();
          m_h_eSA_pt_LargeSpecialplus15in.at(i)->Write();
          m_h_eoff_pt_LargeSpecialminus.at(i)->Write();
          m_h_eL1_pt_LargeSpecialminus.at(i)->Write();
          m_h_eSA_pt_LargeSpecialminus.at(i)->Write();
          m_h_eL1_pt_LargeSpecialminus11.at(i)->Write();
          m_h_eSA_pt_LargeSpecialminus11.at(i)->Write();
          m_h_eL1_pt_LargeSpecialminus15.at(i)->Write();
          m_h_eSA_pt_LargeSpecialminus15.at(i)->Write();
          m_h_eL1_pt_LargeSpecialminus11in.at(i)->Write();
          m_h_eSA_pt_LargeSpecialminus11in.at(i)->Write();
          m_h_eL1_pt_LargeSpecialminus11out.at(i)->Write();
          m_h_eSA_pt_LargeSpecialminus11out.at(i)->Write();
          m_h_eL1_pt_LargeSpecialminus15out.at(i)->Write();
          m_h_eSA_pt_LargeSpecialminus15out.at(i)->Write();
          m_h_eL1_pt_LargeSpecialminus15in.at(i)->Write();
          m_h_eSA_pt_LargeSpecialminus15in.at(i)->Write();
          m_h_eoff_pt_Smallplus.at(i)->Write();
          m_h_eL1_pt_Smallplus.at(i)->Write();
          m_h_eSA_pt_Smallplus.at(i)->Write();
          m_h_eoff_pt_Smallminus.at(i)->Write();
          m_h_eL1_pt_Smallminus.at(i)->Write();
          m_h_eSA_pt_Smallminus.at(i)->Write();
          m_h_eoff_pt_SmallSpecialplus.at(i)->Write();
          m_h_eL1_pt_SmallSpecialplus.at(i)->Write();
          m_h_eSA_pt_SmallSpecialplus.at(i)->Write();
          m_h_eoff_pt_SmallSpecialminus.at(i)->Write();
          m_h_eL1_pt_SmallSpecialminus.at(i)->Write();
          m_h_eSA_pt_SmallSpecialminus.at(i)->Write();
          m_h_eff_poff_etaphi.at(i)->Write();
          m_h_eff_pL1_etaphi.at(i)->Write();
          m_h_eff_pSA_etaphi.at(i)->Write();
          m_h_pSA_respt.at(i)->Write();
          m_h_pCB_respt.at(i)->Write();
          m_h_pEF_respt.at(i)->Write();
          m_h_SA_resptLargeplus.at(i)->Write();
          m_h_SA_resptLargeSpecialplus.at(i)->Write();
          m_h_SA_resptSmallplus.at(i)->Write();
          m_h_SA_resptSmallSpecialplus.at(i)->Write();
          m_h_SA_resptLargeminus.at(i)->Write();
          m_h_SA_resptLargeSpecialminus.at(i)->Write();
          m_h_SA_resptSmallminus.at(i)->Write();
          m_h_SA_resptSmallSpecialminus.at(i)->Write();
          m_h_SA_resptLargeSpecialplus11out.at(i)->Write();
          m_h_SA_resptLargeSpecialplus11in.at(i)->Write();
          m_h_SA_resptLargeSpecialplus15out.at(i)->Write();
          m_h_SA_resptLargeSpecialplus15in.at(i)->Write();
          m_h_SA_resptLargeSpecialminus11out.at(i)->Write();
          m_h_SA_resptLargeSpecialminus11in.at(i)->Write();
          m_h_SA_resptLargeSpecialminus15out.at(i)->Write();
          m_h_SA_resptLargeSpecialminus15in.at(i)->Write();
          m_h_off_ptvsSA_resptLargeplus.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialplus.at(i)->Write();
          m_h_off_ptvsSA_resptSmallplus.at(i)->Write();
          m_h_off_ptvsSA_resptSmallSpecialplus.at(i)->Write();
          m_h_off_ptvsSA_resptLargeminus.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialminus.at(i)->Write();
          m_h_off_ptvsSA_resptSmallminus.at(i)->Write();
          m_h_off_ptvsSA_resptSmallSpecialminus.at(i)->Write();
          m_h_offphivsSA_resptLargeplus.at(i)->Write();
          m_h_offphivsSA_resptLargeSpecialplus.at(i)->Write();
          m_h_offphivsSA_resptSmallplus.at(i)->Write();
          m_h_offphivsSA_resptSmallSpecialplus.at(i)->Write();
          m_h_offphivsSA_resptLargeminus.at(i)->Write();
          m_h_offphivsSA_resptLargeSpecialminus.at(i)->Write();
          m_h_offphivsSA_resptSmallminus.at(i)->Write();
          m_h_offphivsSA_resptSmallSpecialminus.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialplus15out.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialplus15in.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialplus11out.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialplus11in.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialminus15out.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialminus15in.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialminus11out.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialminus11in.at(i)->Write();
          m_h_offetavsSA_resptLargeplus.at(i)->Write();
          m_h_offetavsSA_resptLargeSpecialplus.at(i)->Write();
          m_h_offetavsSA_resptSmallplus.at(i)->Write();
          m_h_offetavsSA_resptSmallSpecialplus.at(i)->Write();
          m_h_offetavsSA_resptLargeminus.at(i)->Write();
          m_h_offetavsSA_resptLargeSpecialminus.at(i)->Write();
          m_h_offetavsSA_resptSmallminus.at(i)->Write();
          m_h_offetavsSA_resptSmallSpecialminus.at(i)->Write();
          m_h_offetavsSA_resptLargeSpecialplus11out.at(i)->Write();
          m_h_offetavsSA_resptLargeSpecialplus11in.at(i)->Write();
          m_h_offetavsSA_resptLargeSpecialplus15out.at(i)->Write();
          m_h_offetavsSA_resptLargeSpecialplus15in.at(i)->Write();
          m_h_offetavsSA_resptLargeSpecialminus11out.at(i)->Write();
          m_h_offetavsSA_resptLargeSpecialminus11in.at(i)->Write();
          m_h_offetavsSA_resptLargeSpecialminus15out.at(i)->Write();
          m_h_offetavsSA_resptLargeSpecialminus15in.at(i)->Write();
          m_h_mdtetavsSA_resptLargeplus.at(i)->Write();
          m_h_mdtetavsSA_resptLargeSpecialplus.at(i)->Write();
          m_h_mdtetavsSA_resptSmallplus.at(i)->Write();
          m_h_mdtetavsSA_resptSmallSpecialplus.at(i)->Write();
          m_h_mdtetavsSA_resptLargeminus.at(i)->Write();
          m_h_mdtetavsSA_resptLargeSpecialminus.at(i)->Write();
          m_h_mdtetavsSA_resptSmallminus.at(i)->Write();
          m_h_mdtetavsSA_resptSmallSpecialminus.at(i)->Write();
          m_h_mdtetavsSA_resptLargeSpecialplus11out.at(i)->Write();
          m_h_mdtetavsSA_resptLargeSpecialplus11in.at(i)->Write();
          m_h_mdtetavsSA_resptLargeSpecialplus15out.at(i)->Write();
          m_h_mdtetavsSA_resptLargeSpecialplus15in.at(i)->Write();
          m_h_mdtetavsSA_resptLargeSpecialminus11out.at(i)->Write();
          m_h_mdtetavsSA_resptLargeSpecialminus11in.at(i)->Write();
          m_h_mdtetavsSA_resptLargeSpecialminus15out.at(i)->Write();
          m_h_mdtetavsSA_resptLargeSpecialminus15in.at(i)->Write();
          m_h_highoffetavsSA_resptLargeSpecialplus11out.at(i)->Write();
          m_h_highoffetavsSA_resptLargeSpecialplus11in.at(i)->Write();
          m_h_highoffetavsSA_resptLargeSpecialplus15out.at(i)->Write();
          m_h_highoffetavsSA_resptLargeSpecialplus15in.at(i)->Write();
          m_h_highoffetavsSA_resptLargeSpecialminus11out.at(i)->Write();
          m_h_highoffetavsSA_resptLargeSpecialminus11in.at(i)->Write();
          m_h_highoffetavsSA_resptLargeSpecialminus15out.at(i)->Write();
          m_h_highoffetavsSA_resptLargeSpecialminus15in.at(i)->Write();
          m_h_etaIndexvsSA_respt.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeplus.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialplus.at(i)->Write();
          m_h_etaIndexvsSA_resptSmallplus.at(i)->Write();
          m_h_etaIndexvsSA_resptSmallSpecialplus.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeminus.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialminus.at(i)->Write();
          m_h_etaIndexvsSA_resptSmallminus.at(i)->Write();
          m_h_etaIndexvsSA_resptSmallSpecialminus.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialplus11out.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialplus11in.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialplus15out.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialplus15in.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialminus11out.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialminus11in.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialminus15out.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialminus15in.at(i)->Write();

         cout<<i*m_thpitch<<"      "<<m_countLarge.at(i)<<"      "<<m_countLargeSpecial.at(i)<<"      "<<m_countSmall.at(i)<<"      "<<m_countSmallSpecial.at(i)<<endl;
     }
     
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
          m_h_poffvsSA_pt.clear();
          m_h_offphivsSA_sAddress.clear();
          m_h_offphivsSAphims.clear();
          m_h_rpchitXY.clear();
          m_h_rpchitXYLargeSpecial.clear();
          m_h_rpchitXYLargeSpecialplus11out.clear();
          m_h_rpchitXYLargeSpecialplus11in.clear();
          m_h_rpchitXYLargeSpecialplus15out.clear();
          m_h_rpchitXYLargeSpecialplus15in.clear();
          m_h_rpchitXYLargeSpecialminus11out.clear();
          m_h_rpchitXYLargeSpecialminus11in.clear();
          m_h_rpchitXYLargeSpecialminus15out.clear();
          m_h_rpchitXYLargeSpecialminus15in.clear();
          m_h_mdthitXY.clear();
          m_h_mdthitXYLargeSpecial.clear();
          m_h_mdthitXYLargeSpecialplus11out.clear();
          m_h_mdthitXYLargeSpecialplus11in.clear();
          m_h_mdthitXYLargeSpecialplus15out.clear();
          m_h_mdthitXYLargeSpecialplus15in.clear();
          m_h_mdthitXYLargeSpecialminus11out.clear();
          m_h_mdthitXYLargeSpecialminus11in.clear();
          m_h_mdthitXYLargeSpecialminus15out.clear();
          m_h_mdthitXYLargeSpecialminus15in.clear();
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
          m_h_eoff_pt_Largeplus.clear();
          m_h_eL1_pt_Largeplus.clear();
          m_h_eSA_pt_Largeplus.clear();
          m_h_eoff_pt_Largeminus.clear();
          m_h_eL1_pt_Largeminus.clear();
          m_h_eSA_pt_Largeminus.clear();
          m_h_eoff_pt_LargeSpecialplus.clear();
          m_h_eL1_pt_LargeSpecialplus.clear();
          m_h_eSA_pt_LargeSpecialplus.clear();
          m_h_eL1_pt_LargeSpecialplus11.clear();
          m_h_eSA_pt_LargeSpecialplus11.clear();
          m_h_eL1_pt_LargeSpecialplus15.clear();
          m_h_eSA_pt_LargeSpecialplus15.clear();
          m_h_eL1_pt_LargeSpecialplus11in.clear();
          m_h_eSA_pt_LargeSpecialplus11in.clear();
          m_h_eL1_pt_LargeSpecialplus11out.clear();
          m_h_eSA_pt_LargeSpecialplus11out.clear();
          m_h_eL1_pt_LargeSpecialplus15out.clear();
          m_h_eSA_pt_LargeSpecialplus15out.clear();
          m_h_eL1_pt_LargeSpecialplus15in.clear();
          m_h_eSA_pt_LargeSpecialplus15in.clear();
          m_h_eoff_pt_LargeSpecialminus.clear();
          m_h_eL1_pt_LargeSpecialminus.clear();
          m_h_eSA_pt_LargeSpecialminus.clear();
          m_h_eL1_pt_LargeSpecialminus11.clear();
          m_h_eSA_pt_LargeSpecialminus11.clear();
          m_h_eL1_pt_LargeSpecialminus15.clear();
          m_h_eSA_pt_LargeSpecialminus15.clear();
          m_h_eL1_pt_LargeSpecialminus11in.clear();
          m_h_eSA_pt_LargeSpecialminus11in.clear();
          m_h_eL1_pt_LargeSpecialminus11out.clear();
          m_h_eSA_pt_LargeSpecialminus11out.clear();
          m_h_eL1_pt_LargeSpecialminus15out.clear();
          m_h_eSA_pt_LargeSpecialminus15out.clear();
          m_h_eL1_pt_LargeSpecialminus15in.clear();
          m_h_eSA_pt_LargeSpecialminus15in.clear();
          m_h_eoff_pt_Smallplus.clear();
          m_h_eL1_pt_Smallplus.clear();
          m_h_eSA_pt_Smallplus.clear();
          m_h_eoff_pt_Smallminus.clear();
          m_h_eL1_pt_Smallminus.clear();
          m_h_eSA_pt_Smallminus.clear();
          m_h_eoff_pt_SmallSpecialplus.clear();
          m_h_eL1_pt_SmallSpecialplus.clear();
          m_h_eSA_pt_SmallSpecialplus.clear();
          m_h_eoff_pt_SmallSpecialminus.clear();
          m_h_eL1_pt_SmallSpecialminus.clear();
          m_h_eSA_pt_SmallSpecialminus.clear();
          m_h_eff_poff_etaphi.clear();
          m_h_eff_pL1_etaphi.clear();
          m_h_eff_pSA_etaphi.clear();
          m_h_pSA_respt.clear();
          m_h_pCB_respt.clear();
          m_h_pEF_respt.clear();
          m_h_SA_resptLargeplus.clear();
          m_h_SA_resptLargeSpecialplus.clear();
          m_h_SA_resptSmallplus.clear();
          m_h_SA_resptSmallSpecialplus.clear();
          m_h_SA_resptLargeminus.clear();
          m_h_SA_resptLargeSpecialminus.clear();
          m_h_SA_resptSmallminus.clear();
          m_h_SA_resptSmallSpecialminus.clear();
          m_h_SA_resptLargeSpecialplus11out.clear();
          m_h_SA_resptLargeSpecialplus11in.clear();
          m_h_SA_resptLargeSpecialplus15out.clear();
          m_h_SA_resptLargeSpecialplus15in.clear();
          m_h_SA_resptLargeSpecialminus11out.clear();
          m_h_SA_resptLargeSpecialminus11in.clear();
          m_h_SA_resptLargeSpecialminus15out.clear();
          m_h_SA_resptLargeSpecialminus15in.clear();
          m_h_off_ptvsSA_resptLargeplus.clear();
          m_h_off_ptvsSA_resptLargeSpecialplus.clear();
          m_h_off_ptvsSA_resptSmallplus.clear();
          m_h_off_ptvsSA_resptSmallSpecialplus.clear();
          m_h_off_ptvsSA_resptLargeminus.clear();
          m_h_off_ptvsSA_resptLargeSpecialminus.clear();
          m_h_off_ptvsSA_resptSmallminus.clear();
          m_h_off_ptvsSA_resptSmallSpecialminus.clear();
          m_h_offphivsSA_resptLargeplus.clear();
          m_h_offphivsSA_resptLargeSpecialplus.clear();
          m_h_offphivsSA_resptSmallplus.clear();
          m_h_offphivsSA_resptSmallSpecialplus.clear();
          m_h_offphivsSA_resptLargeminus.clear();
          m_h_offphivsSA_resptLargeSpecialminus.clear();
          m_h_offphivsSA_resptSmallminus.clear();
          m_h_offphivsSA_resptSmallSpecialminus.clear();
          m_h_off_ptvsSA_resptLargeSpecialplus15out.clear();
          m_h_off_ptvsSA_resptLargeSpecialplus15in.clear();
          m_h_off_ptvsSA_resptLargeSpecialplus11out.clear();
          m_h_off_ptvsSA_resptLargeSpecialplus11in.clear();
          m_h_off_ptvsSA_resptLargeSpecialminus15out.clear();
          m_h_off_ptvsSA_resptLargeSpecialminus15in.clear();
          m_h_off_ptvsSA_resptLargeSpecialminus11out.clear();
          m_h_off_ptvsSA_resptLargeSpecialminus11in.clear();
          m_h_offetavsSA_resptLargeplus.clear();
          m_h_offetavsSA_resptLargeSpecialplus.clear();
          m_h_offetavsSA_resptSmallplus.clear();
          m_h_offetavsSA_resptSmallSpecialplus.clear();
          m_h_offetavsSA_resptLargeminus.clear();
          m_h_offetavsSA_resptLargeSpecialminus.clear();
          m_h_offetavsSA_resptSmallminus.clear();
          m_h_offetavsSA_resptSmallSpecialminus.clear();
          m_h_offetavsSA_resptLargeSpecialplus11out.clear();
          m_h_offetavsSA_resptLargeSpecialplus11in.clear();
          m_h_offetavsSA_resptLargeSpecialplus15out.clear();
          m_h_offetavsSA_resptLargeSpecialplus15in.clear();
          m_h_offetavsSA_resptLargeSpecialminus11out.clear();
          m_h_offetavsSA_resptLargeSpecialminus11in.clear();
          m_h_offetavsSA_resptLargeSpecialminus15out.clear();
          m_h_offetavsSA_resptLargeSpecialminus15in.clear();
          m_h_highoffetavsSA_resptLargeSpecialplus11out.clear();
          m_h_highoffetavsSA_resptLargeSpecialplus11in.clear();
          m_h_highoffetavsSA_resptLargeSpecialplus15out.clear();
          m_h_highoffetavsSA_resptLargeSpecialplus15in.clear();
          m_h_highoffetavsSA_resptLargeSpecialminus11out.clear();
          m_h_highoffetavsSA_resptLargeSpecialminus11in.clear();
          m_h_highoffetavsSA_resptLargeSpecialminus15out.clear();
          m_h_highoffetavsSA_resptLargeSpecialminus15in.clear();
          m_h_etaIndexvsSA_respt.clear();
          m_h_etaIndexvsSA_resptLargeplus.clear();
          m_h_etaIndexvsSA_resptLargeSpecialplus.clear();
          m_h_etaIndexvsSA_resptSmallplus.clear();
          m_h_etaIndexvsSA_resptSmallSpecialplus.clear();
          m_h_etaIndexvsSA_resptLargeminus.clear();
          m_h_etaIndexvsSA_resptLargeSpecialminus.clear();
          m_h_etaIndexvsSA_resptSmallminus.clear();
          m_h_etaIndexvsSA_resptSmallSpecialminus.clear();
          m_h_etaIndexvsSA_resptLargeSpecialplus11out.clear();
          m_h_etaIndexvsSA_resptLargeSpecialplus11in.clear();
          m_h_etaIndexvsSA_resptLargeSpecialplus15out.clear();
          m_h_etaIndexvsSA_resptLargeSpecialplus15in.clear();
          m_h_etaIndexvsSA_resptLargeSpecialminus11out.clear();
          m_h_etaIndexvsSA_resptLargeSpecialminus11in.clear();
          m_h_etaIndexvsSA_resptLargeSpecialminus15out.clear();
          m_h_etaIndexvsSA_resptLargeSpecialminus15in.clear();


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
     m_pSA_rpcX->clear();
     m_pSA_rpcY->clear();
     m_pSA_mdtZ->clear();
     m_pSA_mdtR->clear();
     m_pSA_mdtPhi->clear();
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
