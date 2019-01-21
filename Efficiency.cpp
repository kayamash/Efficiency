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
#include <FillTree.cpp>

void Efficiency::Init(std::string name,Double_t req,Int_t max,const Int_t nh,Int_t proc,string outputfilename){
     	m_method_name = name;
     	m_binmax = max;
     	m_nthreshold = nh;
          m_proc = proc;
    
     	//initialize
          m_reqL1dR = req;
          m_ft.Init(outputfilename);
     }

 bool Efficiency::DicisionBarrel(Double_t eta){
     if(std::fabs(eta) <= 1.05){
          return kTRUE;
     }else{
          return kFALSE;
     }
}

bool Efficiency::CutTagProbe(Int_t pass){
     if(m_sumReqdRL1 < m_tp_extdR && 0.2 < m_tp_extdR && m_sumReqdREF < m_tp_dR && pass > -1 && m_tag_proc == m_proc){
          return kTRUE;
     }else{
          return kFALSE;
     }
}

bool Efficiency::CutL1(Int_t pass){
  if(pass > -1){
	return kTRUE;
     }else{
     return kFALSE;
     }
}

bool Efficiency::CutSA(Int_t pass,Double_t pt,Double_t th){
     if(pass == 1 && std::fabs(pt) > th){
          return kTRUE;
     }else{
          return kFALSE;
     }
}

bool Efficiency::CutCB(Int_t pass){
     if(pass == 1){
          return kTRUE;
     }else{
          return kFALSE;
     }
}

bool Efficiency::CutEF(Int_t pass){
     if(pass == 1){
          return kTRUE;
     }else{
          return kFALSE;
     }
}

int Efficiency::DicisionArea(Double_t roiphi){
     Int_t dicisionarea = 0;
     if(roiphi >= -2.6 && roiphi < -2.4){
          dicisionarea = 1;
     }else if(roiphi >= -2.4 && roiphi < -2.0){
          dicisionarea = 2;
     }else if(roiphi >= -0.8 && roiphi < -0.6){
          dicisionarea = 3;
     }else if(roiphi >= -1.0 && roiphi < -0.8){
          dicisionarea = 4;
     }
     return dicisionarea;
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
          Double_t pSA_sAddress = -1;
          Double_t pSA_phims = -99999;
          Double_t pSA_phibe = -99999;
          float pSA_roiphi = -99999;
          vector<float> *pSA_rpcX = 0;
          vector<float> *pSA_rpcY = 0;
          vector<float> *pSA_rpcZ = 0;
          vector<float> *pSA_rpcR = 0;
          vector<float> *pSA_mdtZ = 0;
          vector<float> *pSA_mdtR = 0;
          vector<float> *pSA_mdtPhi = 0;
          vector<Int_t> *pSA_mdthitChamber = 0;
          Double_t pSA_superpointZ_BI = -99999;
          Double_t pSA_superpointZ_BM = -99999;
          Double_t pSA_superpointZ_BO = -99999;
          Double_t pSA_superpointZ_BME = -99999;
          Double_t pSA_superpointR_BI = -99999;
          Double_t pSA_superpointR_BM = -99999;
          Double_t pSA_superpointR_BO = -99999;
          Double_t pSA_superpointR_BME = -99999;
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
          Double_t segment_parameter[] = {0,0,0,0};
          Int_t areanumber = 0;

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
                    pSA_phibe = m_pSA_phibe->at(method);
                    pSA_roiphi = m_pSA_roiphi->at(method);
                    pSA_rpcX = &(m_pSA_rpcX->at(method));
                    pSA_rpcY = &(m_pSA_rpcY->at(method));
                    pSA_rpcZ = &(m_pSA_rpcZ->at(method));
                    pSA_rpcR = &(m_pSA_rpcR->at(method));
                    pSA_mdtZ = &(m_pSA_mdtZ->at(method));
                    pSA_mdtR = &(m_pSA_mdtR->at(method));
                    pSA_mdtPhi = &(m_pSA_mdtPhi->at(method));
                    pSA_mdthitChamber = &(m_pSA_mdthitChamber->at(method));
                    pSA_superpointZ_BI = m_pSA_superpointZ_BI->at(method);
                    pSA_superpointZ_BM = m_pSA_superpointZ_BM->at(method);
                    pSA_superpointZ_BO = m_pSA_superpointZ_BO->at(method);
                    pSA_superpointZ_BME = m_pSA_superpointZ_BME->at(method);
                    pSA_superpointR_BI = m_pSA_superpointR_BI->at(method);
                    pSA_superpointR_BM = m_pSA_superpointR_BM->at(method);
                    pSA_superpointR_BO = m_pSA_superpointR_BO->at(method);
                    pSA_superpointR_BME = m_pSA_superpointR_BME->at(method);
               }
          }

          tL1_dR = TMath::Sqrt(pow(m_tL1_eta - m_toff_eta,2) + pow(m_tL1_phi - m_toff_phi,2) );
          tEF_dR = TMath::Sqrt(pow(m_tEF_eta - m_toff_eta,2) + pow(m_tEF_phi - m_toff_phi,2) );
          if(std::fabs(m_toff_pt)*0.001 < 10.0)m_reqL1dR = -0.00001*std::fabs(m_toff_pt) + 0.18;
          //Tag and Probe
          if(!CutTagProbe(pEFTAG_pass))continue;
          //L1
          Double_t textL1_dR = TMath::Sqrt(pow(m_tL1_eta - m_toff_exteta,2) + pow(m_tL1_phi - m_toff_extphi,2));
          pextL1_dR = TMath::Sqrt(pow(pL1_eta - m_poff_exteta,2) + pow(pL1_phi - m_poff_extphi,2));
          //SA
          Double_t textSA_dR = TMath::Sqrt(pow(m_tSA_eta - m_toff_exteta,2) + pow(m_tSA_phi - m_toff_extphi,2));
          pextSA_dR = TMath::Sqrt(pow(pSA_eta - m_poff_exteta,2) + pow(pSA_phi - m_poff_extphi,2));
          Double_t resSA_pt = std::fabs(m_poff_pt*0.001)/std::fabs(pSA_pt) - 1.0;
          Double_t buf_pSA_dR = TMath::Sqrt(pow(pSA_eta - m_poff_eta,2) + pow(pSA_phi - m_poff_phi,2));
          areanumber = DicisionArea(pSA_roiphi);
          //CB
          Double_t textCB_dR = TMath::Sqrt(pow(m_tCB_eta - m_toff_exteta,2) + pow(m_tCB_phi - m_toff_extphi,2));
          pextCB_dR = TMath::Sqrt(pow(pCB_eta - m_poff_exteta,2) + pow(pCB_phi - m_poff_extphi,2));
          Double_t resCB_pt = std::fabs(m_poff_pt)/std::fabs(pCB_pt) - 1.0;
          //EF
          Double_t textEF_dR = TMath::Sqrt(pow(m_tEF_eta - m_toff_exteta,2) + pow(m_tEF_phi - m_toff_extphi,2));
          pextEF_dR = TMath::Sqrt(pow(pEF_eta - m_poff_exteta,2) + pow(pEF_phi - m_poff_extphi,2));
          Double_t resEF_pt = std::fabs(m_poff_pt)/std::fabs(pEF_pt) - 1.0;
          DicisionBarrel(m_poff_eta);

          std::vector<Int_t> v_L1pass;
          std::vector<Int_t> v_SApass;
          std::vector<Int_t> v_CBpass;
          std::vector<Int_t> v_EFpass;
          for(Int_t i = 0;i < m_nthreshold;i++){
               v_L1pass.push_back(CutL1(pL1_pass));
               v_SApass.push_back(CutSA(pSA_pass,pSA_pt,i*m_thpitch + m_thmin));
               v_CBpass.push_back(CutCB(pCB_pass));
               v_EFpass.push_back(CutEF(pEF_pass));
          }
          
          vector<Double_t> v_offline_kinematic;//pt,eta,phi,charge
          vector<Double_t> v_L1_kinematic;//pt,eta,phi,dR
          vector<Double_t> v_SA_kinematic;//pt,eta,phi,dR,residual,roiphi
          vector<Double_t> v_CB_kinematic;//pt,eta,phi,dR,residual
          vector<Double_t> v_EF_kinematic;//pt,eta,phi,dR,residual
          vector<Int_t> v_segment_etaIndex;
          vector<Int_t> v_segment_chamberIndex;
          vector<Int_t> v_segment_nPrecisionHits;
          vector<Int_t> v_segment_Sector;
          vector<Double_t> v_segment_x;
          vector<Double_t> v_segment_y;
          vector<Double_t> v_segment_z;
          vector<Double_t> v_SP;//z BI,BM,BO,BME R BI,BM,BO,BME
          vector<float> v_RPC_x;
          vector<float> v_RPC_y;
          vector<float> v_RPC_z;
          vector<float> v_RPC_R;
          vector<float> v_MDT_phi;
          vector<float> v_MDT_z;
          vector<float> v_MDT_R;
          vector<Int_t> v_MDT_chamber;
          vector<Int_t> v_information;//sAddress,areanumber,Qeta/|eta|,the number of stations



          m_ft.Fill();



}//Execute

void Efficiency::Finalize(){
     cout<<"completed!"<<endl;
}
