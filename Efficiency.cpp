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
#include "CalcEff.cpp"

//const Double_t pt_threshold[4] = {3.38,1.25,3.17,3.41};//MU4
const Double_t pt_threshold[4] = {15.87,10.73,12.21,15.87};//MU20

void Efficiency::Init(std::string name,const Int_t np,const Int_t ne,const Double_t mp,const Double_t me,Double_t req,Int_t max,Double_t err,const Int_t nh,Int_t proc){
    	     m_nbin_phi = np;
     	m_nbin_eta = ne;
     	m_phi_max = mp;
     	m_eta_max = me;
     	m_method_name = name;
     	m_binmax = max;
     	m_efficiency_xerr = err;
     	m_nhist = nh;
          m_proc = proc;
    
     	//initialize
          m_reqL1dR = req;
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
        //cout<<"passed L1"<<endl;
	return kTRUE;
     }else{
     return kFALSE;
     }
}

bool Efficiency::CutSA(Int_t pass,Double_t pt,Double_t th){
     if(pass == 1 && std::fabs(pt) > th){
          //cout<<"passed SA"<<endl;
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
     if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) == 1){
          dicisionarea = 1;
     }else if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) == -1){
          dicisionarea = 5;
     }else{
          return 0;
     }
     if(roiphi >= -2.6 && roiphi < -2.4){
          dicisionarea += 0;
     }else if(roiphi >= -2.4 && roiphi < -2.0){
          dicisionarea += 1;
     }else if(roiphi >= -0.8 && roiphi < -0.6){
          dicisionarea += 2;
     }else if(roiphi >= -1.0 && roiphi < -0.8){
          dicisionarea += 3;
     }else{
          return 0;
     }
     return dicisionarea;
}

int Efficiency::EtaDistribution(){
     if(std::fabs(m_poff_eta) < 1.05){
          return 0;
     }else if(std::fabs(m_poff_eta) < 1.5){
          return 1;
     }else if(std::fabs(m_poff_eta) < 2.0){
          return 2;
     }else{
          return 3;
     }
}

void Efficiency::Execute(Int_t ev){
     tChain->GetEntry(ev);
     for(Int_t i = 0;i < m_nhist;i++){
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
          Double_t pSA_ptTGC = -99999;
          Double_t pSA_ptalpha = -99999;
          Double_t pSA_ptbeta = -99999;
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
          Int_t numSP = 0;
          Int_t patternSP = 0;//inner+middle=3,inner+outer=4,middle+outer=5
          Double_t overphi = TMath::Pi();

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
                    pSA_ptTGC = m_pSA_pttgc->at(method);
                    pSA_ptalpha = m_pSA_ptalpha->at(method);
                    pSA_ptbeta = m_pSA_ptbeta->at(method);
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

          Int_t pt_method = -1;//alpha = 0,beta = 1,TGC = 2,alpha & beta = 3
          if(pSA_pt == pSA_ptalpha){
               pt_method = 0;
          }else if(pSA_pt == pSA_ptbeta){
               pt_method = 1;
          }else if(pSA_pt == pSA_ptTGC){
               pt_method = 2;
          }else{
               pt_method = 3;
          }
          //cout<<pt_method<<endl;
          Double_t resptalpha = std::fabs(m_poff_pt*0.001)/std::fabs(pSA_ptalpha) - 1.0;
          Double_t resptbeta = std::fabs(m_poff_pt*0.001)/std::fabs(pSA_ptbeta) - 1.0;
          Double_t respttgc = std::fabs(m_poff_pt*0.001)/std::fabs(pSA_ptTGC) - 1.0;

          for(Int_t index = 0;index < 10;index++){
               if(m_probe_segment_chamberIndex[index] == 1 && (m_probe_segment_sector[index] == 11 || m_probe_segment_sector[index] == 15) )m_h_BIMrvsx.at(i)->Fill(sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)),m_probe_segment_x[index]);
          }

          tL1_dR = TMath::Sqrt(pow(m_tL1_eta - m_toff_eta,2) + pow(m_tL1_phi - m_toff_phi,2) );
          tEF_dR = TMath::Sqrt(pow(m_tEF_eta - m_toff_eta,2) + pow(m_tEF_phi - m_toff_phi,2) );
          if(std::fabs(m_toff_pt)*0.001 < 10.0)m_reqL1dR = -0.00001*std::fabs(m_toff_pt) + 0.18;
          Int_t SPinner = 0;
          Int_t SPmiddle = 0;
          Int_t SPouter = 0;
          if(pSA_superpointZ_BI != 0 && pSA_superpointZ_BI != -99999 && pSA_superpointR_BI != 0 && pSA_superpointR_BI != -99999){
               numSP++;
               patternSP += 1;
               SPinner = 1;
          }
          if(pSA_superpointZ_BM != 0 && pSA_superpointZ_BM != -99999 && pSA_superpointR_BM != 0 && pSA_superpointR_BM != -99999){
               numSP++;
               patternSP += 2;
               SPmiddle = 1;
          }
          if(pSA_superpointZ_BO != 0 && pSA_superpointZ_BO != -99999 && pSA_superpointR_BO != 0 && pSA_superpointR_BO != -99999){
               numSP++;
               patternSP += 3;
               SPouter = 1;
          }
          if(SPinner == 1 && SPmiddle == 1 && SPouter == 0)patternSP = 3;
          if(SPinner == 1 && SPmiddle == 0 && SPouter == 1)patternSP = 4;
          if(SPinner == 0 && SPmiddle == 1 && SPouter == 1)patternSP = 5;
          Int_t decision_noBIM = 0;
          for(Int_t index = 0;index < 10;index++){
               if(m_probe_segment_chamberIndex[index] == 1 && sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)) > 5800 && fabs(m_probe_segment_x[index]) > 4000.)decision_noBIM++;
               //if(m_probe_segment_chamberIndex[index] == 1 && (m_probe_segment_sector[index] == 11 || m_probe_segment_sector[index] == 15) && ((pSA_roiphi > -0.8 && pSA_roiphi < -0.6) || (pSA_roiphi > -2.6 && pSA_roiphi < -2.4)) )decision_noBIM++;
          }

          //offline
          if(!CutTagProbe(pEFTAG_pass))continue;
          if(i == 0 && static_cast<Int_t>(pSA_sAddress) == 1)m_h_offphi_LargeSpecial->Fill(m_poff_phi);
          m_h_poff_pt.at(i)->Fill(m_poff_pt*0.001);
          m_h_eoff_pt.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          //if(std::fabs(m_poff_pt*0.001) > 40){
          if(std::fabs(m_poff_pt*0.001) > 8){
               m_h_eoff_eta.at(i)->Fill(m_poff_eta);
               m_h_eoff_phi.at(i)->Fill(m_poff_phi);
               //m_h_eoff_aipc.at(i)->Fill(m_aipc);
          }
          if(DicisionBarrel(m_poff_eta)){
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
          Int_t nosector9 = 0;
          for(Int_t index = 0;index < 10;index++){
               Int_t nsector = m_probe_segment_sector[index] - 1;
               //cout<<nsector<<endl;
               if(std::fabs(nsector) < 16 && m_probe_segment_chamberIndex[index] >= 0 && m_probe_segment_chamberIndex[index] <= 5)m_h_sectorphi[nsector]->Fill(TMath::ATan2(m_probe_segment_y[index],m_probe_segment_x[index]));
               if(m_probe_segment_sector[index] == 9)nosector9++;
          }

          //L1
          if(!CutL1(pL1_pass))continue;
          Double_t textL1_dR = TMath::Sqrt(pow(m_tL1_eta - m_toff_exteta,2) + pow(m_tL1_phi - m_toff_extphi,2));
          pextL1_dR = TMath::Sqrt(pow(pL1_eta - m_poff_exteta,2) + pow(pL1_phi - m_poff_extphi,2));
          m_h_pL1_pt.at(i)->Fill(std::fabs(pL1_pt*0.001));
          m_h_pL1_dR.at(i)->Fill(pL1_dR);
          m_h_textL1_dR.at(i)->Fill(textL1_dR);
          m_h_pextL1_dR.at(i)->Fill(pextL1_dR);
          m_h_eL1_pt.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          m_h_L1pSA_sAddress->Fill(pSA_sAddress);
          if(DicisionBarrel(m_poff_eta)){
               m_h_eL1_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(decision_noBIM == 0)m_h_eL1_pt_BarrelwithoutBIM.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               m_h_eL1_pt_BarrelincBIM.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(pSA_sAddress == 0 && nosector9 == 0)m_h_eL1_pt_Largenormal.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 1)m_h_eL1SP1_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 2)m_h_eL1SP2_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 3)m_h_eL1SP3_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(patternSP == 3)m_h_eL1innmid_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          }else{
               m_h_eL1_pt_end.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 1)m_h_eL1SP1_pt_endcap.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 2)m_h_eL1SP2_pt_endcap.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 3)m_h_eL1SP3_pt_endcap.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(patternSP == 3)m_h_eL1innmid_pt_endcap.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          }

          //if(std::fabs(m_poff_pt*0.001) > 40){//plateau cut
          if(std::fabs(m_poff_pt*0.001) > 8){
               m_h_eff_pL1_etaphi.at(i)->Fill(m_poff_eta,m_poff_phi);
               m_h_eL1_eta.at(i)->Fill(m_poff_eta);
               m_h_eL1_phi.at(i)->Fill(m_poff_phi);
               //m_h_eL1_aipc.at(i)->Fill(m_aipc);
          }

          areanumber = DicisionArea(pSA_roiphi);

          switch(static_cast<Int_t>(pSA_sAddress)){//switch Large ,LS , Small ,SS
               case 0:
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0)m_h_eL1_pt_Largeplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));//Qeta = +1
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) < 0)m_h_eL1_pt_Largeminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));//Qeta = -1
                    break;
               case 1:
                    if(areanumber > 0 && areanumber < 5)m_h_eL1_pt_LargeSpecialplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));//Qeta = +1
                    if(areanumber > 4 && areanumber < 9)m_h_eL1_pt_LargeSpecialminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));//Qeta = -1
                    if(decision_noBIM == 0)m_h_eL1_pt_LSwithoutBIM.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                    m_h_eL1_pt_LSincBIM.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                    switch(areanumber){
                         case 1://plus11out
                              m_h_eL1_pt_LargeSpecialplus11out.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_eL1_pt_LargeSpecialplus11.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                         case 2://plus11in
                              m_h_eL1_pt_LargeSpecialplus11in.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_eL1_pt_LargeSpecialplus11.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                         case 3://plus15out
                              m_h_eL1_pt_LargeSpecialplus15out.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_eL1_pt_LargeSpecialplus15.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                         case 4://plus15in
                              m_h_eL1_pt_LargeSpecialplus15in.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_eL1_pt_LargeSpecialplus15.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                         case 5://minus11out
                              m_h_eL1_pt_LargeSpecialminus11out.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_eL1_pt_LargeSpecialminus11.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                         case 6://minus11in
                              m_h_eL1_pt_LargeSpecialminus11in.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_eL1_pt_LargeSpecialminus11.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                         case 7://minus15out
                              m_h_eL1_pt_LargeSpecialminus15out.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_eL1_pt_LargeSpecialminus15.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                         case 8://minus15in
                              m_h_eL1_pt_LargeSpecialminus15in.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_eL1_pt_LargeSpecialminus15.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                         default:
                              break;
                    }
                    break;
               case 2:
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0)m_h_eL1_pt_Smallplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));//Qeta = +1
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) < 0)m_h_eL1_pt_Smallminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));//Qeta = -1
                    break;
               case 3:
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0)m_h_eL1_pt_SmallSpecialplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));//Qeta = +1
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) < 0)m_h_eL1_pt_SmallSpecialminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));//Qeta = -1
                    break;
               default:
                    break;
          }

          //SA
          if(!CutSA(pSA_pass,pSA_pt,i*m_thpitch + m_thmin))continue;
          m_h_countSA.at(i)->Fill(m_poff_eta);
          m_h_numSP.at(i)->Fill(numSP);
          Double_t textSA_dR = TMath::Sqrt(pow(m_tSA_eta - m_toff_exteta,2) + pow(m_tSA_phi - m_toff_extphi,2));
          pextSA_dR = TMath::Sqrt(pow(pSA_eta - m_poff_exteta,2) + pow(pSA_phi - m_poff_extphi,2));
          Double_t resSA_pt = std::fabs(m_poff_pt*0.001)/std::fabs(pSA_pt) - 1.0;
          Double_t buf_pSA_dR = TMath::Sqrt(pow(pSA_eta - m_poff_eta,2) + pow(pSA_phi - m_poff_phi,2));
          Double_t buf_eta = 0;
          m_h_SApSA_sAddress->Fill(pSA_sAddress);
          
          m_h_pSA_pt.at(i)->Fill(std::fabs(pSA_pt));
          m_h_pSA_dR.at(i)->Fill(buf_pSA_dR);
          m_h_textSA_dR.at(i)->Fill(textSA_dR);
          m_h_pextSA_dR.at(i)->Fill(pextSA_dR);
          m_h_eSA_pt.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          m_h_pSA_respt.at(i)->Fill(resSA_pt);
          m_h_pSAphivspSAphims.at(i)->Fill(pSA_phi,pSA_phims);
          m_h_pSAphivspSAphibe.at(i)->Fill(pSA_phi,pSA_phibe);
          //cout<<m_poff_eta<<"    "<<pt_method<<endl;
          if(DicisionBarrel(m_poff_eta)){
               m_h_eSA_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               m_h_pSA_respt_barrel.at(i)->Fill(resSA_pt);
               if(decision_noBIM == 0)m_h_eSA_pt_BarrelwithoutBIM.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               m_h_eSA_pt_BarrelincBIM.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               for(Int_t MDTsize = 0;MDTsize < (signed int)pSA_mdthitChamber->size();MDTsize++){
                    m_h_mdtchamber.at(i)->Fill(pSA_mdthitChamber->at(MDTsize));
               }
               if(pSA_sAddress == 0 && nosector9 == 0)m_h_eSA_pt_Largenormal.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 1)m_h_eSASP1_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 2)m_h_eSASP2_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 3)m_h_eSASP3_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(patternSP == 3)m_h_eSAinnmid_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          }else{
               m_h_eSA_pt_end.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               m_h_pSA_respt_endcap.at(i)->Fill(resSA_pt);
               //cout<<pt_method<<endl;
               if(pt_method >= 0 && std::fabs(m_poff_pt*0.001) < pt_threshold[EtaDistribution()])m_h_ptmethod[pt_method]->Fill(std::fabs(m_poff_pt*0.001));
               if(pt_method >= 0 && std::fabs(m_poff_pt*0.001) > pt_threshold[EtaDistribution()])m_h_ptmethodover[pt_method]->Fill(std::fabs(m_poff_pt*0.001));
               if(std::fabs(m_poff_pt*0.001) < pt_threshold[EtaDistribution()]){
                    m_h_ptSA[0]->Fill(std::fabs(pSA_ptalpha));
                    m_h_ptSA[1]->Fill(std::fabs(pSA_ptbeta));
                    m_h_ptSA[2]->Fill(std::fabs(pSA_ptTGC));
                    if(pt_method == 3)m_h_ptSA[3]->Fill(std::fabs(pSA_pt));
                    m_h_resptSA[0]->Fill(resptalpha);
                    m_h_resptSA[1]->Fill(resptbeta);
                    m_h_resptSA[2]->Fill(respttgc);
                    if(pt_method == 3)m_h_resptSA[3]->Fill(resSA_pt);
               }
               if(numSP == 1)m_h_eSASP1_pt_endcap.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 2)m_h_eSASP2_pt_endcap.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 3)m_h_eSASP3_pt_endcap.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(patternSP == 3)m_h_eSAinnmid_pt_endcap.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          }
          for(Int_t size = 0;size < (signed int)pSA_mdtZ->size();size++){
               buf_eta += -TMath::Log((sqrt(pow(pSA_mdtZ->at(size),2) + pow(pSA_mdtR->at(size),2)) - pSA_mdtZ->at(size))/(sqrt(pow(pSA_mdtZ->at(size),2) + pow(pSA_mdtR->at(size),2)) + pow(pSA_mdtZ->at(size),2)))/2.0;
          }
          m_h_numhit.at(i)->Fill(pSA_mdtZ->size());
          for(Int_t index = 0;index < 10;index++){
               m_h_sectorvsphi.at(i)->Fill(m_probe_segment_sector[index],m_poff_phi);
               m_h_indexvseta.at(i)->Fill(m_probe_segment_etaIndex[index],m_poff_eta);
          }
          if(static_cast<Int_t>(pSA_sAddress) == 2 && m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0)m_h_eSA_pt_Smallplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          switch(static_cast<Int_t>(pSA_sAddress)){//switch Large ,LS , Small ,SS
               case 0:
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0){//Qeta = +1
                         m_h_off_ptvsSA_resptLargeplus.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptLargeplus.at(i)->Fill(resSA_pt);
                         m_h_eSA_pt_Largeplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_offphivsSA_resptLargeplus.at(i)->Fill(m_poff_phi,resSA_pt);
                         m_h_offetavsSA_resptLargeplus.at(i)->Fill(m_poff_eta,resSA_pt);

                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeplus.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         }
                    }

                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) < 0){//Qeta = -1
                         m_h_off_ptvsSA_resptLargeminus.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptLargeminus.at(i)->Fill(resSA_pt);
                         m_h_eSA_pt_Largeminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_offphivsSA_resptLargeminus.at(i)->Fill(m_poff_phi,resSA_pt);
                         m_h_offetavsSA_resptLargeminus.at(i)->Fill(m_poff_eta,resSA_pt);
                    }
                    break;
               case 1:
                    if(i == 0){
                         m_h_saphims_LargeSpecial->Fill(pSA_phims);
                         m_h_saroiphi_LargeSpecial->Fill(pSA_roiphi);
                    }
                    if(decision_noBIM == 0)m_h_eSA_pt_LSwithoutBIM.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                    m_h_eSA_pt_LSincBIM.at(i)->Fill(std::fabs(m_poff_pt*0.001));

                    for(Int_t index = 0;index < 10;index++){
                         if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialplus11out.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 2)m_h_etaIndexvsSA_resptLargeSpecialplus11out_2station.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 3)m_h_etaIndexvsSA_resptLargeSpecialplus11out_3station.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0){
                              //m_h_segmentZR_LargeSpecialplus11out.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              //if(m_probe_segment_chamberIndex[index] == 1)m_h_segmentZR_LargeSpecialplus11out_BIL.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         }
                    }

                    if(areanumber > 0 && areanumber < 5){
                         m_h_eSA_pt_LargeSpecialplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_SA_resptLargeSpecialplus.at(i)->Fill(resSA_pt);
                         m_h_off_ptvsSA_resptLargeSpecialplus.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_offphivsSA_resptLargeSpecialplus.at(i)->Fill(m_poff_phi,resSA_pt);
                         m_h_offetavsSA_resptLargeSpecialplus.at(i)->Fill(m_poff_eta,resSA_pt);
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialplus.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              //if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialplus.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         }
                         m_h_mdtSPZR_LargeSpecialplus.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_mdtSPZR_LargeSpecialplus.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_mdtSPZR_LargeSpecialplus.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                    }
                    if(areanumber > 4 && areanumber < 9){
                         m_h_eSA_pt_LargeSpecialminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_SA_resptLargeSpecialminus.at(i)->Fill(resSA_pt);
                         m_h_off_ptvsSA_resptLargeSpecialminus.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_offphivsSA_resptLargeSpecialminus.at(i)->Fill(m_poff_phi,resSA_pt);
                         m_h_offetavsSA_resptLargeSpecialminus.at(i)->Fill(m_poff_eta,resSA_pt);
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialminus.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              //if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialminus.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         }
                         m_h_mdtSPZR_LargeSpecialminus.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_mdtSPZR_LargeSpecialminus.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_mdtSPZR_LargeSpecialminus.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                    }

                    switch(areanumber){
                         case 1://plus 11out
                              m_h_mdtSPX_BI.at(i)->Fill(pSA_superpointR_BI*cos(pSA_roiphi));
                              m_h_mdtSPY_BI.at(i)->Fill(pSA_superpointR_BI*sin(pSA_roiphi));
                              m_h_off_ptvsSA_resptLargeSpecialplus11out.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              if(numSP == 2)m_h_off_ptvsSA_resptLargeSpecialplus11out_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              if(numSP == 3)m_h_off_ptvsSA_resptLargeSpecialplus11out_3station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              m_h_SA_resptLargeSpecialplus11out.at(i)->Fill(resSA_pt);
                              m_h_offetavsSA_resptLargeSpecialplus11out.at(i)->Fill(m_poff_eta,resSA_pt);
                              m_h_eSA_pt_LargeSpecialplus11out.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                   
                              m_h_mdtSPZR_LargeSpecialplus11out.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                              m_h_mdtSPZR_LargeSpecialplus11out.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                              m_h_mdtSPZR_LargeSpecialplus11out.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                              break;
                         case 2://plus 11in
                              m_h_mdtSPX_BI.at(i)->Fill(pSA_superpointR_BI*cos(pSA_roiphi));
                              m_h_mdtSPY_BI.at(i)->Fill(pSA_superpointR_BI*sin(pSA_roiphi));
                              m_h_off_ptvsSA_resptLargeSpecialplus11in.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              if(numSP == 2)m_h_off_ptvsSA_resptLargeSpecialplus11in_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              if(numSP == 3)m_h_off_ptvsSA_resptLargeSpecialplus11in_3station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              m_h_SA_resptLargeSpecialplus11in.at(i)->Fill(resSA_pt);
                              m_h_offetavsSA_resptLargeSpecialplus11in.at(i)->Fill(m_poff_eta,resSA_pt);
                              m_h_eSA_pt_LargeSpecialplus11in.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              for(Int_t index = 0;index < 10;index++){
                                   if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialplus11in.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                   if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 2)m_h_etaIndexvsSA_resptLargeSpecialplus11in_2station.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                   if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 3)m_h_etaIndexvsSA_resptLargeSpecialplus11in_3station.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0){
                                        //m_h_segmentZR_LargeSpecialplus11in.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                                        //if(m_probe_segment_chamberIndex[index] == 1)m_h_segmentZR_LargeSpecialplus11in_BIL.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                                   }
                              }
                              m_h_mdtSPZR_LargeSpecialplus11in.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                              m_h_mdtSPZR_LargeSpecialplus11in.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                              m_h_mdtSPZR_LargeSpecialplus11in.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                              break;
                         case 3://plus 15out
                              m_h_off_ptvsSA_resptLargeSpecialplus15out.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              if(numSP == 2)m_h_off_ptvsSA_resptLargeSpecialplus15out_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              if(numSP == 3)m_h_off_ptvsSA_resptLargeSpecialplus15out_3station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              m_h_SA_resptLargeSpecialplus15out.at(i)->Fill(resSA_pt);
                              m_h_offetavsSA_resptLargeSpecialplus15out.at(i)->Fill(m_poff_eta,resSA_pt);
                              m_h_eSA_pt_LargeSpecialplus15out.at(i)->Fill(std::fabs(m_poff_pt*0.001));

                              for(Int_t index = 0;index < 10;index++){
                                   if(m_probe_segment_etaIndex[index] >= -8 && m_probe_segment_etaIndex[index] <= 8)m_h_etaIndexvsSA_resptLargeSpecialplus15out.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                   if(m_probe_segment_etaIndex[index] >= -8 && m_probe_segment_etaIndex[index] <= 8 && numSP == 2)m_h_etaIndexvsSA_resptLargeSpecialplus15out_2station.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                   if(m_probe_segment_etaIndex[index] >= -8 && m_probe_segment_etaIndex[index] <= 8 && numSP == 3)m_h_etaIndexvsSA_resptLargeSpecialplus15out_3station.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0){
                                        //m_h_segmentZR_LargeSpecialplus15out.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                                        //if(m_probe_segment_chamberIndex[index] == 1)m_h_segmentZR_LargeSpecialplus15out_BIL.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                                   }
                              }
                              m_h_mdtSPZR_LargeSpecialplus15out.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                              m_h_mdtSPZR_LargeSpecialplus15out.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                              m_h_mdtSPZR_LargeSpecialplus15out.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                              break;
                         case 4://plus 15in
                              m_h_off_ptvsSA_resptLargeSpecialplus15in.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              if(numSP == 2)m_h_off_ptvsSA_resptLargeSpecialplus15in_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              if(numSP == 3)m_h_off_ptvsSA_resptLargeSpecialplus15in_3station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              m_h_SA_resptLargeSpecialplus15in.at(i)->Fill(resSA_pt);
                              m_h_offetavsSA_resptLargeSpecialplus15in.at(i)->Fill(m_poff_eta,resSA_pt);
                              m_h_eSA_pt_LargeSpecialplus15in.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              for(Int_t index = 0;index < 10;index++){
                                   if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialplus15in.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                   if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 2)m_h_etaIndexvsSA_resptLargeSpecialplus15in_2station.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                   if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 3)m_h_etaIndexvsSA_resptLargeSpecialplus15in_3station.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0){
                                        //m_h_segmentZR_LargeSpecialplus15in.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                                        //if(m_probe_segment_chamberIndex[index] == 1)m_h_segmentZR_LargeSpecialplus15in_BIL.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                                   }
                              }
                              m_h_mdtSPZR_LargeSpecialplus15in.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                              m_h_mdtSPZR_LargeSpecialplus15in.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                              m_h_mdtSPZR_LargeSpecialplus15in.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                              break;
                         case 5://minus 11out
                              m_h_off_ptvsSA_resptLargeSpecialminus11out.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              if(numSP == 2)m_h_off_ptvsSA_resptLargeSpecialminus11out_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              if(numSP == 3)m_h_off_ptvsSA_resptLargeSpecialminus11out_3station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              m_h_SA_resptLargeSpecialminus11out.at(i)->Fill(resSA_pt);
                              m_h_offetavsSA_resptLargeSpecialminus11out.at(i)->Fill(m_poff_eta,resSA_pt);
                              m_h_eSA_pt_LargeSpecialminus11out.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              for(Int_t index = 0;index < 10;index++){
                                  if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialminus11out.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                  if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 2)m_h_etaIndexvsSA_resptLargeSpecialminus11out_2station.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                  if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 3)m_h_etaIndexvsSA_resptLargeSpecialminus11out_3station.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                  //if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialminus11out.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              }
                              m_h_mdtSPZR_LargeSpecialminus11out.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                              m_h_mdtSPZR_LargeSpecialminus11out.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                              m_h_mdtSPZR_LargeSpecialminus11out.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                              break;
                         case 6://minus 11in
                              m_h_off_ptvsSA_resptLargeSpecialminus11in.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              if(numSP == 2)m_h_off_ptvsSA_resptLargeSpecialminus11in_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              if(numSP == 3)m_h_off_ptvsSA_resptLargeSpecialminus11in_3station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              m_h_SA_resptLargeSpecialminus11in.at(i)->Fill(resSA_pt);
                              m_h_offetavsSA_resptLargeSpecialminus11in.at(i)->Fill(m_poff_eta,resSA_pt);
                              m_h_eSA_pt_LargeSpecialminus11in.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              for(Int_t index = 0;index < 10;index++){
                                   if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialminus11in.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                   if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 2)m_h_etaIndexvsSA_resptLargeSpecialminus11in_2station.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                   if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 3)m_h_etaIndexvsSA_resptLargeSpecialminus11in_3station.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                   //if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialminus11in.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              }
                              m_h_mdtSPZR_LargeSpecialminus11in.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                              m_h_mdtSPZR_LargeSpecialminus11in.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                              m_h_mdtSPZR_LargeSpecialminus11in.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                              break;
                         case 7://minus 15out
                              m_h_off_ptvsSA_resptLargeSpecialminus15out.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              if(numSP == 2)m_h_off_ptvsSA_resptLargeSpecialminus15out_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              if(numSP == 3)m_h_off_ptvsSA_resptLargeSpecialminus15out_3station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              m_h_SA_resptLargeSpecialminus15out.at(i)->Fill(resSA_pt);
                              m_h_offetavsSA_resptLargeSpecialminus15out.at(i)->Fill(m_poff_eta,resSA_pt);
                              m_h_eSA_pt_LargeSpecialminus15out.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              for(Int_t index = 0;index < 10;index++){
                                  if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialminus15out.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                  if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 2)m_h_etaIndexvsSA_resptLargeSpecialminus15out_2station.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                  if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 3)m_h_etaIndexvsSA_resptLargeSpecialminus15out_3station.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                  //if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialminus15out.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              }
                              m_h_mdtSPZR_LargeSpecialminus15out.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                              m_h_mdtSPZR_LargeSpecialminus15out.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                              m_h_mdtSPZR_LargeSpecialminus15out.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                              break;
                         case 8://minus 15in
                              m_h_off_ptvsSA_resptLargeSpecialminus15in.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              if(numSP == 2)m_h_off_ptvsSA_resptLargeSpecialminus15in_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              if(numSP == 3)m_h_off_ptvsSA_resptLargeSpecialminus15in_3station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              m_h_SA_resptLargeSpecialminus15in.at(i)->Fill(resSA_pt);
                              m_h_offetavsSA_resptLargeSpecialminus15in.at(i)->Fill(m_poff_eta,resSA_pt);
                              m_h_eSA_pt_LargeSpecialminus15in.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              for(Int_t index = 0;index < 10;index++){
                                   if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialminus15in.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                   if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 2)m_h_etaIndexvsSA_resptLargeSpecialminus15in_2station.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                   if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 3)m_h_etaIndexvsSA_resptLargeSpecialminus15in_3station.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                   //if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialminus15in.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              }
                              m_h_mdtSPZR_LargeSpecialminus15in.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                              m_h_mdtSPZR_LargeSpecialminus15in.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                              m_h_mdtSPZR_LargeSpecialminus15in.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                              break;
                         default :
                              break;
                    }

               case 2:
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0){//Qeta = +1
                         m_h_off_ptvsSA_resptSmallplus.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptSmallplus.at(i)->Fill(resSA_pt);
                         //m_h_eSA_pt_Smallplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_offphivsSA_resptSmallplus.at(i)->Fill(m_poff_phi,resSA_pt);
                         m_h_offetavsSA_resptSmallplus.at(i)->Fill(m_poff_eta,resSA_pt);
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptSmallplus.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         }
                    }
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) < 0){//Qeta = -1
                         m_h_off_ptvsSA_resptSmallminus.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptSmallminus.at(i)->Fill(resSA_pt);
                         m_h_eSA_pt_Smallminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_offphivsSA_resptSmallminus.at(i)->Fill(m_poff_phi,resSA_pt);
                         m_h_offetavsSA_resptSmallminus.at(i)->Fill(m_poff_eta,resSA_pt);
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptSmallminus.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         }
                    }
                    break;
               case 3:
                    if(i == 0)m_h_saroiphi_SmallSpecial->Fill(pSA_roiphi);
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0){//Qeta = +1
                         m_h_off_ptvsSA_resptSmallSpecialplus.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptSmallSpecialplus.at(i)->Fill(resSA_pt);
                         m_h_eSA_pt_SmallSpecialplus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_offphivsSA_resptSmallSpecialplus.at(i)->Fill(m_poff_phi,resSA_pt);
                         m_h_offetavsSA_resptSmallSpecialplus.at(i)->Fill(m_poff_eta,resSA_pt);
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptSmallSpecialplus.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         }
                    }
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) < 0){//Qeta = -1
                         m_h_off_ptvsSA_resptSmallSpecialminus.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptSmallSpecialminus.at(i)->Fill(resSA_pt);
                         m_h_eSA_pt_SmallSpecialminus.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_offphivsSA_resptSmallSpecialminus.at(i)->Fill(m_poff_phi,resSA_pt);
                         m_h_offetavsSA_resptSmallSpecialminus.at(i)->Fill(m_poff_eta,resSA_pt);
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptSmallSpecialminus.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         }
                    }
                    break;
               default:
                    break;
          }


          for(Int_t mdthit = 0; mdthit < (signed int)pSA_mdtPhi->size(); mdthit++){
               m_h_mdtphi.at(i)->Fill(pSA_mdtPhi->at(mdthit));
          }

          if(numSP == 2){
               if(pSA_sAddress == 0 && std::fabs(m_poff_phi) < TMath::Pi()/8.0)m_h_off_ptvsSA_resptLarge_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
               if(pSA_sAddress == 2 && std::fabs(m_poff_phi) < TMath::Pi()/8.0)m_h_off_ptvsSA_resptSmall_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
          }

          if(pSA_sAddress == 1){
               for(Int_t mdthit = 0; mdthit < (signed int)pSA_mdtZ->size(); mdthit++){
                    m_h_mdtphi_LS.at(i)->Fill(pSA_mdtPhi->at(mdthit));
                    m_h_mdtR.at(i)->Fill(pSA_mdtR->at(mdthit));
               }
          }

          overphi = m_poff_phi;
          if(overphi < 0)overphi += 2.0*TMath::Pi();
          if(pSA_sAddress == 0){
               while(overphi > 3.0*TMath::Pi()/8.0)overphi -= TMath::Pi()/4.0;
               if(overphi < TMath::Pi()/8.0)overphi += TMath::Pi()/4.0;
          }
          if(pSA_sAddress == 2){
               while(overphi > TMath::Pi()/4.0)overphi -= TMath::Pi()/4.0;
          }
          if(pSA_sAddress == 0 || pSA_sAddress == 2)m_h_overphi.at(i)->Fill(overphi);
          Int_t numeta = 0;
          Int_t numphi = 0;
          if(pSA_eta < 0)numeta = 8;
          if(pSA_sAddress == 0)overphi -= TMath::Pi()/8.0;
          numphi = overphi/TMath::Pi()*32.0 - fmod(overphi,TMath::Pi()/32.0);
          numeta += std::fabs(pSA_eta)/0.125 - fmod(std::fabs(pSA_eta),0.125);
          if(pSA_sAddress == 0 && std::fabs(pSA_eta) < 1.0)m_h_divideetaoverphi_resptLarge[numeta][numphi]->Fill(resSA_pt);
          if(pSA_sAddress == 2 && std::fabs(pSA_eta) < 1.0)m_h_divideetaoverphi_resptSmall[numeta][numphi]->Fill(resSA_pt);

          //if(std::fabs(m_poff_pt*0.001) > 40){//plateau cut
          if(std::fabs(m_poff_pt*0.001) > 8){
               m_h_eff_pSA_etaphi.at(i)->Fill(m_poff_eta,m_poff_phi);
               m_h_eSA_eta.at(i)->Fill(m_poff_eta);
               m_h_eSA_phi.at(i)->Fill(m_poff_phi);
               //m_h_eSA_aipc.at(i)->Fill(m_aipc);
          }
          m_h_poffvsSA_pt.at(i)->Fill(std::fabs(m_poff_pt*0.001),std::fabs(pSA_pt));
                         
          if(pSA_sAddress == 1){
               for(Int_t index = 0;index < 10;index++){
                    if(m_probe_segment_chamberIndex[index] == 1 && sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)) > 4000 && fabs(m_probe_segment_x[index]) < 4500.){
                         m_h_nPrecisionHits_normal.at(i)->Fill(m_probe_segment_nPrecisionHits[index]);
                         m_h_sector_normal.at(i)->Fill(m_probe_segment_sector[index]);
                         m_h_etaIndex_normal.at(i)->Fill(m_probe_segment_etaIndex[index]);
                         if(m_probe_segment_sector[index] == 11)m_h_segmentnumberBIR->Fill((Int_t)m_probe_segment_etaIndex[index] - 1,0);
                         if(m_probe_segment_sector[index] == 15)m_h_segmentnumberBIR->Fill((Int_t)m_probe_segment_etaIndex[index] - 1,1);
                    }
                    if(m_probe_segment_chamberIndex[index] == 1 && sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)) > 4000 && fabs(m_probe_segment_x[index]) > 4500.){
                         m_h_nPrecisionHits_special.at(i)->Fill(m_probe_segment_nPrecisionHits[index]);
                         m_h_sector_special.at(i)->Fill(m_probe_segment_sector[index]);
                         m_h_etaIndex_special.at(i)->Fill(m_probe_segment_etaIndex[index]);
                         if(m_probe_segment_sector[index] == 11)m_h_segmentnumberBIM->Fill((Int_t)m_probe_segment_etaIndex[index] - 1,0);
                         if(m_probe_segment_sector[index] == 15)m_h_segmentnumberBIM->Fill((Int_t)m_probe_segment_etaIndex[index] - 1,1);
                    }
               }
          }

          if(pSA_sAddress == 0 && std::fabs(m_poff_phi) < TMath::Pi()/8.0){
               if(numSP == 2)m_h_offptLarge_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 2)m_h_offetaLarge_2station.at(i)->Fill(m_poff_eta);
               if(numSP == 3)m_h_offetaLarge_3station.at(i)->Fill(m_poff_eta);
               switch(patternSP){
                    case 3:
                         if(numSP == 2){
                              m_h_offptLargeinnmid_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaLargeinnmid_2station.at(i)->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptLargeinnmid_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                    case 4:
                         if(numSP == 2){
                              m_h_offptLargeinnout_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaLargeinnout_2station.at(i)->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptLargeinnout_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                    case 5:
                         if(numSP == 2){
                              m_h_offptLargemidout_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaLargemidout_2station.at(i)->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptLargemidout_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                    default :
                         break;
               }
          }
          if(pSA_sAddress == 2 && m_poff_phi > 0 && m_poff_phi < TMath::Pi()/4.0){
               if(numSP == 2)m_h_offptSmall_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 2)m_h_offetaSmall_2station.at(i)->Fill(m_poff_eta);
               if(numSP == 3)m_h_offetaSmall_3station.at(i)->Fill(m_poff_eta);
               switch(patternSP){
                    case 3:
                         if(numSP == 2){
                              m_h_offptSmallinnmid_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaSmallinnmid_2station.at(i)->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptSmallinnmid_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                    case 4:
                         if(numSP == 2){
                              m_h_offptSmallinnout_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaSmallinnout_2station.at(i)->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptSmallinnout_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                    case 5:
                         if(numSP == 2){
                              m_h_offptSmallmidout_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaSmallmidout_2station.at(i)->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptSmallmidout_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                    default :
                         break;
               }
          }
          if(areanumber % 2 == 0 && areanumber != 0){
               if(numSP == 2)m_h_offptBIR_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 2)m_h_offetaBIR_2station.at(i)->Fill(m_poff_eta);
               if(numSP == 3)m_h_offetaBIR_3station.at(i)->Fill(m_poff_eta);
               switch(patternSP){
                    case 3:
                         if(numSP == 2){
                              m_h_offptBIRinnmid_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaBIRinnmid_2station.at(i)->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptBIRinnmid_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                    case 4:
                         if(numSP == 2){
                              m_h_offptBIRinnout_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaBIRinnout_2station.at(i)->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptBIRinnout_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                    case 5:
                         if(numSP == 2){
                              m_h_offptBIRmidout_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaBIRmidout_2station.at(i)->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptBIRmidout_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                    default :
                         break;
               }
          }
          if(areanumber % 2 == 1){
               if(numSP == 2)m_h_offptBIM_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 2)m_h_offetaBIM_2station.at(i)->Fill(m_poff_eta);
               if(numSP == 3)m_h_offetaBIM_3station.at(i)->Fill(m_poff_eta);
               switch(patternSP){
                    case 3:
                         if(numSP == 2){
                              m_h_offptBIMinnmid_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaBIMinnmid_2station.at(i)->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptBIMinnmid_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                    case 4:
                         if(numSP == 2){
                              m_h_offptBIMinnout_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaBIMinnout_2station.at(i)->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptBIMinnout_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                    case 5:
                         if(numSP == 2){
                              m_h_offptBIMmidout_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaBIMmidout_2station.at(i)->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptBIMmidout_2station.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                    default :
                         break;
               }
          }

          Int_t segdisBIR = 0;
          Int_t segdisBIM = 0;
          for(Int_t index = 0;index < 10;index++){
               if(m_probe_segment_chamberIndex[index] == 3 && m_probe_segment_x[index] < 5500. && m_probe_segment_sector[index] == 15)segdisBIR = 1;
               if(m_probe_segment_chamberIndex[index] == 5 && m_probe_segment_x[index] < 7000. && m_probe_segment_sector[index] == 15 && segdisBIR != 0)segdisBIR = 2;
               if(m_probe_segment_chamberIndex[index] == 3 && m_probe_segment_x[index] > 5500. && m_probe_segment_sector[index] == 15)segdisBIM = 1;
               if(m_probe_segment_chamberIndex[index] == 5 && m_probe_segment_x[index] > 7000. && m_probe_segment_sector[index] == 15 && segdisBIM != 0)segdisBIM = 2;
          }
          for(Int_t index = 0;index < 10;index++){
               if(m_probe_segment_chamberIndex[index] == 1 && segdisBIR == 2)m_h_BIRsegment.at(i)->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
               if(m_probe_segment_chamberIndex[index] == 1 && segdisBIM == 2)m_h_BIMsegment.at(i)->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
          }

          for(Int_t chnum = 0; chnum < 10; chnum++){
               if(m_probe_segment_chamberIndex[chnum] == 0)m_h_segmentXY_BIS.at(i)->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
               if(m_probe_segment_chamberIndex[chnum] == 1)m_h_segmentXY_BIL.at(i)->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
               if(m_probe_segment_chamberIndex[chnum] == 2)m_h_segmentXY_BMS.at(i)->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
               if(m_probe_segment_chamberIndex[chnum] == 3)m_h_segmentXY_BML.at(i)->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
               if(m_probe_segment_chamberIndex[chnum] == 4)m_h_segmentXY_BOS.at(i)->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
               if(m_probe_segment_chamberIndex[chnum] == 5)m_h_segmentXY_BOL.at(i)->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
               if(m_probe_segment_chamberIndex[chnum] == 6)m_h_segmentXY_BEE.at(i)->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
               if(m_probe_segment_chamberIndex[chnum] == 7)m_h_segmentXY_EIS.at(i)->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
               if(m_probe_segment_chamberIndex[chnum] == 8)m_h_segmentXY_EIL.at(i)->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
               if(m_probe_segment_chamberIndex[chnum] == 9)m_h_segmentXY_EMS.at(i)->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
               if(m_probe_segment_chamberIndex[chnum] == 10)m_h_segmentXY_EML.at(i)->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
               if(m_probe_segment_chamberIndex[chnum] == 11)m_h_segmentXY_EOS.at(i)->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
               if(m_probe_segment_chamberIndex[chnum] == 12)m_h_segmentXY_EOL.at(i)->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
               if(m_probe_segment_chamberIndex[chnum] == 13)m_h_segmentXY_EES.at(i)->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
               if(m_probe_segment_chamberIndex[chnum] == 14)m_h_segmentXY_EEL.at(i)->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
               if(m_probe_segment_chamberIndex[chnum] == 15)m_h_segmentXY_CSS.at(i)->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
               if(m_probe_segment_chamberIndex[chnum] == 16)m_h_segmentXY_CSL.at(i)->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
               if(m_probe_segment_chamberIndex[chnum] == -1)m_h_segmentXY_unknown.at(i)->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
          }

          if(static_cast<Int_t>(pSA_sAddress) == 0 || static_cast<Int_t>(pSA_sAddress) == 1 || static_cast<Int_t>(pSA_sAddress) == 2 || static_cast<Int_t>(pSA_sAddress) == 3){
               m_h_offphivsSA_sAddress.at(i)->Fill(pSA_phims,pSA_sAddress);
               m_h_mdtSPXY_BI.at(i)->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
               m_h_mdtSPXY_BM.at(i)->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
               m_h_mdtSPXY_BO.at(i)->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
               m_h_mdtSPXY_BME.at(i)->Fill(pSA_superpointR_BME*cos(pSA_roiphi),pSA_superpointR_BME*sin(pSA_roiphi));
               m_h_mdtSPZR.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
               m_h_mdtSPZR.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
               m_h_mdtSPZR.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
               for(Int_t index = 0;index < 10;index++){
                    Int_t buf_index = static_cast<Int_t>(m_probe_segment_etaIndex[index]);
                    if(m_probe_segment_chamberIndex[index] <= 6){
                         switch(buf_index){
                              case -6:
                                   m_h_mdtSPXY_etaIndexminus6.at(i)->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus6.at(i)->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus6.at(i)->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexminus6.at(i)->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                              case -5:
                                   m_h_mdtSPXY_etaIndexminus5.at(i)->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus5.at(i)->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus5.at(i)->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexminus5.at(i)->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                              case -4:
                                   m_h_mdtSPXY_etaIndexminus4.at(i)->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus4.at(i)->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus4.at(i)->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexminus4.at(i)->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                              case -3:
                                   m_h_mdtSPXY_etaIndexminus3.at(i)->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus3.at(i)->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus3.at(i)->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexminus3.at(i)->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                              case -2:
                                   m_h_mdtSPXY_etaIndexminus2.at(i)->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus2.at(i)->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus2.at(i)->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexminus2.at(i)->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                              case -1:
                                   m_h_mdtSPXY_etaIndexminus1.at(i)->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus1.at(i)->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus1.at(i)->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexminus1.at(i)->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                              case 1:
                                   m_h_mdtSPXY_etaIndexplus1.at(i)->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus1.at(i)->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus1.at(i)->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexplus1.at(i)->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                              case 2:
                                   m_h_mdtSPXY_etaIndexplus2.at(i)->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus2.at(i)->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus2.at(i)->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexplus2.at(i)->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                              case 3:
                                   m_h_mdtSPXY_etaIndexplus3.at(i)->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus3.at(i)->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus3.at(i)->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexplus3.at(i)->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                              case 4:
                                   m_h_mdtSPXY_etaIndexplus4.at(i)->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus4.at(i)->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus4.at(i)->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexplus4.at(i)->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                              case 5:
                                   m_h_mdtSPXY_etaIndexplus5.at(i)->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus5.at(i)->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus5.at(i)->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexplus5.at(i)->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                              case 6:
                                   m_h_mdtSPXY_etaIndexplus6.at(i)->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus6.at(i)->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus6.at(i)->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexplus6.at(i)->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                              default:
                                   break;
                         }
                    }
               }
               for(Int_t size = 0;size < (signed int)pSA_mdtZ->size();size++){
                    if(pSA_mdthitChamber->at(size) == 0 || pSA_mdthitChamber->at(size) == 1 || pSA_mdthitChamber->at(size) == 2)m_h_mdthitXY.at(i)->Fill(pSA_mdtR->at(size)*cos(pSA_mdtPhi->at(size)),pSA_mdtR->at(size)*sin(pSA_mdtPhi->at(size)));
                    m_h_mdthitZR.at(i)->Fill(pSA_mdtZ->at(size),pSA_mdtR->at(size));
               }
               for(Int_t size = 0;size < (signed int)pSA_rpcX->size();size++){
                    m_h_rpchitXY.at(i)->Fill(pSA_rpcX->at(size),pSA_rpcY->at(size));
               }

               Int_t buf_numsegment = 0;
               for(Int_t index = 0;index < 10;index++){
                    if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_respt.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                    if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY.at(i)->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                    if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0){
                         //m_h_segmentZR.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         //if(m_probe_segment_chamberIndex[index] == 1)m_h_segmentZR_BIL.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         //if(m_probe_segment_chamberIndex[index] == 3)m_h_segmentZR_BML.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         buf_numsegment++;
                    }
               }
               m_h_numsegment.at(i)->Fill(buf_numsegment);
          }
          m_h_offphivsSAphims.at(i)->Fill(m_poff_phi,pSA_phims);

          //CB
          if(!CutCB(pCB_pass))continue;
          Double_t textCB_dR = TMath::Sqrt(pow(m_tCB_eta - m_toff_exteta,2) + pow(m_tCB_phi - m_toff_extphi,2));
          pextCB_dR = TMath::Sqrt(pow(pCB_eta - m_poff_exteta,2) + pow(pCB_phi - m_poff_extphi,2));
          Double_t resCB_pt = std::fabs(m_poff_pt)/std::fabs(pCB_pt) - 1.0;
          m_h_pCB_pt.at(i)->Fill(std::fabs(pCB_pt*0.001));
          m_h_pCB_dR.at(i)->Fill(pCB_dR);
          m_h_textCB_dR.at(i)->Fill(textCB_dR);
          m_h_pextCB_dR.at(i)->Fill(pextCB_dR);
          m_h_eCB_pt.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          //if(std::fabs(m_poff_pt*0.001) > 40){//plateau cut
          if(std::fabs(m_poff_pt*0.001) > 8){
               m_h_eCB_eta.at(i)->Fill(m_poff_eta);
               m_h_eCB_phi.at(i)->Fill(m_poff_phi);
               //m_h_eCB_aipc.at(i)->Fill(m_aipc);
          }
          m_h_pCB_respt.at(i)->Fill(resCB_pt);
          if(DicisionBarrel(m_poff_eta)){
               m_h_eCB_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               m_h_pCB_respt_barrel.at(i)->Fill(resCB_pt);
          }else{
               m_h_eCB_pt_end.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               m_h_pCB_respt_endcap.at(i)->Fill(resCB_pt);
          }

          //EF
          if(!CutEF(pEF_pass))continue;
          Double_t textEF_dR = TMath::Sqrt(pow(m_tEF_eta - m_toff_exteta,2) + pow(m_tEF_phi - m_toff_extphi,2));
          pextEF_dR = TMath::Sqrt(pow(pEF_eta - m_poff_exteta,2) + pow(pEF_phi - m_poff_extphi,2));
          Double_t resEF_pt = std::fabs(m_poff_pt)/std::fabs(pEF_pt) - 1.0;
          m_h_pEF_pt.at(i)->Fill(std::fabs(pEF_pt*0.001));
          m_h_pEF_dR.at(i)->Fill(pEF_dR);
          m_h_textEF_dR.at(i)->Fill(textEF_dR);
          m_h_pextEF_dR.at(i)->Fill(pextEF_dR);
          m_h_eEF_pt.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          //if(std::fabs(m_poff_pt*0.001) > 40){//plateau cut
          if(std::fabs(m_poff_pt*0.001) > 8){
               m_h_eEF_eta.at(i)->Fill(m_poff_eta);
               m_h_eEF_phi.at(i)->Fill(m_poff_phi);
               //m_h_eEF_aipc.at(i)->Fill(m_aipc);
          }
          m_h_pEF_respt.at(i)->Fill(resEF_pt);
          if(DicisionBarrel(m_poff_eta)){
               m_h_eEF_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               m_h_pEF_respt_barrel.at(i)->Fill(resEF_pt);
          }else{
               m_h_eEF_pt_end.at(i)->Fill(std::fabs(m_poff_pt*0.001));
               m_h_pEF_respt_endcap.at(i)->Fill(resEF_pt);
          }
     }//for

}//Execute

void Efficiency::Finalize(TFile *tf1){
     CalcEff ceff;
     tf1->cd();
     //SetCondition
     //title,file title,yoffset,top margin,bottom margin,left margin,right margin
     m_h_saroiphi_LargeSpecial->Write();
     m_h_saroiphi_SmallSpecial->Write();
     m_h_saphims_LargeSpecial->Write();
     m_h_offphi_LargeSpecial->Write();
     m_h_segmentnumberBIM->Write();
     m_h_segmentnumberBIR->Write();
     for(Int_t i = 0; i < 16;i++){
          m_h_sectorphi[i]->Write();
          for(Int_t j = 0; j < 8;j++){
               m_h_divideetaoverphi_resptLarge[i][j]->Write();
               m_h_divideetaoverphi_resptSmall[i][j]->Write();
         }
     }
     for(Int_t i = 0;i < 4;i++){
          m_h_ptmethod[i]->Write();
          m_h_ptmethodover[i]->Write();
          m_h_ptSA[i]->Write();
          m_h_resptSA[i]->Write();
     }

     for(Int_t i = 0;i < m_nhist;i ++){
          //base,target
          cout<<"efficiency start!"<<endl;
          
          ceff.SetConditionName(Form("L1Efficiency_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt.at(i),m_h_eL1_pt.at(i),m_binmax,200,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiency_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt.at(i),m_h_eSA_pt.at(i),m_binmax,200,m_efficiency_xerr);
          ceff.SetConditionName(Form("CBEfficiency_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("muComb Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eSA_pt.at(i),m_h_eCB_pt.at(i),m_binmax,200,m_efficiency_xerr);
          ceff.SetConditionName(Form("EFEfficiency_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("EventFilter Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eCB_pt.at(i),m_h_eEF_pt.at(i),m_binmax,200,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1Efficiency_eta_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L1 Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyeta(m_h_eoff_eta.at(i),m_h_eL1_eta.at(i));
          ceff.SetConditionName(Form("SAEfficiency_eta_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L2MuonSA Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyeta(m_h_eL1_eta.at(i),m_h_eSA_eta.at(i));
          ceff.SetConditionName(Form("CBEfficiency_eta_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("muComb Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyeta(m_h_eSA_eta.at(i),m_h_eCB_eta.at(i));
          ceff.SetConditionName(Form("EFEfficiency_eta_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("EventFilter Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyeta(m_h_eCB_eta.at(i),m_h_eEF_eta.at(i));
          ceff.SetConditionName(Form("L1Efficiency_phi_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L1 Efficiency;offline phi;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyphi(m_h_eoff_phi.at(i),m_h_eL1_phi.at(i));
          ceff.SetConditionName(Form("SAEfficiency_phi_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L2MuonSA Efficiency;offline phi;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyphi(m_h_eL1_phi.at(i),m_h_eSA_phi.at(i));
          ceff.SetConditionName(Form("CBEfficiency_phi_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("muComb Efficiency;offline phi;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyphi(m_h_eSA_phi.at(i),m_h_eCB_phi.at(i));
          ceff.SetConditionName(Form("EFEfficiency_phi_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("EventFilter Efficiency;offline phi;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyphi(m_h_eCB_phi.at(i),m_h_eEF_phi.at(i));
          /*
          ceff.SetConditionName(Form("L1Efficiency_pileup_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L1 Efficiency;pileup;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencypileup(m_h_eoff_aipc.at(i),m_h_eL1_aipc.at(i));
          ceff.SetConditionName(Form("SAEfficiency_pileup_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L2MuonSA Efficiency;pileup;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencypileup(m_h_eL1_aipc.at(i),m_h_eSA_aipc.at(i));
          ceff.SetConditionName(Form("CBEfficiency_pileup_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("muComb Efficiency;pileup;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencypileup(m_h_eSA_aipc.at(i),m_h_eCB_aipc.at(i));
          ceff.SetConditionName(Form("EFEfficiency_pileup_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("EventFilter Efficiency;pileup;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencypileup(m_h_eCB_aipc.at(i),m_h_eEF_aipc.at(i));
          */
          ceff.SetConditionName(Form("L1Efficiency_barrel_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_barrel.at(i),m_h_eL1_pt_barrel.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiency_barrel_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_barrel.at(i),m_h_eSA_pt_barrel.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("CBEfficiency_barrel_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("muComb Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eSA_pt_barrel.at(i),m_h_eCB_pt_barrel.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("EFEfficiency_barrel_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("EventFilter Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eCB_pt_barrel.at(i),m_h_eEF_pt_barrel.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1Efficiency_end_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_end.at(i),m_h_eL1_pt_end.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiency_end_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_end.at(i),m_h_eSA_pt_end.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("CBEfficiency_end_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("muComb Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eSA_pt_end.at(i),m_h_eCB_pt_end.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("EFEfficiency_end_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("EventFilter Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eCB_pt_end.at(i),m_h_eEF_pt_end.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1EfficiencyLargeplus_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L1 Large Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_Largeplus.at(i),m_h_eL1_pt_Largeplus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeplus_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA Large Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_Largeplus.at(i),m_h_eSA_pt_Largeplus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargenormal_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_Largenormal.at(i),m_h_eSA_pt_Largenormal.at(i),m_binmax,200,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1EfficiencyLargeSpecialplus_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L1 LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_LargeSpecialplus.at(i),m_h_eL1_pt_LargeSpecialplus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialplus_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus.at(i),m_h_eSA_pt_LargeSpecialplus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1EfficiencyLargeminus_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L1 Large Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_Largeminus.at(i),m_h_eL1_pt_Largeminus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeminus_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA Large Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_Largeminus.at(i),m_h_eSA_pt_Largeminus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1EfficiencyLargeSpecialminus_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L1 LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_LargeSpecialminus.at(i),m_h_eL1_pt_LargeSpecialminus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialminus_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus.at(i),m_h_eSA_pt_LargeSpecialminus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialplus11_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus11.at(i),m_h_eSA_pt_LargeSpecialplus11.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialminus11_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus11.at(i),m_h_eSA_pt_LargeSpecialminus11.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialplus15_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus15.at(i),m_h_eSA_pt_LargeSpecialplus15.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialminus15_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus15.at(i),m_h_eSA_pt_LargeSpecialminus15.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialplus11out_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus11out.at(i),m_h_eSA_pt_LargeSpecialplus11out.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialminus11out_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus11out.at(i),m_h_eSA_pt_LargeSpecialminus11out.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialplus15out_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus15out.at(i),m_h_eSA_pt_LargeSpecialplus15out.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialminus15out_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus15out.at(i),m_h_eSA_pt_LargeSpecialminus15out.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialplus11in_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus11in.at(i),m_h_eSA_pt_LargeSpecialplus11in.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialminus11in_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus11in.at(i),m_h_eSA_pt_LargeSpecialminus11in.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialplus15in_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus15in.at(i),m_h_eSA_pt_LargeSpecialplus15in.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeSpecialminus15in_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus15in.at(i),m_h_eSA_pt_LargeSpecialminus15in.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1EfficiencySmallplus_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L1 Small Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_Smallplus.at(i),m_h_eL1_pt_Smallplus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencySmallplus_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA Small Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_Smallplus.at(i),m_h_eSA_pt_Smallplus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1EfficiencySmallminus_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L1 Small Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_Smallminus.at(i),m_h_eL1_pt_Smallminus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencySmallminus_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA Small Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_Smallminus.at(i),m_h_eSA_pt_Smallminus.at(i),m_binmax,300,m_efficiency_xerr);

          ceff.SetConditionName(Form("L1EfficiencySmallSpecialplus_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L1 SmallSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_SmallSpecialplus.at(i),m_h_eL1_pt_SmallSpecialplus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencySmallSpecialplus_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA SmallSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_SmallSpecialplus.at(i),m_h_eSA_pt_SmallSpecialplus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("L1EfficiencySmallSpecialminus_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L1 SmallSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eoff_pt_SmallSpecialminus.at(i),m_h_eL1_pt_SmallSpecialminus.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencySmallSpecialminus_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA SmallSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_SmallSpecialminus.at(i),m_h_eSA_pt_SmallSpecialminus.at(i),m_binmax,300,m_efficiency_xerr);
          
          ceff.SetConditionName(Form("SAEfficiencyBarrelwithoutBIM_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA Barrel without BIM Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_BarrelwithoutBIM.at(i),m_h_eSA_pt_BarrelwithoutBIM.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLSwithoutBIM_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA LS without BIM Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LSwithoutBIM.at(i),m_h_eSA_pt_LSwithoutBIM.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLargeincBIM_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA Barrel inc BIM Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_BarrelincBIM.at(i),m_h_eSA_pt_BarrelincBIM.at(i),m_binmax,300,m_efficiency_xerr);
          ceff.SetConditionName(Form("SAEfficiencyLSincBIM_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("SA LS inc BIM Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1_pt_LSincBIM.at(i),m_h_eSA_pt_LSincBIM.at(i),m_binmax,300,m_efficiency_xerr);

          ceff.SetConditionName(Form("SA2DEfficiency_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L1vsL2MuonSA Efficiency;offline eta;offline phi",1.5,0.1,0.1,0.105,0.165);
          ceff.SetConditionbin(m_nbin_eta,m_nbin_phi,m_eta_max,m_phi_max);
          ceff.DrawEfficiency2D(m_h_eff_pL1_etaphi.at(i),m_h_eff_pSA_etaphi.at(i));
          ceff.SetConditionName(Form("L12DEfficiency_%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L1 Efficiency;offline eta;offline phi",1.5,0.1,0.1,0.105,0.165);
          ceff.SetConditionbin(m_nbin_eta,m_nbin_phi,m_eta_max,m_phi_max);
          ceff.DrawEfficiency2D(m_h_eff_poff_etaphi.at(i),m_h_eff_pL1_etaphi.at(i));

          cout<<"test"<<endl;
          ceff.SetConditionName(Form("SAEfficiency_SP1_barrel%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1SP1_pt_barrel.at(i),m_h_eSASP1_pt_barrel.at(i),m_binmax,200,m_efficiency_xerr); 
          cout<<"test"<<endl;
          ceff.SetConditionName(Form("SAEfficiency_SP2_barrel%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1SP2_pt_barrel.at(i),m_h_eSASP2_pt_barrel.at(i),m_binmax,200,m_efficiency_xerr);
          cout<<"test"<<endl;
          ceff.SetConditionName(Form("SAEfficiency_SP3_barrel%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1SP3_pt_barrel.at(i),m_h_eSASP3_pt_barrel.at(i),m_binmax,200,m_efficiency_xerr);
          cout<<"test"<<endl;
          ceff.SetConditionName(Form("SAEfficiency_innmid_barrel%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1innmid_pt_barrel.at(i),m_h_eSAinnmid_pt_barrel.at(i),m_binmax,200,m_efficiency_xerr);
          cout<<"test"<<endl;
          ceff.SetConditionName(Form("SAEfficiency_SP1_endcap%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1SP1_pt_endcap.at(i),m_h_eSASP1_pt_endcap.at(i),m_binmax,200,m_efficiency_xerr);
          cout<<"test"<<endl;
          ceff.SetConditionName(Form("SAEfficiency_SP2_endcap%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1SP2_pt_endcap.at(i),m_h_eSASP2_pt_endcap.at(i),m_binmax,200,m_efficiency_xerr);
          cout<<"test"<<endl;
          ceff.SetConditionName(Form("SAEfficiency_SP3_endcap%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1SP3_pt_endcap.at(i),m_h_eSASP3_pt_endcap.at(i),m_binmax,200,m_efficiency_xerr);
          cout<<"test"<<endl;
          ceff.SetConditionName(Form("SAEfficiency_innmid_endcap%dGeV",i*m_thpitch + m_thmin));
          ceff.SetCondition("L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiency(m_h_eL1innmid_pt_endcap.at(i),m_h_eSAinnmid_pt_endcap.at(i),m_binmax,200,m_efficiency_xerr);
          cout<<"eff end"<<endl;
          
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
          m_h_pSAphivspSAphims.at(i)->Write();
          m_h_pSAphivspSAphibe.at(i)->Write();
          m_h_countSA.at(i)->Write();
          m_h_ChIndexvshitChamber.at(i)->Write();
          //
          m_h_offptBIM_2station.at(i)->Write();
          m_h_offptBIR_2station.at(i)->Write();
          m_h_offptLarge_2station.at(i)->Write();
          m_h_offptSmall_2station.at(i)->Write();
          m_h_offptBIMinnmid_2station.at(i)->Write();
          m_h_offptBIRinnmid_2station.at(i)->Write();
          m_h_offptLargeinnmid_2station.at(i)->Write();
          m_h_offptSmallinnmid_2station.at(i)->Write();
          m_h_offptBIMinnout_2station.at(i)->Write();
          m_h_offptBIRinnout_2station.at(i)->Write();
          m_h_offptLargeinnout_2station.at(i)->Write();
          m_h_offptSmallinnout_2station.at(i)->Write();
          m_h_offptBIMmidout_2station.at(i)->Write();
          m_h_offptBIRmidout_2station.at(i)->Write();
          m_h_offptLargemidout_2station.at(i)->Write();
          m_h_offptSmallmidout_2station.at(i)->Write();
          //
          m_h_offetaBIM_2station.at(i)->Write();
          m_h_offetaBIR_2station.at(i)->Write();
          m_h_offetaLarge_2station.at(i)->Write();
          m_h_offetaSmall_2station.at(i)->Write();
          m_h_offetaBIMinnmid_2station.at(i)->Write();
          m_h_offetaBIRinnmid_2station.at(i)->Write();
          m_h_offetaLargeinnmid_2station.at(i)->Write();
          m_h_offetaSmallinnmid_2station.at(i)->Write();
          m_h_offetaBIMinnout_2station.at(i)->Write();
          m_h_offetaBIRinnout_2station.at(i)->Write();
          m_h_offetaLargeinnout_2station.at(i)->Write();
          m_h_offetaSmallinnout_2station.at(i)->Write();
          m_h_offetaBIMmidout_2station.at(i)->Write();
          m_h_offetaBIRmidout_2station.at(i)->Write();
          m_h_offetaLargemidout_2station.at(i)->Write();
          m_h_offetaSmallmidout_2station.at(i)->Write();
          //
          m_h_offetaBIM_3station.at(i)->Write();
          m_h_offetaBIR_3station.at(i)->Write();
          m_h_offetaLarge_3station.at(i)->Write();
          m_h_offetaSmall_3station.at(i)->Write();

          m_h_sectorvsphi.at(i)->Write();
          m_h_indexvseta.at(i)->Write();
          m_h_overphi.at(i)->Write();
          m_h_mdtchamber.at(i)->Write();
          m_h_L1pSA_sAddress->Write();
          m_h_SApSA_sAddress->Write();

          m_h_mdthitXY.at(i)->Write();
          m_h_rpchitXY.at(i)->Write();
          m_h_mdthitZR.at(i)->Write();
          m_h_mdtSPXY_BI.at(i)->Write();
          m_h_mdtSPXY_BM.at(i)->Write();
          m_h_mdtSPXY_BO.at(i)->Write();
          m_h_mdtSPXY_BME.at(i)->Write();
          m_h_mdtSPX_BI.at(i)->Write();
          m_h_mdtSPY_BI.at(i)->Write();
          m_h_mdtSPXY_etaIndexplus1.at(i)->Write();
          m_h_mdtSPXY_etaIndexplus2.at(i)->Write();
          m_h_mdtSPXY_etaIndexplus3.at(i)->Write();
          m_h_mdtSPXY_etaIndexplus4.at(i)->Write();
          m_h_mdtSPXY_etaIndexplus5.at(i)->Write();
          m_h_mdtSPXY_etaIndexplus6.at(i)->Write();
          m_h_mdtSPXY_etaIndexminus1.at(i)->Write();
          m_h_mdtSPXY_etaIndexminus2.at(i)->Write();
          m_h_mdtSPXY_etaIndexminus3.at(i)->Write();
          m_h_mdtSPXY_etaIndexminus4.at(i)->Write();
          m_h_mdtSPXY_etaIndexminus5.at(i)->Write();
          m_h_mdtSPXY_etaIndexminus6.at(i)->Write();
          m_h_mdtSPZR.at(i)->Write();
          m_h_mdtSPZR_LargeSpecialplus.at(i)->Write();
          m_h_mdtSPZR_LargeSpecialminus.at(i)->Write();
          m_h_mdtSPZR_LargeSpecialplus11out.at(i)->Write();
          m_h_mdtSPZR_LargeSpecialplus11in.at(i)->Write();
          m_h_mdtSPZR_LargeSpecialplus15out.at(i)->Write();
          m_h_mdtSPZR_LargeSpecialplus15in.at(i)->Write();
          m_h_mdtSPZR_LargeSpecialminus11out.at(i)->Write();
          m_h_mdtSPZR_LargeSpecialminus11in.at(i)->Write();
          m_h_mdtSPZR_LargeSpecialminus15out.at(i)->Write();
          m_h_mdtSPZR_LargeSpecialminus15in.at(i)->Write();
          m_h_segmentXY.at(i)->Write();
          m_h_segmentXYmatch.at(i)->Write();
          m_h_segmentXYnomatch.at(i)->Write();
          m_h_segmentXY_etaIndexplus1.at(i)->Write();
          m_h_segmentXY_etaIndexplus2.at(i)->Write();
          m_h_segmentXY_etaIndexplus3.at(i)->Write();
          m_h_segmentXY_etaIndexplus4.at(i)->Write();
          m_h_segmentXY_etaIndexplus5.at(i)->Write();
          m_h_segmentXY_etaIndexplus6.at(i)->Write();
          m_h_segmentXY_etaIndexminus1.at(i)->Write();
          m_h_segmentXY_etaIndexminus2.at(i)->Write();
          m_h_segmentXY_etaIndexminus3.at(i)->Write();
          m_h_segmentXY_etaIndexminus4.at(i)->Write();
          m_h_segmentXY_etaIndexminus5.at(i)->Write();
          m_h_segmentXY_etaIndexminus6.at(i)->Write();
          m_h_segmentXY_BIL.at(i)->Write();
          m_h_segmentXY_BML.at(i)->Write();
          m_h_segmentXY_BOL.at(i)->Write();
          m_h_segmentXY_BIS.at(i)->Write();
          m_h_segmentXY_BMS.at(i)->Write();
          m_h_segmentXY_BOS.at(i)->Write();
          m_h_segmentXY_BEE.at(i)->Write();
          m_h_segmentXY_EIL.at(i)->Write();
          m_h_segmentXY_EML.at(i)->Write();
          m_h_segmentXY_EOL.at(i)->Write();
          m_h_segmentXY_EIS.at(i)->Write();
          m_h_segmentXY_EMS.at(i)->Write();
          m_h_segmentXY_EOS.at(i)->Write();
          m_h_segmentXY_EES.at(i)->Write();
          m_h_segmentXY_EEL.at(i)->Write();
          m_h_segmentXY_CSS.at(i)->Write();
          m_h_segmentXY_CSL.at(i)->Write();
          m_h_segmentXY_unknown.at(i)->Write();
          m_h_segmentXY_normal.at(i)->Write();
          m_h_segmentXY_special.at(i)->Write();
          /*
          m_h_segmentZR.at(i)->Write();
          m_h_segmentZR_BIL.at(i)->Write();
          m_h_segmentZR_BML.at(i)->Write();
          m_h_segmentZR_LargeSpecialplus.at(i)->Write();
          m_h_segmentZR_LargeSpecialminus.at(i)->Write();
          m_h_segmentZR_LargeSpecialplus11out.at(i)->Write();
          m_h_segmentZR_LargeSpecialplus11in.at(i)->Write();
          m_h_segmentZR_LargeSpecialplus15out.at(i)->Write();
          m_h_segmentZR_LargeSpecialplus15in.at(i)->Write();
          m_h_segmentZR_LargeSpecialminus11out.at(i)->Write();
          m_h_segmentZR_LargeSpecialminus11in.at(i)->Write();
          m_h_segmentZR_LargeSpecialminus15out.at(i)->Write();
          m_h_segmentZR_LargeSpecialminus15in.at(i)->Write();
          m_h_segmentZR_LargeSpecialplus11out_BIL.at(i)->Write();
          m_h_segmentZR_LargeSpecialplus11in_BIL.at(i)->Write();
          m_h_segmentZR_LargeSpecialplus15out_BIL.at(i)->Write();
          m_h_segmentZR_LargeSpecialplus15in_BIL.at(i)->Write();
          */
          m_h_numsegment.at(i)->Write();
          m_h_num_segment_LSBI.at(i)->Write();
          m_h_num_segment_LargeBI.at(i)->Write();
          m_h_numhit.at(i)->Write();
          m_h_nPrecisionHits_normal.at(i)->Write();
          m_h_nPrecisionHits_special.at(i)->Write();
          m_h_sector_normal.at(i)->Write();
          m_h_sector_special.at(i)->Write();
          m_h_etaIndex_normal.at(i)->Write();
          m_h_etaIndex_special.at(i)->Write();
          m_h_mdtphi.at(i)->Write();
          m_h_mdtphi_LS.at(i)->Write();
          m_h_mdtphi_LSBIL.at(i)->Write();
          m_h_mdtR.at(i)->Write();
          m_h_numSP.at(i)->Write();
          m_h_BIRsegment.at(i)->Write();
          m_h_BIMsegment.at(i)->Write();
          m_h_BIMrvsx.at(i)->Write();

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
          m_h_eoff_phi.at(i)->Write();
          m_h_eL1_phi.at(i)->Write();
          m_h_eSA_phi.at(i)->Write();
          m_h_eCB_phi.at(i)->Write();
          m_h_eEF_phi.at(i)->Write();
          //m_h_eoff_aipc.at(i)->Write();
          //m_h_eL1_aipc.at(i)->Write();
          //m_h_eSA_aipc.at(i)->Write();
          //m_h_eCB_aipc.at(i)->Write();
          //m_h_eEF_aipc.at(i)->Write();
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
          m_h_eL1_pt_Largenormal.at(i)->Write();
          m_h_eSA_pt_Largenormal.at(i)->Write();
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
          m_h_eL1_pt_BarrelwithoutBIM.at(i)->Write();
          m_h_eSA_pt_BarrelwithoutBIM.at(i)->Write();
          m_h_eL1_pt_LSwithoutBIM.at(i)->Write();
          m_h_eSA_pt_LSwithoutBIM.at(i)->Write();
          m_h_eL1_pt_BarrelincBIM.at(i)->Write();
          m_h_eSA_pt_BarrelincBIM.at(i)->Write();
          m_h_eL1_pt_LSincBIM.at(i)->Write();
          m_h_eSA_pt_LSincBIM.at(i)->Write();
          m_h_eL1SP1_pt_barrel.at(i)->Write();
          m_h_eSASP1_pt_barrel.at(i)->Write();
          m_h_eL1SP2_pt_barrel.at(i)->Write();
          m_h_eSASP2_pt_barrel.at(i)->Write();
          m_h_eL1SP3_pt_barrel.at(i)->Write();
          m_h_eSASP3_pt_barrel.at(i)->Write();
          m_h_eL1innmid_pt_barrel.at(i)->Write();
          m_h_eSAinnmid_pt_barrel.at(i)->Write();
          m_h_eL1SP1_pt_endcap.at(i)->Write();
          m_h_eSASP1_pt_endcap.at(i)->Write();
          m_h_eL1SP2_pt_endcap.at(i)->Write();
          m_h_eSASP2_pt_endcap.at(i)->Write();
          m_h_eL1SP3_pt_endcap.at(i)->Write();
          m_h_eSASP3_pt_endcap.at(i)->Write();
          m_h_eL1innmid_pt_endcap.at(i)->Write();
          m_h_eSAinnmid_pt_endcap.at(i)->Write();

          m_h_pSA_respt.at(i)->Write();
          m_h_pCB_respt.at(i)->Write();
          m_h_pEF_respt.at(i)->Write();
          m_h_pSA_respt_barrel.at(i)->Write();
          m_h_pCB_respt_barrel.at(i)->Write();
          m_h_pEF_respt_barrel.at(i)->Write();
          m_h_pSA_respt_endcap.at(i)->Write();
          m_h_pCB_respt_endcap.at(i)->Write();
          m_h_pEF_respt_endcap.at(i)->Write();
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
          //
          m_h_off_ptvsSA_resptLarge_2station.at(i)->Write();
          m_h_off_ptvsSA_resptSmall_2station.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialplus15out_2station.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialplus15in_2station.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialplus11out_2station.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialplus11in_2station.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialminus15out_2station.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialminus15in_2station.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialminus11out_2station.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialminus11in_2station.at(i)->Write();
          //
          m_h_off_ptvsSA_resptLargeSpecialplus15out_3station.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialplus15in_3station.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialplus11out_3station.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialplus11in_3station.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialminus15out_3station.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialminus15in_3station.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialminus11out_3station.at(i)->Write();
          m_h_off_ptvsSA_resptLargeSpecialminus11in_3station.at(i)->Write();
          //
          m_h_off_ptvsSA_resptBIMinnmid_2station.at(i)->Write();
          m_h_off_ptvsSA_resptBIRinnmid_2station.at(i)->Write();
          m_h_off_ptvsSA_resptLargeinnmid_2station.at(i)->Write();
          m_h_off_ptvsSA_resptSmallinnmid_2station.at(i)->Write();
          m_h_off_ptvsSA_resptBIMinnout_2station.at(i)->Write();
          m_h_off_ptvsSA_resptBIRinnout_2station.at(i)->Write();
          m_h_off_ptvsSA_resptLargeinnout_2station.at(i)->Write();
          m_h_off_ptvsSA_resptSmallinnout_2station.at(i)->Write();
          m_h_off_ptvsSA_resptBIMmidout_2station.at(i)->Write();
          m_h_off_ptvsSA_resptBIRmidout_2station.at(i)->Write();
          m_h_off_ptvsSA_resptLargemidout_2station.at(i)->Write();
          m_h_off_ptvsSA_resptSmallmidout_2station.at(i)->Write();
          //
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
          //
          m_h_etaIndexvsSA_resptLargeSpecialplus11out_2station.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialplus11in_2station.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialplus15out_2station.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialplus15in_2station.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialminus11out_2station.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialminus11in_2station.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialminus15out_2station.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialminus15in_2station.at(i)->Write();
          //
          m_h_etaIndexvsSA_resptLargeSpecialplus11out_3station.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialplus11in_3station.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialplus15out_3station.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialplus15in_3station.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialminus11out_3station.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialminus11in_3station.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialminus15out_3station.at(i)->Write();
          m_h_etaIndexvsSA_resptLargeSpecialminus15in_3station.at(i)->Write();

          cout<<"residual end"<<endl;
     }

}
