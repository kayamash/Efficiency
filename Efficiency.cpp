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
//const Double_t pt_threshold[4] = {5.17,3.25,4.69,5.14};//MU6
const Double_t pt_threshold[4] = {15.87,10.73,12.21,15.87};//MU20

void Efficiency::Init(std::string name,const Int_t np,const Int_t ne,const Double_t mp,const Double_t me,Double_t req,Int_t max,Double_t err,Int_t proc){
     m_nbin_phi = np;
     m_nbin_eta = ne;
     m_phi_max = mp;
     m_eta_max = me;
     m_method_name = name;
     m_binmax = max;
     m_efficiency_xerr = err;
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
    return kTRUE;
}else{
     return kFALSE;
}
}

bool Efficiency::CutSA(Int_t pass){
     if(pass == 1){
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

bool Efficiency::PlateauCut(Double_t pt){
     Double_t cut_pt = 40.0;
     if(m_method_name == "mu4" || m_method_name == "mu6")cut_pt = 8.0;
     if(pt > cut_pt){
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
     Double_t pSA_superpointZ_BI = 0;
     Double_t pSA_superpointZ_BM = 0;
     Double_t pSA_superpointZ_BO = 0;
     Double_t pSA_superpointZ_BME = 0;
     Double_t pSA_superpointZ_EI = 0;
     Double_t pSA_superpointZ_EM = 0;
     Double_t pSA_superpointZ_EO = 0;
     Double_t pSA_superpointZ_EE = 0;
     Double_t pSA_superpointZ_CSC = 0;
     Double_t pSA_superpointZ_BEE = 0;
     Double_t pSA_superpointR_BI = 0;
     Double_t pSA_superpointR_BM = 0;
     Double_t pSA_superpointR_BO = 0;
     Double_t pSA_superpointR_BME = 0;
     Double_t pSA_superpointR_EI = 0;
     Double_t pSA_superpointR_EM = 0;
     Double_t pSA_superpointR_EO = 0;
     Double_t pSA_superpointR_EE = 0;
     Double_t pSA_superpointR_CSC = 0;
     Double_t pSA_superpointR_BEE = 0;
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
                    pSA_superpointZ_EI = m_pSA_superpointZ_EI->at(method);
                    pSA_superpointZ_EM = m_pSA_superpointZ_EM->at(method);
                    pSA_superpointZ_EO = m_pSA_superpointZ_EO->at(method);
                    pSA_superpointZ_EE = m_pSA_superpointZ_EE->at(method);
                    pSA_superpointZ_CSC = m_pSA_superpointZ_CSC->at(method);
                    pSA_superpointZ_BEE = m_pSA_superpointZ_BEE->at(method);
                    pSA_superpointR_BI = m_pSA_superpointR_BI->at(method);
                    pSA_superpointR_BM = m_pSA_superpointR_BM->at(method);
                    pSA_superpointR_BO = m_pSA_superpointR_BO->at(method);
                    pSA_superpointR_BME = m_pSA_superpointR_BME->at(method);
                    pSA_superpointR_EI = m_pSA_superpointR_EI->at(method);
                    pSA_superpointR_EM = m_pSA_superpointR_EM->at(method);
                    pSA_superpointR_EO = m_pSA_superpointR_EO->at(method);
                    pSA_superpointR_EE = m_pSA_superpointR_EE->at(method);
                    pSA_superpointR_CSC = m_pSA_superpointR_CSC->at(method);
                    pSA_superpointR_BEE = m_pSA_superpointR_BEE->at(method);
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
               if(m_probe_segment_chamberIndex[index] == 1 && (m_probe_segment_sector[index] == 11 || m_probe_segment_sector[index] == 15) )m_h_BIMrvsx->Fill(sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)),m_probe_segment_x[index]);
          }

          tL1_dR = TMath::Sqrt(pow(m_tL1_eta - m_toff_eta,2) + pow(m_tL1_phi - m_toff_phi,2) );
          tEF_dR = TMath::Sqrt(pow(m_tEF_eta - m_toff_eta,2) + pow(m_tEF_phi - m_toff_phi,2) );
          if(std::fabs(m_toff_pt)*0.001 < 10.0)m_reqL1dR = -0.00001*std::fabs(m_toff_pt) + 0.18;
          Int_t SPinner = 0;
          Int_t SPmiddle = 0;
          Int_t SPouter = 0;
          if(pSA_superpointR_BI != 0 && pSA_superpointR_EI != 0){
               numSP++;
               patternSP += 1;
               SPinner = 1;
          }
          if(pSA_superpointR_BM != 0 && pSA_superpointR_EM != 0){
               numSP++;
               patternSP += 2;
               SPmiddle = 1;
          }
          if(pSA_superpointR_BO != 0 && pSA_superpointR_EO != 0){
               numSP++;
               patternSP += 3;
               SPouter = 1;
          }
          if(pSA_superpointR_EE != 0)numSP++;
          if(pSA_superpointR_CSC != 0)numSP++;
          if(pSA_superpointR_BEE != 0)numSP++;
          if(SPinner == 1 && SPmiddle == 1 && SPouter == 0)patternSP = 3;
          if(SPinner == 1 && SPmiddle == 0 && SPouter == 1)patternSP = 4;
          if(SPinner == 0 && SPmiddle == 1 && SPouter == 1)patternSP = 5;
          Int_t decision_noBIM = 0;
          for(Int_t index = 0;index < 10;index++){
               if(m_probe_segment_chamberIndex[index] == 1 && sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)) > 5800 && fabs(m_probe_segment_x[index]) > 4000.)decision_noBIM++;
               //if(m_probe_segment_chamberIndex[index] == 1 && (m_probe_segment_sector[index] == 11 || m_probe_segment_sector[index] == 15) && ((pSA_roiphi > -0.8 && pSA_roiphi < -0.6) || (pSA_roiphi > -2.6 && pSA_roiphi < -2.4)) )decision_noBIM++;
          }

          //offline
          if(!CutTagProbe(pEFTAG_pass))return;
          if(static_cast<Int_t>(pSA_sAddress) == 1)m_h_offphi_LargeSpecial->Fill(m_poff_phi);
          m_h_poff_pt->Fill(m_poff_pt*0.001);
          m_h_eoff_pt->Fill(std::fabs(m_poff_pt*0.001));
          if(PlateauCut(std::fabs(m_poff_pt*0.001))){
               m_h_eoff_eta->Fill(m_poff_eta);
               m_h_eoff_phi->Fill(m_poff_phi);
               //m_h_eoff_aipc->Fill(m_aipc);
          }
          if(DicisionBarrel(m_poff_eta)){
               m_h_eoff_pt_barrel->Fill(std::fabs(m_poff_pt*0.001));
          }else{
               m_h_eoff_pt_end->Fill(std::fabs(m_poff_pt*0.001));
          }
          if(std::fabs(m_poff_pt*0.001) > 40)m_h_eff_poff_etaphi->Fill(m_poff_eta,m_poff_phi);

          switch(static_cast<Int_t>(pSA_sAddress)){
               case 0:
               if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1)m_h_eoff_pt_Largeplus->Fill(std::fabs(m_poff_pt*0.001));
               if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1)m_h_eoff_pt_Largeminus->Fill(std::fabs(m_poff_pt*0.001));
               break;
               case 1:
               if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1)m_h_eoff_pt_LargeSpecialplus->Fill(std::fabs(m_poff_pt*0.001));
               if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1)m_h_eoff_pt_LargeSpecialminus->Fill(std::fabs(m_poff_pt*0.001));
               break;
               case 2:
               if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1)m_h_eoff_pt_Smallplus->Fill(std::fabs(m_poff_pt*0.001));
               if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1)m_h_eoff_pt_Smallminus->Fill(std::fabs(m_poff_pt*0.001));
               break;
               case 3:
               if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1)m_h_eoff_pt_SmallSpecialplus->Fill(std::fabs(m_poff_pt*0.001));
               if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1)m_h_eoff_pt_SmallSpecialminus->Fill(std::fabs(m_poff_pt*0.001));
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
          if(!CutL1(pL1_pass))return;
          Double_t textL1_dR = TMath::Sqrt(pow(m_tL1_eta - m_toff_exteta,2) + pow(m_tL1_phi - m_toff_extphi,2));
          pextL1_dR = TMath::Sqrt(pow(pL1_eta - m_poff_exteta,2) + pow(pL1_phi - m_poff_extphi,2));
          m_h_pL1_pt->Fill(std::fabs(pL1_pt*0.001));
          m_h_pL1_dR->Fill(pL1_dR);
          m_h_textL1_dR->Fill(textL1_dR);
          m_h_pextL1_dR->Fill(pextL1_dR);
          m_h_eL1_pt->Fill(std::fabs(m_poff_pt*0.001));
          m_h_L1pSA_sAddress->Fill(pSA_sAddress);
          if(DicisionBarrel(m_poff_eta)){
               m_h_eL1_pt_barrel->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 1)m_h_eL1SP1_pt_barrel->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 2)m_h_eL1SP2_pt_barrel->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 3)m_h_eL1SP3_pt_barrel->Fill(std::fabs(m_poff_pt*0.001));
               if(patternSP == 3)m_h_eL1innmid_pt_barrel->Fill(std::fabs(m_poff_pt*0.001));
               if(decision_noBIM == 0)m_h_eL1_pt_BarrelwithoutBIM->Fill(std::fabs(m_poff_pt*0.001));
               m_h_eL1_pt_BarrelincBIM->Fill(std::fabs(m_poff_pt*0.001));
               if(pSA_sAddress == 0 && nosector9 == 0)m_h_eL1_pt_Largenormal->Fill(std::fabs(m_poff_pt*0.001));
          }else{
               m_h_eL1_pt_end->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 1)m_h_eL1SP1_pt_endcap->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 2)m_h_eL1SP2_pt_endcap->Fill(std::fabs(m_poff_pt*0.001));
               if(numSP == 3)m_h_eL1SP3_pt_endcap->Fill(std::fabs(m_poff_pt*0.001));
               if(patternSP == 3)m_h_eL1innmid_pt_endcap->Fill(std::fabs(m_poff_pt*0.001));
          }

          if(PlateauCut(std::fabs(m_poff_pt*0.001))){
               m_h_eff_pL1_etaphi->Fill(m_poff_eta,m_poff_phi);
               m_h_eL1_eta->Fill(m_poff_eta);
               m_h_eL1_phi->Fill(m_poff_phi);
               //m_h_eL1_aipc->Fill(m_aipc);
          }

          areanumber = DicisionArea(pSA_roiphi);

          switch(static_cast<Int_t>(pSA_sAddress)){//switch Large ,LS , Small ,SS
               case 0:
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0)m_h_eL1_pt_Largeplus->Fill(std::fabs(m_poff_pt*0.001));//Qeta = +1
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) < 0)m_h_eL1_pt_Largeminus->Fill(std::fabs(m_poff_pt*0.001));//Qeta = -1
                    break;
                    case 1:
                    if(areanumber > 0 && areanumber < 5)m_h_eL1_pt_LargeSpecialplus->Fill(std::fabs(m_poff_pt*0.001));//Qeta = +1
                    if(areanumber > 4 && areanumber < 9)m_h_eL1_pt_LargeSpecialminus->Fill(std::fabs(m_poff_pt*0.001));//Qeta = -1
                    if(decision_noBIM == 0)m_h_eL1_pt_LSwithoutBIM->Fill(std::fabs(m_poff_pt*0.001));
                    m_h_eL1_pt_LSincBIM->Fill(std::fabs(m_poff_pt*0.001));
                    switch(areanumber){
                         case 1://plus11out
                         m_h_eL1_pt_LargeSpecialplus11out->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_eL1_pt_LargeSpecialplus11->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                         case 2://plus11in
                         m_h_eL1_pt_LargeSpecialplus11in->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_eL1_pt_LargeSpecialplus11->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                         case 3://plus15out
                         m_h_eL1_pt_LargeSpecialplus15out->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_eL1_pt_LargeSpecialplus15->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                         case 4://plus15in
                         m_h_eL1_pt_LargeSpecialplus15in->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_eL1_pt_LargeSpecialplus15->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                         case 5://minus11out
                         m_h_eL1_pt_LargeSpecialminus11out->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_eL1_pt_LargeSpecialminus11->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                         case 6://minus11in
                         m_h_eL1_pt_LargeSpecialminus11in->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_eL1_pt_LargeSpecialminus11->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                         case 7://minus15out
                         m_h_eL1_pt_LargeSpecialminus15out->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_eL1_pt_LargeSpecialminus15->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                         case 8://minus15in
                         m_h_eL1_pt_LargeSpecialminus15in->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_eL1_pt_LargeSpecialminus15->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                         default:
                         break;
                    }
                    break;
                    case 2:
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0)m_h_eL1_pt_Smallplus->Fill(std::fabs(m_poff_pt*0.001));//Qeta = +1
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) < 0)m_h_eL1_pt_Smallminus->Fill(std::fabs(m_poff_pt*0.001));//Qeta = -1
                    break;
                    case 3:
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0)m_h_eL1_pt_SmallSpecialplus->Fill(std::fabs(m_poff_pt*0.001));//Qeta = +1
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) < 0)m_h_eL1_pt_SmallSpecialminus->Fill(std::fabs(m_poff_pt*0.001));//Qeta = -1
                    break;
                    default:
                    break;
               }

          //SA
               if(!CutSA(pSA_pass))return;
               m_h_countSA->Fill(m_poff_eta);
               m_h_numSP->Fill(numSP);
               Double_t textSA_dR = TMath::Sqrt(pow(m_tSA_eta - m_toff_exteta,2) + pow(m_tSA_phi - m_toff_extphi,2));
               pextSA_dR = TMath::Sqrt(pow(pSA_eta - m_poff_exteta,2) + pow(pSA_phi - m_poff_extphi,2));
               Double_t resSA_pt = std::fabs(m_poff_pt*0.001)/std::fabs(pSA_pt) - 1.0;
               Double_t buf_pSA_dR = TMath::Sqrt(pow(pSA_eta - m_poff_eta,2) + pow(pSA_phi - m_poff_phi,2));
               Double_t buf_eta = 0;
               m_h_SApSA_sAddress->Fill(pSA_sAddress);

               m_h_pSA_pt->Fill(std::fabs(pSA_pt));
               m_h_pSA_dR->Fill(buf_pSA_dR);
               m_h_textSA_dR->Fill(textSA_dR);
               m_h_pextSA_dR->Fill(pextSA_dR);
               m_h_eSA_pt->Fill(std::fabs(m_poff_pt*0.001));
               m_h_pSA_respt->Fill(resSA_pt);
               m_h_pSAphivspSAphims->Fill(pSA_phi,pSA_phims);
               m_h_pSAphivspSAphibe->Fill(pSA_phi,pSA_phibe);
               if(DicisionBarrel(m_poff_eta)){
                    m_h_eSA_pt_barrel->Fill(std::fabs(m_poff_pt*0.001));
                    m_h_pSA_respt_barrel->Fill(resSA_pt);
                    if(numSP == 1)m_h_eSASP1_pt_barrel->Fill(std::fabs(m_poff_pt*0.001));
                    if(numSP == 2)m_h_eSASP2_pt_barrel->Fill(std::fabs(m_poff_pt*0.001));
                    if(numSP == 3)m_h_eSASP3_pt_barrel->Fill(std::fabs(m_poff_pt*0.001));
                    if(patternSP == 3)m_h_eSAinnmid_pt_barrel->Fill(std::fabs(m_poff_pt*0.001));
                    if(decision_noBIM == 0)m_h_eSA_pt_BarrelwithoutBIM->Fill(std::fabs(m_poff_pt*0.001));
                    m_h_eSA_pt_BarrelincBIM->Fill(std::fabs(m_poff_pt*0.001));
                    for(Int_t MDTsize = 0;MDTsize < (signed int)pSA_mdthitChamber->size();MDTsize++){
                         m_h_mdtchamber->Fill(pSA_mdthitChamber->at(MDTsize));
                    }
                    if(pSA_sAddress == 0 && nosector9 == 0)m_h_eSA_pt_Largenormal->Fill(std::fabs(m_poff_pt*0.001));
               }else{
                    m_h_eSA_pt_end->Fill(std::fabs(m_poff_pt*0.001));
                    m_h_pSA_respt_endcap->Fill(resSA_pt);
                    if(numSP == 1)m_h_eSASP1_pt_endcap->Fill(std::fabs(m_poff_pt*0.001));
                    if(numSP == 2)m_h_eSASP2_pt_endcap->Fill(std::fabs(m_poff_pt*0.001));
                    if(numSP == 3)m_h_eSASP3_pt_endcap->Fill(std::fabs(m_poff_pt*0.001));
                    if(patternSP == 3)m_h_eSAinnmid_pt_endcap->Fill(std::fabs(m_poff_pt*0.001));
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
               }
               for(Int_t size = 0;size < (signed int)pSA_mdtZ->size();size++){
                    buf_eta += -TMath::Log((sqrt(pow(pSA_mdtZ->at(size),2) + pow(pSA_mdtR->at(size),2)) - pSA_mdtZ->at(size))/(sqrt(pow(pSA_mdtZ->at(size),2) + pow(pSA_mdtR->at(size),2)) + pow(pSA_mdtZ->at(size),2)))/2.0;
               }
               m_h_numhit->Fill(pSA_mdtZ->size());
               for(Int_t index = 0;index < 10;index++){
                    m_h_sectorvsphi->Fill(m_probe_segment_sector[index],m_poff_phi);
                    m_h_indexvseta->Fill(m_probe_segment_etaIndex[index],m_poff_eta);
               }
               if(static_cast<Int_t>(pSA_sAddress) == 2 && m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0)m_h_eSA_pt_Smallplus->Fill(std::fabs(m_poff_pt*0.001));
          switch(static_cast<Int_t>(pSA_sAddress)){//switch Large ,LS , Small ,SS
               case 0:
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0){//Qeta = +1
                         m_h_off_ptvsSA_resptLargeplus->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptLargeplus->Fill(resSA_pt);
                         m_h_eSA_pt_Largeplus->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_offphivsSA_resptLargeplus->Fill(m_poff_phi,resSA_pt);
                         m_h_offetavsSA_resptLargeplus->Fill(m_poff_eta,resSA_pt);

                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeplus->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         }
                    }

                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) < 0){//Qeta = -1
                         m_h_off_ptvsSA_resptLargeminus->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptLargeminus->Fill(resSA_pt);
                         m_h_eSA_pt_Largeminus->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_offphivsSA_resptLargeminus->Fill(m_poff_phi,resSA_pt);
                         m_h_offetavsSA_resptLargeminus->Fill(m_poff_eta,resSA_pt);
                    }
                    break;
                    case 1:
                    m_h_saphims_LargeSpecial->Fill(pSA_phims);
                    m_h_saroiphi_LargeSpecial->Fill(pSA_roiphi);
                    if(decision_noBIM == 0)m_h_eSA_pt_LSwithoutBIM->Fill(std::fabs(m_poff_pt*0.001));
                    m_h_eSA_pt_LSincBIM->Fill(std::fabs(m_poff_pt*0.001));

                    for(Int_t index = 0;index < 10;index++){
                         if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialplus11out->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 2)m_h_etaIndexvsSA_resptLargeSpecialplus11out_2station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 3)m_h_etaIndexvsSA_resptLargeSpecialplus11out_3station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0){
                              //m_h_segmentZR_LargeSpecialplus11out->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              //if(m_probe_segment_chamberIndex[index] == 1)m_h_segmentZR_LargeSpecialplus11out_BIL->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         }
                    }

                    if(areanumber > 0 && areanumber < 5){
                         m_h_eSA_pt_LargeSpecialplus->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_SA_resptLargeSpecialplus->Fill(resSA_pt);
                         m_h_off_ptvsSA_resptLargeSpecialplus->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_offphivsSA_resptLargeSpecialplus->Fill(m_poff_phi,resSA_pt);
                         m_h_offetavsSA_resptLargeSpecialplus->Fill(m_poff_eta,resSA_pt);
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialplus->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              //if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialplus->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         }
                         m_h_mdtSPZR_LargeSpecialplus->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_mdtSPZR_LargeSpecialplus->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_mdtSPZR_LargeSpecialplus->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                    }
                    if(areanumber > 4 && areanumber < 9){
                         m_h_eSA_pt_LargeSpecialminus->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_SA_resptLargeSpecialminus->Fill(resSA_pt);
                         m_h_off_ptvsSA_resptLargeSpecialminus->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_offphivsSA_resptLargeSpecialminus->Fill(m_poff_phi,resSA_pt);
                         m_h_offetavsSA_resptLargeSpecialminus->Fill(m_poff_eta,resSA_pt);
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialminus->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              //if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialminus->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         }
                         m_h_mdtSPZR_LargeSpecialminus->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_mdtSPZR_LargeSpecialminus->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_mdtSPZR_LargeSpecialminus->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                    }

                    switch(areanumber){
                         case 1://plus 11out
                         m_h_mdtSPX_BI->Fill(pSA_superpointR_BI*cos(pSA_roiphi));
                         m_h_mdtSPY_BI->Fill(pSA_superpointR_BI*sin(pSA_roiphi));
                         m_h_off_ptvsSA_resptLargeSpecialplus11out->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numSP == 2)m_h_off_ptvsSA_resptLargeSpecialplus11out_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numSP == 3)m_h_off_ptvsSA_resptLargeSpecialplus11out_3station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptLargeSpecialplus11out->Fill(resSA_pt);
                         m_h_offetavsSA_resptLargeSpecialplus11out->Fill(m_poff_eta,resSA_pt);
                         m_h_eSA_pt_LargeSpecialplus11out->Fill(std::fabs(m_poff_pt*0.001));

                         m_h_mdtSPZR_LargeSpecialplus11out->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_mdtSPZR_LargeSpecialplus11out->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_mdtSPZR_LargeSpecialplus11out->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                         break;
                         case 2://plus 11in
                         m_h_mdtSPX_BI->Fill(pSA_superpointR_BI*cos(pSA_roiphi));
                         m_h_mdtSPY_BI->Fill(pSA_superpointR_BI*sin(pSA_roiphi));
                         m_h_off_ptvsSA_resptLargeSpecialplus11in->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numSP == 2)m_h_off_ptvsSA_resptLargeSpecialplus11in_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numSP == 3)m_h_off_ptvsSA_resptLargeSpecialplus11in_3station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptLargeSpecialplus11in->Fill(resSA_pt);
                         m_h_offetavsSA_resptLargeSpecialplus11in->Fill(m_poff_eta,resSA_pt);
                         m_h_eSA_pt_LargeSpecialplus11in->Fill(std::fabs(m_poff_pt*0.001));
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialplus11in->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 2)m_h_etaIndexvsSA_resptLargeSpecialplus11in_2station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 3)m_h_etaIndexvsSA_resptLargeSpecialplus11in_3station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0){
                                        //m_h_segmentZR_LargeSpecialplus11in->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                                        //if(m_probe_segment_chamberIndex[index] == 1)m_h_segmentZR_LargeSpecialplus11in_BIL->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              }
                         }
                         m_h_mdtSPZR_LargeSpecialplus11in->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_mdtSPZR_LargeSpecialplus11in->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_mdtSPZR_LargeSpecialplus11in->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                         break;
                         case 3://plus 15out
                         m_h_off_ptvsSA_resptLargeSpecialplus15out->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numSP == 2)m_h_off_ptvsSA_resptLargeSpecialplus15out_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numSP == 3)m_h_off_ptvsSA_resptLargeSpecialplus15out_3station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptLargeSpecialplus15out->Fill(resSA_pt);
                         m_h_offetavsSA_resptLargeSpecialplus15out->Fill(m_poff_eta,resSA_pt);
                         m_h_eSA_pt_LargeSpecialplus15out->Fill(std::fabs(m_poff_pt*0.001));

                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8 && m_probe_segment_etaIndex[index] <= 8)m_h_etaIndexvsSA_resptLargeSpecialplus15out->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8 && m_probe_segment_etaIndex[index] <= 8 && numSP == 2)m_h_etaIndexvsSA_resptLargeSpecialplus15out_2station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8 && m_probe_segment_etaIndex[index] <= 8 && numSP == 3)m_h_etaIndexvsSA_resptLargeSpecialplus15out_3station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0){
                                        //m_h_segmentZR_LargeSpecialplus15out->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                                        //if(m_probe_segment_chamberIndex[index] == 1)m_h_segmentZR_LargeSpecialplus15out_BIL->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              }
                         }
                         m_h_mdtSPZR_LargeSpecialplus15out->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_mdtSPZR_LargeSpecialplus15out->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_mdtSPZR_LargeSpecialplus15out->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                         break;
                         case 4://plus 15in
                         m_h_off_ptvsSA_resptLargeSpecialplus15in->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numSP == 2)m_h_off_ptvsSA_resptLargeSpecialplus15in_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numSP == 3)m_h_off_ptvsSA_resptLargeSpecialplus15in_3station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptLargeSpecialplus15in->Fill(resSA_pt);
                         m_h_offetavsSA_resptLargeSpecialplus15in->Fill(m_poff_eta,resSA_pt);
                         m_h_eSA_pt_LargeSpecialplus15in->Fill(std::fabs(m_poff_pt*0.001));
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialplus15in->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 2)m_h_etaIndexvsSA_resptLargeSpecialplus15in_2station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 3)m_h_etaIndexvsSA_resptLargeSpecialplus15in_3station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0){
                                        //m_h_segmentZR_LargeSpecialplus15in->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                                        //if(m_probe_segment_chamberIndex[index] == 1)m_h_segmentZR_LargeSpecialplus15in_BIL->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              }
                         }
                         m_h_mdtSPZR_LargeSpecialplus15in->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_mdtSPZR_LargeSpecialplus15in->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_mdtSPZR_LargeSpecialplus15in->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                         break;
                         case 5://minus 11out
                         m_h_off_ptvsSA_resptLargeSpecialminus11out->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numSP == 2)m_h_off_ptvsSA_resptLargeSpecialminus11out_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numSP == 3)m_h_off_ptvsSA_resptLargeSpecialminus11out_3station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptLargeSpecialminus11out->Fill(resSA_pt);
                         m_h_offetavsSA_resptLargeSpecialminus11out->Fill(m_poff_eta,resSA_pt);
                         m_h_eSA_pt_LargeSpecialminus11out->Fill(std::fabs(m_poff_pt*0.001));
                         for(Int_t index = 0;index < 10;index++){
                          if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialminus11out->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                          if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 2)m_h_etaIndexvsSA_resptLargeSpecialminus11out_2station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                          if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 3)m_h_etaIndexvsSA_resptLargeSpecialminus11out_3station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                  //if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialminus11out->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                     }
                     m_h_mdtSPZR_LargeSpecialminus11out->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                     m_h_mdtSPZR_LargeSpecialminus11out->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                     m_h_mdtSPZR_LargeSpecialminus11out->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                     break;
                         case 6://minus 11in
                         m_h_off_ptvsSA_resptLargeSpecialminus11in->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numSP == 2)m_h_off_ptvsSA_resptLargeSpecialminus11in_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numSP == 3)m_h_off_ptvsSA_resptLargeSpecialminus11in_3station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptLargeSpecialminus11in->Fill(resSA_pt);
                         m_h_offetavsSA_resptLargeSpecialminus11in->Fill(m_poff_eta,resSA_pt);
                         m_h_eSA_pt_LargeSpecialminus11in->Fill(std::fabs(m_poff_pt*0.001));
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialminus11in->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 2)m_h_etaIndexvsSA_resptLargeSpecialminus11in_2station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 3)m_h_etaIndexvsSA_resptLargeSpecialminus11in_3station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                   //if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialminus11in->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         }
                         m_h_mdtSPZR_LargeSpecialminus11in->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_mdtSPZR_LargeSpecialminus11in->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_mdtSPZR_LargeSpecialminus11in->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                         break;
                         case 7://minus 15out
                         m_h_off_ptvsSA_resptLargeSpecialminus15out->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numSP == 2)m_h_off_ptvsSA_resptLargeSpecialminus15out_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numSP == 3)m_h_off_ptvsSA_resptLargeSpecialminus15out_3station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptLargeSpecialminus15out->Fill(resSA_pt);
                         m_h_offetavsSA_resptLargeSpecialminus15out->Fill(m_poff_eta,resSA_pt);
                         m_h_eSA_pt_LargeSpecialminus15out->Fill(std::fabs(m_poff_pt*0.001));
                         for(Int_t index = 0;index < 10;index++){
                          if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialminus15out->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                          if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 2)m_h_etaIndexvsSA_resptLargeSpecialminus15out_2station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                          if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 3)m_h_etaIndexvsSA_resptLargeSpecialminus15out_3station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                  //if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialminus15out->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                     }
                     m_h_mdtSPZR_LargeSpecialminus15out->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                     m_h_mdtSPZR_LargeSpecialminus15out->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                     m_h_mdtSPZR_LargeSpecialminus15out->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                     break;
                         case 8://minus 15in
                         m_h_off_ptvsSA_resptLargeSpecialminus15in->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numSP == 2)m_h_off_ptvsSA_resptLargeSpecialminus15in_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numSP == 3)m_h_off_ptvsSA_resptLargeSpecialminus15in_3station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptLargeSpecialminus15in->Fill(resSA_pt);
                         m_h_offetavsSA_resptLargeSpecialminus15in->Fill(m_poff_eta,resSA_pt);
                         m_h_eSA_pt_LargeSpecialminus15in->Fill(std::fabs(m_poff_pt*0.001));
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialminus15in->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 2)m_h_etaIndexvsSA_resptLargeSpecialminus15in_2station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numSP == 3)m_h_etaIndexvsSA_resptLargeSpecialminus15in_3station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                   //if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialminus15in->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         }
                         m_h_mdtSPZR_LargeSpecialminus15in->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_mdtSPZR_LargeSpecialminus15in->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_mdtSPZR_LargeSpecialminus15in->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                         break;
                         default :
                         break;
                    }

                    case 2:
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0){//Qeta = +1
                         m_h_off_ptvsSA_resptSmallplus->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptSmallplus->Fill(resSA_pt);
                         //m_h_eSA_pt_Smallplus->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_offphivsSA_resptSmallplus->Fill(m_poff_phi,resSA_pt);
                         m_h_offetavsSA_resptSmallplus->Fill(m_poff_eta,resSA_pt);
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptSmallplus->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         }
                    }
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) < 0){//Qeta = -1
                         m_h_off_ptvsSA_resptSmallminus->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptSmallminus->Fill(resSA_pt);
                         m_h_eSA_pt_Smallminus->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_offphivsSA_resptSmallminus->Fill(m_poff_phi,resSA_pt);
                         m_h_offetavsSA_resptSmallminus->Fill(m_poff_eta,resSA_pt);
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptSmallminus->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         }
                    }
                    break;
                    case 3:
                    m_h_saroiphi_SmallSpecial->Fill(pSA_roiphi);
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0){//Qeta = +1
                         m_h_off_ptvsSA_resptSmallSpecialplus->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptSmallSpecialplus->Fill(resSA_pt);
                         m_h_eSA_pt_SmallSpecialplus->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_offphivsSA_resptSmallSpecialplus->Fill(m_poff_phi,resSA_pt);
                         m_h_offetavsSA_resptSmallSpecialplus->Fill(m_poff_eta,resSA_pt);
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptSmallSpecialplus->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         }
                    }
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) < 0){//Qeta = -1
                         m_h_off_ptvsSA_resptSmallSpecialminus->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SA_resptSmallSpecialminus->Fill(resSA_pt);
                         m_h_eSA_pt_SmallSpecialminus->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_offphivsSA_resptSmallSpecialminus->Fill(m_poff_phi,resSA_pt);
                         m_h_offetavsSA_resptSmallSpecialminus->Fill(m_poff_eta,resSA_pt);
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptSmallSpecialminus->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         }
                    }
                    break;
                    default:
                    break;
               }


               for(Int_t mdthit = 0; mdthit < (signed int)pSA_mdtPhi->size(); mdthit++){
                    m_h_mdtphi->Fill(pSA_mdtPhi->at(mdthit));
               }

               if(numSP == 2){
                    if(pSA_sAddress == 0 && std::fabs(m_poff_phi) < TMath::Pi()/8.0)m_h_off_ptvsSA_resptLarge_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                    if(pSA_sAddress == 2 && std::fabs(m_poff_phi) < TMath::Pi()/8.0)m_h_off_ptvsSA_resptSmall_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
               }

               if(pSA_sAddress == 1){
                    for(Int_t mdthit = 0; mdthit < (signed int)pSA_mdtZ->size(); mdthit++){
                         m_h_mdtphi_LS->Fill(pSA_mdtPhi->at(mdthit));
                         m_h_mdtR->Fill(pSA_mdtR->at(mdthit));
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
               if(pSA_sAddress == 0 || pSA_sAddress == 2)m_h_overphi->Fill(overphi);
               Int_t numeta = 0;
               Int_t numphi = 0;
               if(pSA_eta < 0)numeta = 8;
               if(pSA_sAddress == 0)overphi -= TMath::Pi()/8.0;
               numphi = overphi/TMath::Pi()*32.0 - fmod(overphi,TMath::Pi()/32.0);
               numeta += std::fabs(pSA_eta)/0.125 - fmod(std::fabs(pSA_eta),0.125);
               if(pSA_sAddress == 0 && std::fabs(pSA_eta) < 1.0)m_h_divideetaoverphi_resptLarge[numeta][numphi]->Fill(resSA_pt);
               if(pSA_sAddress == 2 && std::fabs(pSA_eta) < 1.0)m_h_divideetaoverphi_resptSmall[numeta][numphi]->Fill(resSA_pt);

               if(PlateauCut(std::fabs(m_poff_pt*0.001))){
                    m_h_eff_pSA_etaphi->Fill(m_poff_eta,m_poff_phi);
                    m_h_eSA_eta->Fill(m_poff_eta);
                    m_h_eSA_phi->Fill(m_poff_phi);
               //m_h_eSA_aipc->Fill(m_aipc);
               }
               m_h_poffvsSA_pt->Fill(std::fabs(m_poff_pt*0.001),std::fabs(pSA_pt));

               if(pSA_sAddress == 1){
                    for(Int_t index = 0;index < 10;index++){
                         if(m_probe_segment_chamberIndex[index] == 1 && sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)) > 4000 && fabs(m_probe_segment_x[index]) < 4500.){
                              m_h_nPrecisionHits_normal->Fill(m_probe_segment_nPrecisionHits[index]);
                              m_h_sector_normal->Fill(m_probe_segment_sector[index]);
                              m_h_etaIndex_normal->Fill(m_probe_segment_etaIndex[index]);
                              if(m_probe_segment_sector[index] == 11)m_h_segmentnumberBIR->Fill((Int_t)m_probe_segment_etaIndex[index] - 1,0);
                              if(m_probe_segment_sector[index] == 15)m_h_segmentnumberBIR->Fill((Int_t)m_probe_segment_etaIndex[index] - 1,1);
                         }
                         if(m_probe_segment_chamberIndex[index] == 1 && sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)) > 4000 && fabs(m_probe_segment_x[index]) > 4500.){
                              m_h_nPrecisionHits_special->Fill(m_probe_segment_nPrecisionHits[index]);
                              m_h_sector_special->Fill(m_probe_segment_sector[index]);
                              m_h_etaIndex_special->Fill(m_probe_segment_etaIndex[index]);
                              if(m_probe_segment_sector[index] == 11)m_h_segmentnumberBIM->Fill((Int_t)m_probe_segment_etaIndex[index] - 1,0);
                              if(m_probe_segment_sector[index] == 15)m_h_segmentnumberBIM->Fill((Int_t)m_probe_segment_etaIndex[index] - 1,1);
                         }
                    }
               }

               if(pSA_sAddress == 0 && std::fabs(m_poff_phi) < TMath::Pi()/8.0){
                    if(numSP == 2)m_h_offptLarge_2station->Fill(std::fabs(m_poff_pt*0.001));
                    if(numSP == 2)m_h_offetaLarge_2station->Fill(m_poff_eta);
                    if(numSP == 3)m_h_offetaLarge_3station->Fill(m_poff_eta);
                    switch(patternSP){
                         case 3:
                         if(numSP == 2){
                              m_h_offptLargeinnmid_2station->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaLargeinnmid_2station->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptLargeinnmid_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         case 4:
                         if(numSP == 2){
                              m_h_offptLargeinnout_2station->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaLargeinnout_2station->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptLargeinnout_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         case 5:
                         if(numSP == 2){
                              m_h_offptLargemidout_2station->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaLargemidout_2station->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptLargemidout_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         default :
                         break;
                    }
               }
               if(pSA_sAddress == 2 && m_poff_phi > 0 && m_poff_phi < TMath::Pi()/4.0){
                    if(numSP == 2)m_h_offptSmall_2station->Fill(std::fabs(m_poff_pt*0.001));
                    if(numSP == 2)m_h_offetaSmall_2station->Fill(m_poff_eta);
                    if(numSP == 3)m_h_offetaSmall_3station->Fill(m_poff_eta);
                    switch(patternSP){
                         case 3:
                         if(numSP == 2){
                              m_h_offptSmallinnmid_2station->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaSmallinnmid_2station->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptSmallinnmid_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         case 4:
                         if(numSP == 2){
                              m_h_offptSmallinnout_2station->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaSmallinnout_2station->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptSmallinnout_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         case 5:
                         if(numSP == 2){
                              m_h_offptSmallmidout_2station->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaSmallmidout_2station->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptSmallmidout_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         default :
                         break;
                    }
               }
               if(areanumber % 2 == 0 && areanumber != 0){
                    if(numSP == 2)m_h_offptBIR_2station->Fill(std::fabs(m_poff_pt*0.001));
                    if(numSP == 2)m_h_offetaBIR_2station->Fill(m_poff_eta);
                    if(numSP == 3)m_h_offetaBIR_3station->Fill(m_poff_eta);
                    switch(patternSP){
                         case 3:
                         if(numSP == 2){
                              m_h_offptBIRinnmid_2station->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaBIRinnmid_2station->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptBIRinnmid_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         case 4:
                         if(numSP == 2){
                              m_h_offptBIRinnout_2station->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaBIRinnout_2station->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptBIRinnout_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         case 5:
                         if(numSP == 2){
                              m_h_offptBIRmidout_2station->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaBIRmidout_2station->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptBIRmidout_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         default :
                         break;
                    }
               }
               if(areanumber % 2 == 1){
                    if(numSP == 2)m_h_offptBIM_2station->Fill(std::fabs(m_poff_pt*0.001));
                    if(numSP == 2)m_h_offetaBIM_2station->Fill(m_poff_eta);
                    if(numSP == 3)m_h_offetaBIM_3station->Fill(m_poff_eta);
                    switch(patternSP){
                         case 3:
                         if(numSP == 2){
                              m_h_offptBIMinnmid_2station->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaBIMinnmid_2station->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptBIMinnmid_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         case 4:
                         if(numSP == 2){
                              m_h_offptBIMinnout_2station->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaBIMinnout_2station->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptBIMinnout_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         case 5:
                         if(numSP == 2){
                              m_h_offptBIMmidout_2station->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_offetaBIMmidout_2station->Fill(m_poff_eta);
                              m_h_off_ptvsSA_resptBIMmidout_2station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
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
                    if(m_probe_segment_chamberIndex[index] == 1 && segdisBIR == 2)m_h_BIRsegment->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                    if(m_probe_segment_chamberIndex[index] == 1 && segdisBIM == 2)m_h_BIMsegment->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
               }

               for(Int_t chnum = 0; chnum < 10; chnum++){
                    if(m_probe_segment_chamberIndex[chnum] == 0)m_h_segmentXY_BIS->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 1)m_h_segmentXY_BIL->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 2)m_h_segmentXY_BMS->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 3)m_h_segmentXY_BML->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 4)m_h_segmentXY_BOS->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 5)m_h_segmentXY_BOL->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 6)m_h_segmentXY_BEE->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 7)m_h_segmentXY_EIS->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 8)m_h_segmentXY_EIL->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 9)m_h_segmentXY_EMS->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 10)m_h_segmentXY_EML->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 11)m_h_segmentXY_EOS->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 12)m_h_segmentXY_EOL->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 13)m_h_segmentXY_EES->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 14)m_h_segmentXY_EEL->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 15)m_h_segmentXY_CSS->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 16)m_h_segmentXY_CSL->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == -1)m_h_segmentXY_unknown->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
               }

               if(static_cast<Int_t>(pSA_sAddress) == 0 || static_cast<Int_t>(pSA_sAddress) == 1 || static_cast<Int_t>(pSA_sAddress) == 2 || static_cast<Int_t>(pSA_sAddress) == 3){
                    m_h_offphivsSA_sAddress->Fill(pSA_phims,pSA_sAddress);
                    m_h_mdtSPXY_BI->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                    m_h_mdtSPXY_BM->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                    m_h_mdtSPXY_BO->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                    m_h_mdtSPXY_BME->Fill(pSA_superpointR_BME*cos(pSA_roiphi),pSA_superpointR_BME*sin(pSA_roiphi));
                    m_h_mdtSPZR->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                    m_h_mdtSPZR->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                    m_h_mdtSPZR->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                    for(Int_t index = 0;index < 10;index++){
                         Int_t buf_index = static_cast<Int_t>(m_probe_segment_etaIndex[index]);
                         if(m_probe_segment_chamberIndex[index] <= 6){
                              switch(buf_index){
                                   case -6:
                                   m_h_mdtSPXY_etaIndexminus6->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus6->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus6->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexminus6->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case -5:
                                   m_h_mdtSPXY_etaIndexminus5->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus5->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus5->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexminus5->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case -4:
                                   m_h_mdtSPXY_etaIndexminus4->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus4->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus4->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexminus4->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case -3:
                                   m_h_mdtSPXY_etaIndexminus3->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus3->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus3->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexminus3->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case -2:
                                   m_h_mdtSPXY_etaIndexminus2->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus2->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus2->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexminus2->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case -1:
                                   m_h_mdtSPXY_etaIndexminus1->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus1->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexminus1->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexminus1->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case 1:
                                   m_h_mdtSPXY_etaIndexplus1->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus1->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus1->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexplus1->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case 2:
                                   m_h_mdtSPXY_etaIndexplus2->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus2->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus2->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexplus2->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case 3:
                                   m_h_mdtSPXY_etaIndexplus3->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus3->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus3->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexplus3->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case 4:
                                   m_h_mdtSPXY_etaIndexplus4->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus4->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus4->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexplus4->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case 5:
                                   m_h_mdtSPXY_etaIndexplus5->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus5->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus5->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexplus5->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case 6:
                                   m_h_mdtSPXY_etaIndexplus6->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus6->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_mdtSPXY_etaIndexplus6->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY_etaIndexplus6->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   default:
                                   break;
                              }
                         }
                    }
                    for(Int_t size = 0;size < (signed int)pSA_mdtZ->size();size++){
                         if(pSA_mdthitChamber->at(size) == 0 || pSA_mdthitChamber->at(size) == 1 || pSA_mdthitChamber->at(size) == 2)m_h_mdthitXY->Fill(pSA_mdtR->at(size)*cos(pSA_mdtPhi->at(size)),pSA_mdtR->at(size)*sin(pSA_mdtPhi->at(size)));
                         m_h_mdthitZR->Fill(pSA_mdtZ->at(size),pSA_mdtR->at(size));
                    }
                    for(Int_t size = 0;size < (signed int)pSA_rpcX->size();size++){
                         m_h_rpchitXY->Fill(pSA_rpcX->at(size),pSA_rpcY->at(size));
                    }

                    Int_t buf_numsegment = 0;
                    for(Int_t index = 0;index < 10;index++){
                         if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_respt->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                         if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0){
                         //m_h_segmentZR->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         //if(m_probe_segment_chamberIndex[index] == 1)m_h_segmentZR_BIL->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         //if(m_probe_segment_chamberIndex[index] == 3)m_h_segmentZR_BML->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              buf_numsegment++;
                         }
                    }
                    m_h_numsegment->Fill(buf_numsegment);
               }
               m_h_offphivsSAphims->Fill(m_poff_phi,pSA_phims);

          //CB
               if(!CutCB(pCB_pass))return;
               Double_t textCB_dR = TMath::Sqrt(pow(m_tCB_eta - m_toff_exteta,2) + pow(m_tCB_phi - m_toff_extphi,2));
               pextCB_dR = TMath::Sqrt(pow(pCB_eta - m_poff_exteta,2) + pow(pCB_phi - m_poff_extphi,2));
               Double_t resCB_pt = std::fabs(m_poff_pt)/std::fabs(pCB_pt) - 1.0;
               m_h_pCB_pt->Fill(std::fabs(pCB_pt*0.001));
               m_h_pCB_dR->Fill(pCB_dR);
               m_h_textCB_dR->Fill(textCB_dR);
               m_h_pextCB_dR->Fill(pextCB_dR);
               m_h_eCB_pt->Fill(std::fabs(m_poff_pt*0.001));
               if(PlateauCut(std::fabs(m_poff_pt*0.001))){
                    m_h_eCB_eta->Fill(m_poff_eta);
                    m_h_eCB_phi->Fill(m_poff_phi);
               //m_h_eCB_aipc->Fill(m_aipc);
               }
               m_h_pCB_respt->Fill(resCB_pt);
               if(DicisionBarrel(m_poff_eta)){
                    m_h_eCB_pt_barrel->Fill(std::fabs(m_poff_pt*0.001));
                    m_h_pCB_respt_barrel->Fill(resCB_pt);
               }else{
                    m_h_eCB_pt_end->Fill(std::fabs(m_poff_pt*0.001));
                    m_h_pCB_respt_endcap->Fill(resCB_pt);
               }

          //EF
               if(!CutEF(pEF_pass))return;
               Double_t textEF_dR = TMath::Sqrt(pow(m_tEF_eta - m_toff_exteta,2) + pow(m_tEF_phi - m_toff_extphi,2));
               pextEF_dR = TMath::Sqrt(pow(pEF_eta - m_poff_exteta,2) + pow(pEF_phi - m_poff_extphi,2));
               Double_t resEF_pt = std::fabs(m_poff_pt)/std::fabs(pEF_pt) - 1.0;
               m_h_pEF_pt->Fill(std::fabs(pEF_pt*0.001));
               m_h_pEF_dR->Fill(pEF_dR);
               m_h_textEF_dR->Fill(textEF_dR);
               m_h_pextEF_dR->Fill(pextEF_dR);
               m_h_eEF_pt->Fill(std::fabs(m_poff_pt*0.001));
               if(PlateauCut(std::fabs(m_poff_pt*0.001))){
                    m_h_eEF_eta->Fill(m_poff_eta);
                    m_h_eEF_phi->Fill(m_poff_phi);
               //m_h_eEF_aipc->Fill(m_aipc);
               }
               m_h_pEF_respt->Fill(resEF_pt);
               if(DicisionBarrel(m_poff_eta)){
                    m_h_eEF_pt_barrel->Fill(std::fabs(m_poff_pt*0.001));
                    m_h_pEF_respt_barrel->Fill(resEF_pt);
               }else{
                    m_h_eEF_pt_end->Fill(std::fabs(m_poff_pt*0.001));
                    m_h_pEF_respt_endcap->Fill(resEF_pt);
               }

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

          //base,target
     cout<<"efficiency start!"<<endl;

     ceff.SetCondition("L1Efficiency","L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eoff_pt,m_h_eL1_pt,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiency","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt,m_h_eSA_pt,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("CBEfficiency","muComb Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eSA_pt,m_h_eCB_pt,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("EFEfficiency","EventFilter Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eCB_pt,m_h_eEF_pt,m_binmax,200,m_efficiency_xerr);

     ceff.SetCondition("L1Efficiency_eta","L1 Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyeta(m_h_eoff_eta,m_h_eL1_eta);
     ceff.SetCondition("SAEfficiency_eta","L2MuonSA Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyeta(m_h_eL1_eta,m_h_eSA_eta);
     ceff.SetCondition("CBEfficiency_eta","muComb Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyeta(m_h_eSA_eta,m_h_eCB_eta);
     ceff.SetCondition("EFEfficiency_eta","EventFilter Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyeta(m_h_eCB_eta,m_h_eEF_eta);
     ceff.SetCondition("L1Efficiency_phi","L1 Efficiency;offline phi;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyphi(m_h_eoff_phi,m_h_eL1_phi);
     ceff.SetCondition("SAEfficiency_phi","L2MuonSA Efficiency;offline phi;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyphi(m_h_eL1_phi,m_h_eSA_phi);
     ceff.SetCondition("CBEfficiency_phi","muComb Efficiency;offline phi;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyphi(m_h_eSA_phi,m_h_eCB_phi);
     ceff.SetCondition("EFEfficiency_phi","EventFilter Efficiency;offline phi;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyphi(m_h_eCB_phi,m_h_eEF_phi);
          /*
          ceff.SetCondition("L1Efficiency_pileup","L1 Efficiency;pileup;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencypileup(m_h_eoff_aipc,m_h_eL1_aipc);
          ceff.SetCondition("SAEfficiency_pileup","L2MuonSA Efficiency;pileup;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencypileup(m_h_eL1_aipc,m_h_eSA_aipc);
          ceff.SetCondition("CBEfficiency_pileup","muComb Efficiency;pileup;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencypileup(m_h_eSA_aipc,m_h_eCB_aipc);
          ceff.SetCondition("EFEfficiency_pileup","EventFilter Efficiency;pileup;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencypileup(m_h_eCB_aipc,m_h_eEF_aipc);
          */
     ceff.SetCondition("L1Efficiency_barrel","L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eoff_pt_barrel,m_h_eL1_pt_barrel,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiency_barrel","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_barrel,m_h_eSA_pt_barrel,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("CBEfficiency_barrel","muComb Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eSA_pt_barrel,m_h_eCB_pt_barrel,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("EFEfficiency_barrel","EventFilter Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eCB_pt_barrel,m_h_eEF_pt_barrel,m_binmax,300,m_efficiency_xerr);

     ceff.SetCondition("L1Efficiency_end","L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eoff_pt_end,m_h_eL1_pt_end,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiency_end","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_end,m_h_eSA_pt_end,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("CBEfficiency_end","muComb Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eSA_pt_end,m_h_eCB_pt_end,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("EFEfficiency_end","EventFilter Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eCB_pt_end,m_h_eEF_pt_end,m_binmax,300,m_efficiency_xerr);

     ceff.SetCondition("L1EfficiencyLargeplus","L1 Large Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eoff_pt_Largeplus,m_h_eL1_pt_Largeplus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeplus","SA Large Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_Largeplus,m_h_eSA_pt_Largeplus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargenormal","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_Largenormal,m_h_eSA_pt_Largenormal,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("L1EfficiencyLargeSpecialplus","L1 LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eoff_pt_LargeSpecialplus,m_h_eL1_pt_LargeSpecialplus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeSpecialplus","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus,m_h_eSA_pt_LargeSpecialplus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("L1EfficiencyLargeminus","L1 Large Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eoff_pt_Largeminus,m_h_eL1_pt_Largeminus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeminus","SA Large Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_Largeminus,m_h_eSA_pt_Largeminus,m_binmax,300,m_efficiency_xerr);

     ceff.SetCondition("L1EfficiencyLargeSpecialminus","L1 LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eoff_pt_LargeSpecialminus,m_h_eL1_pt_LargeSpecialminus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeSpecialminus","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus,m_h_eSA_pt_LargeSpecialminus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeSpecialplus11","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus11,m_h_eSA_pt_LargeSpecialplus11,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeSpecialminus11","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus11,m_h_eSA_pt_LargeSpecialminus11,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeSpecialplus15","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus15,m_h_eSA_pt_LargeSpecialplus15,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeSpecialminus15","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus15,m_h_eSA_pt_LargeSpecialminus15,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeSpecialplus11out","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus11out,m_h_eSA_pt_LargeSpecialplus11out,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeSpecialminus11out","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus11out,m_h_eSA_pt_LargeSpecialminus11out,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeSpecialplus15out","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus15out,m_h_eSA_pt_LargeSpecialplus15out,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeSpecialminus15out","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus15out,m_h_eSA_pt_LargeSpecialminus15out,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeSpecialplus11in","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus11in,m_h_eSA_pt_LargeSpecialplus11in,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeSpecialminus11in","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus11in,m_h_eSA_pt_LargeSpecialminus11in,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeSpecialplus15in","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialplus15in,m_h_eSA_pt_LargeSpecialplus15in,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeSpecialminus15in","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_LargeSpecialminus15in,m_h_eSA_pt_LargeSpecialminus15in,m_binmax,300,m_efficiency_xerr);

     ceff.SetCondition("L1EfficiencySmallplus","L1 Small Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eoff_pt_Smallplus,m_h_eL1_pt_Smallplus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencySmallplus","SA Small Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_Smallplus,m_h_eSA_pt_Smallplus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("L1EfficiencySmallminus","L1 Small Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eoff_pt_Smallminus,m_h_eL1_pt_Smallminus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencySmallminus","SA Small Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_Smallminus,m_h_eSA_pt_Smallminus,m_binmax,300,m_efficiency_xerr);

     ceff.SetCondition("L1EfficiencySmallSpecialplus","L1 SmallSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eoff_pt_SmallSpecialplus,m_h_eL1_pt_SmallSpecialplus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencySmallSpecialplus","SA SmallSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_SmallSpecialplus,m_h_eSA_pt_SmallSpecialplus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("L1EfficiencySmallSpecialminus","L1 SmallSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eoff_pt_SmallSpecialminus,m_h_eL1_pt_SmallSpecialminus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencySmallSpecialminus","SA SmallSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_SmallSpecialminus,m_h_eSA_pt_SmallSpecialminus,m_binmax,300,m_efficiency_xerr);

     ceff.SetCondition("SAEfficiencyBarrelwithoutBIM","SA Barrel without BIM Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_BarrelwithoutBIM,m_h_eSA_pt_BarrelwithoutBIM,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLSwithoutBIM","SA LS without BIM Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_LSwithoutBIM,m_h_eSA_pt_LSwithoutBIM,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeincBIM","SA Barrel inc BIM Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_BarrelincBIM,m_h_eSA_pt_BarrelincBIM,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLSincBIM","SA LS inc BIM Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1_pt_LSincBIM,m_h_eSA_pt_LSincBIM,m_binmax,300,m_efficiency_xerr);

     ceff.SetCondition("SA2DEfficiency","L1vsL2MuonSA Efficiency;offline eta;offline phi",1.5,0.1,0.1,0.105,0.165);
     ceff.SetConditionbin(m_nbin_eta,m_nbin_phi,m_eta_max,m_phi_max);
     ceff.DrawEfficiency2D(m_h_eff_pL1_etaphi,m_h_eff_pSA_etaphi);
     ceff.SetCondition("L12DEfficiency","L1 Efficiency;offline eta;offline phi",1.5,0.1,0.1,0.105,0.165);
     ceff.SetConditionbin(m_nbin_eta,m_nbin_phi,m_eta_max,m_phi_max);
     ceff.DrawEfficiency2D(m_h_eff_poff_etaphi,m_h_eff_pL1_etaphi);

     ceff.SetCondition("SAEfficiency_SP1_barrel","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1SP1_pt_barrel,m_h_eSASP1_pt_barrel,m_binmax,200,m_efficiency_xerr); 
     ceff.SetCondition("SAEfficiency_SP2_barrel","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1SP2_pt_barrel,m_h_eSASP2_pt_barrel,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiency_SP3_barrel","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1SP3_pt_barrel,m_h_eSASP3_pt_barrel,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiency_innmid_barrel","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1innmid_pt_barrel,m_h_eSAinnmid_pt_barrel,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiency_SP1_endcap","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1SP1_pt_endcap,m_h_eSASP1_pt_endcap,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiency_SP2_endcap","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1SP2_pt_endcap,m_h_eSASP2_pt_endcap,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiency_SP3_endcap","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1SP3_pt_endcap,m_h_eSASP3_pt_endcap,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiency_innmid_endcap","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1innmid_pt_endcap,m_h_eSAinnmid_pt_endcap,m_binmax,200,m_efficiency_xerr);
     cout<<"eff end"<<endl;

     m_h_poff_pt->Write();
     m_h_pL1_pt->Write();
     m_h_pSA_pt->Write();
     m_h_pCB_pt->Write();
     m_h_pEF_pt->Write();
     m_h_pL1_dR->Write();
     m_h_pSA_dR->Write();
     m_h_pCB_dR->Write();
     m_h_pEF_dR->Write();
     m_h_textL1_dR->Write();
     m_h_textSA_dR->Write();
     m_h_textCB_dR->Write();
     m_h_textEF_dR->Write();
     m_h_pextL1_dR->Write();
     m_h_pextSA_dR->Write();
     m_h_pextCB_dR->Write();
     m_h_pextEF_dR->Write();
     m_h_poffvsSA_pt->Write();
     m_h_offphivsSA_sAddress->Write();
     m_h_offphivsSAphims->Write();
     m_h_pSAphivspSAphims->Write();
     m_h_pSAphivspSAphibe->Write();
     m_h_countSA->Write();
     m_h_ChIndexvshitChamber->Write();
          //
     m_h_offptBIM_2station->Write();
     m_h_offptBIR_2station->Write();
     m_h_offptLarge_2station->Write();
     m_h_offptSmall_2station->Write();
     m_h_offptBIMinnmid_2station->Write();
     m_h_offptBIRinnmid_2station->Write();
     m_h_offptLargeinnmid_2station->Write();
     m_h_offptSmallinnmid_2station->Write();
     m_h_offptBIMinnout_2station->Write();
     m_h_offptBIRinnout_2station->Write();
     m_h_offptLargeinnout_2station->Write();
     m_h_offptSmallinnout_2station->Write();
     m_h_offptBIMmidout_2station->Write();
     m_h_offptBIRmidout_2station->Write();
     m_h_offptLargemidout_2station->Write();
     m_h_offptSmallmidout_2station->Write();
          //
     m_h_offetaBIM_2station->Write();
     m_h_offetaBIR_2station->Write();
     m_h_offetaLarge_2station->Write();
     m_h_offetaSmall_2station->Write();
     m_h_offetaBIMinnmid_2station->Write();
     m_h_offetaBIRinnmid_2station->Write();
     m_h_offetaLargeinnmid_2station->Write();
     m_h_offetaSmallinnmid_2station->Write();
     m_h_offetaBIMinnout_2station->Write();
     m_h_offetaBIRinnout_2station->Write();
     m_h_offetaLargeinnout_2station->Write();
     m_h_offetaSmallinnout_2station->Write();
     m_h_offetaBIMmidout_2station->Write();
     m_h_offetaBIRmidout_2station->Write();
     m_h_offetaLargemidout_2station->Write();
     m_h_offetaSmallmidout_2station->Write();
          //
     m_h_offetaBIM_3station->Write();
     m_h_offetaBIR_3station->Write();
     m_h_offetaLarge_3station->Write();
     m_h_offetaSmall_3station->Write();

     m_h_sectorvsphi->Write();
     m_h_indexvseta->Write();
     m_h_overphi->Write();
     m_h_mdtchamber->Write();
     m_h_L1pSA_sAddress->Write();
     m_h_SApSA_sAddress->Write();

     m_h_mdthitXY->Write();
     m_h_rpchitXY->Write();
     m_h_mdthitZR->Write();
     m_h_mdtSPXY_BI->Write();
     m_h_mdtSPXY_BM->Write();
     m_h_mdtSPXY_BO->Write();
     m_h_mdtSPXY_BME->Write();
     m_h_mdtSPX_BI->Write();
     m_h_mdtSPY_BI->Write();
     m_h_mdtSPXY_etaIndexplus1->Write();
     m_h_mdtSPXY_etaIndexplus2->Write();
     m_h_mdtSPXY_etaIndexplus3->Write();
     m_h_mdtSPXY_etaIndexplus4->Write();
     m_h_mdtSPXY_etaIndexplus5->Write();
     m_h_mdtSPXY_etaIndexplus6->Write();
     m_h_mdtSPXY_etaIndexminus1->Write();
     m_h_mdtSPXY_etaIndexminus2->Write();
     m_h_mdtSPXY_etaIndexminus3->Write();
     m_h_mdtSPXY_etaIndexminus4->Write();
     m_h_mdtSPXY_etaIndexminus5->Write();
     m_h_mdtSPXY_etaIndexminus6->Write();
     m_h_mdtSPZR->Write();
     m_h_mdtSPZR_LargeSpecialplus->Write();
     m_h_mdtSPZR_LargeSpecialminus->Write();
     m_h_mdtSPZR_LargeSpecialplus11out->Write();
     m_h_mdtSPZR_LargeSpecialplus11in->Write();
     m_h_mdtSPZR_LargeSpecialplus15out->Write();
     m_h_mdtSPZR_LargeSpecialplus15in->Write();
     m_h_mdtSPZR_LargeSpecialminus11out->Write();
     m_h_mdtSPZR_LargeSpecialminus11in->Write();
     m_h_mdtSPZR_LargeSpecialminus15out->Write();
     m_h_mdtSPZR_LargeSpecialminus15in->Write();
     m_h_segmentXY->Write();
     m_h_segmentXYmatch->Write();
     m_h_segmentXYnomatch->Write();
     m_h_segmentXY_etaIndexplus1->Write();
     m_h_segmentXY_etaIndexplus2->Write();
     m_h_segmentXY_etaIndexplus3->Write();
     m_h_segmentXY_etaIndexplus4->Write();
     m_h_segmentXY_etaIndexplus5->Write();
     m_h_segmentXY_etaIndexplus6->Write();
     m_h_segmentXY_etaIndexminus1->Write();
     m_h_segmentXY_etaIndexminus2->Write();
     m_h_segmentXY_etaIndexminus3->Write();
     m_h_segmentXY_etaIndexminus4->Write();
     m_h_segmentXY_etaIndexminus5->Write();
     m_h_segmentXY_etaIndexminus6->Write();
     m_h_segmentXY_BIL->Write();
     m_h_segmentXY_BML->Write();
     m_h_segmentXY_BOL->Write();
     m_h_segmentXY_BIS->Write();
     m_h_segmentXY_BMS->Write();
     m_h_segmentXY_BOS->Write();
     m_h_segmentXY_BEE->Write();
     m_h_segmentXY_EIL->Write();
     m_h_segmentXY_EML->Write();
     m_h_segmentXY_EOL->Write();
     m_h_segmentXY_EIS->Write();
     m_h_segmentXY_EMS->Write();
     m_h_segmentXY_EOS->Write();
     m_h_segmentXY_EES->Write();
     m_h_segmentXY_EEL->Write();
     m_h_segmentXY_CSS->Write();
     m_h_segmentXY_CSL->Write();
     m_h_segmentXY_unknown->Write();
     m_h_segmentXY_normal->Write();
     m_h_segmentXY_special->Write();
          /*
          m_h_segmentZR->Write();
          m_h_segmentZR_BIL->Write();
          m_h_segmentZR_BML->Write();
          m_h_segmentZR_LargeSpecialplus->Write();
          m_h_segmentZR_LargeSpecialminus->Write();
          m_h_segmentZR_LargeSpecialplus11out->Write();
          m_h_segmentZR_LargeSpecialplus11in->Write();
          m_h_segmentZR_LargeSpecialplus15out->Write();
          m_h_segmentZR_LargeSpecialplus15in->Write();
          m_h_segmentZR_LargeSpecialminus11out->Write();
          m_h_segmentZR_LargeSpecialminus11in->Write();
          m_h_segmentZR_LargeSpecialminus15out->Write();
          m_h_segmentZR_LargeSpecialminus15in->Write();
          m_h_segmentZR_LargeSpecialplus11out_BIL->Write();
          m_h_segmentZR_LargeSpecialplus11in_BIL->Write();
          m_h_segmentZR_LargeSpecialplus15out_BIL->Write();
          m_h_segmentZR_LargeSpecialplus15in_BIL->Write();
          */
     m_h_numsegment->Write();
     m_h_num_segment_LSBI->Write();
     m_h_num_segment_LargeBI->Write();
     m_h_numhit->Write();
     m_h_nPrecisionHits_normal->Write();
     m_h_nPrecisionHits_special->Write();
     m_h_sector_normal->Write();
     m_h_sector_special->Write();
     m_h_etaIndex_normal->Write();
     m_h_etaIndex_special->Write();
     m_h_mdtphi->Write();
     m_h_mdtphi_LS->Write();
     m_h_mdtphi_LSBIL->Write();
     m_h_mdtR->Write();
     m_h_numSP->Write();
     m_h_BIRsegment->Write();
     m_h_BIMsegment->Write();
     m_h_BIMrvsx->Write();

     m_h_eoff_pt->Write();
     m_h_eL1_pt->Write();
     m_h_eSA_pt->Write();
     m_h_eCB_pt->Write();
     m_h_eEF_pt->Write();
     m_h_eoff_eta->Write();
     m_h_eL1_eta->Write();
     m_h_eSA_eta->Write();
     m_h_eCB_eta->Write();
     m_h_eEF_eta->Write();
     m_h_eoff_phi->Write();
     m_h_eL1_phi->Write();
     m_h_eSA_phi->Write();
     m_h_eCB_phi->Write();
     m_h_eEF_phi->Write();
          //m_h_eoff_aipc->Write();
          //m_h_eL1_aipc->Write();
          //m_h_eSA_aipc->Write();
          //m_h_eCB_aipc->Write();
          //m_h_eEF_aipc->Write();
     m_h_eoff_pt_barrel->Write();
     m_h_eL1_pt_barrel->Write();
     m_h_eSA_pt_barrel->Write();
     m_h_eCB_pt_barrel->Write();
     m_h_eEF_pt_barrel->Write();
     m_h_eoff_pt_end->Write();
     m_h_eL1_pt_end->Write();
     m_h_eSA_pt_end->Write();
     m_h_eCB_pt_end->Write();
     m_h_eEF_pt_end->Write();
     m_h_eoff_pt_Largeplus->Write();
     m_h_eL1_pt_Largeplus->Write();
     m_h_eSA_pt_Largeplus->Write();
     m_h_eL1_pt_Largenormal->Write();
     m_h_eSA_pt_Largenormal->Write();
     m_h_eoff_pt_Largeminus->Write();
     m_h_eL1_pt_Largeminus->Write();
     m_h_eSA_pt_Largeminus->Write();
     m_h_eoff_pt_LargeSpecialplus->Write();
     m_h_eL1_pt_LargeSpecialplus->Write();
     m_h_eSA_pt_LargeSpecialplus->Write();
     m_h_eL1_pt_LargeSpecialplus11->Write();
     m_h_eSA_pt_LargeSpecialplus11->Write();
     m_h_eL1_pt_LargeSpecialplus15->Write();
     m_h_eSA_pt_LargeSpecialplus15->Write();
     m_h_eL1_pt_LargeSpecialplus11in->Write();
     m_h_eSA_pt_LargeSpecialplus11in->Write();
     m_h_eL1_pt_LargeSpecialplus11out->Write();
     m_h_eSA_pt_LargeSpecialplus11out->Write();
     m_h_eL1_pt_LargeSpecialplus15out->Write();
     m_h_eSA_pt_LargeSpecialplus15out->Write();
     m_h_eL1_pt_LargeSpecialplus15in->Write();
     m_h_eSA_pt_LargeSpecialplus15in->Write();
     m_h_eoff_pt_LargeSpecialminus->Write();
     m_h_eL1_pt_LargeSpecialminus->Write();
     m_h_eSA_pt_LargeSpecialminus->Write();
     m_h_eL1_pt_LargeSpecialminus11->Write();
     m_h_eSA_pt_LargeSpecialminus11->Write();
     m_h_eL1_pt_LargeSpecialminus15->Write();
     m_h_eSA_pt_LargeSpecialminus15->Write();
     m_h_eL1_pt_LargeSpecialminus11in->Write();
     m_h_eSA_pt_LargeSpecialminus11in->Write();
     m_h_eL1_pt_LargeSpecialminus11out->Write();
     m_h_eSA_pt_LargeSpecialminus11out->Write();
     m_h_eL1_pt_LargeSpecialminus15out->Write();
     m_h_eSA_pt_LargeSpecialminus15out->Write();
     m_h_eL1_pt_LargeSpecialminus15in->Write();
     m_h_eSA_pt_LargeSpecialminus15in->Write();
     m_h_eoff_pt_Smallplus->Write();
     m_h_eL1_pt_Smallplus->Write();
     m_h_eSA_pt_Smallplus->Write();
     m_h_eoff_pt_Smallminus->Write();
     m_h_eL1_pt_Smallminus->Write();
     m_h_eSA_pt_Smallminus->Write();
     m_h_eoff_pt_SmallSpecialplus->Write();
     m_h_eL1_pt_SmallSpecialplus->Write();
     m_h_eSA_pt_SmallSpecialplus->Write();
     m_h_eoff_pt_SmallSpecialminus->Write();
     m_h_eL1_pt_SmallSpecialminus->Write();
     m_h_eSA_pt_SmallSpecialminus->Write();
     m_h_eff_poff_etaphi->Write();
     m_h_eff_pL1_etaphi->Write();
     m_h_eff_pSA_etaphi->Write();
     m_h_eL1_pt_BarrelwithoutBIM->Write();
     m_h_eSA_pt_BarrelwithoutBIM->Write();
     m_h_eL1_pt_LSwithoutBIM->Write();
     m_h_eSA_pt_LSwithoutBIM->Write();
     m_h_eL1_pt_BarrelincBIM->Write();
     m_h_eSA_pt_BarrelincBIM->Write();
     m_h_eL1_pt_LSincBIM->Write();
     m_h_eSA_pt_LSincBIM->Write();
     m_h_eL1SP1_pt_barrel->Write();
     m_h_eSASP1_pt_barrel->Write();
     m_h_eL1SP2_pt_barrel->Write();
     m_h_eSASP2_pt_barrel->Write();
     m_h_eL1SP3_pt_barrel->Write();
     m_h_eSASP3_pt_barrel->Write();
     m_h_eL1innmid_pt_barrel->Write();
     m_h_eSAinnmid_pt_barrel->Write();
     m_h_eL1SP1_pt_endcap->Write();
     m_h_eSASP1_pt_endcap->Write();
     m_h_eL1SP2_pt_endcap->Write();
     m_h_eSASP2_pt_endcap->Write();
     m_h_eL1SP3_pt_endcap->Write();
     m_h_eSASP3_pt_endcap->Write();
     m_h_eL1innmid_pt_endcap->Write();
     m_h_eSAinnmid_pt_endcap->Write();

     m_h_pSA_respt->Write();
     m_h_pCB_respt->Write();
     m_h_pEF_respt->Write();
     m_h_pSA_respt_barrel->Write();
     m_h_pCB_respt_barrel->Write();
     m_h_pEF_respt_barrel->Write();
     m_h_pSA_respt_endcap->Write();
     m_h_pCB_respt_endcap->Write();
     m_h_pEF_respt_endcap->Write();
     m_h_SA_resptLargeplus->Write();
     m_h_SA_resptLargeSpecialplus->Write();
     m_h_SA_resptSmallplus->Write();
     m_h_SA_resptSmallSpecialplus->Write();
     m_h_SA_resptLargeminus->Write();
     m_h_SA_resptLargeSpecialminus->Write();
     m_h_SA_resptSmallminus->Write();
     m_h_SA_resptSmallSpecialminus->Write();
     m_h_SA_resptLargeSpecialplus11out->Write();
     m_h_SA_resptLargeSpecialplus11in->Write();
     m_h_SA_resptLargeSpecialplus15out->Write();
     m_h_SA_resptLargeSpecialplus15in->Write();
     m_h_SA_resptLargeSpecialminus11out->Write();
     m_h_SA_resptLargeSpecialminus11in->Write();
     m_h_SA_resptLargeSpecialminus15out->Write();
     m_h_SA_resptLargeSpecialminus15in->Write();
     m_h_off_ptvsSA_resptLargeplus->Write();
     m_h_off_ptvsSA_resptLargeSpecialplus->Write();
     m_h_off_ptvsSA_resptSmallplus->Write();
     m_h_off_ptvsSA_resptSmallSpecialplus->Write();
     m_h_off_ptvsSA_resptLargeminus->Write();
     m_h_off_ptvsSA_resptLargeSpecialminus->Write();
     m_h_off_ptvsSA_resptSmallminus->Write();
     m_h_off_ptvsSA_resptSmallSpecialminus->Write();
     m_h_offphivsSA_resptLargeplus->Write();
     m_h_offphivsSA_resptLargeSpecialplus->Write();
     m_h_offphivsSA_resptSmallplus->Write();
     m_h_offphivsSA_resptSmallSpecialplus->Write();
     m_h_offphivsSA_resptLargeminus->Write();
     m_h_offphivsSA_resptLargeSpecialminus->Write();
     m_h_offphivsSA_resptSmallminus->Write();
     m_h_offphivsSA_resptSmallSpecialminus->Write();
     m_h_off_ptvsSA_resptLargeSpecialplus15out->Write();
     m_h_off_ptvsSA_resptLargeSpecialplus15in->Write();
     m_h_off_ptvsSA_resptLargeSpecialplus11out->Write();
     m_h_off_ptvsSA_resptLargeSpecialplus11in->Write();
     m_h_off_ptvsSA_resptLargeSpecialminus15out->Write();
     m_h_off_ptvsSA_resptLargeSpecialminus15in->Write();
     m_h_off_ptvsSA_resptLargeSpecialminus11out->Write();
     m_h_off_ptvsSA_resptLargeSpecialminus11in->Write();
          //
     m_h_off_ptvsSA_resptLarge_2station->Write();
     m_h_off_ptvsSA_resptSmall_2station->Write();
     m_h_off_ptvsSA_resptLargeSpecialplus15out_2station->Write();
     m_h_off_ptvsSA_resptLargeSpecialplus15in_2station->Write();
     m_h_off_ptvsSA_resptLargeSpecialplus11out_2station->Write();
     m_h_off_ptvsSA_resptLargeSpecialplus11in_2station->Write();
     m_h_off_ptvsSA_resptLargeSpecialminus15out_2station->Write();
     m_h_off_ptvsSA_resptLargeSpecialminus15in_2station->Write();
     m_h_off_ptvsSA_resptLargeSpecialminus11out_2station->Write();
     m_h_off_ptvsSA_resptLargeSpecialminus11in_2station->Write();
          //
     m_h_off_ptvsSA_resptLargeSpecialplus15out_3station->Write();
     m_h_off_ptvsSA_resptLargeSpecialplus15in_3station->Write();
     m_h_off_ptvsSA_resptLargeSpecialplus11out_3station->Write();
     m_h_off_ptvsSA_resptLargeSpecialplus11in_3station->Write();
     m_h_off_ptvsSA_resptLargeSpecialminus15out_3station->Write();
     m_h_off_ptvsSA_resptLargeSpecialminus15in_3station->Write();
     m_h_off_ptvsSA_resptLargeSpecialminus11out_3station->Write();
     m_h_off_ptvsSA_resptLargeSpecialminus11in_3station->Write();
          //
     m_h_off_ptvsSA_resptBIMinnmid_2station->Write();
     m_h_off_ptvsSA_resptBIRinnmid_2station->Write();
     m_h_off_ptvsSA_resptLargeinnmid_2station->Write();
     m_h_off_ptvsSA_resptSmallinnmid_2station->Write();
     m_h_off_ptvsSA_resptBIMinnout_2station->Write();
     m_h_off_ptvsSA_resptBIRinnout_2station->Write();
     m_h_off_ptvsSA_resptLargeinnout_2station->Write();
     m_h_off_ptvsSA_resptSmallinnout_2station->Write();
     m_h_off_ptvsSA_resptBIMmidout_2station->Write();
     m_h_off_ptvsSA_resptBIRmidout_2station->Write();
     m_h_off_ptvsSA_resptLargemidout_2station->Write();
     m_h_off_ptvsSA_resptSmallmidout_2station->Write();
          //
     m_h_offetavsSA_resptLargeplus->Write();
     m_h_offetavsSA_resptLargeSpecialplus->Write();
     m_h_offetavsSA_resptSmallplus->Write();
     m_h_offetavsSA_resptSmallSpecialplus->Write();
     m_h_offetavsSA_resptLargeminus->Write();
     m_h_offetavsSA_resptLargeSpecialminus->Write();
     m_h_offetavsSA_resptSmallminus->Write();
     m_h_offetavsSA_resptSmallSpecialminus->Write();
     m_h_offetavsSA_resptLargeSpecialplus11out->Write();
     m_h_offetavsSA_resptLargeSpecialplus11in->Write();
     m_h_offetavsSA_resptLargeSpecialplus15out->Write();
     m_h_offetavsSA_resptLargeSpecialplus15in->Write();
     m_h_offetavsSA_resptLargeSpecialminus11out->Write();
     m_h_offetavsSA_resptLargeSpecialminus11in->Write();
     m_h_offetavsSA_resptLargeSpecialminus15out->Write();
     m_h_offetavsSA_resptLargeSpecialminus15in->Write();
     m_h_etaIndexvsSA_respt->Write();
     m_h_etaIndexvsSA_resptLargeplus->Write();
     m_h_etaIndexvsSA_resptLargeSpecialplus->Write();
     m_h_etaIndexvsSA_resptSmallplus->Write();
     m_h_etaIndexvsSA_resptSmallSpecialplus->Write();
     m_h_etaIndexvsSA_resptLargeminus->Write();
     m_h_etaIndexvsSA_resptLargeSpecialminus->Write();
     m_h_etaIndexvsSA_resptSmallminus->Write();
     m_h_etaIndexvsSA_resptSmallSpecialminus->Write();
     m_h_etaIndexvsSA_resptLargeSpecialplus11out->Write();
     m_h_etaIndexvsSA_resptLargeSpecialplus11in->Write();
     m_h_etaIndexvsSA_resptLargeSpecialplus15out->Write();
     m_h_etaIndexvsSA_resptLargeSpecialplus15in->Write();
     m_h_etaIndexvsSA_resptLargeSpecialminus11out->Write();
     m_h_etaIndexvsSA_resptLargeSpecialminus11in->Write();
     m_h_etaIndexvsSA_resptLargeSpecialminus15out->Write();
     m_h_etaIndexvsSA_resptLargeSpecialminus15in->Write();
          //
     m_h_etaIndexvsSA_resptLargeSpecialplus11out_2station->Write();
     m_h_etaIndexvsSA_resptLargeSpecialplus11in_2station->Write();
     m_h_etaIndexvsSA_resptLargeSpecialplus15out_2station->Write();
     m_h_etaIndexvsSA_resptLargeSpecialplus15in_2station->Write();
     m_h_etaIndexvsSA_resptLargeSpecialminus11out_2station->Write();
     m_h_etaIndexvsSA_resptLargeSpecialminus11in_2station->Write();
     m_h_etaIndexvsSA_resptLargeSpecialminus15out_2station->Write();
     m_h_etaIndexvsSA_resptLargeSpecialminus15in_2station->Write();
          //
     m_h_etaIndexvsSA_resptLargeSpecialplus11out_3station->Write();
     m_h_etaIndexvsSA_resptLargeSpecialplus11in_3station->Write();
     m_h_etaIndexvsSA_resptLargeSpecialplus15out_3station->Write();
     m_h_etaIndexvsSA_resptLargeSpecialplus15in_3station->Write();
     m_h_etaIndexvsSA_resptLargeSpecialminus11out_3station->Write();
     m_h_etaIndexvsSA_resptLargeSpecialminus11in_3station->Write();
     m_h_etaIndexvsSA_resptLargeSpecialminus15out_3station->Write();
     m_h_etaIndexvsSA_resptLargeSpecialminus15in_3station->Write();

     cout<<"residual end"<<endl;

}
