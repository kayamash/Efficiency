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
const Double_t pt_threshold[4] = {5.17,3.25,4.69,5.14};//MU6
//const Double_t pt_threshold[4] = {15.87,10.73,12.21,15.87};//MU20

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

int Efficiency::EtaDistribution(Float_t roieta){
     if(std::fabs(roieta) < 1.05){
          return 0;
     }else if(std::fabs(roieta) < 1.5){
          return 1;
     }else if(std::fabs(roieta) < 2.0){
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
     float pSA_roieta = -99999;
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
                    pSA_roieta = m_pSA_roieta->at(method);
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
          
          /*
          if(pSA_superpointR_BI != 0 || pSA_superpointR_EI != 0 || pSA_superpointR_CSC != 0 || pSA_superpointR_BEE != 0 || pSA_superpointR_EE != 0){
               numSP++;
               patternSP += 1;
               SPinner = 1;
          }
          if(pSA_superpointR_BM != 0 || pSA_superpointR_EM != 0 || pSA_superpointR_BME != 0){
               numSP++;
               patternSP += 2;
               SPmiddle = 1;
          }
          if(pSA_superpointR_BO != 0 || pSA_superpointR_EO != 0){
               numSP++;
               patternSP += 3;
               SPouter = 1;
          }
          */
          if(m_poff_eta < 1.05){
               if(pSA_superpointR_BI != 0 || pSA_superpointR_BEE != 0){
                    numSP++;
                    SPinner = 1;
               }
               if(pSA_superpointR_BM != 0 || pSA_superpointR_EE != 0 || pSA_superpointR_BME != 0){
                    numSP++;
                    SPmiddle = 1;
               }
               if(pSA_superpointR_BO != 0){
                    numSP++;
                    SPouter = 1;
               }
          }else{
               if(pSA_superpointR_EI != 0 || pSA_superpointR_CSC != 0 || pSA_superpointR_EE != 0){
                    numSP++;
                    SPinner = 1;
               }
               if(pSA_superpointR_EM != 0){
                    numSP++;
                    SPmiddle = 1;
               }
               if(pSA_superpointR_EO != 0){
                    numSP++;
                    SPouter = 1;
               }
          }
          
          Int_t numBarrelSP = 0;
          if(pSA_superpointR_BI != 0){
               numBarrelSP++;
          }
          if(pSA_superpointR_BM != 0){
               numBarrelSP++;
          }
          if(pSA_superpointR_BO != 0){
               numBarrelSP++;
          }
          Int_t numEndcapSP = 0;
          if(pSA_superpointR_EI != 0){
               numEndcapSP++;
          }
          if(pSA_superpointR_EM != 0){
               numEndcapSP++;
          }
          if(pSA_superpointR_EO != 0){
               numEndcapSP++;
          }
          if(SPinner == 1 && SPmiddle == 1 && SPouter == 0)patternSP = 3;
          if(SPinner == 1 && SPmiddle == 0 && SPouter == 1)patternSP = 4;
          if(SPinner == 0 && SPmiddle == 1 && SPouter == 1)patternSP = 5;
          Int_t decision_noBIM = 0;
          for(Int_t index = 0;index < 10;index++){
               if(m_probe_segment_chamberIndex[index] == 1 && sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)) > 5800 && fabs(m_probe_segment_x[index]) > 4000.)decision_noBIM++;
               //if(m_probe_segment_chamberIndex[index] == 1 && (m_probe_segment_sector[index] == 11 || m_probe_segment_sector[index] == 15) && ((pSA_roiphi > -0.8 && pSA_roiphi < -0.6) || (pSA_roiphi > -2.6 && pSA_roiphi < -2.4)) )decision_noBIM++;
          }
          Int_t numAllSP = 0;
          if(pSA_superpointR_BI != 0)numAllSP++;
          if(pSA_superpointR_BM != 0)numAllSP++;
          if(pSA_superpointR_BO != 0)numAllSP++;
          if(pSA_superpointR_EI != 0)numAllSP++;
          if(pSA_superpointR_EM != 0)numAllSP++;
          if(pSA_superpointR_EO != 0)numAllSP++;
          if(pSA_superpointR_BME != 0)numAllSP++;
          if(pSA_superpointR_EE != 0)numAllSP++;
          if(pSA_superpointR_CSC != 0)numAllSP++;
          if(pSA_superpointR_BEE != 0)numAllSP++;

          if(numSP != numAllSP)cout<<"the number of SP ="<<numSP<<" the number of all SP ="<<numAllSP<<endl;

          //offline
          if(CutTagProbe(pEFTAG_pass)){
               if(static_cast<Int_t>(pSA_sAddress) == 1)m_h_pOffPhiLS->Fill(m_poff_phi);
               m_h_pOffPt->Fill(m_poff_pt*0.001);
               m_h_eOffPt->Fill(std::fabs(m_poff_pt*0.001));
               if(PlateauCut(std::fabs(m_poff_pt*0.001))){
                    m_h_eOffEta->Fill(m_poff_eta);
                    m_h_eOffPhi->Fill(m_poff_phi);
                    m_h_eOffAipc->Fill(m_aipc);
               }
               switch(EtaDistribution()){
                    case 0:
                    m_h_eOffPtBarrel->Fill(std::fabs(m_poff_pt*0.001));
                    break;
                    case 1:
                    m_h_eOffPtTransition->Fill(std::fabs(m_poff_pt*0.001));
                    break;
                    case 2:
                    m_h_eOffPtEnd->Fill(std::fabs(m_poff_pt*0.001));
                    break;
                    case 3:
                    m_h_eOffPtForward->Fill(std::fabs(m_poff_pt*0.001));
                    break;
                    default:
                    break;
               }
               if(PlateauCut(std::fabs(m_poff_pt*0.001)))m_h_eOffEtaPhi->Fill(m_poff_eta,m_poff_phi);

               switch(static_cast<Int_t>(pSA_sAddress)){
                    case 0:
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1)m_h_eOffPtLargePlus->Fill(std::fabs(m_poff_pt*0.001));
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1)m_h_eOffPtLargeMinus->Fill(std::fabs(m_poff_pt*0.001));
                    break;
                    case 1:
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1)m_h_eOffPtLSPlus->Fill(std::fabs(m_poff_pt*0.001));
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1)m_h_eOffPtLSMinus->Fill(std::fabs(m_poff_pt*0.001));
                    break;
                    case 2:
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1)m_h_eOffPtSmallPlus->Fill(std::fabs(m_poff_pt*0.001));
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1)m_h_eOffPtSmallMinus->Fill(std::fabs(m_poff_pt*0.001));
                    break;
                    case 3:
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==1)m_h_eOffPtSSPlus->Fill(std::fabs(m_poff_pt*0.001));
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta)==-1)m_h_eOffPtSSMinus->Fill(std::fabs(m_poff_pt*0.001));
                    break;
                    default:
                    break;
               }
               Int_t nosector9 = 0;
               for(Int_t index = 0;index < 10;index++){
                    Int_t nsector = m_probe_segment_sector[index] - 1;
                    if(std::fabs(nsector) < 16 && m_probe_segment_chamberIndex[index] >= 0 && m_probe_segment_chamberIndex[index] <= 5)m_h_SectorPhi[nsector]->Fill(TMath::ATan2(m_probe_segment_y[index],m_probe_segment_x[index]));
                    if(std::fabs(nsector) < 16 && m_probe_segment_chamberIndex[index] >= 0 && m_probe_segment_chamberIndex[index] <= 5)m_h_SectorRoIPhi[nsector]->Fill(pSA_roiphi);
                    if(m_probe_segment_sector[index] == 9)nosector9++;
               }

               Int_t numbersector = -1;


          //L1
               if(CutL1(pL1_pass)){
                    Double_t textL1_dR = TMath::Sqrt(pow(m_tL1_eta - m_toff_exteta,2) + pow(m_tL1_phi - m_toff_extphi,2));
                    pextL1_dR = TMath::Sqrt(pow(pL1_eta - m_poff_exteta,2) + pow(pL1_phi - m_poff_extphi,2));
                    m_h_pL1Pt->Fill(std::fabs(pL1_pt*0.001));
                    m_h_pL1dR->Fill(pL1_dR);
                    m_h_tExtL1dR->Fill(textL1_dR);
                    m_h_pExtL1dR->Fill(pextL1_dR);
                    m_h_eL1Pt->Fill(std::fabs(m_poff_pt*0.001));
                    m_h_L1pSAsAddress->Fill(pSA_sAddress);
                    switch(EtaDistribution()){
                         case 0:
                         m_h_eL1PtBarrel->Fill(std::fabs(m_poff_pt*0.001));
                         if(numAllSP == 1)m_h_eL1PtBarrel1SP->Fill(std::fabs(m_poff_pt*0.001));
                         if(numAllSP == 2)m_h_eL1PtBarrel2SP->Fill(std::fabs(m_poff_pt*0.001));
                         if(numAllSP == 3)m_h_eL1PtBarrel3SP->Fill(std::fabs(m_poff_pt*0.001));
                         if(patternSP == 3)m_h_eL1PtBarrelIM->Fill(std::fabs(m_poff_pt*0.001));
                         if(decision_noBIM == 0)m_h_eL1PtBarrelWithoutBIM->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_eL1PtBarrelIncBIM->Fill(std::fabs(m_poff_pt*0.001));
                         if(pSA_sAddress == 0 && nosector9 == 0)m_h_eL1PtLargeNormal->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                         case 1:
                         m_h_eL1PtTransition->Fill(std::fabs(m_poff_pt*0.001));
                         if(numAllSP == 1)m_h_eL1PtTransition1SP->Fill(std::fabs(m_poff_pt*0.001));
                         if(numAllSP == 2)m_h_eL1PtTransition2SP->Fill(std::fabs(m_poff_pt*0.001));
                         if(numAllSP == 3)m_h_eL1PtTransition3SP->Fill(std::fabs(m_poff_pt*0.001));
                         if(patternSP == 3)m_h_eL1PtTransitionIM->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                         case 2:
                         m_h_eL1PtEnd->Fill(std::fabs(m_poff_pt*0.001));
                         if(numAllSP == 1)m_h_eL1PtEndcap1SP->Fill(std::fabs(m_poff_pt*0.001));
                         if(numAllSP == 2)m_h_eL1PtEndcap2SP->Fill(std::fabs(m_poff_pt*0.001));
                         if(numAllSP == 3)m_h_eL1PtEndcap3SP->Fill(std::fabs(m_poff_pt*0.001));
                         if(patternSP == 3)m_h_eL1PtEndcapIM->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                         case 3:
                         m_h_eL1PtForward->Fill(std::fabs(m_poff_pt*0.001));
                         if(numAllSP == 1)m_h_eL1PtForward1SP->Fill(std::fabs(m_poff_pt*0.001));
                         if(numAllSP == 2)m_h_eL1PtForward2SP->Fill(std::fabs(m_poff_pt*0.001));
                         if(numAllSP == 3)m_h_eL1PtForward3SP->Fill(std::fabs(m_poff_pt*0.001));
                         if(patternSP == 3)m_h_eL1PtForwardIM->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                         default:
                         break;
                    }
                    switch(EtaDistribution(pSA_roieta)){
                         case 0:
                         m_h_eL1PtBarrelRoI->Fill(std::fabs(m_poff_pt*0.001));
                         if(numBarrelSP == 1)m_h_eL1PtBarrel1SPRoI->Fill(std::fabs(m_poff_pt*0.001));
                         if(numBarrelSP == 2)m_h_eL1PtBarrel2SPRoI->Fill(std::fabs(m_poff_pt*0.001));
                         if(numBarrelSP == 3)m_h_eL1PtBarrel3SPRoI->Fill(std::fabs(m_poff_pt*0.001));
                         if(pSA_superpointR_BI != 0 && pSA_superpointR_BM != 0 && pSA_superpointR_BO == 0)m_h_eL1PtBarrelIMRoI->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                         default:
                         break;
                    }

                    if(PlateauCut(std::fabs(m_poff_pt*0.001))){
                         m_h_eL1EtaPhi->Fill(m_poff_eta,m_poff_phi);
                         m_h_eL1Eta->Fill(m_poff_eta);
                         m_h_eL1Phi->Fill(m_poff_phi);
                         m_h_eL1Aipc->Fill(m_aipc);
                    }

                    areanumber = DicisionArea(pSA_roiphi);
                    switch(static_cast<Int_t>(pSA_sAddress)){//switch Large ,LS , Small ,SS
                         case 0:
                         if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0)m_h_eL1PtLargePlus->Fill(std::fabs(m_poff_pt*0.001));//Qeta = +1
                         if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) < 0)m_h_eL1PtLargeMinus->Fill(std::fabs(m_poff_pt*0.001));//Qeta = -1
                         break;
                         case 1:
                         if(areanumber > 0 && areanumber < 5)m_h_eL1PtLSPlus->Fill(std::fabs(m_poff_pt*0.001));//Qeta = +1
                         if(areanumber > 4 && areanumber < 9)m_h_eL1PtLSMinus->Fill(std::fabs(m_poff_pt*0.001));//Qeta = -1
                         if(decision_noBIM == 0)m_h_eL1PtLSWithoutBIM->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_eL1PtLSIncBIM->Fill(std::fabs(m_poff_pt*0.001));
                         switch(areanumber){
                              case 1://plus11out
                              m_h_eL1PtLSPlusS11outer->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_eL1PtLSPlusS11->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                              case 2://plus11in
                              m_h_eL1PtLSPlusS11inner->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_eL1PtLSPlusS11->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                              case 3://plus15out
                              m_h_eL1PtLSPlusS15outer->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_eL1PtLSPlusS15->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                              case 4://plus15inner
                              m_h_eL1PtLSPlusS15inner->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_eL1PtLSPlusS15->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                              case 5://minusS11out
                              m_h_eL1PtLSMinusS11outer->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_eL1PtLSMinusS11->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                              case 6://minus11in
                              m_h_eL1PtLSMinusS11inner->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_eL1PtLSMinusS11->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                              case 7://minus15out
                              m_h_eL1PtLSMinusS15outer->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_eL1PtLSMinusS15->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                              case 8://minus15in
                              m_h_eL1PtLSMinusS15inner->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_eL1PtLSMinusS15->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                              default:
                              break;
                         }
                         break;
                         case 2:
                         if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0)m_h_eL1PtSmallPlus->Fill(std::fabs(m_poff_pt*0.001));//Qeta = +1
                         if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) < 0)m_h_eL1PtSmallMinus->Fill(std::fabs(m_poff_pt*0.001));//Qeta = -1
                         break;
                         case 3:
                         if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0)m_h_eL1PtSSPlus->Fill(std::fabs(m_poff_pt*0.001));//Qeta = +1
                         if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) < 0)m_h_eL1PtSSMinus->Fill(std::fabs(m_poff_pt*0.001));//Qeta = -1
                         break;
                         default:
                         break;
                    }

          //SA
                    if(CutSA(pSA_pass)){
                         if(std::fabs(m_poff_pt*0.001) < 2.75 && numSP == 3){
                              std::ofstream ofs;
                              ofs.open("/home/kayamash/LowPtPassed.txt",std::ios::app);
                              ofs<<m_rNumber<<"   "<<m_eNumber<<"   "<<pSA_pt<<"   "<<pSA_eta<<"   "<<pSA_phi<<"   "<<m_poff_pt*0.001<<"   "<<m_poff_eta<<"   "<<m_poff_phi<<std::endl;
                              ofs.close();
                         }
                         if(std::fabs(m_poff_pt*0.001) < 2.75 && numAllSP == 3){
                              std::ofstream ofs;
                              ofs.open("/home/kayamash/LowPtPassedAllSP.txt",std::ios::app);
                              ofs<<m_rNumber<<"   "<<m_eNumber<<"   "<<pSA_pt<<"   "<<pSA_eta<<"   "<<pSA_phi<<"   "<<m_poff_pt*0.001<<"   "<<m_poff_eta<<"   "<<m_poff_phi<<std::endl;
                              ofs.close();
                         }
                         m_h_CountSA->Fill(m_poff_eta);
                         m_h_NumSP->Fill(numAllSP);
                         Double_t textSA_dR = TMath::Sqrt(pow(m_tSA_eta - m_toff_exteta,2) + pow(m_tSA_phi - m_toff_extphi,2));
                         pextSA_dR = TMath::Sqrt(pow(pSA_eta - m_poff_exteta,2) + pow(pSA_phi - m_poff_extphi,2));
                         Double_t resSA_pt = std::fabs(m_poff_pt*0.001)/std::fabs(pSA_pt) - 1.0;
                         Double_t buf_pSA_dR = TMath::Sqrt(pow(pSA_eta - m_poff_eta,2) + pow(pSA_phi - m_poff_phi,2));
                         Double_t buf_eta = 0;
                         m_h_pSAsAddress->Fill(pSA_sAddress);

                         m_h_L2MuonSAPtAlphavsBeta[0]->Fill(std::fabs(pSA_ptalpha),std::fabs(pSA_ptbeta));
                         m_h_L2MuonSAPtAlphavsTGC[0]->Fill(std::fabs(pSA_ptalpha),std::fabs(pSA_ptTGC));
                         m_h_L2MuonSAPtBetavsTGC[0]->Fill(std::fabs(pSA_ptbeta),std::fabs(pSA_ptTGC));
                         if(std::fabs(m_poff_pt*0.001) < pt_threshold[EtaDistribution()]){
                              m_h_L2MuonSAPtAlphavsBeta[1]->Fill(std::fabs(pSA_ptalpha),std::fabs(pSA_ptbeta));
                              m_h_L2MuonSAPtAlphavsTGC[1]->Fill(std::fabs(pSA_ptalpha),std::fabs(pSA_ptTGC));
                              m_h_L2MuonSAPtBetavsTGC[1]->Fill(std::fabs(pSA_ptbeta),std::fabs(pSA_ptTGC));
                         }else{
                              m_h_L2MuonSAPtAlphavsBeta[2]->Fill(std::fabs(pSA_ptalpha),std::fabs(pSA_ptbeta));
                              m_h_L2MuonSAPtAlphavsTGC[2]->Fill(std::fabs(pSA_ptalpha),std::fabs(pSA_ptTGC));
                              m_h_L2MuonSAPtBetavsTGC[2]->Fill(std::fabs(pSA_ptbeta),std::fabs(pSA_ptTGC));
                         }

                         m_h_pSAPt->Fill(std::fabs(pSA_pt));
                         m_h_pSAdR->Fill(buf_pSA_dR);
                         m_h_tExtSAdR->Fill(textSA_dR);
                         m_h_pExtSAdR->Fill(pextSA_dR);
                         m_h_eSAPt->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_pSAResPt->Fill(resSA_pt);
                         m_h_pSAPhivspSAPhims->Fill(pSA_phi,pSA_phims);
                         m_h_pSAPhivspSAPhibe->Fill(pSA_phi,pSA_phibe);
                         switch(EtaDistribution()){
                              case 0:
                              m_h_eSAPtBarrel->Fill(std::fabs(m_poff_pt*0.001));
                              if(numAllSP == 1)m_h_eSAPtBarrel1SP->Fill(std::fabs(m_poff_pt*0.001));
                              if(numAllSP == 2)m_h_eSAPtBarrel2SP->Fill(std::fabs(m_poff_pt*0.001));
                              if(numAllSP == 3)m_h_eSAPtBarrel3SP->Fill(std::fabs(m_poff_pt*0.001));
                              if(patternSP == 3)m_h_eSAPtBarrelIM->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_pSAResPtBarrel->Fill(resSA_pt);
                              if(decision_noBIM == 0)m_h_eSAPtBarrelWithoutBIM->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_eSAPtBarrelIncBIM->Fill(std::fabs(m_poff_pt*0.001));
                              for(Int_t MDTsize = 0;MDTsize < (signed int)pSA_mdthitChamber->size();MDTsize++){
                                   m_h_MDTChamber->Fill(pSA_mdthitChamber->at(MDTsize));
                              }
                              if(pSA_sAddress == 0 && nosector9 == 0)m_h_eSAPtLargeNormal->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                              case 1:
                              m_h_eSAPtTransition->Fill(std::fabs(m_poff_pt*0.001));
                              if(numAllSP == 1)m_h_eSAPtTransition1SP->Fill(std::fabs(m_poff_pt*0.001));
                              if(numAllSP == 2)m_h_eSAPtTransition2SP->Fill(std::fabs(m_poff_pt*0.001));
                              if(numAllSP == 3)m_h_eSAPtTransition3SP->Fill(std::fabs(m_poff_pt*0.001));
                              if(patternSP == 3)m_h_eSAPtTransitionIM->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                              case 2:
                              m_h_eSAPtEnd->Fill(std::fabs(m_poff_pt*0.001));
                              if(numAllSP == 1)m_h_eSAPtEndcap1SP->Fill(std::fabs(m_poff_pt*0.001));
                              if(numAllSP == 2)m_h_eSAPtEndcap2SP->Fill(std::fabs(m_poff_pt*0.001));
                              if(numAllSP == 3)m_h_eSAPtEndcap3SP->Fill(std::fabs(m_poff_pt*0.001));
                              if(patternSP == 3)m_h_eSAPtEndcapIM->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                              case 3:
                              m_h_eSAPtForward->Fill(std::fabs(m_poff_pt*0.001));
                              if(numAllSP == 1)m_h_eSAPtForward1SP->Fill(std::fabs(m_poff_pt*0.001));
                              if(numAllSP == 2)m_h_eSAPtForward2SP->Fill(std::fabs(m_poff_pt*0.001));
                              if(numAllSP == 3)m_h_eSAPtForward3SP->Fill(std::fabs(m_poff_pt*0.001));
                              if(patternSP == 3)m_h_eSAPtForwardIM->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                         }

                         switch(EtaDistribution(pSA_roieta)){
                              case 0:
                              m_h_eSAPtBarrelRoI->Fill(std::fabs(m_poff_pt*0.001));
                              if(numBarrelSP == 1)m_h_eSAPtBarrel1SPRoI->Fill(std::fabs(m_poff_pt*0.001));
                              if(numBarrelSP == 2)m_h_eSAPtBarrel2SPRoI->Fill(std::fabs(m_poff_pt*0.001));
                              if(numBarrelSP == 3)m_h_eSAPtBarrel3SPRoI->Fill(std::fabs(m_poff_pt*0.001));
                              if(pSA_superpointR_BI != 0 && pSA_superpointR_BM != 0 && pSA_superpointR_BO == 0)m_h_eSAPtBarrelIMRoI->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                              default:
                              break;
                         }
                         if(DicisionBarrel(m_poff_eta)){
                         }else{
                              m_h_pSAResPtEndcap->Fill(resSA_pt);
                              m_h_L2MuonSAvsOfflinePt[0]->Fill(std::fabs(pSA_ptalpha),std::fabs(m_poff_pt*0.001));
                              m_h_L2MuonSAvsOfflinePt[1]->Fill(std::fabs(pSA_ptbeta),std::fabs(m_poff_pt*0.001));
                              m_h_L2MuonSAvsOfflinePt[2]->Fill(std::fabs(pSA_ptTGC),std::fabs(m_poff_pt*0.001));
                              if(pt_method >= 0 && std::fabs(m_poff_pt*0.001) < pt_threshold[EtaDistribution()] && EtaDistribution() > 0)m_h_PtMethod[pt_method]->Fill(std::fabs(m_poff_pt*0.001));
                              if(pt_method >= 0 && std::fabs(m_poff_pt*0.001) > pt_threshold[EtaDistribution() && EtaDistribution() > 0])m_h_PtMethodOver[pt_method]->Fill(std::fabs(m_poff_pt*0.001));
                              if(std::fabs(m_poff_pt*0.001) < pt_threshold[EtaDistribution()]){
                                   m_h_PtSA[0]->Fill(std::fabs(pSA_ptalpha));
                                   m_h_PtSA[1]->Fill(std::fabs(pSA_ptbeta));
                                   m_h_PtSA[2]->Fill(std::fabs(pSA_ptTGC));
                                   if(pt_method == 3)m_h_PtSA[3]->Fill(std::fabs(pSA_pt));
                                   m_h_ResPtSA[0]->Fill(resptalpha);
                                   m_h_ResPtSA[1]->Fill(resptbeta);
                                   m_h_ResPtSA[2]->Fill(respttgc);
                                   if(pt_method == 3)m_h_ResPtSA[3]->Fill(resSA_pt);
                              }else{
                                   m_h_PtSAOver[0]->Fill(std::fabs(pSA_ptalpha));
                                   m_h_PtSAOver[1]->Fill(std::fabs(pSA_ptbeta));
                                   m_h_PtSAOver[2]->Fill(std::fabs(pSA_ptTGC));
                                   if(pt_method == 3)m_h_PtSAOver[3]->Fill(std::fabs(pSA_pt));
                                   m_h_ResPtSAOver[0]->Fill(resptalpha);
                                   m_h_ResPtSAOver[1]->Fill(resptbeta);
                                   m_h_ResPtSAOver[2]->Fill(respttgc);
                                   if(pt_method == 3)m_h_ResPtSAOver[3]->Fill(resSA_pt);
                              }
                         }
                         for(Int_t size = 0;size < (signed int)pSA_mdtZ->size();size++){
                              buf_eta += -TMath::Log((sqrt(pow(pSA_mdtZ->at(size),2) + pow(pSA_mdtR->at(size),2)) - pSA_mdtZ->at(size))/(sqrt(pow(pSA_mdtZ->at(size),2) + pow(pSA_mdtR->at(size),2)) + pow(pSA_mdtZ->at(size),2)))/2.0;
                         }
                         m_h_NumHit->Fill(pSA_mdtZ->size());
                         for(Int_t index = 0;index < 10;index++){
                              m_h_SectorvsPhi->Fill(m_probe_segment_sector[index],m_poff_phi);
                              m_h_IndexvsEta->Fill(m_probe_segment_etaIndex[index],m_poff_eta);
                         }
                         if(static_cast<Int_t>(pSA_sAddress) == 2 && m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0)m_h_eSAPtSmallPlus->Fill(std::fabs(m_poff_pt*0.001));
          switch(static_cast<Int_t>(pSA_sAddress)){//switch Large ,LS , Small ,SS
               case 0:
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0){//Qeta = +1
                         m_h_pOffPtvsSAResPtLargePlus->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SAResPtLargePlus->Fill(resSA_pt);
                         m_h_eSAPtLargePlus->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_pOffPhivsSAResPtLargePlus->Fill(m_poff_phi,resSA_pt);
                         m_h_pOffEtavsSAResPtLargePlus->Fill(m_poff_eta,resSA_pt);

                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_EtaIndexvsSAResPtLargePlus->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         }
                    }

                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) < 0){//Qeta = -1
                         m_h_pOffPtvsSAResPtLargeMinus->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SAResPtLargeMinus->Fill(resSA_pt);
                         m_h_eSAPtLargeMinus->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_pOffPhivsSAResPtLargeMinus->Fill(m_poff_phi,resSA_pt);
                         m_h_pOffEtavsSAResPtLargeMinus->Fill(m_poff_eta,resSA_pt);
                    }
                    break;
                    case 1:
                    m_h_SAPhimsLS->Fill(pSA_phims);
                    m_h_SARoIPhiLS->Fill(pSA_roiphi);
                    if(decision_noBIM == 0)m_h_eSAPtLSWithoutBIM->Fill(std::fabs(m_poff_pt*0.001));
                    m_h_eSAPtLSIncBIM->Fill(std::fabs(m_poff_pt*0.001));

                    for(Int_t index = 0;index < 10;index++){
                         if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_EtaIndexvsSAResPtLSPlusS11outer->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numAllSP == 2)m_h_EtaIndexvsSAResPtLSPlusS11outer2Station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numAllSP == 3)m_h_EtaIndexvsSAResPtLSPlusS11outer3Station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0){
                              m_h_SegmentZRLSPlusS11outer->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              if(m_probe_segment_chamberIndex[index] == 1)m_h_SegmentZRLSPlusS11outerBIL->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         }
                    }

                    if(areanumber > 0 && areanumber < 5){
                         m_h_eSAPtLSPlus->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_SAResPtLSPlus->Fill(resSA_pt);
                         m_h_pOffPtvsSAResPtLSPlus->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_pOffPhivsSAResPtLSPlus->Fill(m_poff_phi,resSA_pt);
                         m_h_pOffEtavsSAResPtLSPlus->Fill(m_poff_eta,resSA_pt);
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_EtaIndexvsSAResPtLSPlus->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_SegmentZRLSPlus->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         }
                         m_h_MDTSPZRLSPlus->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_MDTSPZRLSPlus->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_MDTSPZRLSPlus->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                    }
                    if(areanumber > 4 && areanumber < 9){
                         m_h_eSAPtLSMinus->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_SAResPtLSMinus->Fill(resSA_pt);
                         m_h_pOffPtvsSAResPtLSMinus->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_pOffPhivsSAResPtLSMinus->Fill(m_poff_phi,resSA_pt);
                         m_h_pOffEtavsSAResPtLSMinus->Fill(m_poff_eta,resSA_pt);
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_EtaIndexvsSAResPtLSMinus->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_SegmentZRLSMinus->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         }
                         m_h_MDTSPZRLSMinus->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_MDTSPZRLSMinus->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_MDTSPZRLSMinus->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                    }

                    switch(areanumber){
                         case 1://plus 11out
                         m_h_MDTSPXBI->Fill(pSA_superpointR_BI*cos(pSA_roiphi));
                         m_h_MDTSPYBI->Fill(pSA_superpointR_BI*sin(pSA_roiphi));
                         m_h_pOffPtvsSAResPtLSPlusS11outer->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numAllSP == 2)m_h_pOffPtvsSAResPtLSPlusS11outer2Station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numAllSP == 3)m_h_pOffPtvsSAResPtLSPlusS11outer3Station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SAResPtLSPlusS11outer->Fill(resSA_pt);
                         m_h_pOffEtavsSAResPtLSPlusS11outer->Fill(m_poff_eta,resSA_pt);
                         m_h_eSAPtLSPlusS11outer->Fill(std::fabs(m_poff_pt*0.001));

                         m_h_MDTSPZRLSPlusS11outer->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_MDTSPZRLSPlusS11outer->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_MDTSPZRLSPlusS11outer->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                         break;
                         case 2://plus 11in
                         m_h_MDTSPXBI->Fill(pSA_superpointR_BI*cos(pSA_roiphi));
                         m_h_MDTSPYBI->Fill(pSA_superpointR_BI*sin(pSA_roiphi));
                         m_h_pOffPtvsSAResPtLSPlusS11inner->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numAllSP == 2)m_h_pOffPtvsSAResPtLSPlusS11inner2Station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numAllSP == 3)m_h_pOffPtvsSAResPtLSPlusS11inner3Station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SAResPtLSPlusS11inner->Fill(resSA_pt);
                         m_h_pOffEtavsSAResPtLSPlusS11inner->Fill(m_poff_eta,resSA_pt);
                         m_h_eSAPtLSPlusS11inner->Fill(std::fabs(m_poff_pt*0.001));
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_EtaIndexvsSAResPtLSPlusS11inner->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numAllSP == 2)m_h_EtaIndexvsSAResPtLSPlusS11inner2Station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numAllSP == 3)m_h_EtaIndexvsSAResPtLSPlusS11inner3Station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0){
                                   m_h_SegmentZRLSPlusS11inner->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                                   if(m_probe_segment_chamberIndex[index] == 1)m_h_SegmentZRLSPlusS11innerBIL->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              }
                         }
                         m_h_MDTSPZRLSPlusS11inner->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_MDTSPZRLSPlusS11inner->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_MDTSPZRLSPlusS11inner->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                         break;
                         case 3://plus 15out
                         m_h_pOffPtvsSAResPtLSPlusS15outer->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numAllSP == 2)m_h_pOffPtvsSAResPtLSPlusS15outer2Station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numAllSP == 3)m_h_pOffPtvsSAResPtLSPlusS15outer3Station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SAResPtLSPlusS15outer->Fill(resSA_pt);
                         m_h_pOffEtavsSAResPtLSPlusS15outer->Fill(m_poff_eta,resSA_pt);
                         m_h_eSAPtLSPlusS15outer->Fill(std::fabs(m_poff_pt*0.001));

                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8 && m_probe_segment_etaIndex[index] <= 8)m_h_EtaIndexvsSAResPtLSPlusS15outer->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8 && m_probe_segment_etaIndex[index] <= 8 && numAllSP == 2)m_h_EtaIndexvsSAResPtLSPlusS15outer2Station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8 && m_probe_segment_etaIndex[index] <= 8 && numAllSP == 3)m_h_EtaIndexvsSAResPtLSPlusS15outer3Station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0){
                                   m_h_SegmentZRLSPlusS15outer->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                                   if(m_probe_segment_chamberIndex[index] == 1)m_h_SegmentZRLSPlusS15outerBIL->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              }
                         }
                         m_h_MDTSPZRLSPlusS15outer->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_MDTSPZRLSPlusS15outer->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_MDTSPZRLSPlusS15outer->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                         break;
                         case 4://plus 15in
                         m_h_pOffPtvsSAResPtLSPlusS15inner->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numAllSP == 2)m_h_pOffPtvsSAResPtLSPlusS15inner2Station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numAllSP == 3)m_h_pOffPtvsSAResPtLSPlusS15inner3Station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SAResPtLSPlusS15inner->Fill(resSA_pt);
                         m_h_pOffEtavsSAResPtLSPlusS15inner->Fill(m_poff_eta,resSA_pt);
                         m_h_eSAPtLSPlusS15inner->Fill(std::fabs(m_poff_pt*0.001));
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_EtaIndexvsSAResPtLSPlusS15inner->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numAllSP == 2)m_h_EtaIndexvsSAResPtLSPlusS15inner2Station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numAllSP == 3)m_h_EtaIndexvsSAResPtLSPlusS15inner3Station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0){
                                   m_h_SegmentZRLSPlusS15inner->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                                   if(m_probe_segment_chamberIndex[index] == 1)m_h_SegmentZRLSPlusS15innerBIL->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              }
                         }
                         m_h_MDTSPZRLSPlusS15inner->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_MDTSPZRLSPlusS15inner->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_MDTSPZRLSPlusS15inner->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                         break;
                         case 5://minus 11out
                         m_h_pOffPtvsSAResPtLSMinusS11outer->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numAllSP == 2)m_h_pOffPtvsSAResPtLSMinusS11outer2Station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numAllSP == 3)m_h_pOffPtvsSAResPtLSMinusS11outer3Station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SAResPtLSMinusS11outer->Fill(resSA_pt);
                         m_h_pOffEtavsSAResPtLSMinusS11outer->Fill(m_poff_eta,resSA_pt);
                         m_h_eSAPtLSMinusS11outer->Fill(std::fabs(m_poff_pt*0.001));
                         for(Int_t index = 0;index < 10;index++){
                          if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_EtaIndexvsSAResPtLSMinusS11outer->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                          if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numAllSP == 2)m_h_EtaIndexvsSAResPtLSMinusS11outer2Station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                          if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numAllSP == 3)m_h_EtaIndexvsSAResPtLSMinusS11outer3Station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                          if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_SegmentZRLSMinusS11outer->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                     }
                     m_h_MDTSPZRLSMinusS11outer->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                     m_h_MDTSPZRLSMinusS11outer->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                     m_h_MDTSPZRLSMinusS11outer->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                     break;
                         case 6://minus 11in
                         m_h_pOffPtvsSAResPtLSMinusS11inner->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numAllSP == 2)m_h_pOffPtvsSAResPtLSMinusS11inner2Station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numAllSP == 3)m_h_pOffPtvsSAResPtLSMinusS11inner3Station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SAResPtLSMinusS11inner->Fill(resSA_pt);
                         m_h_pOffEtavsSAResPtLSMinusS11inner->Fill(m_poff_eta,resSA_pt);
                         m_h_eSAPtLSMinusS11inner->Fill(std::fabs(m_poff_pt*0.001));
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_EtaIndexvsSAResPtLSMinusS11inner->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numAllSP == 2)m_h_EtaIndexvsSAResPtLSMinusS11inner2Station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numAllSP == 3)m_h_EtaIndexvsSAResPtLSMinusS11inner3Station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_SegmentZRLSMinusS11inner->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         }
                         m_h_MDTSPZRLSMinusS11inner->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_MDTSPZRLSMinusS11inner->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_MDTSPZRLSMinusS11inner->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                         break;
                         case 7://minus 15out
                         m_h_pOffPtvsSAResPtLSMinusS15outer->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numAllSP == 2)m_h_pOffPtvsSAResPtLSMinusS15outer2Station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numAllSP == 3)m_h_pOffPtvsSAResPtLSMinusS15outer3Station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SAResPtLSMinusS15outer->Fill(resSA_pt);
                         m_h_pOffEtavsSAResPtLSMinusS15outer->Fill(m_poff_eta,resSA_pt);
                         m_h_eSAPtLSMinusS15outer->Fill(std::fabs(m_poff_pt*0.001));
                         for(Int_t index = 0;index < 10;index++){
                          if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_EtaIndexvsSAResPtLSMinusS15outer->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                          if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numAllSP == 2)m_h_EtaIndexvsSAResPtLSMinusS15outer2Station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                          if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numAllSP == 3)m_h_EtaIndexvsSAResPtLSMinusS15outer3Station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                          if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_SegmentZRLSMinusS15outer->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                     }
                     m_h_MDTSPZRLSMinusS15outer->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                     m_h_MDTSPZRLSMinusS15outer->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                     m_h_MDTSPZRLSMinusS15outer->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                     break;
                         case 8://minus 15in
                         m_h_pOffPtvsSAResPtLSMinusS15inner->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numAllSP == 2)m_h_pOffPtvsSAResPtLSMinusS15inner2Station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         if(numAllSP == 3)m_h_pOffPtvsSAResPtLSMinusS15inner3Station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SAResPtLSMinusS15inner->Fill(resSA_pt);
                         m_h_pOffEtavsSAResPtLSMinusS15inner->Fill(m_poff_eta,resSA_pt);
                         m_h_eSAPtLSMinusS15inner->Fill(std::fabs(m_poff_pt*0.001));
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_EtaIndexvsSAResPtLSMinusS15inner->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numAllSP == 2)m_h_EtaIndexvsSAResPtLSMinusS15inner2Station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0 && numAllSP == 3)m_h_EtaIndexvsSAResPtLSMinusS15inner3Station->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_SegmentZRLSMinusS15inner->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         }
                         m_h_MDTSPZRLSMinusS15inner->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_MDTSPZRLSMinusS15inner->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_MDTSPZRLSMinusS15inner->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                         break;
                         default :
                         break;
                    }

                    case 2:
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0){//Qeta = +1
                         m_h_pOffPtvsSAResPtSmallPlus->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SAResPtSmallPlus->Fill(resSA_pt);
                         //m_h_eSAPtSmallPlus->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_pOffPhivsSAResPtSmallPlus->Fill(m_poff_phi,resSA_pt);
                         m_h_pOffEtavsSAResPtSmallPlus->Fill(m_poff_eta,resSA_pt);
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_EtaIndexvsSAResPtSmallPlus->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         }
                    }
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) < 0){//Qeta = -1
                         m_h_pOffPtvsSAResPtSmallMinus->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SAResPtSmallMinus->Fill(resSA_pt);
                         m_h_eSAPtSmallMinus->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_pOffPhivsSAResPtSmallMinus->Fill(m_poff_phi,resSA_pt);
                         m_h_pOffEtavsSAResPtSmallMinus->Fill(m_poff_eta,resSA_pt);
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_EtaIndexvsSAResPtSmallMinus->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         }
                    }
                    break;
                    case 3:
                    m_h_SARoIPhiSS->Fill(pSA_roiphi);
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) > 0){//Qeta = +1
                         m_h_pOffPtvsSAResPtSSPlus->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SAResPtSSPlus->Fill(resSA_pt);
                         m_h_eSAPtSSPlus->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_pOffPhivsSAResPtSSPlus->Fill(m_poff_phi,resSA_pt);
                         m_h_pOffEtavsSAResPtSSPlus->Fill(m_poff_eta,resSA_pt);
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_EtaIndexvsSAResPtSSPlus->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         }
                    }
                    if(m_poff_charge*m_poff_eta/std::fabs(m_poff_eta) < 0){//Qeta = -1
                         m_h_pOffPtvsSAResPtSSMinus->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         m_h_SAResPtSSMinus->Fill(resSA_pt);
                         m_h_eSAPtSSMinus->Fill(std::fabs(m_poff_pt*0.001));
                         m_h_pOffPhivsSAResPtSSMinus->Fill(m_poff_phi,resSA_pt);
                         m_h_pOffEtavsSAResPtSSMinus->Fill(m_poff_eta,resSA_pt);
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_EtaIndexvsSAResPtSSMinus->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         }
                    }
                    break;
                    default:
                    break;
               }


               for(Int_t mdthit = 0; mdthit < (signed int)pSA_mdtPhi->size(); mdthit++){
                    m_h_MDTPhi->Fill(pSA_mdtPhi->at(mdthit));
               }

               if(numAllSP == 2){
                    if(pSA_sAddress == 0 && std::fabs(m_poff_phi) < TMath::Pi()/8.0)m_h_pOffPtvsSAResPtLarge2Station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                    if(pSA_sAddress == 2 && std::fabs(m_poff_phi) < TMath::Pi()/8.0)m_h_pOffPtvsSAResPtSmall2Station->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
               }

               if(pSA_sAddress == 1){
                    for(Int_t mdthit = 0; mdthit < (signed int)pSA_mdtZ->size(); mdthit++){
                         m_h_MDTPhiLS->Fill(pSA_mdtPhi->at(mdthit));
                         m_h_MDTR->Fill(pSA_mdtR->at(mdthit));
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
               if(pSA_sAddress == 0 || pSA_sAddress == 2)m_h_OverPhi->Fill(overphi);
               Int_t numeta = 0;
               Int_t numphi = 0;
               if(pSA_eta < 0)numeta = 8;
               if(pSA_sAddress == 0)overphi -= TMath::Pi()/8.0;
               numphi = overphi/TMath::Pi()*32.0 - fmod(overphi,TMath::Pi()/32.0);
               numeta += std::fabs(pSA_eta)/0.125 - fmod(std::fabs(pSA_eta),0.125);
               if(pSA_sAddress == 0 && std::fabs(pSA_eta) < 1.0)m_h_DivideEtaOverPhiResPtLarge[numeta][numphi]->Fill(resSA_pt);
               if(pSA_sAddress == 2 && std::fabs(pSA_eta) < 1.0)m_h_DivideEtaOverPhiResPtSmall[numeta][numphi]->Fill(resSA_pt);

               if(PlateauCut(std::fabs(m_poff_pt*0.001))){
                    m_h_eSAEtaPhi->Fill(m_poff_eta,m_poff_phi);
                    m_h_eSAEta->Fill(m_poff_eta);
                    m_h_eSAPhi->Fill(m_poff_phi);
                    m_h_eSAAipc->Fill(m_aipc);
               }
               m_h_pOffvsSAPt->Fill(std::fabs(m_poff_pt*0.001),std::fabs(pSA_pt));

               if(pSA_sAddress == 1){
                    for(Int_t index = 0;index < 10;index++){
                         if(m_probe_segment_chamberIndex[index] == 1 && sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)) > 4000 && fabs(m_probe_segment_x[index]) < 4500.){
                              m_h_nPrecisionHitsNormal->Fill(m_probe_segment_nPrecisionHits[index]);
                              m_h_SectorNormal->Fill(m_probe_segment_sector[index]);
                              m_h_EtaIndexNormal->Fill(m_probe_segment_etaIndex[index]);
                              if(m_probe_segment_sector[index] == 11)m_h_SegmentNumberBIR->Fill((Int_t)m_probe_segment_etaIndex[index] - 1,0);
                              if(m_probe_segment_sector[index] == 15)m_h_SegmentNumberBIR->Fill((Int_t)m_probe_segment_etaIndex[index] - 1,1);
                         }
                         if(m_probe_segment_chamberIndex[index] == 1 && sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)) > 4000 && fabs(m_probe_segment_x[index]) > 4500.){
                              m_h_nPrecisionHitsSpecial->Fill(m_probe_segment_nPrecisionHits[index]);
                              m_h_SectorSpecial->Fill(m_probe_segment_sector[index]);
                              m_h_EtaIndexSpecial->Fill(m_probe_segment_etaIndex[index]);
                              if(m_probe_segment_sector[index] == 11)m_h_SegmentNumberBIM->Fill((Int_t)m_probe_segment_etaIndex[index] - 1,0);
                              if(m_probe_segment_sector[index] == 15)m_h_SegmentNumberBIM->Fill((Int_t)m_probe_segment_etaIndex[index] - 1,1);
                         }
                    }
               }

               if(pSA_sAddress == 0 && std::fabs(m_poff_phi) < TMath::Pi()/8.0){
                    if(numAllSP == 2)m_h_pOffPtLarge2Station->Fill(std::fabs(m_poff_pt*0.001));
                    if(numAllSP == 2)m_h_pOffEtaLarge2Station->Fill(m_poff_eta);
                    if(numAllSP == 3)m_h_pOffEtaLarge3Station->Fill(m_poff_eta);
                    switch(patternSP){
                         case 3:
                         if(numAllSP == 2){
                              m_h_pOffPtLarge2StationIM->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_pOffEtaLarge2StationIM->Fill(m_poff_eta);
                              m_h_pOffPtvsSAResPtLarge2StationIM->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         case 4:
                         if(numAllSP == 2){
                              m_h_pOffPtLarge2StationIO->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_pOffEtaLarge2StationIO->Fill(m_poff_eta);
                              m_h_pOffPtvsSAResPtLarge2StationIO->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         case 5:
                         if(numAllSP == 2){
                              m_h_pOffPtLarge2StationMO->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_pOffEtaLarge2StationMO->Fill(m_poff_eta);
                              m_h_pOffPtvsSAResPtLarge2StationMO->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         default :
                         break;
                    }
               }
               if(pSA_sAddress == 2 && m_poff_phi > 0 && m_poff_phi < TMath::Pi()/4.0){
                    if(numAllSP == 2)m_h_pOffPtSmall2Station->Fill(std::fabs(m_poff_pt*0.001));
                    if(numAllSP == 2)m_h_pOffEtaSmall2Station->Fill(m_poff_eta);
                    if(numAllSP == 3)m_h_pOffEtaSmall3Station->Fill(m_poff_eta);
                    switch(patternSP){
                         case 3:
                         if(numAllSP == 2){
                              m_h_pOffPtSmall2StationIM->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_pOffEtaSmall2StationIM->Fill(m_poff_eta);
                              m_h_pOffPtvsSAResPtSmall2StationIM->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         case 4:
                         if(numAllSP == 2){
                              m_h_pOffPtSmall2StationIO->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_pOffEtaSmall2StationIO->Fill(m_poff_eta);
                              m_h_pOffPtvsSAResPtSmall2StationIO->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         case 5:
                         if(numAllSP == 2){
                              m_h_pOffPtSmall2StationMO->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_pOffEtaSmall2StationMO->Fill(m_poff_eta);
                              m_h_pOffPtvsSAResPtSmall2StationMO->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         default :
                         break;
                    }
               }
               if(areanumber % 2 == 0 && areanumber != 0){
                    if(numAllSP == 2)m_h_pOffPtBIR2Station->Fill(std::fabs(m_poff_pt*0.001));
                    if(numAllSP == 2)m_h_pOffEtaBIR2Station->Fill(m_poff_eta);
                    if(numAllSP == 3)m_h_pOffEtaBIR3Station->Fill(m_poff_eta);
                    switch(patternSP){
                         case 3:
                         if(numAllSP == 2){
                              m_h_pOffPtBIR2StationIM->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_pOffEtaBIR2StationIM->Fill(m_poff_eta);
                              m_h_pOffPtvsSAResPtBIR2StationIM->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         case 4:
                         if(numAllSP == 2){
                              m_h_pOffPtBIR2StationIO->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_pOffEtaBIR2StationIO->Fill(m_poff_eta);
                              m_h_pOffPtvsSAResPtBIR2StationIO->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         case 5:
                         if(numAllSP == 2){
                              m_h_pOffPtBIR2StationMO->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_pOffEtaBIR2StationMO->Fill(m_poff_eta);
                              m_h_pOffPtvsSAResPtBIR2StationMO->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         default :
                         break;
                    }
               }
               if(areanumber % 2 == 1){
                    if(numAllSP == 2)m_h_pOffPtBIM2Station->Fill(std::fabs(m_poff_pt*0.001));
                    if(numAllSP == 2)m_h_pOffEtaBIM2Station->Fill(m_poff_eta);
                    if(numAllSP == 3)m_h_pOffEtaBIM3Station->Fill(m_poff_eta);
                    switch(patternSP){
                         case 3:
                         if(numAllSP == 2){
                              m_h_pOffPtBIM2StationIM->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_pOffEtaBIM2StationIM->Fill(m_poff_eta);
                              m_h_pOffPtvsSAResPtBIM2StationIM->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         case 4:
                         if(numAllSP == 2){
                              m_h_pOffPtBIM2StationIO->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_pOffEtaBIM2StationIO->Fill(m_poff_eta);
                              m_h_pOffPtvsSAResPtBIM2StationIO->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                         }
                         break;
                         case 5:
                         if(numAllSP == 2){
                              m_h_pOffPtBIM2StationMO->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_pOffEtaBIM2StationMO->Fill(m_poff_eta);
                              m_h_pOffPtvsSAResPtBIM2StationMO->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
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
                    if(m_probe_segment_chamberIndex[index] == 1 && segdisBIR == 2)m_h_BIRSegment->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                    if(m_probe_segment_chamberIndex[index] == 1 && segdisBIM == 2)m_h_BIMSegment->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
               }

               for(Int_t chnum = 0; chnum < 10; chnum++){
                    if(m_probe_segment_chamberIndex[chnum] == 0)m_h_SegmentXYBIS->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 1)m_h_SegmentXYBIL->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 2)m_h_SegmentXYBMS->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 3)m_h_SegmentXYBML->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 4)m_h_SegmentXYBOS->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 5)m_h_SegmentXYBOL->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 6)m_h_SegmentXYBEE->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 7)m_h_SegmentXYEIS->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 8)m_h_SegmentXYEIL->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 9)m_h_SegmentXYEMS->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 10)m_h_SegmentXYEML->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 11)m_h_SegmentXYEOS->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 12)m_h_SegmentXYEOL->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 13)m_h_SegmentXYEES->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 14)m_h_SegmentXYEEL->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 15)m_h_SegmentXYCSS->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == 16)m_h_SegmentXYCSL->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
                    if(m_probe_segment_chamberIndex[chnum] == -1)m_h_SegmentXYUnknown->Fill(m_probe_segment_x[chnum],m_probe_segment_y[chnum]);
               }

               if(static_cast<Int_t>(pSA_sAddress) == 0 || static_cast<Int_t>(pSA_sAddress) == 1 || static_cast<Int_t>(pSA_sAddress) == 2 || static_cast<Int_t>(pSA_sAddress) == 3){
                    m_h_pOffPhivsSAsAddress->Fill(pSA_phims,pSA_sAddress);
                    m_h_MDTSPXYBI->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                    m_h_MDTSPXYBM->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                    m_h_MDTSPXYBO->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                    m_h_MDTSPXYBME->Fill(pSA_superpointR_BME*cos(pSA_roiphi),pSA_superpointR_BME*sin(pSA_roiphi));
                    m_h_MDTSPZR->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                    m_h_MDTSPZR->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                    m_h_MDTSPZR->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                    for(Int_t index = 0;index < 10;index++){
                         Int_t buf_index = static_cast<Int_t>(m_probe_segment_etaIndex[index]);
                         if(m_probe_segment_chamberIndex[index] <= 6){
                              switch(buf_index){
                                   case -6:
                                   m_h_MDTSPXYEtaIndexMinus6->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexMinus6->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexMinus6->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_SegmentXYEtaIndexMinus6->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case -5:
                                   m_h_MDTSPXYEtaIndexMinus5->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexMinus5->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexMinus5->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_SegmentXYEtaIndexMinus5->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case -4:
                                   m_h_MDTSPXYEtaIndexMinus4->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexMinus4->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexMinus4->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_SegmentXYEtaIndexMinus4->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case -3:
                                   m_h_MDTSPXYEtaIndexMinus3->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexMinus3->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexMinus3->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_SegmentXYEtaIndexMinus3->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case -2:
                                   m_h_MDTSPXYEtaIndexMinus2->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexMinus2->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexMinus2->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_SegmentXYEtaIndexMinus2->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case -1:
                                   m_h_MDTSPXYEtaIndexMinus1->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexMinus1->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexMinus1->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_SegmentXYEtaIndexMinus1->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case 1:
                                   m_h_MDTSPXYEtaIndexPlus1->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexPlus1->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexPlus1->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_SegmentXYEtaIndexPlus1->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case 2:
                                   m_h_MDTSPXYEtaIndexPlus2->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexPlus2->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexPlus2->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_SegmentXYEtaIndexPlus2->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case 3:
                                   m_h_MDTSPXYEtaIndexPlus3->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexPlus3->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexPlus3->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_SegmentXYEtaIndexPlus3->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case 4:
                                   m_h_MDTSPXYEtaIndexPlus4->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexPlus4->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexPlus4->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_SegmentXYEtaIndexPlus4->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case 5:
                                   m_h_MDTSPXYEtaIndexPlus5->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexPlus5->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexPlus5->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_SegmentXYEtaIndexPlus5->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   case 6:
                                   m_h_MDTSPXYEtaIndexPlus6->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexPlus6->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
                                   m_h_MDTSPXYEtaIndexPlus6->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_SegmentXYEtaIndexPlus6->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                                   break;
                                   default:
                                   break;
                              }
                         }
                    }
                    for(Int_t size = 0;size < (signed int)pSA_mdtZ->size();size++){
                         if(pSA_mdthitChamber->at(size) == 0 || pSA_mdthitChamber->at(size) == 1 || pSA_mdthitChamber->at(size) == 2)m_h_MDTHitXY->Fill(pSA_mdtR->at(size)*cos(pSA_mdtPhi->at(size)),pSA_mdtR->at(size)*sin(pSA_mdtPhi->at(size)));
                         m_h_MDTHitZR->Fill(pSA_mdtZ->at(size),pSA_mdtR->at(size));
                    }
                    for(Int_t size = 0;size < (signed int)pSA_rpcX->size();size++){
                         m_h_RPCHitXY->Fill(pSA_rpcX->at(size),pSA_rpcY->at(size));
                    }

                    Int_t buf_numsegment = 0;
                    for(Int_t index = 0;index < 10;index++){
                         if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_EtaIndexvsSAResPt->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                         if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_SegmentXY->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                         if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0){
                              m_h_SegmentZR->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              if(m_probe_segment_chamberIndex[index] == 1)m_h_SegmentZRBIL->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              if(m_probe_segment_chamberIndex[index] == 3)m_h_SegmentZRBML->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              buf_numsegment++;
                         }
                    }
                    m_h_NumSegment->Fill(buf_numsegment);
               }
               m_h_pOffPhivsSAPhims->Fill(m_poff_phi,pSA_phims);

          //CB
               if(CutCB(pCB_pass)){
                    Double_t textCB_dR = TMath::Sqrt(pow(m_tCB_eta - m_toff_exteta,2) + pow(m_tCB_phi - m_toff_extphi,2));
                    pextCB_dR = TMath::Sqrt(pow(pCB_eta - m_poff_exteta,2) + pow(pCB_phi - m_poff_extphi,2));
                    Double_t resCB_pt = std::fabs(m_poff_pt)/std::fabs(pCB_pt) - 1.0;
                    m_h_pCBPt->Fill(std::fabs(pCB_pt*0.001));
                    m_h_pCBdR->Fill(pCB_dR);
                    m_h_tExtCBdR->Fill(textCB_dR);
                    m_h_pExtCBdR->Fill(pextCB_dR);
                    m_h_eCBPt->Fill(std::fabs(m_poff_pt*0.001));
                    if(PlateauCut(std::fabs(m_poff_pt*0.001))){
                         m_h_eCBEta->Fill(m_poff_eta);
                         m_h_eCBPhi->Fill(m_poff_phi);
                         m_h_eCBAipc->Fill(m_aipc);
                    }
                    m_h_pCBResPt->Fill(resCB_pt);
                    switch(EtaDistribution()){
                         case 0:
                         m_h_eCBPtBarrel->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                         case 1:
                         m_h_eCBPtTransition->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                         case 2:
                         m_h_eCBPtEnd->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                         case 3:
                         m_h_eCBPtForward->Fill(std::fabs(m_poff_pt*0.001));
                         break;
                         default:
                         break;
                    }
                    if(DicisionBarrel(m_poff_eta)){
                         m_h_pCBResPtBarrel->Fill(resCB_pt);
                    }else{
                         m_h_pCBResPtEndcap->Fill(resCB_pt);
                    }

          //EF
                    if(CutEF(pEF_pass)){
                         Double_t textEF_dR = TMath::Sqrt(pow(m_tEF_eta - m_toff_exteta,2) + pow(m_tEF_phi - m_toff_extphi,2));
                         pextEF_dR = TMath::Sqrt(pow(pEF_eta - m_poff_exteta,2) + pow(pEF_phi - m_poff_extphi,2));
                         Double_t resEF_pt = std::fabs(m_poff_pt)/std::fabs(pEF_pt) - 1.0;
                         m_h_pEFPt->Fill(std::fabs(pEF_pt*0.001));
                         m_h_pEFdR->Fill(pEF_dR);
                         m_h_tExtEFdR->Fill(textEF_dR);
                         m_h_pExtEFdR->Fill(pextEF_dR);
                         m_h_eEFPt->Fill(std::fabs(m_poff_pt*0.001));
                         if(PlateauCut(std::fabs(m_poff_pt*0.001))){
                              m_h_eEFEta->Fill(m_poff_eta);
                              m_h_eEFPhi->Fill(m_poff_phi);
                              m_h_eEFAipc->Fill(m_aipc);
                         }
                         m_h_pEFResPt->Fill(resEF_pt);
                         switch(EtaDistribution()){
                              case 0:
                              m_h_eEFPtBarrel->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                              case 1:
                              m_h_eEFPtTransition->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                              case 2:
                              m_h_eEFPtEnd->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                              case 3:
                              m_h_eEFPtForward->Fill(std::fabs(m_poff_pt*0.001));
                              break;
                              default:
                              break;
                         }
                         if(DicisionBarrel(m_poff_eta)){
                              m_h_eEFPtBarrel->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_pEFResPtBarrel->Fill(resEF_pt);
                         }else{
                              m_h_eEFPtEnd->Fill(std::fabs(m_poff_pt*0.001));
                              m_h_pEFResPtEndcap->Fill(resEF_pt);
                         }
                    }
               }
          }
     }
}

}//Execute

void Efficiency::Finalize(TFile *tf1){
     CalcEff ceff;
     tf1->cd();
     //SetCondition
     //title,file title,yoffset,top margin,bottom margin,left margin,right margin
     m_h_SARoIPhiLS->Write();
     m_h_SARoIPhiSS->Write();
     m_h_SAPhimsLS->Write();
     m_h_pOffPhiLS->Write();
     m_h_SegmentNumberBIM->Write();
     m_h_SegmentNumberBIR->Write();
     for(Int_t i = 0; i < 16;i++){
          m_h_SectorPhi[i]->Write();
          m_h_SectorRoIPhi[i]->Write();
          for(Int_t j = 0; j < 8;j++){
               m_h_DivideEtaOverPhiResPtLarge[i][j]->Write();
               m_h_DivideEtaOverPhiResPtSmall[i][j]->Write();
          }
     }
     for(Int_t i = 0;i < 4;i++){
          m_h_PtMethod[i]->Write();
          m_h_PtMethodOver[i]->Write();
          m_h_PtSA[i]->Write();
          m_h_ResPtSA[i]->Write();
          m_h_PtSAOver[i]->Write();
          m_h_ResPtSAOver[i]->Write();
     }

     for(Int_t i = 0;i < 3;i++){
          m_h_L2MuonSAvsOfflinePt[i]->Write();
          m_h_L2MuonSAPtAlphavsBeta[i]->Write();
          m_h_L2MuonSAPtAlphavsTGC[i]->Write();
          m_h_L2MuonSAPtBetavsTGC[i]->Write();
     }

          //base,target
     cout<<"efficiency start!"<<endl;

     ceff.SetCondition("L1Efficiency","L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eOffPt,m_h_eL1Pt,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiency","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1Pt,m_h_eSAPt,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("CBEfficiency","muComb Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eSAPt,m_h_eCBPt,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("EFEfficiency","EventFilter Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eCBPt,m_h_eEFPt,m_binmax,200,m_efficiency_xerr);

     ceff.SetCondition("L1EfficiencyEta","L1 Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyeta(m_h_eOffEta,m_h_eL1Eta);
     ceff.SetCondition("SAEfficiencyEta","L2MuonSA Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyeta(m_h_eL1Eta,m_h_eSAEta);
     ceff.SetCondition("CBEfficiencyEta","muComb Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyeta(m_h_eSAEta,m_h_eCBEta);
     ceff.SetCondition("EFEfficiencyEta","EventFilter Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyeta(m_h_eCBEta,m_h_eEFEta);
     ceff.SetCondition("L1EfficiencyPhi","L1 Efficiency;offline phi;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyphi(m_h_eOffPhi,m_h_eL1Phi);
     ceff.SetCondition("SAEfficiencyPhi","L2MuonSA Efficiency;offline phi;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyphi(m_h_eL1Phi,m_h_eSAPhi);
     ceff.SetCondition("CBEfficiencyPhi","muComb Efficiency;offline phi;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyphi(m_h_eSAPhi,m_h_eCBPhi);
     ceff.SetCondition("EFEfficiencyPhi","EventFilter Efficiency;offline phi;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyphi(m_h_eCBPhi,m_h_eEFPhi);

     ceff.SetCondition("L1EfficiencyPileup","L1 Efficiency;pileup;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencypileup(m_h_eOffAipc,m_h_eL1Aipc);
     ceff.SetCondition("SAEfficiencyPileup","L2MuonSA Efficiency;pileup;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencypileup(m_h_eL1Aipc,m_h_eSAAipc);
     ceff.SetCondition("CBEfficiencyPileup","muComb Efficiency;pileup;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencypileup(m_h_eSAAipc,m_h_eCBAipc);
     ceff.SetCondition("EFEfficiencyPileup","EventFilter Efficiency;pileup;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencypileup(m_h_eCBAipc,m_h_eEFAipc);

     ceff.SetCondition("L1EfficiencyBarrel","L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eOffPtBarrel,m_h_eL1PtBarrel,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyBarrel","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrel,m_h_eSAPtBarrel,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("CBEfficiencyBarrel","muComb Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eSAPtBarrel,m_h_eCBPtBarrel,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("EFEfficiencyBarrel","EventFilter Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eCBPtBarrel,m_h_eEFPtBarrel,m_binmax,300,m_efficiency_xerr);

     ceff.SetCondition("L1EfficiencyEnd","L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eOffPtEnd,m_h_eL1PtEnd,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyEnd","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtEnd,m_h_eSAPtEnd,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("CBEfficiencyEnd","muComb Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eSAPtEnd,m_h_eCBPtEnd,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("EFEfficiencyEnd","EventFilter Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eCBPtEnd,m_h_eEFPtEnd,m_binmax,300,m_efficiency_xerr);

     ceff.SetCondition("L1EfficiencyForward","L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eOffPtForward,m_h_eL1PtForward,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyForward","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtForward,m_h_eSAPtForward,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("CBEfficiencyForward","muComb Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eSAPtForward,m_h_eCBPtForward,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("EFEfficiencyForward","EventFilter Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eCBPtForward,m_h_eEFPtForward,m_binmax,300,m_efficiency_xerr);

     ceff.SetCondition("L1EfficiencyTransition","L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eOffPtTransition,m_h_eL1PtTransition,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyTransition","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtTransition,m_h_eSAPtTransition,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("CBEfficiencyTransition","muComb Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eSAPtTransition,m_h_eCBPtTransition,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("EFEfficiencyTransition","EventFilter Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eCBPtTransition,m_h_eEFPtTransition,m_binmax,300,m_efficiency_xerr);

     ceff.SetCondition("L1EfficiencyLargePlus","L1 Large Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eOffPtLargePlus,m_h_eL1PtLargePlus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargePlus","SA Large Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtLargePlus,m_h_eSAPtLargePlus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeNormal","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtLargeNormal,m_h_eSAPtLargeNormal,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("L1EfficiencyLSPlus","L1 LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eOffPtLSPlus,m_h_eL1PtLSPlus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLSPlus","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtLSPlus,m_h_eSAPtLSPlus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("L1EfficiencyLargeMinus","L1 Large Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eOffPtLargeMinus,m_h_eL1PtLargeMinus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeMinus","SA Large Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtLargeMinus,m_h_eSAPtLargeMinus,m_binmax,300,m_efficiency_xerr);

     ceff.SetCondition("L1EfficiencyLSMinus","L1 LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eOffPtLSMinus,m_h_eL1PtLSMinus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLSMinus","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtLSMinus,m_h_eSAPtLSMinus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLSPlusS11","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtLSPlusS11,m_h_eSAPtLSPlusS11,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLSMinusS11","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtLSMinusS11,m_h_eSAPtLSMinusS11,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLSPlusS15","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtLSPlusS15,m_h_eSAPtLSPlusS15,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLSMinusS15","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtLSMinusS15,m_h_eSAPtLSMinusS15,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLSPlusS11outer","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtLSPlusS11outer,m_h_eSAPtLSPlusS11outer,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLSMinusS11outer","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtLSMinusS11outer,m_h_eSAPtLSMinusS11outer,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLSPlusS15outer","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtLSPlusS15outer,m_h_eSAPtLSPlusS15outer,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLSMinusS15outer","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtLSMinusS15outer,m_h_eSAPtLSMinusS15outer,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLSPlusS11inner","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtLSPlusS11inner,m_h_eSAPtLSPlusS11inner,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLSMinusS11inner","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtLSMinusS11inner,m_h_eSAPtLSMinusS11inner,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLSPlusS15inner","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtLSPlusS15inner,m_h_eSAPtLSPlusS15inner,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLSMinusS15inner","SA LargeSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtLSMinusS15inner,m_h_eSAPtLSMinusS15inner,m_binmax,300,m_efficiency_xerr);

     ceff.SetCondition("L1EfficiencySmallPlus","L1 Small Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eOffPtSmallPlus,m_h_eL1PtSmallPlus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencySmallPlus","SA Small Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtSmallPlus,m_h_eSAPtSmallPlus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("L1EfficiencySmallMinus","L1 Small Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eOffPtSmallMinus,m_h_eL1PtSmallMinus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencySmallMinus","SA Small Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtSmallMinus,m_h_eSAPtSmallMinus,m_binmax,300,m_efficiency_xerr);

     ceff.SetCondition("L1EfficiencySSPlus","L1 SmallSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eOffPtSSPlus,m_h_eL1PtSSPlus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencySSPlus","SA SmallSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtSSPlus,m_h_eSAPtSSPlus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("L1EfficiencySSMinus","L1 SmallSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eOffPtSSMinus,m_h_eL1PtSSMinus,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencySSMinus","SA SmallSpecial Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtSSMinus,m_h_eSAPtSSMinus,m_binmax,300,m_efficiency_xerr);

     ceff.SetCondition("SAEfficiencyBarrelWithoutBIM","SA Barrel without BIM Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrelWithoutBIM,m_h_eSAPtBarrelWithoutBIM,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLSWithoutBIM","SA LS without BIM Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtLSWithoutBIM,m_h_eSAPtLSWithoutBIM,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLargeIncBIM","SA Barrel inc BIM Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrelIncBIM,m_h_eSAPtBarrelIncBIM,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyLSIncBIM","SA LS inc BIM Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtLSIncBIM,m_h_eSAPtLSIncBIM,m_binmax,300,m_efficiency_xerr);

     ceff.SetCondition("SA2DEfficiency","L1vsL2MuonSA Efficiency;offline eta;offline phi",1.5,0.1,0.1,0.105,0.165);
     ceff.SetConditionbin(m_nbin_eta,m_nbin_phi,m_eta_max,m_phi_max);
     ceff.DrawEfficiency2D(m_h_eL1EtaPhi,m_h_eSAEtaPhi);
     ceff.SetCondition("L12DEfficiency","L1 Efficiency;offline eta;offline phi",1.5,0.1,0.1,0.105,0.165);
     ceff.SetConditionbin(m_nbin_eta,m_nbin_phi,m_eta_max,m_phi_max);
     ceff.DrawEfficiency2D(m_h_eOffEtaPhi,m_h_eL1EtaPhi);

     ceff.SetCondition("SAEfficiencyBarrel1SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrel1SP,m_h_eSAPtBarrel1SP,m_binmax,200,m_efficiency_xerr); 
     ceff.SetCondition("SAEfficiencyBarrel2SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrel2SP,m_h_eSAPtBarrel2SP,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyBarrel3SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrel3SP,m_h_eSAPtBarrel3SP,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyBarrelIM","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrelIM,m_h_eSAPtBarrelIM,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyEndcap1SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtEndcap1SP,m_h_eSAPtEndcap1SP,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyEndcap2SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtEndcap2SP,m_h_eSAPtEndcap2SP,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyEndcap3SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtEndcap3SP,m_h_eSAPtEndcap3SP,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyEndcapIM","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtEndcapIM,m_h_eSAPtEndcapIM,m_binmax,200,m_efficiency_xerr);

     ceff.SetCondition("SAEfficiencyForward1SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtForward1SP,m_h_eSAPtForward1SP,m_binmax,200,m_efficiency_xerr); 
     ceff.SetCondition("SAEfficiencyForward2SP","L2MP Fficienuoncy;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtForward2SP,m_h_eSAPtForward2SP,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyForward3SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtForward3SP,m_h_eSAPtForward3SP,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyForwardIM","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtForwardIM,m_h_eSAPtForwardIM,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyTransition1SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtTransition1SP,m_h_eSAPtTransition1SP,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyTransition2SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtTransition2SP,m_h_eSAPtTransition2SP,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyTransition3SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtTransition3SP,m_h_eSAPtTransition3SP,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyTransitionIM","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtTransitionIM,m_h_eSAPtTransitionIM,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyBarrelRoI","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrelRoI,m_h_eSAPtBarrelRoI,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyBarrel1SPRoI","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrel1SPRoI,m_h_eSAPtBarrel1SPRoI,m_binmax,200,m_efficiency_xerr); 
     ceff.SetCondition("SAEfficiencyBarrel2SPRoI","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrel2SPRoI,m_h_eSAPtBarrel2SPRoI,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyBarrel3SPRoI","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrel3SPRoI,m_h_eSAPtBarrel3SPRoI,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyBarrelIMRoI","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrelIMRoI,m_h_eSAPtBarrelIMRoI,m_binmax,200,m_efficiency_xerr);
     cout<<"eff end"<<endl;

     m_h_pOffPt->Write();
     m_h_pL1Pt->Write();
     m_h_pSAPt->Write();
     m_h_pCBPt->Write();
     m_h_pEFPt->Write();
     m_h_pL1dR->Write();
     m_h_pSAdR->Write();
     m_h_pCBdR->Write();
     m_h_pEFdR->Write();
     m_h_tExtL1dR->Write();
     m_h_tExtSAdR->Write();
     m_h_tExtCBdR->Write();
     m_h_tExtEFdR->Write();
     m_h_pExtL1dR->Write();
     m_h_pExtSAdR->Write();
     m_h_pExtCBdR->Write();
     m_h_pExtEFdR->Write();
     m_h_pOffvsSAPt->Write();
     m_h_pOffPhivsSAsAddress->Write();
     m_h_pOffPhivsSAPhims->Write();
     m_h_pSAPhivspSAPhims->Write();
     m_h_pSAPhivspSAPhibe->Write();
     m_h_CountSA->Write();
          //
     m_h_pOffPtBIM2Station->Write();
     m_h_pOffPtBIR2Station->Write();
     m_h_pOffPtLarge2Station->Write();
     m_h_pOffPtSmall2Station->Write();
     m_h_pOffPtBIM2StationIM->Write();
     m_h_pOffPtBIR2StationIM->Write();
     m_h_pOffPtLarge2StationIM->Write();
     m_h_pOffPtSmall2StationIM->Write();
     m_h_pOffPtBIM2StationIO->Write();
     m_h_pOffPtBIR2StationIO->Write();
     m_h_pOffPtLarge2StationIO->Write();
     m_h_pOffPtSmall2StationIO->Write();
     m_h_pOffPtBIM2StationMO->Write();
     m_h_pOffPtBIR2StationMO->Write();
     m_h_pOffPtLarge2StationMO->Write();
     m_h_pOffPtSmall2StationMO->Write();
          //
     m_h_pOffEtaBIM2Station->Write();
     m_h_pOffEtaBIR2Station->Write();
     m_h_pOffEtaLarge2Station->Write();
     m_h_pOffEtaSmall2Station->Write();
     m_h_pOffEtaBIM2StationIM->Write();
     m_h_pOffEtaBIR2StationIM->Write();
     m_h_pOffEtaLarge2StationIM->Write();
     m_h_pOffEtaSmall2StationIM->Write();
     m_h_pOffEtaBIM2StationIO->Write();
     m_h_pOffEtaBIR2StationIO->Write();
     m_h_pOffEtaLarge2StationIO->Write();
     m_h_pOffEtaSmall2StationIO->Write();
     m_h_pOffEtaBIM2StationMO->Write();
     m_h_pOffEtaBIR2StationMO->Write();
     m_h_pOffEtaLarge2StationMO->Write();
     m_h_pOffEtaSmall2StationMO->Write();
          //
     m_h_pOffEtaBIM3Station->Write();
     m_h_pOffEtaBIR3Station->Write();
     m_h_pOffEtaLarge3Station->Write();
     m_h_pOffEtaSmall3Station->Write();

     m_h_SectorvsPhi->Write();
     m_h_IndexvsEta->Write();
     m_h_OverPhi->Write();
     m_h_MDTChamber->Write();
     m_h_L1pSAsAddress->Write();
     m_h_pSAsAddress->Write();

     m_h_MDTHitXY->Write();
     m_h_RPCHitXY->Write();
     m_h_MDTHitZR->Write();
     m_h_MDTSPXYBI->Write();
     m_h_MDTSPXYBM->Write();
     m_h_MDTSPXYBO->Write();
     m_h_MDTSPXYBME->Write();
     m_h_MDTSPXBI->Write();
     m_h_MDTSPYBI->Write();
     m_h_MDTSPXYEtaIndexPlus1->Write();
     m_h_MDTSPXYEtaIndexPlus2->Write();
     m_h_MDTSPXYEtaIndexPlus3->Write();
     m_h_MDTSPXYEtaIndexPlus4->Write();
     m_h_MDTSPXYEtaIndexPlus5->Write();
     m_h_MDTSPXYEtaIndexPlus6->Write();
     m_h_MDTSPXYEtaIndexMinus1->Write();
     m_h_MDTSPXYEtaIndexMinus2->Write();
     m_h_MDTSPXYEtaIndexMinus3->Write();
     m_h_MDTSPXYEtaIndexMinus4->Write();
     m_h_MDTSPXYEtaIndexMinus5->Write();
     m_h_MDTSPXYEtaIndexMinus6->Write();
     m_h_MDTSPZR->Write();
     m_h_MDTSPZRLSPlus->Write();
     m_h_MDTSPZRLSMinus->Write();
     m_h_MDTSPZRLSPlusS11outer->Write();
     m_h_MDTSPZRLSPlusS11inner->Write();
     m_h_MDTSPZRLSPlusS15outer->Write();
     m_h_MDTSPZRLSPlusS15inner->Write();
     m_h_MDTSPZRLSMinusS11outer->Write();
     m_h_MDTSPZRLSMinusS11inner->Write();
     m_h_MDTSPZRLSMinusS15outer->Write();
     m_h_MDTSPZRLSMinusS15inner->Write();
     m_h_SegmentXY->Write();
     m_h_SegmentXYEtaIndexPlus1->Write();
     m_h_SegmentXYEtaIndexPlus2->Write();
     m_h_SegmentXYEtaIndexPlus3->Write();
     m_h_SegmentXYEtaIndexPlus4->Write();
     m_h_SegmentXYEtaIndexPlus5->Write();
     m_h_SegmentXYEtaIndexPlus6->Write();
     m_h_SegmentXYEtaIndexMinus1->Write();
     m_h_SegmentXYEtaIndexMinus2->Write();
     m_h_SegmentXYEtaIndexMinus3->Write();
     m_h_SegmentXYEtaIndexMinus4->Write();
     m_h_SegmentXYEtaIndexMinus5->Write();
     m_h_SegmentXYEtaIndexMinus6->Write();
     m_h_SegmentXYBIL->Write();
     m_h_SegmentXYBML->Write();
     m_h_SegmentXYBOL->Write();
     m_h_SegmentXYBIS->Write();
     m_h_SegmentXYBMS->Write();
     m_h_SegmentXYBOS->Write();
     m_h_SegmentXYBEE->Write();
     m_h_SegmentXYEIL->Write();
     m_h_SegmentXYEML->Write();
     m_h_SegmentXYEOL->Write();
     m_h_SegmentXYEIS->Write();
     m_h_SegmentXYEMS->Write();
     m_h_SegmentXYEOS->Write();
     m_h_SegmentXYEES->Write();
     m_h_SegmentXYEEL->Write();
     m_h_SegmentXYCSS->Write();
     m_h_SegmentXYCSL->Write();
     m_h_SegmentXYUnknown->Write();
     m_h_SegmentZR->Write();
     m_h_SegmentZRBIL->Write();
     m_h_SegmentZRBML->Write();
     m_h_SegmentZRLSPlus->Write();
     m_h_SegmentZRLSMinus->Write();
     m_h_SegmentZRLSPlusS11outer->Write();
     m_h_SegmentZRLSPlusS11inner->Write();
     m_h_SegmentZRLSPlusS15outer->Write();
     m_h_SegmentZRLSPlusS15inner->Write();
     m_h_SegmentZRLSMinusS11outer->Write();
     m_h_SegmentZRLSMinusS11inner->Write();
     m_h_SegmentZRLSMinusS15outer->Write();
     m_h_SegmentZRLSMinusS15inner->Write();
     m_h_SegmentZRLSPlusS11outerBIL->Write();
     m_h_SegmentZRLSPlusS11innerBIL->Write();
     m_h_SegmentZRLSPlusS15outerBIL->Write();
     m_h_SegmentZRLSPlusS15innerBIL->Write();
     m_h_NumSegment->Write();
     m_h_NumSegmentLSBI->Write();
     m_h_NumSegmentLargeBI->Write();
     m_h_NumHit->Write();
     m_h_nPrecisionHitsNormal->Write();
     m_h_nPrecisionHitsSpecial->Write();
     m_h_SectorNormal->Write();
     m_h_SectorSpecial->Write();
     m_h_EtaIndexNormal->Write();
     m_h_EtaIndexSpecial->Write();
     m_h_MDTPhi->Write();
     m_h_MDTPhiLS->Write();
     m_h_MDTPhiLSBIL->Write();
     m_h_MDTR->Write();
     m_h_NumSP->Write();
     m_h_BIRSegment->Write();
     m_h_BIMSegment->Write();
     m_h_BIMrvsx->Write();

     m_h_eOffPt->Write();
     m_h_eL1Pt->Write();
     m_h_eSAPt->Write();
     m_h_eCBPt->Write();
     m_h_eEFPt->Write();
     m_h_eOffEta->Write();
     m_h_eL1Eta->Write();
     m_h_eSAEta->Write();
     m_h_eCBEta->Write();
     m_h_eEFEta->Write();
     m_h_eOffPhi->Write();
     m_h_eL1Phi->Write();
     m_h_eSAPhi->Write();
     m_h_eCBPhi->Write();
     m_h_eEFPhi->Write();
     m_h_eOffAipc->Write();
     m_h_eL1Aipc->Write();
     m_h_eSAAipc->Write();
     m_h_eCBAipc->Write();
     m_h_eEFAipc->Write();
     m_h_eOffPtBarrel->Write();
     m_h_eL1PtBarrel->Write();
     m_h_eSAPtBarrel->Write();
     m_h_eCBPtBarrel->Write();
     m_h_eEFPtBarrel->Write();
     m_h_eOffPtEnd->Write();
     m_h_eL1PtEnd->Write();
     m_h_eSAPtEnd->Write();
     m_h_eCBPtEnd->Write();
     m_h_eEFPtEnd->Write();
     m_h_eOffPtForward->Write();
     m_h_eL1PtForward->Write();
     m_h_eSAPtForward->Write();
     m_h_eCBPtForward->Write();
     m_h_eEFPtForward->Write();
     m_h_eOffPtTransition->Write();
     m_h_eL1PtTransition->Write();
     m_h_eSAPtTransition->Write();
     m_h_eCBPtTransition->Write();
     m_h_eEFPtTransition->Write();
     m_h_eOffPtLargePlus->Write();
     m_h_eL1PtLargePlus->Write();
     m_h_eSAPtLargePlus->Write();
     m_h_eL1PtLargeNormal->Write();
     m_h_eSAPtLargeNormal->Write();
     m_h_eOffPtLargeMinus->Write();
     m_h_eL1PtLargeMinus->Write();
     m_h_eSAPtLargeMinus->Write();
     m_h_eOffPtLSPlus->Write();
     m_h_eL1PtLSPlus->Write();
     m_h_eSAPtLSPlus->Write();
     m_h_eL1PtLSPlusS11->Write();
     m_h_eSAPtLSPlusS11->Write();
     m_h_eL1PtLSPlusS15->Write();
     m_h_eSAPtLSPlusS15->Write();
     m_h_eL1PtLSPlusS11inner->Write();
     m_h_eSAPtLSPlusS11inner->Write();
     m_h_eL1PtLSPlusS11outer->Write();
     m_h_eSAPtLSPlusS11outer->Write();
     m_h_eL1PtLSPlusS15outer->Write();
     m_h_eSAPtLSPlusS15outer->Write();
     m_h_eL1PtLSPlusS15inner->Write();
     m_h_eSAPtLSPlusS15inner->Write();
     m_h_eOffPtLSMinus->Write();
     m_h_eL1PtLSMinus->Write();
     m_h_eSAPtLSMinus->Write();
     m_h_eL1PtLSMinusS11->Write();
     m_h_eSAPtLSMinusS11->Write();
     m_h_eL1PtLSMinusS15->Write();
     m_h_eSAPtLSMinusS15->Write();
     m_h_eL1PtLSMinusS11inner->Write();
     m_h_eSAPtLSMinusS11inner->Write();
     m_h_eL1PtLSMinusS11outer->Write();
     m_h_eSAPtLSMinusS11outer->Write();
     m_h_eL1PtLSMinusS15outer->Write();
     m_h_eSAPtLSMinusS15outer->Write();
     m_h_eL1PtLSMinusS15inner->Write();
     m_h_eSAPtLSMinusS15inner->Write();
     m_h_eOffPtSmallPlus->Write();
     m_h_eL1PtSmallPlus->Write();
     m_h_eSAPtSmallPlus->Write();
     m_h_eOffPtSmallMinus->Write();
     m_h_eL1PtSmallMinus->Write();
     m_h_eSAPtSmallMinus->Write();
     m_h_eOffPtSSPlus->Write();
     m_h_eL1PtSSPlus->Write();
     m_h_eSAPtSSPlus->Write();
     m_h_eOffPtSSMinus->Write();
     m_h_eL1PtSSMinus->Write();
     m_h_eSAPtSSMinus->Write();
     m_h_eOffEtaPhi->Write();
     m_h_eL1EtaPhi->Write();
     m_h_eSAEtaPhi->Write();
     m_h_eL1PtBarrelWithoutBIM->Write();
     m_h_eSAPtBarrelWithoutBIM->Write();
     m_h_eL1PtLSWithoutBIM->Write();
     m_h_eSAPtLSWithoutBIM->Write();
     m_h_eL1PtBarrelIncBIM->Write();
     m_h_eSAPtBarrelIncBIM->Write();
     m_h_eL1PtLSIncBIM->Write();
     m_h_eSAPtLSIncBIM->Write();
     m_h_eL1PtBarrel1SP->Write();
     m_h_eSAPtBarrel1SP->Write();
     m_h_eL1PtBarrel2SP->Write();
     m_h_eSAPtBarrel2SP->Write();
     m_h_eL1PtBarrel3SP->Write();
     m_h_eSAPtBarrel3SP->Write();
     m_h_eL1PtBarrelIM->Write();
     m_h_eSAPtBarrelIM->Write();
     m_h_eL1PtEndcap1SP->Write();
     m_h_eSAPtEndcap1SP->Write();
     m_h_eL1PtEndcap2SP->Write();
     m_h_eSAPtEndcap2SP->Write();
     m_h_eL1PtEndcap3SP->Write();
     m_h_eSAPtEndcap3SP->Write();
     m_h_eL1PtEndcapIM->Write();
     m_h_eSAPtEndcapIM->Write();
     m_h_eL1PtForward1SP->Write();
     m_h_eSAPtForward1SP->Write();
     m_h_eL1PtForward2SP->Write();
     m_h_eSAPtForward2SP->Write();
     m_h_eL1PtForward3SP->Write();
     m_h_eSAPtForward3SP->Write();
     m_h_eL1PtForwardIM->Write();
     m_h_eSAPtForwardIM->Write();
     m_h_eL1PtTransition1SP->Write();
     m_h_eSAPtTransition1SP->Write();
     m_h_eL1PtTransition2SP->Write();
     m_h_eSAPtTransition2SP->Write();
     m_h_eL1PtTransition3SP->Write();
     m_h_eSAPtTransition3SP->Write();
     m_h_eL1PtTransitionIM->Write();
     m_h_eSAPtTransitionIM->Write();
     m_h_eL1PtBarrelRoI->Write();
     m_h_eSAPtBarrelRoI->Write();
     m_h_eL1PtBarrel1SPRoI->Write();
     m_h_eSAPtBarrel1SPRoI->Write();
     m_h_eL1PtBarrel2SPRoI->Write();
     m_h_eSAPtBarrel2SPRoI->Write();
     m_h_eL1PtBarrel3SPRoI->Write();
     m_h_eSAPtBarrel3SPRoI->Write();
     m_h_eL1PtBarrelIMRoI->Write();
     m_h_eSAPtBarrelIMRoI->Write();

     m_h_pSAResPt->Write();
     m_h_pCBResPt->Write();
     m_h_pEFResPt->Write();
     m_h_pSAResPtBarrel->Write();
     m_h_pCBResPtBarrel->Write();
     m_h_pEFResPtBarrel->Write();
     m_h_pSAResPtEndcap->Write();
     m_h_pCBResPtEndcap->Write();
     m_h_pEFResPtEndcap->Write();
     m_h_SAResPtLargePlus->Write();
     m_h_SAResPtLSPlus->Write();
     m_h_SAResPtSmallPlus->Write();
     m_h_SAResPtSSPlus->Write();
     m_h_SAResPtLargeMinus->Write();
     m_h_SAResPtLSMinus->Write();
     m_h_SAResPtSmallMinus->Write();
     m_h_SAResPtSSMinus->Write();
     m_h_SAResPtLSPlusS11outer->Write();
     m_h_SAResPtLSPlusS11inner->Write();
     m_h_SAResPtLSPlusS15outer->Write();
     m_h_SAResPtLSPlusS15inner->Write();
     m_h_SAResPtLSMinusS11outer->Write();
     m_h_SAResPtLSMinusS11inner->Write();
     m_h_SAResPtLSMinusS15outer->Write();
     m_h_SAResPtLSMinusS15inner->Write();
     m_h_pOffPtvsSAResPtLargePlus->Write();
     m_h_pOffPtvsSAResPtLSPlus->Write();
     m_h_pOffPtvsSAResPtSmallPlus->Write();
     m_h_pOffPtvsSAResPtSSPlus->Write();
     m_h_pOffPtvsSAResPtLargeMinus->Write();
     m_h_pOffPtvsSAResPtLSMinus->Write();
     m_h_pOffPtvsSAResPtSmallMinus->Write();
     m_h_pOffPtvsSAResPtSSMinus->Write();
     m_h_pOffPhivsSAResPtLargePlus->Write();
     m_h_pOffPhivsSAResPtLSPlus->Write();
     m_h_pOffPhivsSAResPtSmallPlus->Write();
     m_h_pOffPhivsSAResPtSSPlus->Write();
     m_h_pOffPhivsSAResPtLargeMinus->Write();
     m_h_pOffPhivsSAResPtLSMinus->Write();
     m_h_pOffPhivsSAResPtSmallMinus->Write();
     m_h_pOffPhivsSAResPtSSMinus->Write();
     m_h_pOffPtvsSAResPtLSPlusS15outer->Write();
     m_h_pOffPtvsSAResPtLSPlusS15inner->Write();
     m_h_pOffPtvsSAResPtLSPlusS11outer->Write();
     m_h_pOffPtvsSAResPtLSPlusS11inner->Write();
     m_h_pOffPtvsSAResPtLSMinusS15outer->Write();
     m_h_pOffPtvsSAResPtLSMinusS15inner->Write();
     m_h_pOffPtvsSAResPtLSMinusS11outer->Write();
     m_h_pOffPtvsSAResPtLSMinusS11inner->Write();
          //
     m_h_pOffPtvsSAResPtLarge2Station->Write();
     m_h_pOffPtvsSAResPtSmall2Station->Write();
     m_h_pOffPtvsSAResPtLSPlusS15outer2Station->Write();
     m_h_pOffPtvsSAResPtLSPlusS15inner2Station->Write();
     m_h_pOffPtvsSAResPtLSPlusS11outer2Station->Write();
     m_h_pOffPtvsSAResPtLSPlusS11inner2Station->Write();
     m_h_pOffPtvsSAResPtLSMinusS15outer2Station->Write();
     m_h_pOffPtvsSAResPtLSMinusS15inner2Station->Write();
     m_h_pOffPtvsSAResPtLSMinusS11outer2Station->Write();
     m_h_pOffPtvsSAResPtLSMinusS11inner2Station->Write();
          //
     m_h_pOffPtvsSAResPtLSPlusS15outer3Station->Write();
     m_h_pOffPtvsSAResPtLSPlusS15inner3Station->Write();
     m_h_pOffPtvsSAResPtLSPlusS11outer3Station->Write();
     m_h_pOffPtvsSAResPtLSPlusS11inner3Station->Write();
     m_h_pOffPtvsSAResPtLSMinusS15outer3Station->Write();
     m_h_pOffPtvsSAResPtLSMinusS15inner3Station->Write();
     m_h_pOffPtvsSAResPtLSMinusS11outer3Station->Write();
     m_h_pOffPtvsSAResPtLSMinusS11inner3Station->Write();
          //
     m_h_pOffPtvsSAResPtBIM2StationIM->Write();
     m_h_pOffPtvsSAResPtBIR2StationIM->Write();
     m_h_pOffPtvsSAResPtLarge2StationIM->Write();
     m_h_pOffPtvsSAResPtSmall2StationIM->Write();
     m_h_pOffPtvsSAResPtBIM2StationIO->Write();
     m_h_pOffPtvsSAResPtBIR2StationIO->Write();
     m_h_pOffPtvsSAResPtLarge2StationIO->Write();
     m_h_pOffPtvsSAResPtSmall2StationIO->Write();
     m_h_pOffPtvsSAResPtBIM2StationMO->Write();
     m_h_pOffPtvsSAResPtBIR2StationMO->Write();
     m_h_pOffPtvsSAResPtLarge2StationMO->Write();
     m_h_pOffPtvsSAResPtSmall2StationMO->Write();
          //
     m_h_pOffEtavsSAResPtLargePlus->Write();
     m_h_pOffEtavsSAResPtLSPlus->Write();
     m_h_pOffEtavsSAResPtSmallPlus->Write();
     m_h_pOffEtavsSAResPtSSPlus->Write();
     m_h_pOffEtavsSAResPtLargeMinus->Write();
     m_h_pOffEtavsSAResPtLSMinus->Write();
     m_h_pOffEtavsSAResPtSmallMinus->Write();
     m_h_pOffEtavsSAResPtSSMinus->Write();
     m_h_pOffEtavsSAResPtLSPlusS11outer->Write();
     m_h_pOffEtavsSAResPtLSPlusS11inner->Write();
     m_h_pOffEtavsSAResPtLSPlusS15outer->Write();
     m_h_pOffEtavsSAResPtLSPlusS15inner->Write();
     m_h_pOffEtavsSAResPtLSMinusS11outer->Write();
     m_h_pOffEtavsSAResPtLSMinusS11inner->Write();
     m_h_pOffEtavsSAResPtLSMinusS15outer->Write();
     m_h_pOffEtavsSAResPtLSMinusS15inner->Write();
     m_h_EtaIndexvsSAResPt->Write();
     m_h_EtaIndexvsSAResPtLargePlus->Write();
     m_h_EtaIndexvsSAResPtLSPlus->Write();
     m_h_EtaIndexvsSAResPtSmallPlus->Write();
     m_h_EtaIndexvsSAResPtSSPlus->Write();
     m_h_EtaIndexvsSAResPtLargeMinus->Write();
     m_h_EtaIndexvsSAResPtLSMinus->Write();
     m_h_EtaIndexvsSAResPtSmallMinus->Write();
     m_h_EtaIndexvsSAResPtSSMinus->Write();
     m_h_EtaIndexvsSAResPtLSPlusS11outer->Write();
     m_h_EtaIndexvsSAResPtLSPlusS11inner->Write();
     m_h_EtaIndexvsSAResPtLSPlusS15outer->Write();
     m_h_EtaIndexvsSAResPtLSPlusS15inner->Write();
     m_h_EtaIndexvsSAResPtLSMinusS11outer->Write();
     m_h_EtaIndexvsSAResPtLSMinusS11inner->Write();
     m_h_EtaIndexvsSAResPtLSMinusS15outer->Write();
     m_h_EtaIndexvsSAResPtLSMinusS15inner->Write();
          //
     m_h_EtaIndexvsSAResPtLSPlusS11outer2Station->Write();
     m_h_EtaIndexvsSAResPtLSPlusS11inner2Station->Write();
     m_h_EtaIndexvsSAResPtLSPlusS15outer2Station->Write();
     m_h_EtaIndexvsSAResPtLSPlusS15inner2Station->Write();
     m_h_EtaIndexvsSAResPtLSMinusS11outer2Station->Write();
     m_h_EtaIndexvsSAResPtLSMinusS11inner2Station->Write();
     m_h_EtaIndexvsSAResPtLSMinusS15outer2Station->Write();
     m_h_EtaIndexvsSAResPtLSMinusS15inner2Station->Write();
          //
     m_h_EtaIndexvsSAResPtLSPlusS11outer3Station->Write();
     m_h_EtaIndexvsSAResPtLSPlusS11inner3Station->Write();
     m_h_EtaIndexvsSAResPtLSPlusS15outer3Station->Write();
     m_h_EtaIndexvsSAResPtLSPlusS15inner3Station->Write();
     m_h_EtaIndexvsSAResPtLSMinusS11outer3Station->Write();
     m_h_EtaIndexvsSAResPtLSMinusS11inner3Station->Write();
     m_h_EtaIndexvsSAResPtLSMinusS15outer3Station->Write();
     m_h_EtaIndexvsSAResPtLSMinusS15inner3Station->Write();

     cout<<"residual end"<<endl;

}
