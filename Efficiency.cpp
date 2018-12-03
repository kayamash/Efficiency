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
          //offline
          if(!Cut_tagprobe(pEFTAG_pass))continue;
          if(i == 0 && static_cast<Int_t>(pSA_sAddress) == 1)m_h_offphi_LargeSpecial->Fill(m_poff_phi);
          m_h_poff_pt.at(i)->Fill(m_poff_pt*0.001);
          m_h_eoff_pt.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          if(std::fabs(m_poff_pt*0.001) > 40){
          //if(std::fabs(m_poff_pt*0.001) > 8){
               m_h_eoff_eta.at(i)->Fill(m_poff_eta);
               m_h_eoff_phi.at(i)->Fill(m_poff_phi);
               m_h_eoff_aipc.at(i)->Fill(m_aipc);
          }
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
          if(!Cut_L1(pL1_pass))continue;
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
          //if(std::fabs(m_poff_pt*0.001) > 8){
               m_h_eff_pL1_etaphi.at(i)->Fill(m_poff_eta,m_poff_phi);
               m_h_eL1_eta.at(i)->Fill(m_poff_eta);
               m_h_eL1_phi.at(i)->Fill(m_poff_phi);
               m_h_eL1_aipc.at(i)->Fill(m_aipc);
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
          if(!Cut_SA(pSA_pass,pSA_pt,i*m_thpitch))continue;
          Double_t textSA_dR = TMath::Sqrt(pow(m_tSA_eta - m_toff_exteta,2) + pow(m_tSA_phi - m_toff_extphi,2));
          pextSA_dR = TMath::Sqrt(pow(pSA_eta - m_poff_exteta,2) + pow(pSA_phi - m_poff_extphi,2));
          Double_t resSA_pt = std::fabs(m_poff_pt*0.001)/std::fabs(pSA_pt) - 1.0;
          Double_t buf_pSA_dR = TMath::Sqrt(pow(pSA_eta - m_poff_eta,2) + pow(pSA_phi - m_poff_phi,2));
          Double_t buf_eta = 0;
          for(Int_t size = 0;size < (signed int)pSA_mdtZ->size();size++){
               buf_eta += -TMath::Log((sqrt(pow(pSA_mdtZ->at(size),2) + pow(pSA_mdtR->at(size),2)) - pSA_mdtZ->at(size))/(sqrt(pow(pSA_mdtZ->at(size),2) + pow(pSA_mdtR->at(size),2)) + pow(pSA_mdtZ->at(size),2)))/2.0;
          }
          Double_t ave_mdteta = buf_eta/static_cast<Double_t>(pSA_mdtZ->size());
          m_h_avemdteta.at(i)->Fill(ave_mdteta);

          m_h_pSA_pt.at(i)->Fill(std::fabs(pSA_pt));
          m_h_pSA_dR.at(i)->Fill(buf_pSA_dR);
          m_h_textSA_dR.at(i)->Fill(textSA_dR);
          m_h_pextSA_dR.at(i)->Fill(pextSA_dR);
          m_h_eSA_pt.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          m_h_pSA_respt.at(i)->Fill(resSA_pt);
          m_h_pSAphivspSAphims.at(i)->Fill(pSA_phi,pSA_phims);
          m_h_pSAphivspSAphibe.at(i)->Fill(pSA_phi,pSA_phibe);
          if(Dicision_barrel(m_poff_eta)){
               m_h_eSA_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          }else{
               m_h_eSA_pt_end.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          }
          if(std::fabs(m_poff_pt*0.001) > 40){
          //if(std::fabs(m_poff_pt*0.001) > 8){
               m_h_eff_pSA_etaphi.at(i)->Fill(m_poff_eta,m_poff_phi);
               m_h_eSA_eta.at(i)->Fill(m_poff_eta);
               m_h_eSA_phi.at(i)->Fill(m_poff_phi);
               m_h_eSA_aipc.at(i)->Fill(m_aipc);
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
                         Double_t buf_BIsegmentR = 0;
                         Double_t buf_resR = 999999;
                         for(Int_t index = 0;index < 10;index++){
                              if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeminus.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                              if(sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)) <= 5500){
                                   if(buf_resR >= sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2))){
                                        buf_BIsegmentR = sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2));
                                        buf_resR = buf_BIsegmentR - pSA_superpointR_BI;
                                   }
                              }
                         }
                         if(buf_BIsegmentR != 0){
                              m_h_segSP_diffR_LargeBI.at(i)->Fill(buf_BIsegmentR - pSA_superpointR_BI);
                              m_h_segSP_resR_LargeBI.at(i)->Fill(1.0/pSA_superpointR_BI - 1.0/buf_BIsegmentR);
                              Int_t buf_index = fabs(static_cast<Int_t>(m_probe_segment_etaIndex[index]));
                                        switch(buf_index){
                                             case 1:
                                                  m_h_segSP_resR_LargeBIetaindex1.at(i)->Fill(1.0/pSA_superpointR_BI - 1.0/buf_BIsegmentR);
                                                  break;
                                             case 2:
                                                  m_h_segSP_resR_LargeBIetaindex2.at(i)->Fill(1.0/pSA_superpointR_BI - 1.0/buf_BIsegmentR);
                                                  break;
                                             case 3:
                                                  m_h_segSP_resR_LargeBIetaindex3.at(i)->Fill(1.0/pSA_superpointR_BI - 1.0/buf_BIsegmentR);
                                                  break;
                                             case 4:
                                                  m_h_segSP_resR_LargeBIetaindex4.at(i)->Fill(1.0/pSA_superpointR_BI - 1.0/buf_BIsegmentR);
                                                  break;
                                             case 5:
                                                  m_h_segSP_resR_LargeBIetaindex5.at(i)->Fill(1.0/pSA_superpointR_BI - 1.0/buf_BIsegmentR);
                                                  break;
                                             case 6:
                                                  m_h_segSP_resR_LargeBIetaindex6.at(i)->Fill(1.0/pSA_superpointR_BI - 1.0/buf_BIsegmentR);
                                                  break;
                                             default :
                                                  break;
                                        }
                         }
                    }
                    m_countLarge.at(i)++;
                    break;
               case 1:
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
                              if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialplus.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         }
                         m_h_mdtSPZR_LargeSpecialplus.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_mdtSPZR_LargeSpecialplus.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_mdtSPZR_LargeSpecialplus.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);

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
                                        if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialplus15out.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                                   }
                                   m_h_mdtSPZR_LargeSpecialplus15out.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                                   m_h_mdtSPZR_LargeSpecialplus15out.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                                   m_h_mdtSPZR_LargeSpecialplus15out.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);

                              }else{
                                   m_h_off_ptvsSA_resptLargeSpecialplus15in.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                   m_h_SA_resptLargeSpecialplus15in.at(i)->Fill(resSA_pt);
                                   m_h_offetavsSA_resptLargeSpecialplus15in.at(i)->Fill(m_poff_eta,resSA_pt);
                                   m_h_mdtetavsSA_resptLargeSpecialplus15in.at(i)->Fill(ave_mdteta,resSA_pt);
                                   m_h_eSA_pt_LargeSpecialplus15in.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                   for(Int_t index = 0;index < 10;index++){
                                       if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialplus15in.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                       if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialplus15in.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                                   }
                                   if(std::fabs(m_poff_pt*0.001) > 30.0 && std::fabs(m_poff_pt*0.001) < 50.0){
                                        m_h_highoffetavsSA_resptLargeSpecialplus15in.at(i)->Fill(m_poff_eta,resSA_pt);
                                   }
                                   m_h_mdtSPZR_LargeSpecialplus15in.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                                   m_h_mdtSPZR_LargeSpecialplus15in.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                                   m_h_mdtSPZR_LargeSpecialplus15in.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                              }
                         }
                         if(pSA_roiphi < -2.0 && pSA_roiphi > -2.6){//11
                              m_h_mdtSPX_BI.at(i)->Fill(pSA_superpointR_BI*cos(pSA_roiphi));
                              m_h_mdtSPY_BI.at(i)->Fill(pSA_superpointR_BI*sin(pSA_roiphi));
                              if(pSA_roiphi > -2.4){
                                   m_h_off_ptvsSA_resptLargeSpecialplus11in.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                   m_h_SA_resptLargeSpecialplus11in.at(i)->Fill(resSA_pt);
                                   m_h_offetavsSA_resptLargeSpecialplus11in.at(i)->Fill(m_poff_eta,resSA_pt);
                                   m_h_mdtetavsSA_resptLargeSpecialplus11in.at(i)->Fill(ave_mdteta,resSA_pt);
                                   m_h_eSA_pt_LargeSpecialplus11in.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                   for(Int_t index = 0;index < 10;index++){
                                        if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialplus11in.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                        if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialplus11in.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                                   }
                                   if(std::fabs(m_poff_pt*0.001) > 30.0 && std::fabs(m_poff_pt*0.001) < 50.0){
                                        m_h_highoffetavsSA_resptLargeSpecialplus11in.at(i)->Fill(m_poff_eta,resSA_pt);
                                   }
                                   m_h_mdtSPZR_LargeSpecialplus11in.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                                   m_h_mdtSPZR_LargeSpecialplus11in.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                                   m_h_mdtSPZR_LargeSpecialplus11in.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);

                              }else{
                                   m_h_off_ptvsSA_resptLargeSpecialplus11out.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                                   m_h_SA_resptLargeSpecialplus11out.at(i)->Fill(resSA_pt);
                                   m_h_offetavsSA_resptLargeSpecialplus11out.at(i)->Fill(m_poff_eta,resSA_pt);
                                   m_h_mdtetavsSA_resptLargeSpecialplus11out.at(i)->Fill(ave_mdteta,resSA_pt);
                                   m_h_eSA_pt_LargeSpecialplus11out.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                                   Double_t buf_BIsegmentR = 0;
                                   Double_t buf_resR = 99999;
                                   for(Int_t index = 0;index < 10;index++){
                                        if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialplus11out.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                        if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialplus11out.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                                        if(m_probe_segment_x[index] <= -3000.0 && m_probe_segment_x[index] >= -6000.0 && m_probe_segment_y[index] <= -2000.0 && m_probe_segment_y[index] >= -6000.0 && sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)) <= 5500 && sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)) >= 4900){
                                             if(buf_resR >= sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)) - pSA_superpointR_BI){
                                                  buf_BIsegmentR = sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2));
                                                  buf_resR = sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)) - pSA_superpointR_BI;
                                             }
                                        }
                                   }
                                   if(buf_BIsegmentR != 0){
                                        m_h_segSP_diffR_LSBI.at(i)->Fill(buf_resR);
                                        m_h_segSP_resR_LSBI.at(i)->Fill(1.0/pSA_superpointR_BI - 1.0/buf_BIsegmentR);
                                        Int_t buf_index = fabs(static_cast<Int_t>(m_probe_segment_etaIndex[index]));
                                        switch(buf_index){
                                             case 1:
                                                  m_h_segSP_resR_LSBIetaindex1.at(i)->Fill(1.0/pSA_superpointR_BI - 1.0/buf_BIsegmentR);
                                                  break;
                                             case 2:
                                                  m_h_segSP_resR_LSBIetaindex2.at(i)->Fill(1.0/pSA_superpointR_BI - 1.0/buf_BIsegmentR);
                                                  break;
                                             case 3:
                                                  m_h_segSP_resR_LSBIetaindex3.at(i)->Fill(1.0/pSA_superpointR_BI - 1.0/buf_BIsegmentR);
                                                  break;
                                             case 4:
                                                  m_h_segSP_resR_LSBIetaindex4.at(i)->Fill(1.0/pSA_superpointR_BI - 1.0/buf_BIsegmentR);
                                                  break;
                                             case 5:
                                                  m_h_segSP_resR_LSBIetaindex5.at(i)->Fill(1.0/pSA_superpointR_BI - 1.0/buf_BIsegmentR);
                                                  break;
                                             case 6:
                                                  m_h_segSP_resR_LSBIetaindex6.at(i)->Fill(1.0/pSA_superpointR_BI - 1.0/buf_BIsegmentR);
                                                  break;
                                             default :
                                                  break;
                                        }
                                   }
                                   
                                   if(std::fabs(m_poff_pt*0.001) > 30.0 && std::fabs(m_poff_pt*0.001) < 50.0){
                                        m_h_highoffetavsSA_resptLargeSpecialplus11out.at(i)->Fill(m_poff_eta,resSA_pt);
                                   }
                                   m_h_mdtSPZR_LargeSpecialplus11out.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                                   m_h_mdtSPZR_LargeSpecialplus11out.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                                   m_h_mdtSPZR_LargeSpecialplus11out.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
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
                              if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialminus.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         }
                         m_h_mdtSPZR_LargeSpecialminus.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                         m_h_mdtSPZR_LargeSpecialminus.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                         m_h_mdtSPZR_LargeSpecialminus.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);

                         if(pSA_roiphi < -0.6 && pSA_roiphi > -1.0){//15
                              if(pSA_roiphi > -0.8){
                              m_h_off_ptvsSA_resptLargeSpecialminus15out.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              m_h_SA_resptLargeSpecialminus15out.at(i)->Fill(resSA_pt);
                              m_h_offetavsSA_resptLargeSpecialminus15out.at(i)->Fill(m_poff_eta,resSA_pt);
                              m_h_mdtetavsSA_resptLargeSpecialminus15out.at(i)->Fill(ave_mdteta,resSA_pt);
                              m_h_eSA_pt_LargeSpecialminus15out.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              for(Int_t index = 0;index < 10;index++){
                                  if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialminus15out.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                  if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialminus15out.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              }
                              if(std::fabs(m_poff_pt*0.001) > 30.0 && std::fabs(m_poff_pt*0.001) < 50.0){
                                   m_h_highoffetavsSA_resptLargeSpecialminus15out.at(i)->Fill(m_poff_eta,resSA_pt);
                              }
                              m_h_mdtSPZR_LargeSpecialminus15out.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                              m_h_mdtSPZR_LargeSpecialminus15out.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                              m_h_mdtSPZR_LargeSpecialminus15out.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                         }else{
                              m_h_off_ptvsSA_resptLargeSpecialminus15in.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              m_h_SA_resptLargeSpecialminus15in.at(i)->Fill(resSA_pt);
                              m_h_offetavsSA_resptLargeSpecialminus15in.at(i)->Fill(m_poff_eta,resSA_pt);
                              m_h_mdtetavsSA_resptLargeSpecialminus15in.at(i)->Fill(ave_mdteta,resSA_pt);
                              m_h_eSA_pt_LargeSpecialminus15in.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              for(Int_t index = 0;index < 10;index++){
                                   if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialminus15in.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                   if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialminus15in.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              }
                              if(std::fabs(m_poff_pt*0.001) > 30.0 && std::fabs(m_poff_pt*0.001) < 50.0){
                                   m_h_highoffetavsSA_resptLargeSpecialminus15in.at(i)->Fill(m_poff_eta,resSA_pt);
                              }
                              m_h_mdtSPZR_LargeSpecialminus15in.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                              m_h_mdtSPZR_LargeSpecialminus15in.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                              m_h_mdtSPZR_LargeSpecialminus15in.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
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
                                  if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialminus11in.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              }
                              if(std::fabs(m_poff_pt*0.001) > 30.0 && std::fabs(m_poff_pt*0.001) < 50.0){
                                   m_h_highoffetavsSA_resptLargeSpecialminus11in.at(i)->Fill(m_poff_eta,resSA_pt);
                              }
                              m_h_mdtSPZR_LargeSpecialminus11in.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                              m_h_mdtSPZR_LargeSpecialminus11in.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                              m_h_mdtSPZR_LargeSpecialminus11in.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
                         }else{
                              m_h_off_ptvsSA_resptLargeSpecialminus11out.at(i)->Fill(std::fabs(m_poff_pt*0.001),resSA_pt);
                              m_h_SA_resptLargeSpecialminus11out.at(i)->Fill(resSA_pt);
                              m_h_offetavsSA_resptLargeSpecialminus11out.at(i)->Fill(m_poff_eta,resSA_pt);
                              m_h_mdtetavsSA_resptLargeSpecialminus11out.at(i)->Fill(ave_mdteta,resSA_pt);
                              m_h_eSA_pt_LargeSpecialminus11out.at(i)->Fill(std::fabs(m_poff_pt*0.001));
                              for(Int_t index = 0;index < 10;index++){
                                  if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_resptLargeSpecialminus11out.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                                  if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0)m_h_segmentZR_LargeSpecialminus11out.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                              }
                              if(std::fabs(m_poff_pt*0.001) > 30.0 && std::fabs(m_poff_pt*0.001) < 50.0){
                                   m_h_highoffetavsSA_resptLargeSpecialminus11out.at(i)->Fill(m_poff_eta,resSA_pt);
                              }
                              m_h_mdtSPZR_LargeSpecialminus11out.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
                              m_h_mdtSPZR_LargeSpecialminus11out.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
                              m_h_mdtSPZR_LargeSpecialminus11out.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
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
               m_h_mdtSPXY_BI.at(i)->Fill(pSA_superpointR_BI*cos(pSA_roiphi),pSA_superpointR_BI*sin(pSA_roiphi));
               m_h_mdtSPXY_BM.at(i)->Fill(pSA_superpointR_BM*cos(pSA_roiphi),pSA_superpointR_BM*sin(pSA_roiphi));
               m_h_mdtSPXY_BO.at(i)->Fill(pSA_superpointR_BO*cos(pSA_roiphi),pSA_superpointR_BO*sin(pSA_roiphi));
               m_h_mdtSPXY_BME.at(i)->Fill(pSA_superpointR_BME*cos(pSA_roiphi),pSA_superpointR_BME*sin(pSA_roiphi));
               m_h_mdtSPZR.at(i)->Fill(pSA_superpointZ_BI,pSA_superpointR_BI);
               m_h_mdtSPZR.at(i)->Fill(pSA_superpointZ_BM,pSA_superpointR_BM);
               m_h_mdtSPZR.at(i)->Fill(pSA_superpointZ_BO,pSA_superpointR_BO);
               //if(pSA_superpointR_BI*cos(pSA_roiphi) < -4500 && 5200 pSA_superpointR_BI*cos(pSA_roiphi) pSA_superpointR_BI*sin(pSA_roiphi) pSA_superpointR_BI*sin(pSA_roiphi))m_h_etaIndexout->Fill(m_probe_segment_etaIndex);
               //if(pSA_superpointR_BI*cos(pSA_roiphi) pSA_superpointR_BI*cos(pSA_roiphi) pSA_superpointR_BI*sin(pSA_roiphi) pSA_superpointR_BI*sin(pSA_roiphi))m_h_etaIndexin->Fill(m_probe_segment_etaIndex);
               for(Int_t index = 0;index < 10;index++){
                    Int_t buf_index = static_cast<Int_t>(m_probe_segment_etaIndex[index]);
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
               for(Int_t size = 0;size < (signed int)pSA_rpcX->size();size++){
                    m_h_rpchitXY.at(i)->Fill(pSA_rpcX->at(size),pSA_rpcY->at(size));
                    if(m_g_rpchitXY.at(i)->GetN() <= 1000000)m_g_rpchitXY.at(i)->SetPoint(m_g_rpchitXY.at(i)->GetN(),pSA_rpcX->at(size),pSA_rpcY->at(size));
               }
               for(Int_t size = 0;size < (signed int)pSA_rpcZ->size();size++){
                    m_h_rpchitZR.at(i)->Fill(pSA_rpcZ->at(size),pSA_rpcR->at(size));
               }
               for(Int_t size = 0;size < (signed int)pSA_mdtZ->size();size++){
                    m_h_mdthitXY.at(i)->Fill(pSA_mdtR->at(size)*cos(pSA_mdtPhi->at(size)),pSA_mdtR->at(size)*sin(pSA_mdtPhi->at(size)));
                    m_h_mdthitZR.at(i)->Fill(pSA_mdtZ->at(size),pSA_mdtR->at(size));
                    if(m_g_mdthitXY.at(i)->GetN() <= 1000000)m_g_mdthitXY.at(i)->SetPoint(m_g_mdthitXY.at(i)->GetN(),pSA_mdtR->at(size)*cos(pSA_mdtPhi->at(size)),pSA_mdtR->at(size)*sin(pSA_mdtPhi->at(size)));
               }
               Int_t buf_numsegment = 0;
               for(Int_t index = 0;index < 10;index++){
                    if(m_probe_segment_etaIndex[index] >= -8.0 && m_probe_segment_etaIndex[index] <= 8.0)m_h_etaIndexvsSA_respt.at(i)->Fill(m_probe_segment_etaIndex[index],resSA_pt);
                    if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0)m_h_segmentXY.at(i)->Fill(m_probe_segment_x[index],m_probe_segment_y[index]);
                    if(m_probe_segment_x[index] != -77777.0 && m_probe_segment_y[index] != -77777.0 && m_probe_segment_z[index] != -77777.0){
                         m_h_segmentZR.at(i)->Fill(m_probe_segment_z[index],sqrt(pow(m_probe_segment_x[index],2) + pow(m_probe_segment_y[index],2)));
                         buf_numsegment++;
                    }
               }
               m_h_numsegment.at(i)->Fill(buf_numsegment);
          }
          m_h_offphivsSAphims.at(i)->Fill(m_poff_phi,pSA_phims);

          //CB
          if(!Cut_CB(pCB_pass))continue;
          Double_t textCB_dR = TMath::Sqrt(pow(m_tCB_eta - m_toff_exteta,2) + pow(m_tCB_phi - m_toff_extphi,2));
          pextCB_dR = TMath::Sqrt(pow(pCB_eta - m_poff_exteta,2) + pow(pCB_phi - m_poff_extphi,2));
          Double_t resCB_pt = std::fabs(m_poff_pt)/std::fabs(pCB_pt) - 1.0;
          m_h_pCB_pt.at(i)->Fill(std::fabs(pCB_pt*0.001));
          m_h_pCB_dR.at(i)->Fill(pCB_dR);
          m_h_textCB_dR.at(i)->Fill(textCB_dR);
          m_h_pextCB_dR.at(i)->Fill(pextCB_dR);
          m_h_eCB_pt.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          if(std::fabs(m_poff_pt*0.001) > 40){
          //if(std::fabs(m_poff_pt*0.001) > 8){
               m_h_eCB_eta.at(i)->Fill(m_poff_eta);
               m_h_eCB_phi.at(i)->Fill(m_poff_phi);
               m_h_eCB_aipc.at(i)->Fill(m_aipc);
          }
          m_h_pCB_respt.at(i)->Fill(resCB_pt);
          if(Dicision_barrel(m_poff_eta)){
               m_h_eCB_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          }else{
               m_h_eCB_pt_end.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          }

          //EF
          if(!Cut_EF(pEF_pass))continue;
          Double_t textEF_dR = TMath::Sqrt(pow(m_tEF_eta - m_toff_exteta,2) + pow(m_tEF_phi - m_toff_extphi,2));
          pextEF_dR = TMath::Sqrt(pow(pEF_eta - m_poff_exteta,2) + pow(pEF_phi - m_poff_extphi,2));
          Double_t resEF_pt = std::fabs(m_poff_pt)/std::fabs(pEF_pt) - 1.0;
          m_h_pEF_pt.at(i)->Fill(std::fabs(pEF_pt*0.001));
          m_h_pEF_dR.at(i)->Fill(pEF_dR);
          m_h_textEF_dR.at(i)->Fill(textEF_dR);
          m_h_pextEF_dR.at(i)->Fill(pextEF_dR);
          m_h_eEF_pt.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          if(std::fabs(m_poff_pt*0.001) > 40){
          //if(std::fabs(m_poff_pt*0.001) > 8){
               m_h_eEF_eta.at(i)->Fill(m_poff_eta);
               m_h_eEF_phi.at(i)->Fill(m_poff_phi);
               m_h_eEF_aipc.at(i)->Fill(m_aipc);
          }
          m_h_pEF_respt.at(i)->Fill(resEF_pt);
          if(Dicision_barrel(m_poff_eta)){
               m_h_eEF_pt_barrel.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          }else{
               m_h_eEF_pt_end.at(i)->Fill(std::fabs(m_poff_pt*0.001));
          }
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

     for(Int_t i = 0;i < m_nhist;i++){

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

          ceff.SetConditionName(Form("L1Efficiency_phi_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 Efficiency;offline phi;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyphi(m_h_eoff_phi.at(i),m_h_eL1_phi.at(i));
          ceff.SetConditionName(Form("SAEfficiency_phi_%dGeV",i*m_thpitch));
          ceff.SetCondition("L2MuonSA Efficiency;offline phi;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyphi(m_h_eL1_phi.at(i),m_h_eSA_phi.at(i));
          ceff.SetConditionName(Form("CBEfficiency_phi_%dGeV",i*m_thpitch));
          ceff.SetCondition("muComb Efficiency;offline phi;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyphi(m_h_eSA_phi.at(i),m_h_eCB_phi.at(i));
          ceff.SetConditionName(Form("EFEfficiency_phi_%dGeV",i*m_thpitch));
          ceff.SetCondition("EventFilter Efficiency;offline phi;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencyphi(m_h_eCB_phi.at(i),m_h_eEF_phi.at(i));

          ceff.SetConditionName(Form("L1Efficiency_pileup_%dGeV",i*m_thpitch));
          ceff.SetCondition("L1 Efficiency;pileup;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencypileup(m_h_eoff_aipc.at(i),m_h_eL1_aipc.at(i));
          ceff.SetConditionName(Form("SAEfficiency_pileup_%dGeV",i*m_thpitch));
          ceff.SetCondition("L2MuonSA Efficiency;pileup;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencypileup(m_h_eL1_aipc.at(i),m_h_eSA_aipc.at(i));
          ceff.SetConditionName(Form("CBEfficiency_pileup_%dGeV",i*m_thpitch));
          ceff.SetCondition("muComb Efficiency;pileup;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencypileup(m_h_eSA_aipc.at(i),m_h_eCB_aipc.at(i));
          ceff.SetConditionName(Form("EFEfficiency_pileup_%dGeV",i*m_thpitch));
          ceff.SetCondition("EventFilter Efficiency;pileup;Efficiency",1.0,0.1,0.1,0.105,0.165);
          ceff.DrawEfficiencypileup(m_h_eCB_aipc.at(i),m_h_eEF_aipc.at(i));

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

          m_g_rpchitXY.at(i)->SetName(Form("g_rpchitxy_%dGeV",i*m_thpitch));
          m_g_mdthitXY.at(i)->SetName(Form("g_mdthitxy_%dGeV",i*m_thpitch));

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
          m_h_rpchitZR.at(i)->Write();
          m_h_mdthitXY.at(i)->Write();
          m_g_rpchitXY.at(i)->Write();
          m_g_mdthitXY.at(i)->Write();
          m_h_mdthitZR.at(i)->Write();
          m_h_pSAphivspSAphims.at(i)->Write();
          m_h_pSAphivspSAphibe.at(i)->Write();
          m_h_mdtSPXY_BI.at(i)->Write();
          m_h_mdtSPXY_BM.at(i)->Write();
          m_h_mdtSPXY_BO.at(i)->Write();
          m_h_mdtSPXY_BME.at(i)->Write();
          m_h_avemdteta.at(i)->Write();
          //m_h_etaIndexout.at(i)->Write();
          //m_h_etaIndexin.at(i)->Write();
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
          m_h_segmentZR.at(i)->Write();
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
          m_h_numsegment.at(i)->Write();

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
          m_h_eoff_aipc.at(i)->Write();
          m_h_eL1_aipc.at(i)->Write();
          m_h_eSA_aipc.at(i)->Write();
          m_h_eCB_aipc.at(i)->Write();
          m_h_eEF_aipc.at(i)->Write();
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

          m_h_segSP_diffR_LSBI.at(i)->Write();
          m_h_segSP_resR_LSBI.at(i)->Write();
          m_h_segSP_diffR_LargeBI.at(i)->Write();
          m_h_segSP_resR_LargeBI.at(i)->Write();
          m_h_segSP_resR_LSBIetaindex1.at(i)->Write();
          m_h_segSP_resR_LargeBIetaindex1.at(i)->Write();
          m_h_segSP_resR_LSBIetaindex2.at(i)->Write();
          m_h_segSP_resR_LargeBIetaindex2.at(i)->Write();
          m_h_segSP_resR_LSBIetaindex3.at(i)->Write();
          m_h_segSP_resR_LargeBIetaindex3.at(i)->Write();
          m_h_segSP_resR_LSBIetaindex4.at(i)->Write();
          m_h_segSP_resR_LargeBIetaindex4.at(i)->Write();
          m_h_segSP_resR_LSBIetaindex5.at(i)->Write();
          m_h_segSP_resR_LargeBIetaindex5.at(i)->Write();
          m_h_segSP_resR_LSBIetaindex6.at(i)->Write();
          m_h_segSP_resR_LargeBIetaindex6.at(i)->Write();

          if(m_countLarge.size() != 0 && m_countLargeSpecial.size() != 0 && m_countSmall.size() != 0 && m_countSmallSpecial.size() != 0)cout<<i*m_thpitch<<"      "<<m_countLarge.at(i)<<"      "<<m_countLargeSpecial.at(i)<<"      "<<m_countSmall.at(i)<<"      "<<m_countSmallSpecial.at(i)<<endl;
     }

}
