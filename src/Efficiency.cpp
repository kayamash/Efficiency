#include "../Efficiency/Efficiency.chh"
#include "/home/kayamash/working/PtNewMethod/src/kayamashForLUT.cpp"
#include "CalcEff.cpp"
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
#include <TProfile.h>

const Double_t pt_threshold[4] = {3.38,1.25,3.17,3.41};//MU4
//const Double_t pt_threshold[4] = {5.17,3.25,4.69,5.14};//MU6
//const Double_t pt_threshold[4] = {15.87,10.73,12.21,15.87};//MU20
const Double_t LargeRegion[9] = {0.25,0.55,1.025,1.325,1.80,2.125,2.575,2.90,TMath::Pi()};
//const Double_t correctionAlpha[3] = {0,0.0283101,-0.00993997};//1SP 2SP 3SP
//const Double_t correctionBeta[3] = {0,0.024428,-0.00257345};//1SP 2SP 3SP
const Double_t correctionAlpha[3] = {0.,0.,0.};//1SP 2SP 3SP
const Double_t correctionBeta[3] = {0.,0.,0.};//1SP 2SP 3SP

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

     kayamashForLUT LUT(0.,0.);
     //LUT.ReadLUT("NewMethodAlphaJPZ.LUT",m_LUTAlphaSectorChargeEtaPhi);
     //LUT.ReadLUT("NewMethodBetaJPZ.LUT",m_LUTBetaSectorChargeEtaPhi);
     LUT.ReadLUT("/gpfs/fs7001/kayamash/Mywork/LUT/test/NewMethodAlphaJPZ.LUT",m_LUTAlphaSectorChargeEtaPhi);
     LUT.ReadLUT("/gpfs/fs7001/kayamash/Mywork/LUT/test/NewMethodBetaJPZ.LUT",m_LUTBetaSectorChargeEtaPhi);
     for(Int_t sector = 0; sector < 5; ++sector){
          for(Int_t charge = 0; charge < 2; ++charge){
               for(Int_t eta = 0; eta < 30; ++eta){
                    for(Int_t phi = 0; phi < 30; ++phi){
                         m_h_lutAvsBAlpha->Fill(m_LUTAlphaSectorChargeEtaPhi[sector][charge][eta][phi][0],m_LUTAlphaSectorChargeEtaPhi[sector][charge][eta][phi][1]);
                         m_h_lutAvsBBeta->Fill(m_LUTBetaSectorChargeEtaPhi[sector][charge][eta][phi][0],m_LUTBetaSectorChargeEtaPhi[sector][charge][eta][phi][1]);
                    }
               }
          }
     }
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

bool Efficiency::PlateauCut(Double_t pt){
     Double_t cut_pt = 40.0;
     if(m_method_name == "mu4" || m_method_name == "mu6")cut_pt = 8.0;
     if(pt > cut_pt){
          return kTRUE;
     }else{
          return kFALSE;
     }
}

bool Efficiency::CutSAMyLUT(Double_t pt,Int_t L1Number,Int_t L1Sector,Int_t SANumber,Int_t SASector){
     if(L1Number == SANumber && L1Sector == SASector && std::fabs(pt) >= pt_threshold[0])return kTRUE;
     //if(std::fabs(pt) >= pt_threshold[0])return kTRUE;
     return kFALSE;
}

void Efficiency::CalcPtByAlpha(Double_t parA,Double_t parB,Double_t angle,Double_t charge,Int_t nSP,Double_t &pt){
     Double_t discriminant = 0;
     Double_t sqrtDiscriminant = 0;
     Double_t Angle = 0;
     Double_t correction = 0;
     if(nSP == 2)correction = correctionAlpha[1];
     if(nSP == 3)correction = correctionAlpha[2];
     if(charge == 1.){
          Angle = - std::fabs(angle);
     }else{
          Angle = std::fabs(angle);
     }
     discriminant = (parA*parA) + 4*parB*(correction + Angle);
     if(discriminant > 0){
          sqrtDiscriminant = std::sqrt(discriminant);
     }
     if(charge == 1.){
          if(-parA - sqrtDiscriminant == 0){
               pt = 0;
          }else{
               pt = 2*parB/(-parA - sqrtDiscriminant);
          }
     }else{
          if(-parA + sqrtDiscriminant == 0){
               pt = 0;
          }else{
               pt = 2*parB/(-parA + sqrtDiscriminant);
          }
     }
}

void Efficiency::CalcPtByBeta(Double_t parA,Double_t parB,Double_t angle,Double_t charge,Int_t nSP,Double_t &pt){
     Double_t discriminant = 0;
     Double_t sqrtDiscriminant = 0;
     Double_t Angle = 0;
     Double_t correction = 0;
     if(nSP == 2)correction = correctionBeta[1];
     if(nSP == 3)correction = correctionBeta[2];
     if(charge == 1.){
          Angle = std::fabs(angle);
     }else{
          Angle = - std::fabs(angle);
     }
     discriminant = (parA*parA) + 4*parB*(correction + Angle);
     if(discriminant > 0){
          sqrtDiscriminant = std::sqrt(discriminant);
     }
     if(charge == 1.){
          if(-parA - sqrtDiscriminant == 0){
               pt = 0;
          }else{
               pt = 2*parB/(-parA + sqrtDiscriminant);
          }
     }else{
          if(-parA + sqrtDiscriminant == 0){
               pt = 0;
          }else{
               pt = 2*parB/(-parA - sqrtDiscriminant);
          }
     }
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

int Efficiency::getLSSector(Double_t roiphi,Int_t address){
     if(address != 1)goto exit;
     if(roiphi >= -2.6 && roiphi < -2.4){
          return 0;//15 outer
     }else if(roiphi >= -2.4 && roiphi < -2.0){
          return 1;//15 inner
     }else if(roiphi >= -1.0 && roiphi < -0.8){
          return 2;//11 inner 
     }else if(roiphi >= -0.8 && roiphi < -0.6){
          return 3;//11outer
     }
     exit:
          return -1;

     return -1;
}

bool Efficiency::getBarrelMuFastRes(Double_t ptGeV,Double_t &res){
     Double_t par[6] = {0.042, -0.00046, 3.5, -1.8, 0.35, -0.017};
     Double_t pt = ptGeV*1000.;
     Double_t AbsPtInv = std::fabs(1./pt);
     if(pt == 0){
          res = 1.0e33;
          return kTRUE;
     }else{
          if(AbsPtInv < 0.000186){
               res = par[0]*AbsPtInv + par[1]/1000.;
               return kTRUE;
          }else{
               res = par[2]*pow(AbsPtInv,3)/(1000.*1000.) + par[3]*pow(AbsPtInv,2)/(1000.) + par[4]*AbsPtInv + par[5]/1000.;
               return kTRUE;
          }
     }
     return kFALSE;
}

bool getBarrelIDSCANRes(Double_t ptMeV,Double_t &res){
     Double_t par[2] = {0.017, 0.000000418};
     if(pt == 0){
          res = 1.0e33;
          return kTRUE;
     }else{
          Double_t AbsPtInv = std::fabs(1./pt);
          res = par[0]*AbsPtInv + par[1]/1000.;
          return kTRUE;
     }
     return kFALSE;
}

void Efficiency::Execute(Int_t ev){
     tChain->GetEntry(ev);
     Double_t pL1_pt = -99999;
     Double_t pL1_eta = 0;
     Double_t pL1_phi = 0;
     Double_t pL1_dR = 1;
     Int_t pL1_pass = 0;
     Int_t pL1_roiNumber = -1;
     Int_t pL1_roiSector = -1;
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
     Int_t pSA_roiNumber = -1;
     Int_t pSA_roiSector = -1;
     vector<float> *pSA_rpcX = 0;
     vector<float> *pSA_rpcY = 0;
     vector<float> *pSA_rpcZ = 0;
     vector<float> *pSA_rpcR = 0;
     vector<float> *pSA_mdtZ = 0;
     vector<float> *pSA_mdtR = 0;
     vector<float> *pSA_mdtPhi = 0;
     vector<Int_t> *pSA_mdthitChamber = 0;
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
     Double_t pSA_superpointSlope_BI = 0;
     Double_t pSA_superpointSlope_BM = 0;
     Double_t pSA_superpointSlope_BME = 0;
     Double_t pSA_superpointSlope_BEE = 0;
     Double_t pCB_pt = -99999;
     Double_t pCB_eta = 0;
     Double_t pCB_phi = 0;
     Double_t pCB_dR = 1;
     Int_t pCB_pass = 0;
     Int_t pEFTAG_pass = -1;

     for(Int_t method = 0;method < 25;method++){
          if(m_mes_name->at(method) == m_method_name){
               pL1_pt = m_pL1_pt->at(method);
               pSA_pt = m_pSA_pt->at(method);
               pCB_pt = m_pCB_pt->at(method);
               pL1_eta = m_pL1_eta->at(method);
               pSA_eta = m_pSA_eta->at(method);
               pCB_eta = m_pCB_eta->at(method);
               pL1_phi = m_pL1_phi->at(method);
               pSA_phi = m_pSA_phi->at(method);
               pCB_phi = m_pCB_phi->at(method);
               pL1_pass = m_pL1_pass->at(method);
               pSA_pass = m_pSA_pass->at(method);
               pCB_pass = m_pCB_pass->at(method);
               pL1_dR = m_pL1_dR->at(method);
               pSA_dR = m_pSA_dR->at(method);
               pCB_dR = m_pCB_dR->at(method);
               pEFTAG_pass = m_pEFTAG_pass->at(method);
               pL1_roiNumber = m_pL1_roiNumber->at(method);
               pL1_roiSector = m_pL1_roiSector->at(method);
               pSA_sAddress = m_pSA_sAddress->at(method);
               pSA_phims = m_pSA_phims->at(method);
               pSA_phibe = m_pSA_phibe->at(method);
               pSA_roieta = m_pSA_roieta->at(method);
               pSA_roiphi = m_pSA_roiphi->at(method);
               pSA_roiNumber = m_pSA_roiNumber->at(method);
               pSA_roiSector = m_pSA_roiSector->at(method);
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
               pSA_superpointSlope_BI = m_pSA_superpointSlope_BI->at(method);
               pSA_superpointSlope_BM = m_pSA_superpointSlope_BM->at(method);
               pSA_superpointSlope_BME = m_pSA_superpointSlope_BME->at(method);
               pSA_superpointSlope_BEE = m_pSA_superpointSlope_BEE->at(method);
          }
     }

     Int_t numBarrelSP = 0;
     Int_t numSP = 0;
     Double_t SPZInner = 0;
     Double_t SPRInner = 0;
     Double_t SPSlopeInner = 0;
     Double_t SPZMiddle = 0;
     Double_t SPRMiddle = 0;
     Double_t SPSlopeMiddle = 0;
     Double_t SPROuter = 0;
     Double_t SPZOuter = 0;
     Double_t SPSlopeOuter = 0;
     Double_t phiInteg = 0;
     Double_t barrelalpha = -99999;
     Double_t barrelbeta = -99999;
     Double_t truthBeta = -99999;
     Double_t segmentBISlope = 0;
     Double_t segmentBMSlope = 0;
     Double_t AlphaPt = 0;
     Double_t BetaPt = 0;
     Double_t deltaTheta = 0;
     Double_t newMethodSAPt = 0;
     Double_t NewMethodResPt = 0;

     if(std::fabs(m_toff_pt)*0.001 < 10.0)m_reqL1dR = -0.00001*std::fabs(m_toff_pt) + 0.18;
     if(pSA_superpointR_BI != 0){
          numBarrelSP++;
     }
     if(pSA_superpointR_BM != 0){
          numBarrelSP++;
     }
     if(pSA_superpointR_BO != 0){
          numBarrelSP++;
     }

     if(pSA_superpointR_BI > 0)numSP++;
     if(pSA_superpointR_BM > 0)numSP++;
     if(pSA_superpointR_BO > 0)numSP++;
     if(pSA_superpointR_EI > 0)numSP++;
     if(pSA_superpointR_EM > 0)numSP++;
     if(pSA_superpointR_EO > 0)numSP++;
     if(pSA_superpointR_CSC > 0)numSP++;
     if(pSA_superpointR_BME > 0)numSP++;
     if(pSA_superpointR_BEE > 0)numSP++;
     if(pSA_superpointR_EE > 0)numSP++;

          //offline
     if(!CutTagProbe(pEFTAG_pass))return;
     m_h_pOffPt->Fill(m_poff_pt*0.001);
     m_h_eOffPt->Fill(std::fabs(m_poff_pt*0.001));
     if(PlateauCut(std::fabs(m_poff_pt*0.001))){//Plateau Efficiency
          m_h_eOffEta->Fill(m_poff_eta);
          m_h_eOffPhi->Fill(m_poff_phi);
     }
     switch(EtaDistribution()){
          case 0://Barrel
          m_h_eOffPtBarrel->Fill(std::fabs(m_poff_pt*0.001));
          break;
          case 1://Transition
          m_h_eOffPtTransition->Fill(std::fabs(m_poff_pt*0.001));
          break;
          case 2://Endcap
          m_h_eOffPtEndcap->Fill(std::fabs(m_poff_pt*0.001));
          break;
          case 3://Forward
          m_h_eOffPtForward->Fill(std::fabs(m_poff_pt*0.001));
          break;
          default:
          break;
     }

          //L1
     if(!CutL1(pL1_pass))return;
     m_h_pL1Pt->Fill(std::fabs(pL1_pt*0.001));
     m_h_eL1Pt->Fill(std::fabs(m_poff_pt*0.001));
     switch(EtaDistribution(pSA_roieta)){
          case 0://Barrel
          m_h_eL1PtBarrel->Fill(std::fabs(m_poff_pt*0.001));
          m_h_eL1PtBarrelSP[numSP]->Fill(std::fabs(m_poff_pt*0.001));
          m_h_eL1PtBarrelMyLUTAlpha->Fill(std::fabs(m_poff_pt*0.001));
          break;
          case 1://Transition
          m_h_eL1PtTransition->Fill(std::fabs(m_poff_pt*0.001));
          m_h_eL1PtTransitionSP[numSP]->Fill(std::fabs(m_poff_pt*0.001));
          break;
          case 2://Endcap
          m_h_eL1PtEndcap->Fill(std::fabs(m_poff_pt*0.001));
          m_h_eL1PtEndcapSP[numSP]->Fill(std::fabs(m_poff_pt*0.001));
          break;
          case 3://Forward
          m_h_eL1PtForward->Fill(std::fabs(m_poff_pt*0.001));
          m_h_eL1PtForwardSP[numSP]->Fill(std::fabs(m_poff_pt*0.001));
          break;
          default:
          break;
     }
     if(PlateauCut(std::fabs(m_poff_pt*0.001))){//Plateau Efficiency
          m_h_eL1Eta->Fill(m_poff_eta);
          m_h_eL1Phi->Fill(m_poff_phi);
     }

     kayamashForLUT LUT(0.,0.);
     if(pSA_superpointR_BI != 0){
          SPZInner = pSA_superpointZ_BI;
          SPRInner = pSA_superpointR_BI;
          SPSlopeInner = pSA_superpointSlope_BI;
     }
     if(pSA_superpointR_BEE != 0){
          SPZInner = pSA_superpointZ_BEE;
          SPRInner = pSA_superpointR_BEE;
          SPSlopeInner = pSA_superpointSlope_BEE;
     }
     if(pSA_superpointR_BM != 0){
          SPZMiddle = pSA_superpointZ_BM;
          SPRMiddle = pSA_superpointR_BM;
          SPSlopeMiddle = pSA_superpointSlope_BM;
     }
     if(pSA_superpointR_BME != 0){
          SPZMiddle = pSA_superpointZ_BME;
          SPRMiddle = pSA_superpointR_BME;
          SPSlopeMiddle = pSA_superpointSlope_BME;
     }
     if(pSA_superpointR_BO != 0)SPROuter = pSA_superpointR_BO;

     if(SPRMiddle != 0)barrelalpha = atan(SPZMiddle/SPRMiddle) - atan(SPSlopeMiddle);//Reciprocal number?;
     if(SPRInner != 0 && SPRMiddle != 0)barrelbeta = atan(1.0/SPSlopeInner) - atan(1.0/SPSlopeMiddle);//Reciprocal number?
     bool BIsegmentcheck = kFALSE;
     bool BMsegmentcheck = kFALSE;
     for(Int_t segmentNumber = 0; segmentNumber < 10; ++segmentNumber){
          Double_t tmp_segmentR = sqrt(m_probe_segment_x[segmentNumber]*m_probe_segment_x[segmentNumber] + m_probe_segment_y[segmentNumber]*m_probe_segment_y[segmentNumber]);
          Double_t tmp_segmentPR = sqrt(m_probe_segment_px[segmentNumber]*m_probe_segment_px[segmentNumber] + m_probe_segment_py[segmentNumber]*m_probe_segment_py[segmentNumber]);
          if(tmp_segmentR > 4000 && tmp_segmentR < 6500 && BIsegmentcheck == kFALSE){
               BIsegmentcheck = kTRUE;
               segmentBISlope = tmp_segmentPR/m_probe_segment_pz[segmentNumber];
          }
          if(tmp_segmentR > 6000 && tmp_segmentR < 9000 && BMsegmentcheck == kFALSE){
               BMsegmentcheck = kTRUE;
               segmentBMSlope = tmp_segmentPR/m_probe_segment_pz[segmentNumber];
          }
     }
     if(BIsegmentcheck)truthBeta = atan(segmentBISlope) - atan(segmentBMSlope);
     if(BIsegmentcheck && BMsegmentcheck && EtaDistribution(pSA_roieta) == 0){
          m_h_BarrelBetavsTruth->Fill(barrelbeta,truthBeta);
     }
     Int_t tmp_LUTpar[5] = {0,0,0,0,0};
     bool LUTcheck = LUT.getLUTparameter(pSA_sAddress,m_poff_charge,pSA_eta,pSA_phi,tmp_LUTpar,phiInteg);
     Int_t LUTparameter[4] = {tmp_LUTpar[0],tmp_LUTpar[1],tmp_LUTpar[2],tmp_LUTpar[3]};
     if(barrelalpha != -99999)m_h_ChargevsAlpha->Fill(LUTparameter[1],barrelalpha);
     if(barrelbeta != -99999)m_h_ChargevsBeta->Fill(LUTparameter[1],barrelbeta);
     if(LUTcheck && barrelalpha != -99999){
          Double_t parA = m_LUTAlphaSectorChargeEtaPhi[LUTparameter[0]][LUTparameter[1]][LUTparameter[2]][LUTparameter[3]][0];
          Double_t parB = m_LUTAlphaSectorChargeEtaPhi[LUTparameter[0]][LUTparameter[1]][LUTparameter[2]][LUTparameter[3]][1];
          CalcPtByAlpha(parA,parB,barrelalpha,m_poff_charge,numBarrelSP,AlphaPt);
     }
     if(LUTcheck && barrelbeta != -99999){
          Double_t parA = m_LUTBetaSectorChargeEtaPhi[LUTparameter[0]][LUTparameter[1]][LUTparameter[2]][LUTparameter[3]][0];
          Double_t parB = m_LUTBetaSectorChargeEtaPhi[LUTparameter[0]][LUTparameter[1]][LUTparameter[2]][LUTparameter[3]][1];
          CalcPtByBeta(parA,parB,barrelbeta,m_poff_charge,numBarrelSP,BetaPt);
     }

     if(SPRMiddle != 0 && EtaDistribution(pSA_roieta) == 0){//barrel alpha
          m_h_BarrelAlpha->Fill(barrelalpha);
          m_h_PtvsBarrelAlpha->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelalpha);
          m_h_pSAResPtBarrelAlpha->Fill(std::fabs(m_poff_pt*0.001)/std::fabs(AlphaPt) - 1.0);
          m_h_pSAResPtBarrelAlphaSector[LUTparameter[0]]->Fill(std::fabs(m_poff_pt*0.001)/std::fabs(AlphaPt) - 1.0);
     }//barrel alpha end

     if(SPRInner != 0 && SPRMiddle != 0 && EtaDistribution(pSA_roieta) == 0){//barrel beta
          m_h_BarrelBeta->Fill(barrelbeta);
          m_h_PtvsBarrelBeta->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelbeta);
          m_h_eL1PtBarrelMyLUTBeta->Fill(std::fabs(m_poff_pt*0.001));
          m_h_pSAResPtBarrelBeta->Fill(std::fabs(m_poff_pt*0.001)/std::fabs(BetaPt) - 1.0);
          m_h_pSAResPtBarrelBetaSector[LUTparameter[0]]->Fill(std::fabs(m_poff_pt*0.001)/std::fabs(BetaPt) - 1.0);
     }//barrel beta end

     if(SPRMiddle != 0 && EtaDistribution(pSA_roieta) == 0 && AlphaPt != 0){
          if(numBarrelSP == 1)m_h_pSAResPtBarrelAlpha1SP->Fill(std::fabs(m_poff_pt*0.001)/std::fabs(AlphaPt) - 1.0);
          if(numBarrelSP == 2)m_h_pSAResPtBarrelAlpha2SP->Fill(std::fabs(m_poff_pt*0.001)/std::fabs(AlphaPt) - 1.0);
          if(numBarrelSP == 3)m_h_pSAResPtBarrelAlpha3SP->Fill(std::fabs(m_poff_pt*0.001)/std::fabs(AlphaPt) - 1.0);
     }
     if(SPRInner != 0 && SPRMiddle != 0 && EtaDistribution(pSA_roieta) == 0 && BetaPt != 0){
          if(numBarrelSP == 2)m_h_pSAResPtBarrelBeta2SP->Fill(std::fabs(m_poff_pt*0.001)/std::fabs(BetaPt) - 1.0);
          if(numBarrelSP == 3)m_h_pSAResPtBarrelBeta3SP->Fill(std::fabs(m_poff_pt*0.001)/std::fabs(BetaPt) - 1.0);
     }

     deltaTheta = atan(SPRInner/SPZInner) - atan(1.0/SPSlopeInner);
     if(EtaDistribution(pSA_roieta) == 0){
          if(SPRMiddle != 0 && numBarrelSP <= 2){
               if(CutSAMyLUT(AlphaPt,pL1_roiNumber,pL1_roiSector,pSA_roiNumber,pSA_roiSector))m_h_eSAPtBarrelMyLUTAlpha->Fill(std::fabs(m_poff_pt*0.001));
               newMethodSAPt = AlphaPt;
          }else{
               if(CutSA(pSA_pass))m_h_eSAPtBarrelMyLUTAlpha->Fill(std::fabs(m_poff_pt*0.001));
               newMethodSAPt = pSA_pt;
          }
          if(SPRInner != 0 && SPRMiddle != 0 && SPROuter == 0){
               if(CutSAMyLUT(BetaPt,pL1_roiNumber,pL1_roiSector,pSA_roiNumber,pSA_roiSector))m_h_eSAPtBarrelMyLUTBeta->Fill(std::fabs(m_poff_pt*0.001));
               newMethodSAPt = BetaPt;
          }else{
               if(CutSA(pSA_pass))m_h_eSAPtBarrelMyLUTBeta->Fill(std::fabs(m_poff_pt*0.001));
               newMethodSAPt = pSA_pt;
          }

          if(numBarrelSP == 3){
               if(CutSA(pSA_pass))m_h_eSAPtBarrelkayamashMethod->Fill(std::fabs(m_poff_pt*0.001));
          }else if(SPRInner != 0 && SPRMiddle != 0 && SPROuter == 0 && std::fabs(deltaTheta) <= 0.02){
               if(CutSAMyLUT(BetaPt,pL1_roiNumber,pL1_roiSector,pSA_roiNumber,pSA_roiSector))m_h_eSAPtBarrelkayamashMethod->Fill(std::fabs(m_poff_pt*0.001));
          }else if(SPRMiddle != 0){
               if(CutSAMyLUT(AlphaPt,pL1_roiNumber,pL1_roiSector,pSA_roiNumber,pSA_roiSector))m_h_eSAPtBarrelkayamashMethod->Fill(std::fabs(m_poff_pt*0.001));
          }else{
               if(CutSA(pSA_pass))m_h_eSAPtBarrelkayamashMethod->Fill(std::fabs(m_poff_pt*0.001));
          }
     }
     NewMethodResPt = std::fabs(m_poff_pt*0.001)/std::fabs(newMethodSAPt) - 1.0;

     Double_t resSA_pt = std::fabs(m_poff_pt*0.001)/std::fabs(pSA_pt) - 1.0;
     if(SPRInner != 0 && EtaDistribution(pSA_roieta) == 0){
          m_h_pSAResPtvsDeltaTheta->Fill(deltaTheta,std::fabs(m_poff_pt*0.001)/std::fabs(BetaPt) - 1.0);
          if(barrelbeta != -99999)m_h_pSAResPtvsDeltaThetaForProf->Fill(std::fabs(deltaTheta),std::fabs(std::fabs(m_poff_pt*0.001)/std::fabs(BetaPt) - 1.0));
          if(barrelalpha != -99999)m_h_pAlphaResPtvsDeltaThetaForProf->Fill(std::fabs(deltaTheta),std::fabs(std::fabs(m_poff_pt*0.001)/std::fabs(AlphaPt) - 1.0));
          m_h_pOffPtvsDeltaTheta->Fill(std::fabs(m_poff_pt*0.001),deltaTheta);
          m_h_pOffPtvsDeltaThetaEtaIndex[static_cast<Int_t>(std::fabs(pSA_roieta)*6./1.05)]->Fill(std::fabs(m_poff_pt*0.001),deltaTheta);
          m_h_pOffPtvsDeltaThetaForProf->Fill(std::fabs(m_poff_pt*0.001),std::fabs(deltaTheta));
          if(deltaTheta < 0.05){
               m_h_pASResPtBarrelBetaSmallDeltaTheta->Fill(std::fabs(m_poff_pt*0.001)/std::fabs(BetaPt) - 1.0);
          }else{
               m_h_pASResPtBarrelBetaLargeDeltaTheta->Fill(std::fabs(m_poff_pt*0.001)/std::fabs(BetaPt) - 1.0);
          }
          for(Int_t seg = 0; seg < 10; ++seg){
               if(m_probe_segment_chamberIndex[seg] == 0 || m_probe_segment_chamberIndex[seg] == 1)m_h_DeltaThetavsHits->Fill(deltaTheta,m_probe_segment_nPrecisionHits[seg]);
          }
     }

     if(pSA_sAddress == 1)m_h_pSAResPtLargeSpecial->Fill(NewMethodResPt);
     if(getLSSector(pSA_roiphi,pSA_sAddress) >= 0 && SPRMiddle != 0)m_h_pSAResPtLS[getLSSector(pSA_roiphi,pSA_sAddress)]->Fill(std::fabs(m_poff_pt*0.001)/std::fabs(AlphaPt) - 1.0);
     if(getLSSector(pSA_roiphi,pSA_sAddress) >= 0 && SPRMiddle != 0)m_h_pSAResPtLSRadius[getLSSector(pSA_roiphi,pSA_sAddress)]->Fill(resSA_pt);

          //SA
     if(!CutSA(pSA_pass))return;
     m_h_pSAPt->Fill(std::fabs(pSA_pt));
     m_h_eSAPt->Fill(std::fabs(m_poff_pt*0.001));
     m_h_pSAResPt->Fill(resSA_pt);
     switch(EtaDistribution(pSA_roieta)){
          case 0://Barrel
          m_h_eSAPtBarrel->Fill(std::fabs(m_poff_pt*0.001));
          m_h_eSAPtBarrelSP[numBarrelSP]->Fill(std::fabs(m_poff_pt*0.001));
          if(numBarrelSP == 1)m_h_pSAResPtBarrel1SP->Fill(resSA_pt);
          if(numBarrelSP == 2)m_h_pSAResPtBarrel2SP->Fill(resSA_pt);
          if(numBarrelSP == 3)m_h_pSAResPtBarrel3SP->Fill(resSA_pt);
          break;
          case 1://Transition
          m_h_eSAPtTransition->Fill(std::fabs(m_poff_pt*0.001));
          m_h_eSAPtTransitionSP[numSP]->Fill(std::fabs(m_poff_pt*0.001));
          break;
          case 2://Endcap
          m_h_eSAPtEndcap->Fill(std::fabs(m_poff_pt*0.001));
          m_h_eSAPtEndcapSP[numSP]->Fill(std::fabs(m_poff_pt*0.001));
          break;
          case 3://Forward
          m_h_eSAPtForward->Fill(std::fabs(m_poff_pt*0.001));
          m_h_eSAPtForwardSP[numSP]->Fill(std::fabs(m_poff_pt*0.001));
          break;
     }
     if(PlateauCut(std::fabs(m_poff_pt*0.001))){//Plateau Efficiency
          m_h_eSAEta->Fill(m_poff_eta);
          m_h_eSAPhi->Fill(m_poff_phi);
     }

     //CB
     if(!CutCB(pCB_pass))return;
     switch(EtaDistribution(pSA_roieta)){
          case 0:
          m_h_eCBPtBarrel->Fill(std::fabs(m_poff_pt*0.001));
          break;
          default:
          break;
     }

}//Execute

void Efficiency::Finalize(TFile *tf1){
     CalcEff ceff;
     tf1->cd();
     //SetCondition
     //title,file title,yoffset,top margin,bottom margin,left margin,right margin

          //base,target
     cout<<"efficiency start!"<<endl;

     ceff.SetCondition("L1Efficiency","L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eOffPt,m_h_eL1Pt,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiency","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1Pt,m_h_eSAPt,m_binmax,200,m_efficiency_xerr);

     ceff.SetCondition("L1EfficiencyEta","L1 Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyeta(m_h_eOffEta,m_h_eL1Eta);
     ceff.SetCondition("SAEfficiencyEta","L2MuonSA Efficiency;offline eta;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyeta(m_h_eL1Eta,m_h_eSAEta);
     ceff.SetCondition("L1EfficiencyPhi","L1 Efficiency;offline phi;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyphi(m_h_eOffPhi,m_h_eL1Phi);
     ceff.SetCondition("SAEfficiencyPhi","L2MuonSA Efficiency;offline phi;Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiencyphi(m_h_eL1Phi,m_h_eSAPhi);

     ceff.SetCondition("L1EfficiencyBarrel","L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eOffPtBarrel,m_h_eL1PtBarrel,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyBarrel","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrel,m_h_eSAPtBarrel,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("CBEfficiencyBarrel","L2MuonComb Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eSAPtBarrel,m_h_eCBPtBarrel,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyBarrel0SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrelSP[0],m_h_eSAPtBarrelSP[0],m_binmax,200,m_efficiency_xerr); 
     ceff.SetCondition("SAEfficiencyBarrel1SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrelSP[1],m_h_eSAPtBarrelSP[1],m_binmax,200,m_efficiency_xerr); 
     ceff.SetCondition("SAEfficiencyBarrel2SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrelSP[2],m_h_eSAPtBarrelSP[2],m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyBarrel3SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrelSP[3],m_h_eSAPtBarrelSP[3],m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("L1EfficiencyTransition","L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eOffPtTransition,m_h_eL1PtTransition,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyTransition","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtTransition,m_h_eSAPtTransition,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyTransition0SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtTransitionSP[0],m_h_eSAPtTransitionSP[0],m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyTransition1SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtTransitionSP[1],m_h_eSAPtTransitionSP[1],m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyTransition2SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtTransitionSP[2],m_h_eSAPtTransitionSP[2],m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyTransition3SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtTransitionSP[3],m_h_eSAPtTransitionSP[3],m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("L1EfficiencyEnd","L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eOffPtEndcap,m_h_eL1PtEndcap,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyEnd","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtEndcap,m_h_eSAPtEndcap,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyEndcap0SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtEndcapSP[0],m_h_eSAPtEndcapSP[0],m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyEndcap1SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtEndcapSP[1],m_h_eSAPtEndcapSP[1],m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyEndcap2SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtEndcapSP[2],m_h_eSAPtEndcapSP[2],m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyEndcap3SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtEndcapSP[3],m_h_eSAPtEndcapSP[3],m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("L1EfficiencyForward","L1 Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eOffPtForward,m_h_eL1PtForward,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyForward","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtForward,m_h_eSAPtForward,m_binmax,300,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyForward0SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtForwardSP[0],m_h_eSAPtForwardSP[0],m_binmax,200,m_efficiency_xerr); 
     ceff.SetCondition("SAEfficiencyForward1SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtForwardSP[1],m_h_eSAPtForwardSP[1],m_binmax,200,m_efficiency_xerr); 
     ceff.SetCondition("SAEfficiencyForward2SP","L2MP Fficienuoncy;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtForwardSP[2],m_h_eSAPtForwardSP[2],m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyForward3SP","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtForwardSP[3],m_h_eSAPtForwardSP[3],m_binmax,200,m_efficiency_xerr);

     ceff.SetCondition("SAEfficiencyBarrelMyLUTAlpha","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrel,m_h_eSAPtBarrelMyLUTAlpha,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyBarrelMyLUTBeta","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrel,m_h_eSAPtBarrelMyLUTBeta,m_binmax,200,m_efficiency_xerr);
     ceff.SetCondition("SAEfficiencyBarrelkayamashMethod","L2MuonSA Efficiency;offline pt[GeV];Efficiency",1.0,0.1,0.1,0.105,0.165);
     ceff.DrawEfficiency(m_h_eL1PtBarrel,m_h_eSAPtBarrelkayamashMethod,m_binmax,200,m_efficiency_xerr);
     cout<<"efficiency end"<<endl;

     m_h_pOffPt->Write();
     m_h_pL1Pt->Write();
     m_h_pSAPt->Write();

     m_h_BarrelAlpha->Write();
     m_h_BarrelBeta->Write();
     m_h_PtvsBarrelAlpha->Write();
     m_h_PtvsBarrelBeta->Write();
     m_h_ChargevsAlpha->Write();
     m_h_ChargevsBeta->Write();
     m_h_lutAvsBAlpha->Write();
     m_h_lutAvsBBeta->Write();
     m_h_pOffPtvsDeltaTheta->Write();
     m_h_pOffPtvsDeltaThetaForProf->Write();
     m_prof_pOffPtvsDeltaTheta = m_h_pOffPtvsDeltaThetaForProf->ProfileX();
     m_prof_pOffPtvsDeltaTheta->Write();
     for(Int_t i = 0; i < 6; ++i){
          m_h_pOffPtvsDeltaThetaEtaIndex[i]->Write();
     }
     m_h_DeltaThetavsHits->Write();
     m_h_BarrelBetavsTruth->Write();

     m_h_eOffPt->Write();
     m_h_eL1Pt->Write();
     m_h_eSAPt->Write();
     m_h_eOffEta->Write();
     m_h_eL1Eta->Write();
     m_h_eSAEta->Write();
     m_h_eOffPhi->Write();
     m_h_eL1Phi->Write();
     m_h_eSAPhi->Write();
     m_h_eOffPtBarrel->Write();
     m_h_eL1PtBarrel->Write();
     m_h_eSAPtBarrel->Write();
     m_h_eCBPtBarrel->Write();
     m_h_eOffPtTransition->Write();
     m_h_eL1PtTransition->Write();
     m_h_eSAPtTransition->Write();
     m_h_eOffPtEndcap->Write();
     m_h_eL1PtEndcap->Write();
     m_h_eSAPtEndcap->Write();
     m_h_eOffPtForward->Write();
     m_h_eL1PtForward->Write();
     m_h_eSAPtForward->Write();
     for(Int_t i = 0;i < 4;++i){
          m_h_eL1PtBarrelSP[i]->Write();
          m_h_eSAPtBarrelSP[i]->Write();
          m_h_eL1PtTransitionSP[i]->Write();
          m_h_eSAPtTransitionSP[i]->Write();
          m_h_eL1PtEndcapSP[i]->Write();
          m_h_eSAPtEndcapSP[i]->Write();
          m_h_eL1PtForwardSP[i]->Write();
          m_h_eSAPtForwardSP[i]->Write();
     }
     m_h_eL1PtBarrelMyLUTAlpha->Write();
     m_h_eSAPtBarrelMyLUTAlpha->Write();
     m_h_eL1PtBarrelMyLUTBeta->Write();
     m_h_eSAPtBarrelMyLUTBeta->Write();
     m_h_eSAPtBarrelkayamashMethod->Write();

     m_h_pSAResPt->Write();
     for(Int_t i = 0; i < 4; ++i){
          m_h_pSAResPtLS[i]->Write();
          m_h_pSAResPtLSRadius[i]->Write();
     }
     m_h_pSAResPtLargeSpecial->Write();
     m_h_pSAResPtBarrelAlpha->Write();
     m_h_pSAResPtBarrelBeta->Write();
     for(Int_t i = 0; i < 5;++i){
          m_h_pSAResPtBarrelAlphaSector[i]->Write();
          m_h_pSAResPtBarrelBetaSector[i]->Write();
     }
     m_h_pSAResPtBarrelAlpha1SP->Write();
     m_h_pSAResPtBarrelAlpha2SP->Write();
     m_h_pSAResPtBarrelBeta2SP->Write();
     m_h_pSAResPtBarrelAlpha3SP->Write();
     m_h_pSAResPtBarrelBeta3SP->Write();
     m_h_pSAResPtBarrel1SP->Write();
     m_h_pSAResPtBarrel2SP->Write();
     m_h_pSAResPtBarrel3SP->Write();
     m_h_pSAResPtvsDeltaTheta->Write();
     m_h_pSAResPtvsDeltaThetaForProf->Write();
     m_h_pAlphaResPtvsDeltaThetaForProf->Write();
     m_prof_pAlphaResPtvsDeltaTheta = m_h_pAlphaResPtvsDeltaThetaForProf->ProfileX();
     m_prof_pAlphaResPtvsDeltaTheta->Write();
     m_prof_pSAResPtvsDeltaTheta = m_h_pSAResPtvsDeltaThetaForProf->ProfileX();
     m_prof_pSAResPtvsDeltaTheta->Write();
     m_h_pASResPtBarrelBetaSmallDeltaTheta->Write();//!
     m_h_pASResPtBarrelBetaLargeDeltaTheta->Write();//!

     TF1 *SAResolution = new TF1("SAResolution","([0]*x*x*x/(1000.*1000.)+[1]*x*x/1000.+[2]*x+[3]/1000.)/1000.",0.05,0.5);
     SAResolution->SetParameter(0,3.5);
     SAResolution->SetParameter(1,-1.8);
     SAResolution->SetParameter(2,0.35);
     SAResolution->SetParameter(3,-0.017);
     SAResolution->Draw();
     SAResolution->Write();
     TF1 *IDResolution = new TF1("IDResolution","([0]*x+[1]/1000.)/1000.",0.05,0.5);
     IDResolution->SetParameter(0,0.017);
     IDResolution->SetParameter(1,0.000000418);
     IDResolution->Draw();
     IDResolution->Write();
     cout<<"finish!"<<endl;
}