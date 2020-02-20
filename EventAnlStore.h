// $Id: EventAnlStore.h 755 2011-05-20 08:04:11Z linev $
//-----------------------------------------------------------------------
//       The GSI Online Offline Object Oriented (Go4) Project
//         Experiment Data Processing at EE department, GSI
//-----------------------------------------------------------------------
// Copyright (C) 2000- GSI Helmholtzzentrum f�r Schwerionenforschung GmbH
//                     Planckstr. 1, 64291 Darmstadt, Germany
// Contact:            http://go4.gsi.de
//-----------------------------------------------------------------------
// This software can be used under the license agreements as stated
// in Go4License.txt file which is part of the distribution.
//-----------------------------------------------------------------------

#ifndef TSCNCALEVENT_H
#define TSCNCALEVENT_H

#include "TGo4EventElement.h"
//#include "TSCNUnpackEvent.h"
#include "EventUnpackStore.h"
#include "AIDA_Event.h"

  
struct AidaAnlData {
      std::vector<AidaHit> Implants;
      std::vector<AidaHit> Decays;
      
      
      };
class EventAnlStore : public TGo4EventElement {
   public:
      EventAnlStore() : TGo4EventElement() {}
      EventAnlStore(const char* name) : TGo4EventElement(name) {}
      virtual ~EventAnlStore() {}

      virtual void  Clear(Option_t *t="");

      ///General Outputs 
      Int_t pEvent_Number;
      Int_t pPrcID_Conv[7];
      Int_t pUsed_Systems[7];
      
      ///White Rabbit Outputs
      Long64_t pFRS_WR;
      Long64_t pbPLAS_WR;
      Long64_t pFAT_WR;
      Long64_t pAIDA_WR;
      Long64_t pGAL_WR;
      ///FRS Outputs
    
      Float_t pFRS_AoQ;
      Float_t pFRS_ID_x2;
      Float_t pFRS_ID_x4;
      Float_t pFRS_z;
      Float_t pFRS_z2;
      Int_t   pSci_num;
      Float_t pFRS_sci_l[12];
      Float_t pFRS_sci_r[12];
      Float_t pFRS_sci_e[12];
      Float_t pFRS_sci_tx[12];
      Float_t pFRS_sci_x[12];
      
      Bool_t pFRS_ZAoQ_pass;
      Bool_t pFRS_x2AoQ_pass;
      Bool_t pFRS_x4AoQ_pass;
      Bool_t pFRS_Z_Z2_pass;
      
     ///AIDA output 
       AidaAnlData pAida;
   

      ///Plastic AnlProc Outputs
      Int_t    pbPlas_QDCFired; 
      Int_t    pbPlas_QDCID[32];
      Double_t pbPlas_QDCGainMatch_i[32];
      Int_t    pbPlas_TDCFired;
      Int_t    pbPlas_TDC_Multiplicity[32]; 
      Int_t    pbPlas_TDCID[50];
      Double_t pbPlasTDC_T[32];
//       Double_t pbPlas_SC41_dT[32];
//       Double_t pbPlas_SiPM_dT_Calib[32];
     
     
      Int_t    pFat_QDCFired;
      Int_t    pFat_QDCID[50];
      Double_t pFat_QDCGainMatch[50];
      Int_t    pFat_TDCFired;
      Int_t    pFat_TDCID[50];
      Long64_t pFat_TDC_T[50];
      //Double_t pFat_SC41_dT_Calib[50];
      Int_t    pFat_TDC_Multipl_perCh[50];
      //Double_t pFat_Ch_dT[50];
    //  Double_t pFat_Ch0_TDC;
      
      
      
      Double_t pFat_ToTCalib[50];
      Double_t pFat_LeadT[50][10];
      Double_t pFat_TrailT[50][10];
      Int_t    pFat_LeadHits;
      Int_t    pFat_TrailHits;
      ///Disabled 06.01.20 AKM 
//       Int_t    pFat_firedQDC_Comb;
//       Double_t pFat_QDC_E_Comb[50];
//       Int_t    pFat_QDC_ID_Comb[50];
      
      Double_t pbPlas_ToTCalib[48];
      Int_t    pbPlas_PMT_Lead_N[48];
      Double_t pbPlas_LeadT[48][10];
      Double_t pbPlas_LeadT_Avg;
      Double_t pbPlas_TrailT[48][10];
      Int_t    pbPlas_LeadHits;
      Int_t    pbPlas_TrailHits;
      
     
      
      //Int_t    pGalFired;
      //Int_t    pGalID[32];
      //Long64_t pGalT[32];
      //Double_t pGalE_Cal_i[32];
     // Double_t pGal_Energy;
      //Long64_t pGal_dT;
      //Double_t pGalE_Addback;

      // New GALILEO Outputs NH 19.02.2020
      Long64_t pGal_T[GALILEO_MAX_DETS][GALILEO_CRYSTALS];
      double   pGal_E[GALILEO_MAX_DETS][GALILEO_CRYSTALS];
      double   pGal_EAddback[GALILEO_MAX_DETS][GALILEO_CRYSTALS];
      
//       Int_t pFing_firedTamex; 
//       Int_t pFing_iterator[4];
//       Int_t pFing_LeadChan[4][32];
//       Double_t pFing_leadT[4][32];
//       Double_t pFing_LeadDiff[100];
//       Double_t pFing_LeadPlus[100];
//       Double_t pFing_SC41_diff[100];
//       Double_t pFing_TOT[4][32]; 
//       Double_t pFing_pos_ToT[4][32];
//       Double_t pFing_downData;
//       Double_t pFing_total_time;
      
      Int_t pFing_tot;
      Double_t pFing_lead_lead[52];
      Double_t pFing_Lead_Up[52][100];
      Double_t pFing_Lead_Down[52][100];
      Double_t pFing_Trail_Up[52][100];
      Double_t pFing_Trail_Down[52][100];   
      Int_t pFing_stripID;
      Double_t pFing_maxtot;
      Int_t pFing_maxtotchan;

   ClassDef(EventAnlStore,1)
};
#endif //TSCNCALEVENT_H
