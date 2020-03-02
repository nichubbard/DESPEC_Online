// $Id: TSCNCalEvent.cxx 755 2011-05-20 08:04:11Z linev $
//-----------------------------------------------------------------------
//       The GSI Online Offline Object Oriented (Go4) Project
//         Experiment Data Processing at EE department, GSI
//-----------------------------------------------------------------------
// Copyright (C) 2000- GSI Helmholtzzentrum für Schwerionenforschung GmbH
//                     Planckstr. 1, 64291 Darmstadt, Germany
// Contact:            http://go4.gsi.de
//-----------------------------------------------------------------------
// This software can be used under the license agreements as stated
// in Go4License.txt file which is part of the distribution.
//-----------------------------------------------------------------------

#include "EventAnlStore.h"
void EventAnlStore::Clear(Option_t *t)
{
    
      pFRS_WR = 0;
      pbPLAS_WR = 0;
      pFAT_WR = 0;
    //  pAIDA_WR = 0;
      pGAL_WR = 0;
      
      pAida.Implants.clear();
      pAida.Decays.clear();
      pEvent_Number=0;
    for(int i=0; i<7; i++){
     pPrcID_Conv[i] = 0;
     pUsed_Systems[i] = 0;   
    }
       pFRS_AoQ = 0;
       pFRS_ID_x2 = 0 ;
       pFRS_ID_x4 = 0 ;
       pFRS_z = 0 ;
       pFRS_z2 = 0 ;
       pSci_num = 0;
       
       
       for(int l=0; l<12;l++){
   
        pFRS_sci_l[l] = 0;
        pFRS_sci_r[l] = 0;
        pFRS_sci_e[l] = 0;
        pFRS_sci_tx[l] = 0;
        pFRS_sci_x[l] = 0;
     
     }
   
       pFRS_ZAoQ_pass = false;
       pFRS_x2AoQ_pass = false;
       pFRS_x4AoQ_pass = false;
       pFRS_Z_Z2_pass = false;

    for(int i=0; i< 50; i++){
        //Fatima QDC ID
    pFat_QDCID[i] = -1;
    //Fatima Gainmatched energy
    pFat_QDCGainMatch[i] = 0;
    //Fatimsa TDC ID
    pFat_TDCID[i] = -1;
    //Fatima gainmatched time 
    pFat_TDC_T[i] = 0;

    //Fatima TDC hits/channel
    pFat_TDC_Multipl_perCh[i] = 0;

    pFat_LeadHits = 0;
    

    for(int j=0; j<10; j++){
      pFat_ToTCalib[i][j] =0;   
        }
    }
    pbPlas_LeadHits = 0;
    pbPlas_TrailHits = 0;
    pbPlas_LeadT_Avg = 0;
    for(int i =0; i<3;i++){
           for(int j =0; j<16; j++){
                pbPlas_PMT_Lead_N[i][j] = 0;
                pbPlas_PMT_Trail_N[i][j] = 0;
                
                for(int k =0; k<10; k++){
                    pbPlas_ToTCalib[i][j][k] = 0;
                    pbPlas_LeadT[i][j][k] = 0;
                    pbPlas_TrailT[i][j][k] = 0;
            }
        }
    }
///      pFat_LeadT[i][j] = 0;
///         pFat_TrailT[i][j] = 0;
    for (int i = 0; i < GALILEO_MAX_DETS; i++)
        {
          for (int j = 0; j < GALILEO_CRYSTALS; j++)
          {
            pGal_T[i][j] = 0;
            pGal_E[i][j] = 0;
            pGal_EAddback[i][j] = 0;
          }
        }
    
    
    
    //FINGER
    
    pFing_tot = 0;
    pFing_stripID = 0;
    pFing_maxtot= 0;
    pFing_maxtotchan= 0;
    for(int i=0; i<52; i++){
        pFing_lead_lead[i] = 0;
    
        for(int j=0; j<100; j++){
            pFing_Lead_Up[i][j] = 0;
            pFing_Lead_Down[i][j] = 0;
            pFing_Trail_Up[i][j] = 0;
            pFing_Trail_Down[i][j] = 0;
    
        }
    }
//     pFing_firedTamex = -1;
//     pFing_downData = -1;
//     pFing_total_time = -1;
//     
//     for(int i =0; i<4; i++){
//         pFing_iterator[i]=-1;
//         for(int j =0; j<32; j++){
//             pFing_leadT[i][j] = 0;
//             pFing_TOT[i][j] = 0;
//             pFing_pos_ToT[i][j] = 0;
//             pFing_LeadChan[4][32]= 0;
//         }
//     }
//     for(int k=0; k<100; k++){
//         pFing_LeadDiff[k] = 0;
//         pFing_LeadPlus[k] = 0;
//         pFing_SC41_diff[k] = 0;
//     }
}
