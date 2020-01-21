//                 DESPEC analysis software AM 07.03.19
//
//---------------------------------------------------------------
//       The GSI Online Offline Object Oriented (Go4) Project
//       Experiment Data Processing at EE department, GSI
//---------------------------------------------------------------
//
//Copyright (C) 2000- Gesellschaft f. Schwerionenforschung, GSI
//                    Planckstr. 1, 64291 Darmstadt, Germany
//Contact:            http://go4.gsi.de
//----------------------------------------------------------------
//This software can be used under the license agreements as stated
//in Go4License.txt file which is part of the distribution.
//----------------------------------------------------------------

#include "EventUnpackStore.h"

#include "Riostream.h"

//***********************************************************
EventUnpackStore::EventUnpackStore() :
        TGo4EventElement()
{
     //   cout << "**** EventUnpackStore: Create instance" << endl;
}
//***********************************************************
EventUnpackStore::EventUnpackStore(const char* name) :
        TGo4EventElement(name)
{
       // cout << "**** EventUnpackStore: Create instance " << name << endl;
}
//***********************************************************
EventUnpackStore::~EventUnpackStore()
{
      // cout << "**** EventUnpackStore: Delete instance " << endl;
}

//-----------------------------------------------------------
void  EventUnpackStore::Clear(Option_t *t)
{
//     for (int i =0; i<3; ++i){
//         fFRS_Music_dE[i] = 0; 
//         fFRS_Music_dE[i] = 0;
//         fFRS_ID_brho[i] = 0;
//         fFRS_ID_rho[i] = 0;
//     }
//     for(int l=0;l<12;++l){
//         fFRS_sci_l[l] = 0;
//         fFRS_sci_r[l] = 0;
//         fFRS_sci_e[l] = 0;
//         fFRS_sci_tx[l] = 0;
//         fFRS_sci_x[l] = 0;
//     }
    for(int l=0;l<12;++l) fFRS_sci_e[l] = 0;
    
//    fFRS_sci_tofll2 = 0;
//    fFRS_sci_tofll3  = 0;
//    fFRS_sci_tof2 = 0;
//    fFRS_sci_tofrr2 = 0;
//    fFRS_sci_tofrr3 = 0;
//    fFRS_sci_tof3 = 0;
   fFRS_ID_x2 = 0;
//    fFRS_ID_y2 = 0; 
//    fFRS_ID_a2 = 0;
//    fFRS_ID_b2 = 0;
   fFRS_ID_x4 = 0;
//    fFRS_ID_y4 = 0;
//    fFRS_ID_a4 = 0;
//    fFRS_ID_b4 = 0;
//    fFRS_sci_dt_21l_21r = 0;
//    fFRS_sci_dt_41l_41r = 0;
//    fFRS_sci_dt_42l_42r = 0;
//    fFRS_sci_dt_43l_43r = 0;
//    fFRS_sci_dt_21l_41l = 0;
//    fFRS_sci_dt_21r_41r = 0;
//    fFRS_sci_dt_21l_42l = 0;
//    fFRS_sci_dt_21r_42r = 0;

//    fFRS_beta = 0;
//    fFRS_beta3 = 0;
//    fFRS_gamma = 0;
   fFRS_AoQ = 0;
   fFRS_AoQ_corr = 0;
   fFRS_z = 0;
   fFRS_z2 = 0;
//   fFRS_z3 = 0;
   /*fFRS_timestamp = 0;
   fFRS_ts = 0;
   fFRS_ts2 = 0; */  
   fFRS_WR = 0;
   fAIDAHits = 0;
    
    AIDATime = 0;
    AIDAHits = 0;
    Aida.clear();
    //fAIDA_WR =0;
    
        fevent_number = 0;
        fArray_count = 0;
//         fVME_TAMEX_bPlas = false;
//         fVME_TAMEX_Fatima = false;
//         fVME_AND_TAMEX_Fatima = false;
        //fVME_AND_TAMEX_Fatima = false;
        for (int i=0; i<7;i++){
         fProcID[i] = -1;
         fUsed_Systems[i] = 0;
         }
        
        fbPlas_VME_firedQDC = -1;    
        fbPlas_VME_firedTDC = -1;
        fbPlas_VME_QDC_Multiplicity = 0;

        fScalar_fired = -1;
        fScalar_ID = -1;
      
  
        for (int i=0; i<32;i++){         
         fbPlas_VME_QDC_E[i] = 0;
         fbPlas_VME_QDC_ID[i] = -1;
         fbPlas_VME_QDC_E_AIDA[i] = 0;
         }
        for (int i=0; i<50; i++){         
         fbPlas_VME_TDC_ID[i] = -1;
        
         for (int j=0; j<32; j++){
          fbPlas_VME_TDC_TS[i][j] = 0;
          fbPlas_VME_TDC_Multiplicity[i] = 0;
         }
        }
         
        //fFat_QDC_ID = -1;
        for (int i=0; i<50; i++){
            fFat_TDC_ID[i] = -1;
            fFat_QDC_ID[i] = -1;
            fFat_QDC_E[i] = 0;
            fFat_QDC_T[i] = 0;
            fSC41[i] = -1;
            fFat_TDC_Multiplicity[i] = 0;
            
          for (int j=0; j<50; j++){
            fFat_TDC_TS[i][j] = 0;        
                }
            }
                  
            fFat_firedQDC = -1;
            fFat_firedTDC = -1;
            fFat_WR = 0;
            fFat_QDC_Multiplicity = 0;
      
            fGal_Pileup = -1;
            fGal_fired = -1;
            fGal_WR = 0;
            
                for (int i=0; i<32; i++){
                    fGal_ID[i] = 0;
                    fGal_E[i] = 0;
                    fGal_T[i] = 0;
       }        
       
        //FINGER 
//         ffing_tamexhits = -1;
//        
//         
//         for(int i =0; i<4; i++){
//             ffing_leadHits[i] = -1;
//             ffing_trailHits[i] = -1;
//             ffing_iterator[i] = -1;
//             ffing_Trig[i] = 0;
//             for(int j=0; j<32; j++){
//                 ffing_tamexCh[i][j] = 0;
//                 ffing_Lead_Phys_Chan[i][j] = 0;
//                 ffing_Trail_Phys_Chan[i][j] = 0;
//                 ffing_chID[i][j] = 0;
//                     }
//             for(int k=0; k<100; k++){    
//                 ffing_Lead_T[i][k] = 0;
//                 ffing_Trail_T[i][k] = 0;
//                 ffing_TOT_added[i][k] =0;
//                 ffing_TOT[i][k] = 0;
//                 ffing_lead_coarse[i][k] = 0;
//                 ffing_lead_fine[i][k]= 0;
//                 ffing_trail_coarse[i][k]= 0;
//                 ffing_trail_fine[i][k]= 0;              
//                     }
//         }

            fFinger_WR = 0;
         for(int i=0; i<2; i++){
            for(int j=0; j<50; j++){
                fFing_SC41_lead[i][j] = 0;
                fFing_SC41_trail[i][j] = 0;
                }
        }
        for (int i = 0; i < 52; i++)
        {
          fFing_Strip_N[i] = 0;
          fFing_Strip_N_LU[i] = 0;
          fFing_Strip_N_TU[i] = 0;
          fFing_Strip_N_LD[i] = 0;
          fFing_Strip_N_TD[i] = 0;
          fFing_PMT_Lead_N[i] = 0;
          fFing_PMT_Trail_N[i] = 0;
          for(int j=0;j<100; j++){
          fFing_Lead_Up[i][j] = 0;
          fFing_Lead_Down[i][j] = 0;
        }
        }
          for (int i = 0; i < 50; i++){
          fFat_Strip_N[i] = 0;
          fFat_PMT_Lead_N[i] = 0;
          fFat_PMT_Trail_N[i] = 0;
          fFat_PMT_Lead_N[i] = 0;
          fFat_PMT_Trail_N[i] = 0;
          for(int j =0; j<10;j++){
          fFat_Lead_PMT[i][j] = 0;
          fFat_Trail_PMT[i][j] = 0;
         
            }
         }
         
         fbPlas_WR = 0;
        for(int i =0; i<48;i++){
          fbPlas_TAMEX_ID[i] = 0;
          fbPlas_Strip_N[i] = 0;
          fbPlas_PMT_Lead_N[i] = 0;
          fbPlas_PMT_Trail_N[i] = 0;
          fbPlas_PMT_Lead_N[i] = 0;
          fbPlas_PMT_Trail_N[i] = 0;
            for(int j =0; j<10;j++){
         fbPlas_Lead_PMT[i][j] = 0;   
         fbPlas_Trail_PMT[i][j] = 0;
        }
        }
        //      fPatternUnit = 0;
}
void  EventUnpackStore::ClearEvent()

{
//         for (int i =0; i<3; ++i){
//         fFRS_Music_dE[i] = 0; 
//         fFRS_Music_dE_corr[i] = 0;
//         
//         }
//         for (int i =0; i<2; ++i){
//         fFRS_ID_brho[i] = 0;
//         fFRS_ID_rho[i] = 0;
//     }
//     for(int l=0;l<12;++l){
//         fFRS_sci_l[l] = 0;
//         fFRS_sci_r[l] = 0;
//         fFRS_sci_e[l] = 0;
//         fFRS_sci_tx[l] = 0;
//         fFRS_sci_x[l] = 0;
//     }
    for(int l=0;l<12;++l){
       
        fFRS_sci_e[l] = 0;
        
    }
    
//    fFRS_sci_tofll2 = 0;
//    fFRS_sci_tofll3  = 0;
//    fFRS_sci_tof2 = 0;
//    fFRS_sci_tofrr2 = 0;
//    fFRS_sci_tofrr3 = 0;
//    fFRS_sci_tof3 = 0;
   fFRS_ID_x2 = 0;
//    fFRS_ID_y2 = 0; 
//    fFRS_ID_a2 = 0;
//    fFRS_ID_b2 = 0;
   fFRS_ID_x4 = 0;
//    fFRS_ID_y4 = 0;
//    fFRS_ID_a4 = 0;
//    fFRS_ID_b4 = 0;
//    fFRS_sci_dt_21l_21r = 0;
//    fFRS_sci_dt_41l_41r = 0;
//    fFRS_sci_dt_42l_42r = 0;
//    fFRS_sci_dt_43l_43r = 0;
//    fFRS_sci_dt_21l_41l = 0;
//    fFRS_sci_dt_21r_41r = 0;
//    fFRS_sci_dt_21l_42l = 0;
//    fFRS_sci_dt_21r_42r = 0;

//    fFRS_beta = 0;
//    fFRS_beta3 = 0;
//    fFRS_gamma = 0;
   fFRS_AoQ = 0;
   fFRS_AoQ_corr = 0;
   fFRS_z = 0;
   fFRS_z2 = 0;
  // fFRS_z3 = 0;
  /* fFRS_timestamp = 0;
   fFRS_ts = 0;
   fFRS_ts2 = 0; */  
   fFRS_WR = 0;
   fAIDAHits = 0;
     
    AIDAHits = 0;
    AIDATime = 0;
    Aida.clear();
   
    fevent_number = 0;
        for (int i=0; i<7;i++){
         fProcID[i] = -1;
         fUsed_Systems[i] = 0;
         }
        
        fbPlas_VME_firedQDC = -1;    
        fbPlas_VME_firedTDC = -1;
        fbPlas_VME_QDC_Multiplicity = 0;
        
        fScalar_fired = -1;
        fScalar_ID = -1;
        
        
        for (int i=0; i<32;i++){         
         fbPlas_VME_QDC_E[i] = 0;
         fbPlas_VME_QDC_ID[i] = -1;
         fbPlas_VME_TDC_Multiplicity[i] = 0;
         }
         for (int i=0; i<50; i++){         
         fbPlas_VME_TDC_ID[i] = -1;
        
         for (int j=0; j<32; j++){
          fbPlas_VME_TDC_TS[i][j] = 0;
         }
        }
         
      //  fFat_QDC_ID = -1;
        for (int i=0; i<50; i++){
         fFat_TDC_ID[i] = -1;
         fFat_QDC_ID[i] = -1;
         fFat_QDC_E[i] = 0;
         fFat_QDC_T[i] = 0;
         fSC41[i] = -1;
         fFat_TDC_Multiplicity[i] = 0;
         
        for (int j=0; j<50; j++){
         fFat_TDC_TS[i][j] = 0;        
                }
            }
                  
        fFat_firedQDC = -1;
        fFat_firedTDC = -1;
        fFat_QDC_Multiplicity = 0;
        fFat_WR = 0;
      
        fGal_Pileup = -1;
        fGal_fired = -1;
        fGal_WR = 0;
       for (int i=0; i<32; i++){
        fGal_ID[i] = 0;
        fGal_E[i] = 0;
        fGal_T[i] = 0;
       }
         //FINGER 
         
         for (int i = 0; i < 52; i++)
        {
          fFing_Strip_N[i] = 0;
          fFing_Strip_N_LU[i] = 0;
          fFing_Strip_N_TU[i] = 0;
          fFing_Strip_N_LD[i] = 0;
          fFing_Strip_N_TD[i] = 0;
          fFing_PMT_Lead_N[i] = 0;
          fFing_PMT_Trail_N[i] = 0;
          for(int j=0;j<100; j++){
          fFing_Lead_Up[i][j] = 0;
          fFing_Lead_Down[i][j] = 0;
        }
        }
        fbPlas_WR =0;
         for (int i = 0; i < 50; i++){
          fFat_Strip_N[i] = 0;
          fFat_PMT_Lead_N[i] = 0;
          fFat_PMT_Trail_N[i] = 0;
          fFat_PMT_Lead_N[i] = 0;
          fFat_PMT_Trail_N[i] = 0;
          for(int j =0; j<10;j++){
          fFat_Lead_PMT[i][j] = 0;
          fFat_Trail_PMT[i][j] = 0;
         
            }
         }
        
        for(int i =0; i<48;i++){
          fbPlas_TAMEX_ID[i] = 0;
          fbPlas_Strip_N[i] = 0;
          fbPlas_PMT_Lead_N[i] = 0;
          fbPlas_PMT_Trail_N[i] = 0;
          fbPlas_PMT_Lead_N[i] = 0;
          fbPlas_PMT_Trail_N[i] = 0;
            for(int j =0; j<10;j++){
         fbPlas_Lead_PMT[i][j] = 0;   
         fbPlas_Trail_PMT[i][j] = 0;
            }
        }
      //  ffing_tamexhits = -1;
       
        
//         for(int i =0; i<4; i++){
//             ffing_leadHits[i] = -1;
//             ffing_trailHits[i] = -1;
//             ffing_iterator[i] = -1;
//             ffing_Trig[i] = 0;
//             for(int j=0; j<32; j++){
//                 ffing_tamexCh[i][j] = 0;
//                 ffing_Lead_Phys_Chan[i][j] = 0;
//                 ffing_Trail_Phys_Chan[i][j] = 0;
//                 ffing_chID[i][j] = 0;
//                     }
//                     
//             for(int k=0; k<100; k++){    
//                 ffing_Lead_T[i][k] = 0;
//                 ffing_Trail_T[i][k] = 0;
//                 ffing_TOT_added[i][k] = 0;
//                 ffing_TOT[i][k] = 0;
//                 ffing_lead_coarse[i][k]= 0;
//                 ffing_lead_fine[i][k]= 0;
//                 ffing_trail_coarse[i][k]= 0;
//                 ffing_trail_fine[i][k]= 0;
//        
//              }
//         }
         
}
//----------------------------END OF GO4 SOURCE FILE ---------------------
