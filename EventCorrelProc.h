// $Id: EventCorrelProc.h 755 2011-05-20 08:04:11Z linev $
//Mistry 10.04.19
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

#ifndef EVENTCORRELPROCESSOR_H
#define EVENTCORRELPROCESSOR_H

#include "TGo4EventProcessor.h"

#include "CorrelParameter.h"
#include "EventUnpackStore.h"
#include "EventAnlStore.h"
#include "AIDA_Headers.h"
#include "AIDA_Event.h"
#include "AIDA_Data_Types.h"
#include "Go4ConditionsBase/TGo4WinCond.h"
#include "Go4ConditionsBase/TGo4PolyCond.h"
#include "DESPEC_Array_Sizes.h"

#include "TTree.h"
#include "TFile.h"

//Max lengths
#define  maxdef(a,b) ( ((a) > (b)) ? (a) : (b) )

class EventAnlStore;
class EventCorrelStore;
class TSCNParameter;

class EventCorrelProc : public TGo4EventProcessor {
   public:
      EventCorrelProc();
      EventCorrelProc(const char * name);
      virtual ~EventCorrelProc();
       CorrelParameter *fCorrel;
 
      virtual Bool_t BuildEvent(TGo4EventElement* dest);
        std::vector<AidaCluster> EventsToClusters(std::vector<AidaEvent> const&);
        AidaHit ClusterPairToHit(std::pair<AidaCluster, AidaCluster> const&);

            
    void Make_FRS_AIDA_Histos();
    void Make_FRS_Prompt_AIDA_FATIMA_Ge_Histos();
    void Make_FRS_AIDA_bPlast_Histos();
    void Make_FRS_Delayed_AIDA_Gamma_Histos();
    void Make_FRS_AIDA_bPlas_FATIMA_Histos();
    
     
    void Process_FRS_AIDA(EventAnlStore* cInput, EventCorrelStore* cOutput);
    void Process_FRS_Prompt_AIDA_FATIMA_Ge(EventAnlStore* cInput, EventCorrelStore* cOutput);
    void Process_FRS_AIDA_bPlast(EventAnlStore* cInput, EventCorrelStore* cOutput);
    void Process_FRS_Delayed_AIDA_Gamma(EventAnlStore* cInput, EventCorrelStore* cOutput);
    void Process_FRS_AIDA_FATIMA(EventAnlStore* cInput, EventCorrelStore* cOutput);
    
    TGo4WinCond* MakeWindowCond(const char* foldername,
                  const char* condname,
                  float left = 0.,
                  float right = 4096.,
                  const char* HistoName = 0);
                                 
        TGo4PolyCond* MakePolyCond(const char* foldername,
                 const char* condname,
                 Int_t size,
                 Float_t (*points)[2],
                 const char* HistoName = 0);
     
      std::vector<TH2*> hA_FRS_Z1Z2_implants_strip_xy;
      std::vector<TH2*> hA_FRS_Z1Z2_implants_pos_xy;
      std::vector<TH1*> hA_FRS_Z1Z2_implants_e;
      std::vector<TH2*> hA_FRS_Z1Z2_implants_e_xy;
      
      std::vector<TH2*> hA_FRS_Z1Z2_x2x4AoQ_implants_strip_xy;
      std::vector<TH2*> hA_FRS_Z1Z2_x2x4AoQ_implants_pos_xy;
      std::vector<TH1*> hA_FRS_Z1Z2_x2x4AoQ_implants_e;
      std::vector<TH2*> hA_FRS_Z1Z2_x2x4AoQ_implants_e_xy;
      std::vector<TH1*> hA_FRS_Z1Z2_x2x4AoQ_implants_time_delta;
      std::vector<TH1*> hA_FRS_Z1Z2_x2x4AoQ_implants_strip_1d;
      std::vector<TH1*> hA_FRS_Z1Z2_x2x4AoQ_implants_per_event;
      std::vector<TH2*> hA_FRS_Z1Z2_x2x4AoQ_implants_strip_xy_dssdg;
      TH1 *hA_FRS_dT;
      TH1 *hA_bPlas_dT;
      TH1 *hADec_bPlas_dT;
      
      TH1 *hFRS_bPlas_dT;
      
      TH1 *hA_impdec_dT;
      TH1 *hA_impdec_dT_FRS_gated;
      
      
      TGo4WinCond  *fconHis1;
      TGo4WinCond  *fWinCon1;
      
      TH1 *hA_implant_FatdT;
      TH1 *hA_implant_FatE;
      TH2 *hA_implant_FatdT_FatE;
      
      TH1 *hA_implant_GaldT;
      TH1 *hA_implant_GalE;
      TH2 *hA_implant_GaldT_GalE;
      
      TH1 *hA_implantBplas_GalE;
      TH1 *hA_implant_bPlasdT;
      
      std::vector<TH2*> hA_bPlas_decays_strip_xy;
      std::vector<TH1*> hA_bPlas_decays_e;
      TH1 *hA_bPlas_Energy_AidaImpGated[2][16];
      TH1 *hA_bPlas_Energy_AidaDecayGated[2][16];
      TH1 *hA_bPlas_E;
      TH1 *hADecay_bPlas_E;
      TH1 *hA_FRSgated_bPlas_E;
      TH1 *hA_gatedE_bPlas_E;
      
      TH1 *hA_Ge_WRdT;
      TH1 *hA_dT_GeE;
      TH2 *hA_dT_imp_decay_vs_GeE;
      TH1 *hA_dT_FRS_Gated_GeE;
      TH2 *hA_dT_FRS_Gated_imp_decay_vs_GeE;
      
      TH1 *hA_Fat_WRdT;
      TH1 *hA_dT_FatE;
      TH2 *hA_dT_imp_decay_vs_FatE;
      TH1 *hA_dT_FRS_Gated_FatE;
      TH2 *hA_dT_FRS_Gated_imp_decay_vs_FatE;
     
      TH1 *hbPlas_FRS[3][16];
      TH1 *dT_bPlas_FRS;
      TH1 *hA_bPlas_E_Ch[3][16];
      TH1 *hA_gatedE_bPlas_E_Ch[3][16];
      
    //  bool AIDA_implantation_gate;
      TGo4PolyCond  *cAIDA_IMPgate_DSSD1;
      TGo4PolyCond  *cAIDA_IMPgate_DSSD2;
      TGo4PolyCond  *cAIDA_IMPgate_DSSD3;

      int event_number;
      Long64_t AIDA_WR;
      Long64_t FRS_WR;
      Long64_t bPLAS_WR;
      Long64_t FAT_WR;
      Long64_t GAL_WR;
      Long64_t dT_AIDA_FRS;
      Long64_t dT_AIDA_bPlas;
      Long64_t dT_FRS_bPlas;
      Long64_t dT_AIDA_Gal;
      Long64_t dT_AIDA_Fat;
      TSCNParameter *fParam1;
      
      //Long64_t*** aidatime;

       Long64_t aida_imptime[3][128][128];
       Long64_t aida_imptime_FRS_gated[3][128][128];
      
      ClassDef(EventCorrelProc, 1)
	};
#endif //TSCNANLPROCESSOR_H
