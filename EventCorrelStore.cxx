// $Id: TSCNAnlEvent.cxx 755 2011-05-20 08:04:11Z linev $
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

#include "EventCorrelStore.h"
#include "TAidaConfiguration.h"
void EventCorrelStore::Clear(Option_t *t)
{
       cFRS_AoQ = 0;
       cFRS_ID_x2 = 0;
       cFRS_ID_x4 = 0;
       cFRS_z = 0;
       cFRS_z2 = 0; 
       
    TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
   
         cdT_AIDA_FRS = 0;
         AIDA_implantation_gate_DSSD1 = false;
         AIDA_implantation_gate_DSSD2 = false;
         AIDA_implantation_gate_DSSD3 = false;
         cImplantIterator = -1;
         cDecayIterator = -1;
         for(int i=0; i<512;i++){
           AIDA_implantation_gate_DSSD1_StripX[i] = 0;
           AIDA_implantation_gate_DSSD1_StripY[i] = 0;
         }
         cAIDA_dT_imp_decay_hits=0;
         
         cGalE =-1;
         for(int i =0; i<50; i++) cFatE[i] =-1;
         for(int i =0; i<48; i++) cbPlasE[i] =-1;
            cAIDA_WR = -1;
            cbPlas_WR= -1;
            cFAT_WR= -1;
            cGAL_WR= -1;
            cFRS_WR= -1;
        for(int i =0; i<50; i++){
          cAIDA_dT_imp_decay[i] =0;
          cAIDA_dT_imp_decay_FRS_gated[i]=0;        
            }
            for(int i=0; i<200; i++) {               
                cAIDAImplantE[i]=0;
                         
                }
                 for(int i=0; i<1000; i++) {   
                 cAIDADecayE[i] = 0; 
                 }
         //cAIDA_dT_imp_decay = 0;
       //  for(int i =0; i<128; i++){
//            AIDA_implantation_gate_DSSD1_StripX[i] = 0;
//            AIDA_implantation_gate_DSSD1_StripY[i] = 0;
        // }
}
