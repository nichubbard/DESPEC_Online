// $Id: TSCNAnalysis.h 524 2009-11-11 09:53:34Z adamczew $
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

#ifndef DESPECANALYSIS_H
#define DESPECANALYSIS_H

#include "TGo4Analysis.h"
#include "EventUnpackFact.h"

class TH1D;
class TGo4MbsEvent;
class EventUnpackStore;
class TSCNParameter;
class EventAnlStore;
class EventCorrelStore;


class DESPECAnalysis : public TGo4Analysis  {
   public:
      DESPECAnalysis();
      DESPECAnalysis(int argc, char** argv);
      virtual ~DESPECAnalysis() ;
      virtual Int_t UserPreLoop();
      virtual Int_t UserEventFunc();
      virtual Int_t UserPostLoop();
      
   private:
      TGo4MbsEvent       *fMbsEvent;
      EventUnpackStore   *fRawEvent;
     // EventUnpackStore   *fOutput;
      AIDA_Event          *Aida_inp;
      EventAnlStore       *fAnlEvent;
      EventCorrelStore       *fCorrelEvent;
      TSCNParameter      *fPar;


      Int_t               fEvents;
      Int_t               fLastEvent;

   ClassDef(DESPECAnalysis,1)
};

#endif //TSCNANALYSIS_H



