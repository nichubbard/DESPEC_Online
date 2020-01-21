// $Id: TSCNUnpackEvent.h 755 2011-05-20 08:04:11Z linev $
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

#ifndef TSCNEVENT_H
#define TSCNEVENT_H

//#define SCN_NUM_CHAN 32
#define RAW_DATA_LENGTH 1024
#define E_DATA_LENGTH   512
#define FADC_CHAN 8
#define ADC_NUM_CHAN 32
#include "AIDA_Decay_Event_Store.h"

#include "AIDA_Headers.h"
#include "AIDA_Event.h"
#include "AIDA_Data_Types.h"
#include "TGo4MbsEvent.h"
#include "Data_Stream.cxx"
#include "EventBuilder.cxx"
//#include "TSCNUnpackEvent.h"
#include "EventUnpackStore.h"
#include "AIDA_Decay_Event_Store.h"
//#include "AIDA_Processor.h"

#include "Detector_System.cxx"
#include "AIDA_Headers.h"
#include "AIDA_Event.h"
#include "AIDA_Data_Types.h"
#include "TGo4MbsEvent.h"
#include "Data_Stream.cxx"
#include "EventBuilder.cxx"
//#include "TSCNUnpackEvent.h"
#include "EventUnpackStore.h"
#include "AIDA_Decay_Event_Store.h"
//#include "AIDA_Processor.h"


#include "AIDA_Event.h"
//#include "AIDA_Processor.h"

#include "Detector_System.cxx"

#include "TGo4EventElement.h"

class AidaMover : public TGo4EventElement 
{
 public:
 AidaMover() : TGo4EventElement() {}
 AidaMover(const char* name) : TGo4EventElement(name) {}
  virtual ~AidaMover() {}
  
  /**
   * Method called by the event owner (analysis step) to clear the
   * event element.
   */
  virtual void Clear(Option_t *t="");
    
   // EventAida Povver;

  
  ClassDef(AidaMover,1)
    };
#endif //TSCNEVENT_H



