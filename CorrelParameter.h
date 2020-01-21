// parameter holding data needed for correlation routine

#ifndef CORRELPAR_H
#define CORRELPAR_H

#include "Go4StatusBase/TGo4Parameter.h"
#include <iomanip>
#include <fstream>
class CorrelParameter : public TGo4Parameter {

public:
          CorrelParameter();
          CorrelParameter(const Text_t* name);
  virtual ~CorrelParameter();
  Int_t    PrintParameter(Text_t *buf, Int_t);
  Bool_t   UpdateFrom(TGo4Parameter *);
  int      IsData(std::ifstream &f);

 Int_t GFat_Egate_low, GFat_Egate_high;
 Int_t GFat_TRawgate_low, GFat_TRawgate_high;
 Int_t GbPlas_Egate_low, GbPlas_Egate_high;
 Int_t GFRS_AIDA_TLow, GFRS_AIDA_THigh;
 Int_t GAIDA_bPlas_TLow, GAIDA_bPlas_THigh;
 Bool_t GFRS_AIDA_DSSD1,  GFRS_AIDA_DSSD2, GFRS_AIDA_DSSD3;
//  Double_t Aff[60], Bff[60];
//  Double_t Ab[32], Bb[32];
//  Double_t Ag[4], Bg[4];

ClassDef(CorrelParameter,6)
};

#endif
