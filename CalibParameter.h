#ifndef CALIBPAR_H
#define CALIBPAR_H

#include "Go4StatusBase/TGo4Parameter.h"
#include <iomanip>

class CalibParameter : public TGo4Parameter {

public:
          CalibParameter();
          CalibParameter(const Text_t *name);
  virtual ~CalibParameter();
  Int_t   PrintParameter(Text_t *buf, Int_t);
  Bool_t  UpdateFrom(TGo4Parameter *);
  int      IsData(std::ifstream &f);

 Double_t AplasQDC[32], BplasQDC[32];
 Double_t AplasTDC_Chref_dT[32],BplasTDC_SC41dT[32];
 Double_t Abplas_TAMEX[32], Bbplas_TAMEX[32];
 Int_t DetIDPlas_TAMEX;
 Double_t Afat[50], Bfat[50], Cfat[50], Dfat[50], Extra1fat[50], TFatTDC_Chref_dT[50], TFatTDC_SC41dT[50];
 Double_t Afat_TAMEX[50], Bfat_TAMEX[50], Cfat_TAMEX[50], Dfat_TAMEX[50], Extra1fat_TAMEX[50];
 Int_t DetIDFAT_TAMEX;
 Int_t DetIDFAT, TDCfatID;
 Int_t DetIDGal;
 Double_t AGal[36], BGal[36], CGal[36];
//  Double_t Aff[60], Bff[60];
//  Double_t Ab[32], Bb[32];
//  Double_t Ag[4], Bg[4];

ClassDef(CalibParameter,5)
};

#endif
