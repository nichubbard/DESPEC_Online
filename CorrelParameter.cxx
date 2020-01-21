#include "CorrelParameter.h"
#include <iostream>

using namespace std;

CorrelParameter::CorrelParameter()
: TGo4Parameter("AnaOnlineCorPar")
{
    GFat_Egate_low = 0;
    GFat_Egate_high = 0;
    GFat_TRawgate_low = 0;
    GFat_TRawgate_high = 0;
    GbPlas_Egate_low = 0;
    GbPlas_Egate_high = 0;
    
    GFRS_AIDA_TLow = 0;
    GFRS_AIDA_THigh = 0;
    GAIDA_bPlas_TLow = 0;
    GAIDA_bPlas_THigh = 0;
    
    GFRS_AIDA_DSSD1 = false;
    GFRS_AIDA_DSSD2 = false;
    GFRS_AIDA_DSSD3 = false;
    
}
//-------------------------------------------------------------------------------------------//
// constructor
// reads correlation parameters from the file Correlations.dat
CorrelParameter::CorrelParameter(const Text_t* name)
: TGo4Parameter(name)
{
  ifstream    file;
  
      file.open("Configuration_Files/Correlations.dat");
  if (file.fail()) {
        cout << "ERROR: CorrelParameter - Could not open file: Configuration_Files/Correlations.dat ! Params set to nominal\n";    
         GFRS_AIDA_TLow = -1000;
         GFRS_AIDA_THigh = 1000;
         GAIDA_bPlas_TLow = -1000;
         GAIDA_bPlas_THigh = 1000;
         GFRS_AIDA_DSSD1 = true;
         GFRS_AIDA_DSSD2 = true;
         GFRS_AIDA_DSSD3 = true;
         
         GFat_Egate_low = 0;
         GFat_Egate_high = 10000;
         GFat_TRawgate_low = 0;
         GFat_TRawgate_high = 10000;
         GbPlas_Egate_low = 0;
         GbPlas_Egate_high = 5000;
  }
  
else {
cout << "CorrelParameter - reading from Configuration_Files/Correlations.dat";
       //AIDA FRS WR Time Gate
       if(IsData(file)) file >> GFRS_AIDA_TLow >>GFRS_AIDA_THigh;
       //AIDA bPlas WR Time Gate
       if(IsData(file)) file >> GAIDA_bPlas_TLow >>GAIDA_bPlas_THigh;
       //Use AIDA  DSSD Gate
       if(IsData(file)) file >> GFRS_AIDA_DSSD1 >>GFRS_AIDA_DSSD2 >>GFRS_AIDA_DSSD3;
       //Fatima energy gates
       if(IsData(file)) file >> GFat_Egate_low>> GFat_Egate_high;
       //Fatima Raw TDC gate
       if(IsData(file)) file >> GFat_TRawgate_low>> GFat_TRawgate_high;
       // Plastic energy gates 
       if(IsData(file)) file >> GbPlas_Egate_low>> GbPlas_Egate_high;
       
       
       if (file.fail()) cout << " !!!There is a problem in the correlation file !!!";
    
  }
  file.close();
}
//--------------------------------------------------------------------------//
CorrelParameter::~CorrelParameter()
{}
//---------------------------------------------------------------------------//
Int_t CorrelParameter::PrintParameter(Text_t *buf, Int_t)
{
  cout << "\n Online Analysis Correlation Parameters: " << GetName() << endl;
  cout << "//////AIDA FRS WR Time Gate\n";
  cout << "////Lower: " << GFRS_AIDA_TLow << "  \t Upper: = " <<  GFRS_AIDA_THigh <<  endl;
  
  cout << "//////AIDA bPlas WR Time Gate\n";
  cout << "////Lower: " << GAIDA_bPlas_TLow << "  \t Upper: = " <<  GAIDA_bPlas_THigh <<  endl;
  
  cout << "//////Use AIDA DSSD 2D position Gate\n";
  cout << "////DSSD1: " << GFRS_AIDA_DSSD1 << "  \t DSSD2: " <<  GFRS_AIDA_DSSD2 << "  \t DSSD3: "<<GFRS_AIDA_DSSD3 <<endl;

  cout << "////FATIMA Gamma Gate limits\n";
  cout << "////Lower: " << GFat_Egate_low << "  \t Upper: = " <<  GFat_Egate_high <<  endl;
  
  cout << "////FATIMA Raw TDC Gate limits\n";
  cout << "////Lower: " << GFat_TRawgate_low << "  \t Upper: = " << GFat_TRawgate_high <<  endl;
  
  cout << "////bPlastic Energy Gate limits\n";
  cout << "////Lower: " << GbPlas_Egate_low << "  \t Upper: = " << GbPlas_Egate_high <<  endl;
  return 1;
}
//---------------------------------------------------------------------------//
Bool_t CorrelParameter::UpdateFrom(TGo4Parameter *pp)
{
  if(pp->InheritsFrom("CorrelParameter"))
  {
    CorrelParameter *from = (CorrelParameter *) pp;
      GFRS_AIDA_TLow = from->GFRS_AIDA_TLow;
      GFRS_AIDA_THigh = from->GFRS_AIDA_THigh;
      GAIDA_bPlas_TLow = from->GAIDA_bPlas_TLow;
      GAIDA_bPlas_THigh = from->GAIDA_bPlas_THigh;
      
      GFRS_AIDA_DSSD1 = from->GFRS_AIDA_DSSD1;
      GFRS_AIDA_DSSD2 = from->GFRS_AIDA_DSSD2;
      GFRS_AIDA_DSSD3 = from->GFRS_AIDA_DSSD3;
      
      GFat_Egate_low  = from->GFat_Egate_low;
      GFat_Egate_high = from->GFat_Egate_high;
      GFat_TRawgate_low = from->GFat_TRawgate_low;
      GFat_TRawgate_high = from->GFat_TRawgate_high;
      GbPlas_Egate_low  = from->GbPlas_Egate_low;
      GbPlas_Egate_high = from->GbPlas_Egate_high;

       cout << "CorrelParameter - Parameter : " << GetName() << " UPDATED\n";
  }
  else {
      cout << "ERROR: CorrelParameter - Wrong parameter object: " << pp->ClassName() << endl;}
        return kTRUE;
  }
  //---------------------------------------------------------------------------//
  int CorrelParameter::IsData(ifstream &f) {
  char dum;
  char dumstr[300];
  int retval = 0;

  /* 'operator >>' does not read End-of-line, therefore check if read 
      character is not EOL (10) */
  do {
    dum = f.get();
    if (dum == '#')    // comment line => read whole line and throw it away
      f.getline(dumstr,300);
  }
  while ((dum == '#') || ((int)dum == 10)); 

  f.unget();   // go one character back
  retval = 1;
  return retval;
}
ClassImp(CorrelParameter)