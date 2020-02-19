// $Id: EventAnlProc.cxx 754 2011-05-18 11:04:52Z adamczew $
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

// Uncomment this to align the AIDA ASICs with a pulser
//  Only needed if the ASICs didn't align properly
//#define AIDA_PULSER_ALIGN

#include "EventAnlProc.h"

#include <cstdlib>
#include <math.h>

#include "TH1.h"
#include "TH2.h"
#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"

#include "TGo4WinCond.h"
#include "TGo4Analysis.h"

#include "TGo4Picture.h"

#include "EventAnlStore.h"
//#include "TSCNUnpackEvent.h"
#include "EventUnpackStore.h"
#include "TSCNParameter.h"

#include "TAidaConfiguration.h"

#define ABS(x)  ((x)>=0 ? (x):-(x))  // absolute_value(x)
//-----------------------------------------------------------
EventAnlProc::EventAnlProc() :
TGo4EventProcessor(),
fParam1(0)
{
}
//-----------------------------------------------------------
EventAnlProc::EventAnlProc(const char* name) :
TGo4EventProcessor(name)
{
  //Clear up for AIDA
  implantEvents = 0;
  decayEvents = 0;
  pulserEvents = 0;
  nonsenseEvents = 0;

  cout << "**** EventAnlProc: Create" << endl;

  checkTAMEXorVME();

  TFile* f;
  f = new TFile("Parameters.root");
  if(f->IsOpen()) //The File Parameter.root exist
  {
    cout<<"Reading initial parameters from file Parameters.root"<<endl;
    fParam1=(TSCNParameter*)f->Get("SCNParameter");
    f->Close();
  }
  else //Data file with parameters
  {
    ifstream myfile;
    myfile.open ("Parameters.dat", ios::in);
    if(myfile.is_open()) //I have the file
    {
      fParam1 = (TSCNParameter*)  GetParameter("SCNParameter");
      string intro;
      myfile>>intro;
      cout<<intro<<endl;
      cout<<"Reading data file Parameters.dat"<<endl;
      int a;
      double b;
      for(int i=0;i<SCN_NUM_CHAN;i++)
      {
        myfile>>a;
        fParam1->SetPedestal(i,a);
      }
      myfile>>intro;
      for(int i=0;i<SCN_NUM_CHAN;i++)
      {
        myfile>>b;
        fParam1->SetFactor(i,b);
      }
      myfile.close();
    }
    else  //No Parameters at all. Go to default constructor
    {
      cout<<"No parameters file. Using default constructor"<<endl;
      fParam1 = (TSCNParameter*)  GetParameter("SCNParameter");
    }
  }
  fCal = new CalibParameter("CalibPar");
  AddParameter(fCal);
  // fCal = (MoDSSCalibParameter*) GetParameter("CalibPar");
  if (fCal) fCal->PrintParameter(0,0);
  else cout << "**** ERRR - CalibPar doesn't exist - program will crash.\n";

    fCorrel = new CorrelParameter("CorrelPar");
  AddParameter(fCorrel);
  if (fCorrel) fCorrel->PrintParameter(0,0);
  else cout << "**** ERRR - CorrelPar doesn't exist - program will crash.\n";

  ///FINGER 2D Polygon ToT vs strip gate window test
    Double_t ToTvalues[15]={17,24,30,36,42,48,52,52,48,42,36,30,24,17,17};
    Double_t Stripvalues[15]={0,0,0,0,0,0,0,1,1,1,1,1,1,1,0};
    TCutG* ToTvsStripcut = new TCutG("initialcut",15,ToTvalues,Stripvalues);
    fCond_FingToTvsStrip = new TGo4PolyCond("FING_TOTvsStrip");
    fCond_FingToTvsStrip -> SetValues(ToTvsStripcut);
    AddAnalysisCondition(fCond_FingToTvsStrip);
    fCond_FingToTvsStrip -> Enable();
    delete ToTvsStripcut;

    read_setup_parameters();
    get_used_Systems();

}
//-----------------------------------------------------------
EventAnlProc::~EventAnlProc()
{
  cout << "**** EventAnlProc: Delete" << endl;

}
//-----------------------------------------------------------


Bool_t EventAnlProc::BuildEvent(TGo4EventElement* dest)
{
     for(int i=0; i<7; i++){
        PrcID_Conv[i]=0;
        }

  Bool_t isValid=kFALSE; // validity of output event

  EventUnpackStore* pInput  = (EventUnpackStore*) GetInputEvent();
  EventAnlStore* pOutput = (EventAnlStore*) dest;

  pAida.Implants.clear();
  pAida.Decays.clear();


  if((pInput==0) || !pInput->IsValid()){ // input invalid
    pOutput->SetValid(isValid); // invalid
    return isValid; // must be same is for SetValid
  }
  isValid=kTRUE;

   ///general inputs from the unpacker
    event_number = pInput->fevent_number;
    pOutput->pEvent_Number = event_number;
    VMEorTAMEX_bPlas = VME_TAMEX_bPlas;
    VMEorTAMEX_fatima = VME_TAMEX_Fatima;
    //cout<<"pInput->fVME_TAMEX_Fatima " <<pInput->fVME_TAMEX_Fatima<<endl;
    VMEandTAMEX_fatima = VME_AND_TAMEX_Fatima;
  //  cout<<"VMEorTAMEX_bPlas " << VMEorTAMEX_bPlas <<endl;

  for (int i = 0; i<7; i++){
      if(pInput->fProcID[i]>-1){
      PrcID_Conv[i] = pInput->fProcID[i];
      pOutput->pPrcID_Conv[i] = pInput->fProcID[i];
      Used_Systems[i] = pInput->fUsed_Systems[i];
//  cout<<"event_number " << event_number<< " PrcID_Conv " << PrcID_Conv[i] <<" i " << i << "  Used_Systems[i] " <<  Used_Systems[i] << endl;
    }
  
    pOutput->pUsed_Systems[i] = Used_Systems[i];

  }
 

   static bool create =false;
  //Create histograms
  if (!create)
  {
    Make_WR_Histos();
    if(Used_Systems[0])  Make_FRS_Histos();
    if (Used_Systems[1]) Make_Aida_Histos();
   // if (Used_Systems[2] && VMEorTAMEX_bPlas==true) Make_Plastic_VME_Histos();
    if (Used_Systems[2] && VMEorTAMEX_bPlas==false) Make_Plastic_Tamex_Histos();
    if (Used_Systems[3] && VMEorTAMEX_fatima==true) Make_Fatima_Histos();
   // if (Used_Systems[3] && VMEorTAMEX_fatima==false && VMEandTAMEX_fatima==false) Make_Fatima_Tamex_Histos();

   // if ((Used_Systems[3] || Used_Systems[4]) && VMEorTAMEX_fatima==false && VMEandTAMEX_fatima==true) Make_Fatima_VME_Tamex_Histos();


    if (Used_Systems[5]) Make_Galileo_Histos();
    if (Used_Systems[6]) Make_Finger_Histos();
//     if (Used_Systems[2] && Used_Systems[3]) Make_Fat_Plas_Histos();
//     if (Used_Systems[6] && Used_Systems[2]) Make_Fing_Plas_Histos();
        create = true;
        }
 
        Do_WR_Histos(pInput);
                 /** Now extract the data from the stored Unpacker array (root tree)**/
    ///--------------------------------------/**FRS Input**/------------------------------------------///

      if (PrcID_Conv[0]==0){
          pOutput->pFRS_WR = pInput->fFRS_WR;
        ///MUSIC
//        for(int i =0; i<2; ++i){
//             FRS_dE[i] = pInput->fFRS_Music_dE[i];
//             FRS_dE_cor[i] = pInput->fFRS_Music_dE_corr[i];
//            }
        ///SCI
          
        for(int l=0;l<12;++l){
          
           FRS_sci_l[l] = pInput->fFRS_sci_l[l];
           FRS_sci_r[l] = pInput->fFRS_sci_r[l];
           FRS_sci_e[l] = pInput->fFRS_sci_e[l];
           FRS_sci_tx[l] = pInput->fFRS_sci_tx[l];
           FRS_sci_x[l] = pInput->fFRS_sci_x[l];
           
        }
            
        ///SCI TOF
       // FRS_sci_tofll2 = pInput->fFRS_sci_tofll2;
       // FRS_sci_tofll3 = pInput->fFRS_sci_tofll3;
        //FRS_sci_tof2 = pInput->fFRS_sci_tof2;
       // FRS_sci_tofrr2 = pInput->fFRS_sci_tofrr2;
        //FRS_sci_tofrr3 = pInput->fFRS_sci_tofrr3;
      //  FRS_sci_tof3 = pInput->fFRS_sci_tof3;
        ///ID 2 4
        FRS_ID_x2 = pInput->fFRS_ID_x2;
       
       // FRS_ID_y2 = pInput->fFRS_ID_y2;
       // FRS_ID_a2 = pInput->fFRS_ID_a2;
        //FRS_ID_b2 = pInput->fFRS_ID_b2;

        FRS_ID_x4 = pInput->fFRS_ID_x4;
       // FRS_ID_y4 = pInput->fFRS_ID_y4;
       // FRS_ID_a4 = pInput->fFRS_ID_a4;
        //FRS_ID_b4 = pInput->fFRS_ID_b4;
            ///SCI dT
//         FRS_sci_dt_21l_21r = pInput->fFRS_sci_dt_21l_21r;
//         FRS_sci_dt_41l_41r = pInput->fFRS_sci_dt_41l_41r;
//         FRS_sci_dt_42l_42r = pInput->fFRS_sci_dt_42l_42r;
//         FRS_sci_dt_43l_43r = pInput->fFRS_sci_dt_43l_43r;
//
//         FRS_sci_dt_21l_41l = pInput->fFRS_sci_dt_21l_41l;
//         FRS_sci_dt_21r_41r = pInput->fFRS_sci_dt_21r_41r;
//
//         FRS_sci_dt_21l_42l = pInput->fFRS_sci_dt_21l_42l;
//         FRS_sci_dt_21r_42r = pInput->fFRS_sci_dt_21r_42r;
            ///ID Beta Rho
//         for(int i =0; i<2; ++i){
//        // FRS_ID_brho[i] = pInput->fFRS_ID_brho[i];
//       //  FRS_ID_rho[i] = pInput->fFRS_ID_rho[i];
//         }
//         FRS_beta = pInput->fFRS_beta;
//         FRS_beta3 = pInput->fFRS_beta3;
      //  FRS_gamma  = pInput->fFRS_gamma;
            ///ID Z AoQ
        FRS_AoQ = pInput->fFRS_AoQ;
        FRS_AoQ_corr = pInput->fFRS_AoQ_corr;
//         if(FRS_AoQ_corr>0){
// //cout<<"ANL STAGE FRS_AoQ_corr " <<pInput->fevent_number <<" AoQ CORR " << FRS_AoQ_corr << endl;
//            // cout<<" " << endl;
//         }
        FRS_z = pInput->fFRS_z;
        FRS_z2 = pInput->fFRS_z2;
       // FRS_z3 = pInput->fFRS_z3;
            ///ID Timestamp
       // FRS_timestamp = pInput->fFRS_timestamp;
//         FRS_ts = pInput->fFRS_ts;
//         FRS_ts2 = pInput->fFRS_ts2;
         Do_FRS_Histos(pOutput);
            }


   ///-------------------------------- /**AIDA Input**/ --------------------------------///
      //  if (Used_Systems[1]&&  PrcID_Conv[1]==1) {ProcessAida(pInput);}
        ProcessAida(pInput, pOutput);
        Aida_Fired = 0;
        Aida_Fired = pInput->fAIDAHits;

   ///-------------------------------- /**bPlastic VME Input**/ --------------------------------///
        ///Switched off A.M. 11.12.19


 ///--------------------------------------/**Scalar Input**/------------------------------------------///

//     ScalarFired = pInput->fScalar_fired;
//     for (int i = 0; i<ScalarFired; i++){
//       ScalarID = pInput->fScalar_ID;
//     }
//   }
        
///--------------------------------------/**bPlastic TAMEX Input**/------------------------------------------///
       for(int j=0; j<100;j++){
//         bPlas_TAM_SC41L_ANA[j] = 0;
//         bPlas_TAM_SC41R_ANA[j] = 0;
//         bPlas_TAM_SC41L_DIG[j] = 0;
//         bPlas_TAM_SC41R_DIG[j] = 0;
// 	bPlas_RefCh[j] = 0;
// 	//Fat_RefCh[j] = 0;
//         bPlas_AND_Coinc[j] = 0;
        }
        
       
   if (PrcID_Conv[2] ==2){
      
        pOutput->pbPLAS_WR = pInput->fbPlas_WR;
         
//      cout<< "PLAS " << pOutput->pEvent_Number<<" pOutput->pbPLAS_WR " << pOutput->pbPLAS_WR<< endl;
       
        for (int i = 0; i < 48; i++){
        for(int j=0; j<pInput->fbPlas_PMT_Lead_N[i];j++){
        bPlas_TAM_SC41L_ANA[j] = pInput->fbPlas_Lead_PMT[12][j];

        bPlas_TAM_SC41R_ANA[j] = pInput->fbPlas_Lead_PMT[11][j];
        bPlas_TAM_SC41L_DIG[j] = pInput->fbPlas_Lead_PMT[13][j];
        bPlas_TAM_SC41R_DIG[j] = pInput->fbPlas_Lead_PMT[14][j];

       // bPlas_RefCh[j] = pInput->fbPlas_Lead_PMT[16][j];
   
      //  Fat_RefCh[j] = pInput->fbPlas_Lead_PMT[0][j];
//
        bPlas_AND_Coinc[j] = pInput->fbPlas_Lead_PMT[9][j];
       
        }
    }
      Do_Plastic_Tamex_Histos(pInput,pOutput);
     
      
   }
  
  

  ///--------------------------------------/**Fatima Input**/------------------------------------------///
        FatQDCFired = 0;
        FatTDCFired = 0;
        SC41 = 0;
        SC41_ns =0;
        Fat_CHA_0_TDC = 0;
        Fat_WR = 0;
  for (int i=0; i<50; i++){
         FatQDC[i] = 0;
         FatQDC_T[i] =0;
  }
  for (int j=0; j<50; j++){
          FatTDCID[j] = -1;
          FatQDCID[j] = -1;
          FatTDC_Multipl[j] = 0;

  for (int k=0; k<50; k++){
            FatTDC_TS[j][k] = 0;

    }
  }
  
 // cout<<" event " << pOutput->pEvent_Number << " Used_Systems[3] "<< Used_Systems[3] <<" pOutput->pUsed_Systems[3] " << pOutput->pUsed_Systems[3] << endl;
 

  if (( PrcID_Conv[3]==3) &&  VMEorTAMEX_fatima==true){
    FatQDCFired =  pInput->fFat_firedQDC;
    

    Fat_WR = pInput->fFat_WR;
    //QDC
    for (int i = 0; i<FatQDCFired; i++){
      FatQDCID[i] = pInput->fFat_QDC_ID[i];
      FatQDC[i] = pInput->fFat_QDC_E[i];
      FatQDC_T[i] = pInput->fFat_QDC_T[i];
                                       }

    //TDC
    FatTDCFired =  pInput->fFat_firedTDC;

    for (int i = 0; i<FatTDCFired; i++){
      FatTDCID[i] = pInput->fFat_TDC_ID[i];
            FatTDC_TS[i][FatTDCID[i]] = (pInput->fFat_TDC_TS[i][FatTDCID[i]])*25; //in ps;

            FatTDC_Multipl[FatTDCID[i]] = pInput-> fFat_TDC_Multiplicity[FatTDCID[i]];
            ///Reference channel
             if( FatTDC_TS[i][0]>0  && FatTDCID[i] ==0){
                           Fat_CHA_0_TDC =   FatTDC_TS[i][0];
                          // pOutput->pFat_Ch0_TDC = Fat_CHA_0_TDC;
                         }

            //SC41 Trigger
      if (FatTDCID[i]==40){
            SC41 = pInput ->fSC41[i]*25; //in 25ps
            SC41_ns = SC41*0.025; //SC41 signal in ns
      }
    }

        Do_Fatima_Histos(pOutput);
  }

   ///--------------------------------------/**Fatima TAMEX Input**/------------------------------------------///
//    if (Used_Systems[4]&& PrcID_Conv[4]==4 && VMEorTAMEX_fatima==false && VMEandTAMEX_fatima==false){
//         Do_Fatima_Tamex_Histos(pInput,pOutput);
//
//    }

//      if (((Used_Systems[4]&& PrcID_Conv[4]==4) ||(Used_Systems[3]&& PrcID_Conv[3]==3)) && VMEorTAMEX_fatima==false && VMEandTAMEX_fatima==true){
//          Do_Fatima_VME_Tamex_Histos(pInput,pOutput);
//      }


 ///--------------------------------------/**Galileo Input**/------------------------------------------///
   GalFired = -1;
   //Gal_WR = 0;
   for(int g = 0; g<GALILEO_MAX_HITS; g++){
      GalDet[g] = -1;
      GalCrys[g] = -1;
      GalE[g] = -1;
      GalE_Cal[g] = -1;
      GalT[g] =-1;
      GalPileUp[g] = false;
      GalOverFlow[g] = false;
   }
//    cout<<" event " << pOutput->pEvent_Number << " Used_Systems[5] "<< Used_Systems[5] <<" pOutput->pUsed_Systems[5] " << pOutput->pUsed_Systems[5] << endl;

   if ( PrcID_Conv[5]==5){


    GalFired =  pInput->fGal_fired;
//    GalPileup = pInput->fGal_Pileup;
    Gal_WR = pInput->fGal_WR;
   
    //for(int f=0;f<1000;f++){

    //}

    for (int i = 0; i<GalFired; i++)
    {
      GalDet[i] = pInput->fGal_Detector[i];
      GalCrys[i] = pInput->fGal_Crystal[i];
      GalPileUp[i] = pInput->fGal_Pileup[i];
      GalOverFlow[i] = pInput->fGal_Overflow[i];
      GalE[i] = pInput->fGal_E[i];
      //TODO: NH - this is dumb for now to avoid cahgning calib file
      int id = GalDet[i] * 3 + GalCrys[i];
      GalE_Cal[i] = fCal->AGal[id]* pow( GalE[i],2) + fCal->BGal[id]*  GalE[i] + fCal->CGal[id];
    }

        Do_Galileo_Histos(pOutput);
  }

 ///--------------------------------------/**Finger Input**/------------------------------------------///
        Fing_firedTamex = -1;
        for (int i = 0; i<4; i++)
        {
            Fing_leadHits[i] = -1;
            Fing_trailHits[i] = -1;
            Fing_iterator[i] = -1;
            Fing_trig[i] = 0;

        for (int j = 0; j<32; j++){
            Fing_tamex_ch[i][j] = -1;
            Fing_leadChan[i][j] = -1;
            Fing_leadT[i][j] = 0;
            Fing_trailChan[i][j] = -1;
            Fing_trailT[i][j] = 0;
            Fing_TOT[i][j] = 0;
            Fing_TOT_added[i][j] = 0;
            Fing_chID[i][j] = -1;
            Fing_lead_coarse[i][j] =  -1;
            Fing_lead_fine[i][j]  =-1;
            Fing_trail_coarse[i][j] = -1;
            Fing_trail_fine[i][j] = -1;
            }
        }
    if(Used_Systems[6] && PrcID_Conv[6]==6){
//           Fing_firedTamex = pInput->ffing_tamexhits;
//           maxToT = Fing_TOT[0][0];
//           maxToT_added = Fing_TOT_added[0][0];
//
//           for (int i=0; i< Fing_firedTamex; i++){
//             Fing_leadHits[i] = pInput-> ffing_leadHits[i];
//             Fing_trailHits[i] = pInput-> ffing_trailHits[i];
//             Fing_iterator[i] = pInput-> ffing_iterator[i];
//             Fing_trig[i] = pInput->ffing_Trig[i];
//
//                 for (int j =0; j<Fing_iterator[i]; j++){
//                     Fing_tamex_ch[i][j] =  pInput->ffing_tamexCh[i][j];
//                     Fing_lead_coarse[i][j] =  pInput->ffing_lead_coarse[i][j];
//                     Fing_lead_fine[i][j]  = pInput->ffing_lead_fine[i][j];
//                     Fing_trail_coarse[i][j] =  pInput->ffing_trail_coarse[i][j];
//                     Fing_trail_fine[i][j] =  pInput->ffing_trail_fine[i][j];
//
//                     Fing_chID[i][j] = pInput->ffing_chID[i][j];
//                   //  cout << "  Fing_chID[i][j] " << Fing_chID[i][j]<< " i " << i << " j " << j << " ffing_Lead_T[i][j] " << pInput->ffing_Lead_T[i][j] <<endl;
//                     if(Fing_chID[i][j] % 2 == 1){
//                         Fing_leadChan[i][j] = pInput->ffing_Lead_Phys_Chan[i][j];
//                         Fing_leadT[i][j] = pInput->ffing_Lead_T[i][j];
//
//                         }
//
//                 else{
//                      Fing_trailChan[i][j] = pInput->ffing_Trail_Phys_Chan[i][j];
//                      Fing_trailT[i][j] = pInput->ffing_Trail_T[i][j];
//
//                 }
//
//                 ///Note: ToT value here is only for the 'up' PMTs
//                     Fing_TOT[i][j] = pInput->ffing_TOT[i][j];
//                    // cout <<"1) ev "<< event_number << " Fing_TOT[i][j] " << Fing_TOT[i][j] <<"  Fing_leadChan[i][j] " <<  Fing_leadChan[i][j] << endl;
//                     /// Up PMT ToT + Down PMT ToT
//                     Fing_TOT_added[i][j] =   pInput->ffing_TOT_added[i][j];
//
//
//                         if(Fing_chID[i][j] % 2 == 1){
//
//                           if(maxToT<Fing_TOT[i][j] && Fing_leadChan[i][j]>0 ){
//
//                                 maxToT = Fing_TOT[i][j];
//                                 maxToTChan = Fing_leadChan[i][j];
//                         }
//
//                         //Get max ToT for PMT pairings
//                         if(maxToT_added<Fing_TOT_added[i][j]&& Fing_leadChan[i][j]>0 ){
//                             maxToT_added = Fing_TOT_added[i][j];
//                             maxToT_added_Chan = Fing_leadChan[i][j];
//
//                         }
//
//
//                 }
//             }
//           }
   Do_Finger_Histos(pInput,pOutput);

  }

///------------------------/**Setup for some temporary correlations**/----------------------------------///
//   if (Used_Systems[2]&& PrcID_Conv[2]==2 && Used_Systems[3]&& PrcID_Conv[3]==3){
//             Do_Fat_Plas_Histos(pOutput);
//             }
//   if (Used_Systems[6]&& PrcID_Conv[6]==6 && Used_Systems[2]&& PrcID_Conv[2]==2){
//     Do_Fing_Plas_Histos(pOutput);
//             }
                        /** End of Unpack Tree input**/

  BuildEvent2(dest);
  pOutput->SetValid(isValid);
  return isValid;

 }  //End of BuildEvent

    Bool_t EventAnlProc::BuildEvent2(TGo4EventElement* dest){
//        EventUnpackStore* pInput  = (EventUnpackStore*) GetInputEvent();
//         int event_number1 = pInput -> fevent_number;
//         int event_counter = pInput -> fArray_count;

       //cout <<"event " << event_number << " WR " <<  pInput ->fWR_main_array[event_counter] << " event_counter " << event_counter << endl;
      return(1);
        }
///End of Input from Unpacker ///
 ///-----------------------------------------------------------------------------------------------------------------///
   void EventAnlProc::get_used_Systems(){
    for(int i = 0;i < 7;i++) Used_Systems[i] = false;

    ifstream data("Configuration_Files/Used_Systems.txt");
    if(data.fail()){
        cerr << "Could not find Used_Systems config file!" << endl;
        exit(0);
    }
    int i = 0;
    int id = 0;
    string line;
    char s_tmp[100];
    while(data.good()){
        getline(data,line,'\n');
        if(line[0] == '#') continue;
        sscanf(line.c_str(),"%s %d",s_tmp,&id);
        Used_Systems[i] = (id == 1);
        i++;
    }
    string DET_NAME[7] = {"FRS","AIDA","PLASTIC","FATIMA_VME","FATIMA_TAMEX","GALILEO","FINGER"};

//     cout << "\n=====================================================" << endl;
//     cout << "USED SYSTEMS" << endl;
//     cout << "-----------------------------------------------------" << endl;
//     for(int j = 0;j < 6;++j){
//         if(Used_Systems[j]) cout << DET_NAME[j] << endl;
//     }
//     cout << "=====================================================" << endl;


}
///-----------------------------------------------------------------------------------------------------------------///
  void EventAnlProc::read_setup_parameters(){

    // unused // const char* format = "%s %d";

    ifstream file("Configuration_Files/Detector_System_Setup_File.txt");

    if(file.fail()){
        cerr << "Could not find File for setup parameters!" << endl;
        exit(0);
    }

    string line;
    string var_name;
    // unused //int dummy_var;
    //file.ignore(256,'GENERAL_CONFIGURATION');

    file.ignore(256,':');
    file >> FAT_exclusion_dist;//dummy_var;

    file.ignore(256,':');
    file >> FAT_nearest_neighbour_exclusion;//dummy_var;

    file.ignore(256,':');
    file >> same_ring_exclusion;//dummy_var;

    file.ignore(256,':');
    file >> output_position_matrix;//dummy_var;

    cout<<endl;
    cout<<endl;
    cout<<"////////////////////////////////////////////////////////////////////////"<<endl;
    cout<<"Setup Parameters List Analysis Proc: "<<endl;
    if(FAT_exclusion_dist > 0) cout<<"FATIMA Detectors Excluded if Linear Difference Exceeds "<<FAT_exclusion_dist<<" mm"<<endl;
    else if(FAT_exclusion_dist == 0) cout<<"'Nearest Neighbour Exclusion': Disabled (Distance set to 0)"<<endl;
    cout<<"////////////////////////////////////////////////////////////////////////"<<endl;
    cout<<endl;
    cout<<endl;
    /*while(file.good()){
        getline(file,line,'\n');
        if(line[0] == '#') continue;
        sscanf(line.c_str(),format,&var_name,&dummy_var);

        cout<<"Hello Again?"<<endl;

        if (var_name == "White_Rabbit_Enabled:" && dummy_var == 1)  WHITE_RABBIT_USED = true;
        else if (var_name == "White_Rabbit_Enabled:" && dummy_var == 0)  WHITE_RABBIT_USED = false;

        if (var_name == "FATIMA_Gain_Match_Enabled:" && dummy_var == 1)  FAT_gain_match_used = true;
        else if (var_name == "FATIMA_Gain_Match_Enabled:" && dummy_var == 0) FAT_gain_match_used  = false;

    }*/

}
///-----------------------------------------------------------------------------------------------------------------///
void EventAnlProc::load_GalileoMap_File(){

}
/**----------------------------------------------------------------------------------------------**/
 /**--------------------------------------    White Rabbit   ---------------------------------------------**/
 /**----------------------------------------------------------------------------------------------**/
 void EventAnlProc::Make_WR_Histos(){
   hAida_Fat_WRdT = MakeTH1('I',"WR/Aida_Fatima_dT","White Rabbit Aida-Fatima",1000,-1000,1000);
   hAida_Gal_WRdT = MakeTH1('I',"WR/Aida_Galileo_dT","White Rabbit Aida-Galileo",1000,-1000,1000);
   hAida_bPlas_WRdT = MakeTH1('I',"WR/Aida_bPlas_dT","White Rabbit Aida-bPlas",1000,-1000,1000);
   hbPlas_Fat_WRdT = MakeTH1('I',"WR/bPlas_Fatima_dT","White Rabbit bPlas-Fatima",1000,-1000,1000);
   hbPlas_Gal_WRdT = MakeTH1('I',"WR/bPlas_Galileo_dT","White Rabbit bPlas-Galileo",1000,-1000,1000);
   hFat_Gal_WRdT = MakeTH1('I',"WR/Fatima_Galileo_dT","White Rabbit Fatima_Galileo",1000,-1000,1000);



}
 void EventAnlProc::Do_WR_Histos(EventUnpackStore* pInputMain){
   for(AidaUnpackData& pInputD : pInputMain->Aida)
  {
AidaUnpackData* pInput = &pInputD;
  hAida_Fat_WRdT->Fill(pInputMain->fAIDA_WR -pInputMain->fFat_WR );
  hAida_Gal_WRdT->Fill(pInputMain->fAIDA_WR -pInputMain->fGal_WR );
  hAida_bPlas_WRdT->Fill(pInputMain->fAIDA_WR -pInputMain->fbPlas_WR);
  hbPlas_Fat_WRdT->Fill(pInputMain->fbPlas_WR -pInputMain->fFat_WR);
  hbPlas_Gal_WRdT->Fill(pInputMain->fbPlas_WR - pInputMain->fGal_WR);
  hFat_Gal_WRdT->Fill(pInputMain->fFat_WR - pInputMain->fGal_WR);

 // cout<<"pInputMain->fAIDA_WR " << pInputMain->fAIDA_WR<<endl;
  }
 }


/**----------------------------------------------------------------------------------------------**/
 /**--------------------------------------    FRS   ---------------------------------------------**/
 /**----------------------------------------------------------------------------------------------**/
void EventAnlProc::Make_FRS_Histos(){
    // hFRS_z1_z2 = MakeTH2('D',"FRS/z1vsz2","FRS Z1 vs Z2", 500, 0, 5000,  500, 0, 5000);
    hID_x2AoQ = MakeH2I("FRS/ID","ID_x2AoQ", 300,1.4,2.5, 200,-100.,100.,"A/Q s2-s4", "X at S2 [mm]", 2);
    hID_x4AoQ = MakeH2I("FRS/ID","ID_x4AoQ", 300,1.4,2.5, 200,-100.,100.,"A/Q s2-s4", "X at S4 [mm]", 2);

    hID_Z_AoQ = MakeH2I("FRS/ID","ID_Z_AoQ", 1500,1.3,2.8, 1500,35.,90.,"A/Q s2-s4", "Z s2-s4", 2);
    //hID_Z_AoQ_zsame = MakeH2I("FRS/ID","ID_Z_AoQ_zsame", 600,1.4,2.5, 600,1.,30.,"Z1==Z2 A/Q s2-s4", "Z s2-s4", 2);
    hID_Z_AoQ_corr = MakeH2I("FRS/ID","ID_Z_AoQ_S2-S4corr", 1500,1.3,2.8, 1500,35.,90., "A/Q s2-s4", "Z s2-s4", 2);

    hID_Z_Z2 = MakeH2I("FRS/ID","ID_Z_Z2", 2000,35,100, 2000,35.,100.,"Z", "Z2", 2);

   int num_ID_x2AoQ[3] = {5, 5, 5};
  Float_t init_ID_x2AoQ[3][5][2] =
     {
       {{  2.114 ,   102.2},
    {  2.111,  -98.6},
    {  2.49673,  -98.3},
    {  2.51517,  103.09},
    {  2.114,  102.19}
       },
    {{  2.009 ,   87.121 },
    {  2.0113,  -82.0475},
    {  2.3872,  -82.3590},
    {  2.3872,  -82.3590},
    {  2.3886,   85.8756}
       },
       {{  2.009 ,   87.121 },
    {  2.0113,  -82.0475},
    {  2.3872,  -82.3590},
    {  2.3872,  -82.3590},
    {  2.3886,   85.8756}
       }
     };

       int num_ID_x4AoQ[3] = {5, 5, 5};
  Float_t init_ID_x4AoQ[3][5][2] =
     {
       {{  2.118 ,   99.2 },
    {  2.12051,  -69.777},
    {  2.5,  -70.079},
    {  2.5,  97.98},
    {  2.118,   99.2}
       },
    {{  2.009 ,   87.121 },
    {  2.0113,  -82.0475},
    {  2.3872,  -82.3590},
    {  2.3872,  -82.3590},
    {  2.3886,   85.8756}
       },
       {{  2.009 ,   87.121 },
    {  2.0113,  -82.0475},
    {  2.3872,  -82.3590},
    {  2.3872,  -82.3590},
    {  2.3886,   85.8756}
       }
     };

       int num_ID_Z_Z2[3] = {5, 5, 5};
  Float_t init_ID_Z_Z2[3][5][2] =
     {
       {{  12 ,   4 },
    {  13.8,  8},
    {  14.6,  8},
    {  14.6,  4},
    {  14.4,   4}
       },
    {{  2.009 ,   87.121 },
    {  2.0113,  -82.0475},
    {  2.3872,  -82.3590},
    {  2.3872,  -82.3590},
    {  2.3886,   85.8756}
       },
       {{  2.009 ,   87.121 },
    {  2.0113,  -82.0475},
    {  2.3872,  -82.3590},
    {  2.3872,  -82.3590},
    {  2.3886,   85.8756}
       },
     };


  char name[50], title[100];

  for(int i=0;i<3;++i)
    {
      sprintf(name,"cID_x2AoQ%d",i);
      cID_x2AoQ[i] = MakePolyCond("FRS_ID_Gated",name,num_ID_x2AoQ[i],init_ID_x2AoQ[i], hID_x2AoQ->GetName());

      sprintf(name,"cID_x4AoQ%d",i);
      cID_x4AoQ[i] = MakePolyCond("FRS_ID_Gated",name,num_ID_x4AoQ[i],init_ID_x4AoQ[i], hID_x4AoQ->GetName());

      sprintf(name,"cID_Z_Z2%d",i);
      cID_Z_Z2gate[i] = MakePolyCond("FRS_Z1_Z2_Gated",name,num_ID_Z_Z2[i],init_ID_Z_Z2[i], hID_Z_Z2 ->GetName());
      //cID_Z_AoQ[i] = MakePolyCond("ID", name, num_ID_Z_AoQ[i], init_ID_Z_AoQ[i], hID_Z_AoQ->GetName());

      sprintf(name,"hID_x2AoQ_x2AoQgate%d",i);
      hID_x2AoQ_x2AoQgate[i] = MakeH2I("FRS/ID_Gated", name, 300,1.,2.4, 200,-100.,100.,"A/Q s2-s4", "gate on FRS AoQ, ID X2: X at S2 [mm]", 2);

      sprintf(name,"hID_x4AoQ_x2AoQgate%d",i);
      hID_x4AoQ_x2AoQgate[i] = MakeH2I("FRS/ID_Gated", name, 300,1.,2.4, 200,-100.,100.,"A/Q s2-s4", "gate on FRS AoQ, ID X2: X at S4 [mm]", 2);

      sprintf(name,"hID_x2AoQ_x4AoQgate%d",i);
      hID_x2AoQ_x4AoQgate[i] = MakeH2I("FRS/ID_Gated", name, 300,1.,2.4, 200,-100.,100.,"A/Q s2-s4", "gate on FRS AoQ, ID X4: X at S2 [mm]", 2);

      sprintf(name,"hID_x4AoQ_x4AoQgate%d",i);
      hID_x4AoQ_x4AoQgate[i] = MakeH2I("FRS/ID_Gated", name, 300,1.,2.4, 200,-100.,100.,"A/Q s2-s4", "gate on FRS AoQ, ID X4: X at S4 [mm]", 2);

      sprintf(name,"ID_ZAoQ_x2AoQgate%d",i);
      hID_ZAoQ_x2AoQgate[i] = MakeH2I("FRS/ID_Gated", name, 300,1.,2.4, 400,30.,90.,"A/Q s2-s4", " Z music2", 2);

      sprintf(name,"ID_ZAoQ_x4AoQgate%d",i);
      hID_ZAoQ_x4AoQgate[i] = MakeH2I("FRS/ID_Gated", name, 300,1.,2.4, 400,30.,90.,"A/Q s2-s4", " Z music2", 2);

      ///////////////////////////////////////////////////////
      sprintf(name,"ID_x2AoQ_Z1Z2gate%d",i);
      hID_x2AoQ_Z1Z2gate[i] = MakeH2I("FRS/ID_Gated", name, 300,1.,2.4, 200,-100.,100.,"A/Q s2-s4", "gate on Z    X at S2 [mm]", 2);

      sprintf(name,"ID_x4AoQ_Z1Z2gate%d",i);
      hID_x4AoQ_Z1Z2gate[i] = MakeH2I("FRS/ID_Gated", name, 300,1.,2.4, 200,-100.,100.,"A/Q s2-s4", "gate on Z    X at S4 [mm]", 2);

      sprintf(name,"ID_ZAoQ_Z1Z2gate%d",i);
      hID_ZAoQ_Z1Z2gate[i] = MakeH2I("FRS/ID_Gated", name, 300,1.4,2.5, 400,1.,20.,"A/Q s2-s4", " Z music2", 2);

      sprintf(name,"ID_SC43Z1_Z1Z2gate%d",i);
      hID_SC43Z1_Z1Z2gate[i] = MakeH2I("FRS/ID_Gated", name, 300,1.,2.4, 400,30.,90.,"SC41 dE", " Z music1", 2);


    }
      //------------------------------------------------------------------------------------------------//

      int num_ID_Z_AoQ[3] = {5, 5, 5};
  Float_t init_ID_Z_AoQ[3][5][2] =
     {
       // ID_Z_AOQ(1)
       /* 34Si setting */
//        {{2.30,        13.36},
//         {2.33,        14.22},
//         {2.239,        15.10},
//         {2.18471,        14.1637},
//         {2.3011,        13.3654}},
// /* 168Re setting */
        {{2.236,        77.6},
        {2.23389,        77.45},
        {2.23322,        76.4},
        {2.247,        76.5},
        {2.247,        77.5}},
	
// 	/* 168W setting */
//         {{2.248,        76.5},
//         {2.236,        76.6},
//         {2.236,        75.5},
//         {2.25,        75.6},
//         {2.25,        76.5}},
	
 	/* 168Ta setting */
//         {{2.248,        75.5},
//         {2.236,        75.6},
//         {2.236,        74.5},
//         {2.25,        74.5},
//         {2.25,        75.5}},
	
       /* 223Th setting */
       //      {{2.48025,        89.3017},
       //      {2.49183,        89.2808},
       //      {2.49122,        90.7808},
       //      {2.48005,        90.8122},
       //      {2.48025,        89.3017}},
       /* 70 Ni @ 70Ni setting shifted by 10*/
       {{2.482,        38.0000},
        {2.508,        38.5000},
        {2.522,        38.0000},
        {2.489,        37.5000},
        {2.488,        37.5000}},
       // ID_Z_AOQ(3)
       {{2.25029,      51.22740},
        {2.27401,      50.88124},
        {2.30980,      51.80895},
        {2.28475,      51.93358},
        {2.25208,      52.07204}}};
     for(int i=0;i<3;i++){
       sprintf(name,"cID_Z_AoQ%d",i);
       cID_Z_AoQ[i] = MakePolyCond("FRS_ID_Gated", name, num_ID_Z_AoQ[i], init_ID_Z_AoQ[i], hID_Z_AoQ->GetName());

       sprintf(name,"ID_Z_AoQgate%d",i);
       hID_Z_AoQgate[i] = MakeH2I("FRS/ID_Gated",name,  600,1.3,2.8, 600,35.,90.,"A/Q s2-s4", "Z s2-s4", 2);

//        sprintf(name,"ID_Z3_gate%d",i);
//        char name_x_title[256];
//        sprintf(name_x_title,"Z3 gate%d",i);
//        hID_Z3_gate[i] = MakeH1I("FRS/ID_Gated",name,  500,1.,20.,name_x_title, 2, 6);
    }
}

void EventAnlProc::Do_FRS_Histos(EventAnlStore* pOutput){
  ////WARNING!!!!!!!!!!!!!!!!!!
    FRS_AoQ = FRS_AoQ_corr;
    
     pOutput->pFRS_AoQ = FRS_AoQ;
     pOutput->pFRS_ID_x2 = FRS_ID_x2;
     pOutput->pFRS_ID_x4 = FRS_ID_x4;
     pOutput->pFRS_z = FRS_z;
     pOutput->pFRS_z2 = FRS_z2;
     for(int l=0; l<12;l++){
     pOutput->pSci_num = l;
     pOutput->pFRS_sci_l[l] = FRS_sci_l[l];
     pOutput->pFRS_sci_r[l] = FRS_sci_r[l];
     pOutput->pFRS_sci_e[l] = FRS_sci_e[l];
     pOutput->pFRS_sci_tx[l] = FRS_sci_tx[l];
     pOutput->pFRS_sci_x[l] = FRS_sci_x[l];
     
     }
     
    ///AoQ vs X
    
     hID_x2AoQ->Fill(FRS_AoQ, FRS_ID_x2);
     hID_x4AoQ->Fill(FRS_AoQ, FRS_ID_x4);

     ///AoQ vs Z
     hID_Z_AoQ->Fill(FRS_AoQ, FRS_z);
     hID_Z_AoQ_corr->Fill(FRS_AoQ_corr, FRS_z);  //S2-S4 correlated
     hID_Z_Z2->Fill(FRS_z,FRS_z2);
//        if(TMath::Abs(FRS_z-FRS_z2-0.3)<0.6)
//             {
//               hID_Z_AoQ_zsame->Fill(FRS_AoQ, FRS_z);
//             }

    for(int i=0;i<3;++i)
    {
        ///GATE: ID vs x2AoQ
      if(cID_x2AoQ[i]->Test(FRS_AoQ, FRS_ID_x2)==true)
        {
          pOutput->pFRS_x2AoQ_pass = true;
           if(FRS_AoQ>2 && FRS_AoQ<2.4 &&   FRS_ID_x2 > -100 && FRS_ID_x2<100){

             //  cout<<"FRS_AoQ " << FRS_AoQ << " FRS_ID_x2 " << FRS_ID_x2 <<endl;
          hID_x2AoQ_x2AoQgate[i]->Fill(FRS_AoQ, FRS_ID_x2);


        }
          if(FRS_AoQ>2 && FRS_AoQ<2.4 &&   FRS_ID_x4 > -100 && FRS_ID_x4<100){
          hID_x4AoQ_x2AoQgate[i]->Fill(FRS_AoQ, FRS_ID_x4);}
          if (FRS_z2) hID_ZAoQ_x2AoQgate[i]->Fill(FRS_AoQ, FRS_z2);
            }

        ///GATE: ID vs x4AoQ
      if(cID_x4AoQ[i]->Test(FRS_AoQ, FRS_ID_x4)==true)
        {
                pOutput->pFRS_x4AoQ_pass = true;
                hID_x2AoQ_x4AoQgate[i]->Fill(FRS_AoQ, FRS_ID_x2);
                hID_x4AoQ_x4AoQgate[i]->Fill(FRS_AoQ, FRS_ID_x4);
          if (FRS_z2) hID_ZAoQ_x4AoQgate[i]->Fill(FRS_AoQ, FRS_z2);
            }

          ///GATE: Z1 vs Z2
            if(cID_Z_Z2gate[i]->Test(FRS_z, FRS_z2)==true)
            {
                 pOutput->pFRS_Z_Z2_pass = true;
              // cout<<"111event " <<  pOutput->pEvent_Number << " ANL pOutput->pFRS_Z_Z2_pass " << pOutput->pFRS_Z_Z2_pass<< endl;

                 if(FRS_AoQ>2 && FRS_AoQ<2.4 &&   FRS_ID_x2 > -100 && FRS_ID_x2<100){
                    // cout<<"FRS_AoQ " << FRS_AoQ << " FRS_ID_x2 " << FRS_ID_x2 <<endl;
                hID_x2AoQ_Z1Z2gate[i]->Fill(FRS_AoQ, FRS_ID_x2);
                }
                   if(FRS_AoQ>2 && FRS_AoQ<2.4   && FRS_ID_x4 > -100 && FRS_ID_x4<100){
                        //cout<<"FRS_AoQ " << FRS_AoQ << " FRS_ID_x4 " << FRS_ID_x4 <<endl;
                hID_x4AoQ_Z1Z2gate[i]->Fill(FRS_AoQ, FRS_ID_x4);
                   }
                   if(FRS_AoQ>2 && FRS_AoQ<2.4){
                   hID_ZAoQ_Z1Z2gate[i] ->Fill(FRS_AoQ, FRS_z);
                   }
                 //if(FRS_sci_e[7]>1 && FRS_sci_e[7]<2.4 &&  && FRS_z>30 && FRS_z<10){
                hID_SC43Z1_Z1Z2gate[i]->Fill(FRS_sci_e[7], FRS_z);
          if (FRS_z2) hID_ZAoQ_x2AoQgate[i]->Fill(FRS_AoQ, FRS_z2);
            }
          //  else pOutput->pFRS_Z_Z2_pass = false;

         //    cout<<"222event " << pOutput->pEvent_Number << " ANL pOutput->pFRS_Z_Z2_pass " << pOutput->pFRS_Z_Z2_pass<< endl;
        }
     for(int i =0;i<3;i++){
        ///GATE: AoQ vs Z
       if( cID_Z_AoQ[i]->Test(FRS_AoQ, FRS_z)==true){
          pOutput->pFRS_ZAoQ_pass =true;
          hID_Z_AoQgate[i]->Fill(FRS_AoQ, FRS_z);
        }
       }
    }


/**----------------------------------------------------------------------------------------------**/
/**--------------------------------------    AIDA   ---------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/

void EventAnlProc::Make_Aida_Histos(){
  TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
  implants_strip_xy.resize(conf->DSSDs());
  implants_pos_xy.resize(conf->DSSDs());
  implants_e.resize(conf->DSSDs());
  implants_e_xy.resize(conf->DSSDs());
  implants_time_delta.resize(conf->DSSDs());
  implants_strip_1d.resize(conf->DSSDs());
  implants_per_event.resize(conf->DSSDs());
  decays_strip_xy.resize(conf->DSSDs());
  decays_pos_xy.resize(conf->DSSDs());
  decays_e.resize(conf->DSSDs());
  decays_e_xy.resize(conf->DSSDs());
  decays_time_delta.resize(conf->DSSDs());
  decays_strip_1d.resize(conf->DSSDs());
  decays_per_event.resize(conf->DSSDs());
  implants_channels.resize(conf->DSSDs());
  decays_channels.resize(conf->DSSDs());
  implants_x_ex.resize(conf->DSSDs());
  implants_y_ey.resize(conf->DSSDs());
#ifdef AIDA_PULSER_ALIGN
  aida_pulser_time = MakeTH2('I', "AIDA/Pulser_Time", "AIDA Pulser Time Comparison", 768, 0, 768, 2000, -4000, 4000);
#endif
  for (int i = 0; i < conf->DSSDs(); ++i)
  {
    implants_strip_xy[i] = MakeTH2('I', Form("AIDA/Implants/DSSD%d_implants_strip_XY", i+1), Form("DSSD %d implant hit pattern", i+1), 128, 0, 128, 128, 0, 128, "X strip", "Y strip");
    //implants_pos_xy[i] = MakeTH2('D', Form("AIDA/Implants/DSSD%d_implants_pos_XY", i+1), Form("DSSD %d implant position", i+1), 128, -37.8, 37.8, 128, -37.8, 37.8, "X position/mm", "Y position/mm");
    implants_e[i] = MakeTH1('F', Form("AIDA/Implants/DSSD%d_implants_energy", i+1), Form("DSSD %d implant energy", i+1), 1000, 0, 10000, "Implant Energy/MeV");
    //implants_e_xy[i] = MakeTH2('F', Form("AIDA/Implants/DSSD%d_implants_energy_XY", i+1), Form("DSSD %d implant front energy vs back energy", i+1), 1000, 0, 10000, 1000, 0, 10000, "X Energy", "Y Energy");
    implants_time_delta[i] = MakeTH1('F', Form("AIDA/Implants/DSSD%d_implants_time_delta", i+1), Form("DSSD %d implant front vs back time", i+1), 1000, -10000, 10000, "Time Difference/ns");
    implants_strip_1d[i] = MakeTH1('I', Form("AIDA/Implants/DSSD%d_implants_strip_1d", i+1), Form("DSSD %d implant 1D hit pattern", i+1), 256, 0, 256, "Strip number");
    implants_per_event[i] = MakeTH1('I', Form("AIDA/Implants/DSSD%d_implants_per_event", i+1), Form("DSSD %d implants per event", i+1), 100, 0, 100, "Number of implants");
    implants_x_ex[i] = MakeTH2('F', Form("AIDA/Implants/DSSD%d_implants_x_ex", i+1), Form("DSSD %d Ex vs X position", i+1), 128, 0, 128, 1000, 0, 10000, "X Strip", "X Energy");
    implants_y_ey[i] = MakeTH2('F', Form("AIDA/Implants/DSSD%d_implants_y_ey", i+1), Form("DSSD %d Ey vs Y position", i+1), 128, 0, 128, 1000, 0, 10000, "Y Strip", "Y Energy");

    decays_strip_xy[i] = MakeTH2('I', Form("AIDA/Decays/DSSD%d_decays_strip_XY", i+1), Form("DSSD %d decay hit pattern", i+1), 128, 0, 128, 128, 0, 128, "X strip", "Y strip");
   // decays_pos_xy[i] = MakeTH2('D', Form("AIDA/Decays/DSSD%d_decays_pos_XY", i+1), Form("DSSD %d decay position", i+1), 128, -37.8, 37.8, 128, -37.8, 37.8, "X position/mm", "Y position/mm");
    decays_e[i] = MakeTH1('F', Form("AIDA/Decays/DSSD%d_decays_energy", i+1), Form("DSSD %d decay energy", i+1), 1000, 0, 20000, "Decay Energy/keV");
    //decays_e_xy[i] = MakeTH2('F', Form("AIDA/Decays/DSSD%d_decays_energy_XY", i+1), Form("DSSD %d decay front energy vs back energy", i+1), 1000, 0, 10000, 1000, 0, 20000, "X Energy", "Y Energy");
    decays_time_delta[i] = MakeTH1('F', Form("AIDA/Decays/DSSD%d_decays_time_delta", i+1), Form("DSSD %d decay front vs back time", i+1), 1000, -10000, 10000, "Time Difference/ns");
    decays_strip_1d[i] = MakeTH1('I', Form("AIDA/Decays/DSSD%d_decays_strip_1d", i+1), Form("DSSD %d decay 1D hit pattern", i+1), 256, 0, 256, "Strip number");
    decays_per_event[i] = MakeTH1('I', Form("AIDA/Decays/DSSD%d_decays_per_event", i+1), Form("DSSD %d decays per event", i+1), 100, 0, 100, "Number of decays");

    implants_channels[i] = MakeTH1('I', Form("AIDA/DSSD%d_implants_channels", i+1), Form("DSSD %d number of implant channels", i+1), 769, 0, 769);
    decays_channels[i] = MakeTH1('I', Form("AIDA/DSSD%d_decays_channels", i+1), Form("DSSD %d number of decay channels", i+1), 769, 0, 769);
  }
}

///-----------------------------------------------------------------------------------------------------------------------------------------------------------------------///

void EventAnlProc::ProcessAida(EventUnpackStore* pInputMain, EventAnlStore* pOutput) {
 // int Aida_hits =0;
//       double bPlasQDCGainMatch_AIDA[32] ={0};
  TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();


  for(AidaUnpackData& pInputD : pInputMain->Aida)
  {

    AidaUnpackData* pInput = &pInputD;
    //pOutput->pAida = pAida;
    pOutput->pAIDA_WR = pInputMain->fAIDA_WR;
    if (pInput->ImplantEvents.size() > 1)
    {

      //         cout << " pInput->ImplantEvents.size() " << pInput->ImplantEvents.size() <<  endl;

      implantEvents++;


      // Cluster events on adjecent strips into one
      std::vector<AidaCluster> clusters = EventsToClusters(pInput->ImplantEvents);
      //
      //     // Match front-back clusters which define physical hits on the detector
      std::vector<std::pair<AidaCluster, AidaCluster>> hits;
      //
      //
      // Track the number of implants in every DSSD
      std::vector<int> counts(conf->DSSDs(), 0);
      int max_dssd = 0;

      for (auto& i : clusters)
      {
        if (i.DSSD == -1) continue;

        counts[i.DSSD - 1]++;
        if (i.DSSD > max_dssd) max_dssd = i.DSSD;

        if(i.Side != conf->DSSD(i.DSSD -1).XSide) continue;

        for (auto& j : clusters)
        {
          if(j.DSSD != i.DSSD || j.Side != conf->DSSD(j.DSSD -1).YSide) continue;
          /// Gates (set in TAidaConfiguration)
          if (abs(i.Energy - j.Energy) < conf->FrontBackEnergyH() && i.IsGoodTime(j, conf->FrontBackWindow()))
          {
            hits.push_back({i, j});
          }
        }
      }

      int channels[768] = {0};
      for (auto& i : pInput->ImplantEvents)
      {
        channels[i.Module * 64 + i.Channel]++;
      }
      int channelM = 0;
      for (int i = 0; i < 768; ++i){
        if (channels[i]) ++channelM;
      }
      implants_channels[0]->Fill(channelM);

      //     // Generate stored data for hits and plot the histograms

      for (auto& i : hits)
      {
        AidaHit hit = ClusterPairToHit(i);

        hit.Stopped = (hit.DSSD == max_dssd);

        // Check that every DSSD before has at least one implant event
        for(int j = hit.DSSD - 1; j > 0; j--)
        {
            if (counts[j - 1] == 0) hit.Stopped = false;
        }
        //
        pInput->Implants.push_back(hit);

        pOutput->pAida.Implants.push_back(hit);

        implants_strip_xy[hit.DSSD - 1]->Fill(hit.StripX, hit.StripY);
       // implants_pos_xy[hit.DSSD - 1]->Fill(hit.PosX, hit.PosY);
        implants_e[hit.DSSD - 1]->Fill(hit.Energy);
    //   cout<<"ANALYSIS AIDA " <<pOutput->pEvent_Number<< " Energy " <<  hit.Energy << endl;

      //  implants_e_xy[hit.DSSD - 1]->Fill(hit.EnergyFront, hit.EnergyBack);
        //implants_time_delta[hit.DSSD - 1]->Fill(hit.TimeFront - hit.TimeBack);
        implants_time_delta[hit.DSSD - 1]->Fill(hit.FastTimeFront - hit.FastTimeBack);

        implants_x_ex[hit.DSSD - 1]->Fill(hit.StripX, hit.EnergyFront);
        implants_y_ey[hit.DSSD - 1]->Fill(hit.StripY, hit.EnergyBack);

        int channel = i.first.Strip;
        implants_strip_1d[hit.DSSD - 1]->Fill(channel);
        channel = i.second.Strip + 128;
        implants_strip_1d[hit.DSSD - 1]->Fill(channel);
      }
      //
      if (hits.size() > 0)
      {
        implants_per_event[0]->Fill(hits.size());
      }
    }
    else if (pInput->DecayEvents.size() > 1)
    {
      decayEvents++;

      int channels[768] = {0};
#ifdef AIDA_PULSER_ALIGN
      int64_t wr_base = 0;
#endif
      for (auto& i : pInput->DecayEvents)
      {
        channels[i.Module * 64 + i.Channel]++;
#ifdef AIDA_PULSER_ALIGN
        if(i.Module == 0 && i.Channel == 0) wr_base = i.Time;
#endif
      }
      int channelM = 0;
      for (int i = 0; i < 768; ++i)
      if (channels[i]) ++channelM;
      // cout  << " channelM " << channelM<<"  channels[i] " << channels[i]<< endl;}
      decays_channels[0]->Fill(channelM);

      if (channelM > 400)
      {
        decayEvents--;
        pulserEvents++;
#ifdef AIDA_PULSER_ALIGN
        if(channelM > 750)
        {
          std::cout << "Identified a pulser event!" << std::endl;
          std::vector<int> offsets;
          offsets.resize(conf->FEEs());
          for (int i = 0; i < offsets.size(); i++) offsets[i] = 0;

          for (auto& i : pInput->DecayEvents)
          {
            if (i.Energy < 1000) continue;
            std::cout << i.Module << ", " << i.Channel << ", " << i.Energy << ", " << std::hex << i.Time << "," << i.FastTime << std::dec << std::endl;
            int offset = (i.Time - wr_base) % 2000;
            if (offsets[i.Module] == 0)
            {
              offsets[i.Module] = offset;
            }
            else if (offset > offsets[i.Module])
            {
              // confirm the offset is 2us out
              if (abs(offset - offsets[i.Module]) % 2000 != 0)
                std::cout << "LOGICAL MISTAKE IN ALIGNMENT" << std::endl;
            }
            else if (offset < offsets[i.Module])
            {
              if (abs(offset - offsets[i.Module]) % 2000 != 0)
                std::cout << "LOGICAL MISTAKE IN ALIGNMENT" << std::endl;
              offsets[i.Module] = offset;
            }

            aida_pulser_time->Fill(i.Module * 64 + i.Channel, offsets[i.Module]);
          }

          std::cout << std::endl;

          std::cout << "Put this into AIDA_Times.txt" << std::endl;
          for (int i  = 0; i < offsets.size(); i++)
          {
            std::cout << i << " " << offsets[i] << std::endl;
          }

          throw TGo4UserException(3,"");
        }
#endif
        return;
      }

      // Clean up huge event buffers - for now we just destroy them
      if (pInput->DecayEvents.size() > 400)
      {
        decayEvents--;
        nonsenseEvents++;
//        pInput->SetValid(kFALSE); ///NEEDED!?
        return;
      }

      std::vector<AidaCluster> clusters = EventsToClusters(pInput->DecayEvents);

      std::vector<std::pair<AidaCluster, AidaCluster>> hits;

      for (auto& i : clusters)
      {

        if(i.DSSD == -1 || i.Side != conf->DSSD(i.DSSD -1).XSide) continue;

        //if(i.Energy < 100) continue;
        for (auto& j : clusters)
        {
          if(j.DSSD != i.DSSD || j.Side != conf->DSSD(j.DSSD -1).YSide) continue;
          //if(j.Energy < 100) continue;
          // Gates
          if (abs(i.Energy - j.Energy) < conf->FrontBackEnergyL() && i.IsGoodTime(j, conf->FrontBackWindow()))
          hits.push_back({i, j});
        }
      }

      for (auto& i : hits)
      {
        AidaHit hit = ClusterPairToHit(i);

        pInput->Decays.push_back(hit);
        //pOutput->pAida.push_back(hit); ///TEST

        pOutput->pAida.Decays.push_back(hit);
        decays_strip_xy[hit.DSSD - 1]->Fill(hit.StripX, hit.StripY);
        //decays_pos_xy[hit.DSSD - 1]->Fill(hit.PosX, hit.PosY);
        decays_e[hit.DSSD - 1]->Fill(hit.Energy);
       // decays_e_xy[hit.DSSD - 1]->Fill(hit.EnergyFront, hit.EnergyBack);
        //decays_time_delta[hit.DSSD - 1]->Fill(hit.TimeFront - hit.TimeBack);
        decays_time_delta[hit.DSSD - 1]->Fill(hit.FastTimeFront - hit.FastTimeBack);

        int channel = i.first.Strip;
        decays_strip_1d[hit.DSSD - 1]->Fill(channel);
        channel = i.second.Strip + 128;
        decays_strip_1d[hit.DSSD - 1]->Fill(channel);
      }

      if (clusters.size() > 0)
      {
        decays_per_event[0]->Fill(clusters.size());
      }
      //   cout <<" decayEvents " << decayEvents << endl;
    }
    else
    {
      nonsenseEvents++;
      //pInput->SetValid(kFALSE);
    }

  }
}



std::vector<AidaCluster> EventAnlProc::EventsToClusters(std::vector<AidaEvent> const& events)
{
  std::vector<AidaCluster> clusters;
  for (auto& i : events)
  {
    // Don't cluster invalid events
    if (i.DSSD == -1) continue;

    bool added = false;

    // Try to add the event to an existing cluster
    for (auto& j : clusters)
    {
      if(j.IsAdjacent(i) && j.IsGoodTime(i))
      {
        j.AddEvent(i);
        added = true;
        break;
      }
    }

    // Otherwise make a new cluster for the event
    if (!added)
    {
      AidaCluster c_test;
      c_test.AddEvent(i);
      clusters.push_back(c_test);
    }
  }
  return clusters;
}

AidaHit EventAnlProc::ClusterPairToHit(std::pair<AidaCluster, AidaCluster> const& i)
{
  AidaHit hit;
  hit.DSSD = i.first.DSSD;

  hit.StripX = i.first.Strip;
  hit.StripY = i.second.Strip;
  hit.PosX = 75.6 * i.first.Strip / 128. - 37.75;
  hit.PosY = 75.6 * i.second.Strip / 128. - 37.75;

  hit.StripXMin = i.first.StripMin;
  hit.StripXMax = i.first.StripMax;
  hit.StripYMin = i.second.StripMin;
  hit.StripYMax = i.second.StripMax;

  hit.Energy = (i.first.Energy + i.second.Energy) / 2;
  hit.EnergyFront = i.first.Energy;
  hit.EnergyBack = i.second.Energy;

  hit.Time = std::min(i.first.Time, i.second.Time);
  hit.TimeFront = i.first.Time;
  hit.TimeBack = i.second.Time;
  hit.FastTime = std::min(i.first.FastTime, i.second.FastTime);
  hit.FastTimeFront = i.first.FastTime;
  hit.FastTimeBack = i.second.FastTime;

  return hit;
}
                                                    ///End of Aida ///
/**----------------------------------------------------------------------------------------------**/
/**--------------------------------------  bPLASTIC VME (+Scalar)  ----------------------------------------**/
/**----------------------------------------------------------------------------------------------**/

// void EventAnlProc::Make_Plastic_VME_Histos(){
//   for (int i=0; i<32; i++){
//     hPLAS_QDCCalib1[i] =  MakeTH1('D', Form("bPlastic/Energy/QDC1Calib/QDC1Calib_Ch.%2d",i), Form("QDC1 Calib Ch. %2d",i), 20000, 0., 20000.);
//    // hPLAS_TDCCalib1[i] =  MakeTH1('D', Form("bPlastic/Timing/TDC1Calib/TDC1Calib_Ch.%2d",i), Form("TDC1 Ch. %2d",i), 2E4, 0, 2E5);
//     hPLAS_TimeDiffSiPM_Ch_Raw[i] = MakeTH1('D',Form("bPlastic/Timing/TDCdt_SiPM1-SiPMRaw/TDCSiPM_dT_Ch.%02d",i),Form("SiPM dT(ns) Ch1 - Ch%0d",i),250,-2E4,2E4);
//     hPLAS_TimeDiffSiPM_Ch_Calib[i] = MakeTH1('D',Form("bPlastic/Timing/TDCdt_SiPM1-SiPMCalib/TDCSiPM1_dT_Ch.%02d",i),Form("TimeDiff Gainmatched Ch1 - Ch%0d",i),250,-2E4,2E4);
//     hPLAS_TimeDiffSiPM_Ch_Calib_Egated[i] = MakeTH1('D',Form("bPlastic/Timing/EGated/TDCdt_SiPM1-SiPM/TDCSiPM1_dT_EGated_Ch.%02d",i),Form("TimeDiff Gainmatched Ch1 - Ch%0d",i),250,-2E4,2E4);
//
//     //hPLAS_TimeDiff_Ch_Raw[i] = MakeTH1('D',Form("bPlastic/Timing/Raw/TDCdt_ref-plas/TDCSiPM1_dt_Raw_Ch.%02d",i),Form("TimeDiff Ch0 - Ch%0d",i),1E5,-2E5,2E5);
//    // hPLAS_TimeDiffSiPM_Ch[i] = MakeTH1('D',Form("bPlastic/Timing/TDCdt_SiPM1-SiPm/TDCSiPMdt_Ch.%02d",i),Form("TimeDiff Ch1 - Ch%0d",i),1E5,-2E5,2E5);
//     ///new
//     hPLAS_TDC_FiredRatio[i] = MakeTH1('D', Form("bPlastic/Timing/TDC1_FiredRatio/TDC_FiredRatio_Ch.%2d",i), Form("TDC1 Ratio CalibTDC/FiredTDC Ch. %2d",i), 4000, 0, 4000);
//     hPLAS_TimeDiff_SC41_Raw[i] = MakeTH1('D',Form("bPlastic/Timing/TDCdt_SC41-SiPMRaw/TDCdT_SC41_Plas_RawCh.%02d",i),Form("SC41-SiPM dT(ns) Raw Ch%0d",i),5E3,0,5E5);
//     hPLAS_TimeDiff_SC41_Calib[i] = MakeTH1('D',Form("bPlastic/Timing/TDCdt_SC41-SiPMCalib/TDCdT_SC41_Plas_Ch.%02d",i),Form("SC41-SiPM dT Calib  Ch%0d",i),5E3,0,5E5);
//     hPLAS_TimeDiff_SC41_Calib_Egated[i] = MakeTH1('D',Form("bPlastic/Timing/EGated/TDCdt_SC41-SiPM/TDCdT_SC41_Plas_EGated_Ch.%02d",i),Form("SC41-SiPM dT Calib  Ch%0d",i),5E4,0,5E5);
//
//     hPLAS_TDC_multich[i] = MakeTH1('D', Form("bPlastic/Stats/TDC_MultiCh/TDCMch%2d",i), Form("TDC channel Multi %2d",i), 50, 0, 50);
//
// //     hPLAS_EvsdT[i] = MakeTH2('D',Form("bPlastic/Timing/E_vs_dT/E_vs_dT_Ch.%02d",i), Form("E_vs_dT_Ch.%0d",i), 2500, -5000, 0,  2500, 0, 5000);
//    // hPLAS_CoincE1E2[i] = MakeTH2('D',Form("bPlastic/Energy/CoincEnergy/Coinc_Energy_Energy_Ch.%02d",i), Form("Coinc_Energy_Energy_Sum_Ch.%0d",i), 500, 0, 5000,  500, 0, 5000);
//     hPLAS_CoincE_dTSiPM_Ch[i] =   MakeTH2('D',Form("bPlastic/Timing/Coinc_Energy_SiPM_dT/Coinc_E_dTCh1-Ch.%02d",i), Form("Energy vs. SiPM1-SiPMCh.%0d",i), 2500, 0, 5000, 2500,-5000,5000);
//
//
//   }
//     //hPLAS_E1E2TimeGated[i] = MakeTH2('D',Form("bPlastic/CoincEEGated/CoincEEGated%02d",i), Form("CoincGated%0d",i), 400, 0, 4000.0,  400, 0, 4000.0);
//
//
//     //Sum spectra
//     hPLAS_QDCCalib1Sum = MakeTH1('D',"bPlastic/Energy/QDCSum","QDC Calibrated Sum",2000,0,20000);
//     hPLAS_TDCCalib1Sum = MakeTH1('D',"bPlastic/Timing/TDCSum","TDC Calibrated Sum (ns)",2E4, 0, 2E5);
//     hPLAS_TimeDiffSiPM_Ch_Sum = MakeTH1('D',"bPlastic/Timing/TDCdT_SiPM1-SiPM","SiPM1 - SiPM Ch.x dT(ns) (calibrated)",1000,-2E6,2E6);
//     hPLAS_TimeDiffSiPM_Ch_Sum_M1 = MakeTH1('D',"bPlastic/Timing/Multiplicity/TDCdT_SiPM1-SiPM_Multi1","SiPM1 - SiPM Ch.x dT(ns) Multiplicity1 (calibrated)",250,-2E4,2E4);
//     hPLAS_TimeDiffSiPM_Ch_Sum_M2 = MakeTH1('D',"bPlastic/Timing/Multiplicity/TDCdT_SiPM1-SiPM_Multi2","SiPM1 - SiPM Ch.x dT(ns) Multiplicity2 (calibrated)",250,-2E4,2E4);
//     hPLAS_TimeDiffSiPM_Ch_Sum_M3P = MakeTH1('D',"bPlastic/Timing/Multiplicity/TDCdT_SiPM1-SiPM_Multi3","SiPM1 - SiPM Ch.x dT(ns) Multiplicity>3 (calibrated)",250,-2E4,2E4);
//
//     hPLAS_TimeDiffSiPM_Ch_Sum_Egated = MakeTH1('D',"bPlastic/Energy/EGated/TDCdT_SiPM1-SiPM_EGated","Energy gated SiPM1 - SiPM  dT(ns) Ch.x (calibrated)",250,-2E4,2E4);
//
//     hPLAS_TimeDiff_SC41_Sum = MakeTH1('D',"bPlastic/Timing/TDCdT_SC41-SiPM","SC41 - SiPM dT(ns) (calibrated)",4E4,0,2E6);
//     hPLAS_TimeDiff_SC41_Sum_M1 = MakeTH1('D',"bPlastic/Timing/Multiplicity/TDCdT_SC41-SiPM_Multi1","SC41 - SiPM  Ch.x dT(ns) Multiplicity1 (calibrated)",4E4,0,2E6);
//     hPLAS_TimeDiff_SC41_Sum_M2 = MakeTH1('D',"bPlastic/Timing/Multiplicity/TDCdT_SC41-SiPM_Multi2","SC41 - SiPM  Ch.x dT(ns) Multiplicity2 (calibrated)",4E4,0,2E6);
//     hPLAS_TimeDiff_SC41_Sum_M3P = MakeTH1('D',"bPlastic/Timing/Multiplicity/TDCdT_SC41-SiPM_Multi3+","SC41 - SiPM  Ch.x dT(ns) Multiplicity >3 (calibrated)",4E4,0,2E6);
//
//
//     hPLAS_TimeDiff_SC41_Sum_Egated  = MakeTH1('D',"bPlastic/Timing/EGated/TDCdT_SC41-SiPM_EGated","Energy gated SC41 - SiPM dT(ns) Ch.x (calibrated)",4E4,0,2E6);
//     hPLAS_TDC_FiredRatio_Sum = MakeTH1('D',"bPlastic/Timing/TDC_FiredRatio","TDC1 Ratio CalibTDC/FiredTDC",4000, 0, 4000);
//
//     hPLAS_QDC1_hits  = MakeTH1('D',"bPlastic/Stats/QDC1_hits","bPlastic hit pattern QDC1",32,0,32);
//     hPLAS_TDC_hits  = MakeTH1('D',"bPlastic/Stats/TDC_hits","bPlastic hit pattern TDC",32,0,32);
//     hPLAS_TDC_multi = MakeTH1('D',"bPlastic/Stats/TDC_Multi","bPlastic TDC Multiplicity",50,0,50);
//
//     hPLAS_CoincE1E2_Sum = MakeTH2('D',"bPlastic/Energy/Coinc_Energy_Energy_Sum","bPlastic Energy-Energy Sum", 500, 0, 5000,  500, 0, 5000);
//     hPLAS_CoincE_dTSiPM_Sum = MakeTH2('D',"bPlastic/Energy/CoincEnergy_SiPM1-SiPMCh.x","Energy vs. SiPMCh.1-SiPM_all", 5000, 0, 5000, 100,-100,100);
//
//     hScalar_hit_pattern = MakeTH1('D',"Scalar/HitPat","Scalar Hit pattern",32,0,32);
// }
///-----------------------------------------------------------------------------------------------------------------------------------------------------------------------///
// void EventAnlProc::Do_Plastic_VME_Histos(EventAnlStore* pOutput){
//
//   int bPlasTDCIDMain;
//     int bPlasQDCID_i, bPlasQDCID_j;
//     double bPlasTDC_TS_Raw[32],    bPlasTDC_T_Calib[32] ;
//     double bPlasTDC_TS_Firedratio[32];
//     double bPlas_SiPM_dT_Raw[32], bPlas_SiPM_dT_Calib[32];
//     double bPlas_SC41_dT_Raw[32], bPlas_SC41_dT_Calib[32];
//     double bPlasQDCGainMatch_i[32],bPlasQDCGainMatch_j[32];
//     double bPlas_TDC_Cha1;
//   bPlasTDCIDMain = -1;
//     bPlas_TDC_Cha1 = 0;
//     bPlasQDCID_i=0;
//     bPlasQDCID_j=0;
//   for (int i=0; i< 32; i++){
//
//         bPlasTDC_TS_Raw[i]=0;
//         bPlasTDC_T_Calib[i]=0;
//         bPlasQDCGainMatch_i[i] = 0;
//         bPlasQDCGainMatch_j[i] = 0;
//         bPlas_SiPM_dT_Raw[i] = 0;
//         bPlas_SiPM_dT_Calib[i] = 0;
//         bPlas_SC41_dT_Raw[i] = 0;
//         bPlas_SC41_dT_Calib[i] = 0;
//   }
//   /**------------------bPlastic Energy -----------------------------------------**/
//          pOutput->pbPlas_QDCFired = bPlasQDCFired;
//     for (int i=0; i<bPlasQDCFired; i++){
//          bPlasQDCID_i = bPlasQDCID[i];
//          pOutput->pbPlas_QDCID[i] = bPlasQDCID[i];
//          bPlasQDCGainMatch_i[bPlasQDCID_i] = fCal->AplasQDC[bPlasQDCID_i]*bPlasQDC[bPlasQDCID_i] + fCal -> BplasQDC[bPlasQDCID_i]; //gain matching
//           hPLAS_QDC1_hits->Fill(bPlasQDCID_i);
//
//           if(bPlasQDCGainMatch_i[bPlasQDCID_i]>30){
//                if(bPlasQDCID_i<16 ){//temporarily remove the noise for chans>15
//
//              pOutput->pbPlas_QDCGainMatch_i[bPlasQDCID_i] = bPlasQDCGainMatch_i[bPlasQDCID_i];
//
//              hPLAS_QDCCalib1[bPlasQDCID_i]->Fill(bPlasQDCGainMatch_i[bPlasQDCID_i]);
//              hPLAS_QDCCalib1Sum->Fill(bPlasQDCGainMatch_i[bPlasQDCID_i]);
//                              }
//           }
//                ///Energy-Energy Matrix
//                 for(int j=0; j<bPlasQDCFired; j++){
//                   bPlasQDCID_j = bPlasQDCID[j];
//                    ///Dont loop on the first hit again: (check)
//                   if(bPlasQDCID_i<bPlasQDCID_j){
//                   bPlasQDCGainMatch_j[bPlasQDCID_j] = fCal->AplasQDC[bPlasQDCID_j]*bPlasQDC[bPlasQDCID_j] + fCal -> BplasQDC[bPlasQDCID_j]; //gain matching
//                //   hPLAS_CoincE1E2[bPlasQDCID_i] ->Fill(bPlasQDCGainMatch_i[bPlasQDCID_i],bPlasQDCGainMatch_j[bPlasQDCID_j]);
//                   hPLAS_CoincE1E2_Sum->Fill(bPlasQDCGainMatch_i[bPlasQDCID_i],bPlasQDCGainMatch_j[bPlasQDCID_j]);
//
//                }
//             }
//           }
//      /**----------------------------bPlastic Timing -----------------------------------------**/
//      ///Channel 0 is SC41 trigger
//         pOutput->pbPlas_TDCFired = bPlasTDCFired;
//   for (int i=0; i<bPlasTDCFired;i++){
//
//           bPlasTDCIDMain = bPlasTDCID[i];
//           pOutput->pbPlas_TDCID[i] = bPlasTDCID[i];
//           bPlasTDC_TS_Raw[bPlasTDCIDMain] = bPlasTDC_TS[i][bPlasTDCIDMain];
//           bPlas_TDC_Multiplicity[bPlasTDCIDMain]++;
//           pOutput-> pbPlas_TDC_Multiplicity[bPlasTDCIDMain] = bPlas_TDC_Multiplicity[bPlasTDCIDMain];
//           hPLAS_TDC_hits->Fill(bPlasTDCIDMain);
//           hPLAS_TDC_multi->Fill(bPlas_TDC_Multiplicity[bPlasTDCIDMain]);
//        //  cout <<"FAT Event " << event_number <<" bPlasTDC_TS_Raw[bPlasTDCIDMain]" <<bPlasTDC_TS_Raw[bPlasTDCIDMain]  << " bPlasTDCID[i] " << bPlasTDCID[i] << endl;
//
//           //cout <<"ev " << event_number << " fired " <<  bPlasTDCFired << " i " << i <<" bPlasTDCID[i] " << bPlasTDCID[i]  <<endl;
//
//           ///Calibrate raw TDC
//         //  bPlasTDC_T_Calib[bPlasTDCIDMain] =  bPlasTDC_TS_Raw[bPlasTDCIDMain] + fCal->AplasTDC_Raw[bPlasTDCIDMain];
//           pOutput->pbPlasTDC_T[bPlasTDCIDMain] =   bPlasTDC_TS_Raw[bPlasTDCIDMain]; //Output bPlas Raw TDC
//             ///Get the first hit of the reference channel (Cha.1)
//            if(bPlasTDC_TS[0][1]>0 && bPlas_TDC_Multiplicity[bPlasTDCIDMain]==1 &&bPlasTDCIDMain ==1){
//                            bPlas_TDC_Cha1 =   bPlasTDC_TS[0][1];
//                                    }
//
//           //hPLAS_TDCCalib1[bPlasTDCIDMain] -> Fill(bPlasTDC_T_Calib[bPlasTDCIDMain]);
//           hPLAS_TDCCalib1Sum -> Fill(bPlasTDC_T_Calib[bPlasTDCIDMain]);
//
//           bPlasTDC_TS_Firedratio[bPlasTDCIDMain] = bPlasTDC_T_Calib[bPlasTDCIDMain]/bPlasTDCFired;
//           hPLAS_TDC_FiredRatio[bPlasTDCIDMain] -> Fill(bPlasTDC_TS_Firedratio[bPlasTDCIDMain]);
//           hPLAS_TDC_FiredRatio_Sum -> Fill( bPlasTDC_TS_Firedratio[bPlasTDCIDMain]);
//
//           ///Ref SiPMCh.1 - SiPM Ch.x
//        //   if(bPlas_TDC_Multiplicity[bPlasTDCIDMain]>1 && bPlasTDCIDMain==1){
//
//        //   }
//        //   if(bPlas_TDC_Cha1>0 &&bPlasTDC_TS_Raw[bPlasTDCIDMain]>0){
//         //  cout<<"event " << event_number<<" bPlas_TDC_Cha1 " <<bPlas_TDC_Cha1 << " bPlasTDC_TS_Raw[bPlasTDCIDMain] " << bPlasTDC_TS_Raw[bPlasTDCIDMain]<< endl;
//            bPlas_SiPM_dT_Raw[bPlasTDCIDMain] = (bPlasTDC_TS[0][1] - bPlasTDC_TS_Raw[bPlasTDCIDMain]);
//            bPlas_SiPM_dT_Calib[bPlasTDCIDMain] = ((bPlasTDC_TS[0][1] - bPlasTDC_TS_Raw[bPlasTDCIDMain] )+ fCal->AplasTDC_Chref_dT[bPlasTDCIDMain]);
//            pOutput ->  pbPlas_SiPM_dT_Calib[bPlasTDCIDMain] =   bPlas_SiPM_dT_Calib[bPlasTDCIDMain];
//
//
//            if(bPlasTDCIDMain!=1){
//             hPLAS_TimeDiffSiPM_Ch_Raw[bPlasTDCIDMain] -> Fill(bPlas_SiPM_dT_Raw[bPlasTDCIDMain]);
//             hPLAS_TimeDiffSiPM_Ch_Calib[bPlasTDCIDMain] -> Fill(bPlas_SiPM_dT_Calib[bPlasTDCIDMain]);
//             hPLAS_TimeDiffSiPM_Ch_Sum -> Fill(bPlas_SiPM_dT_Calib[bPlasTDCIDMain]);
//
//                       ///Fill for Multiplicites
// if (bPlas_TDC_Multiplicity[bPlasTDCIDMain] ==1) hPLAS_TimeDiffSiPM_Ch_Sum_M1 ->Fill(bPlas_SiPM_dT_Calib[bPlasTDCIDMain]);
// if (bPlas_TDC_Multiplicity[bPlasTDCIDMain] ==2) hPLAS_TimeDiffSiPM_Ch_Sum_M2 ->Fill(bPlas_SiPM_dT_Calib[bPlasTDCIDMain]);
// if (bPlas_TDC_Multiplicity[bPlasTDCIDMain] >2)  hPLAS_TimeDiffSiPM_Ch_Sum_M3P ->Fill(bPlas_SiPM_dT_Calib[bPlasTDCIDMain]);
//         }
//         //  }
//           ///SC41 - SiPM Ch.x
//           if(bPlasTDC_TS_Raw[0]>0 && bPlasTDC_TS_Raw[bPlasTDCIDMain]>0){
//
//            bPlas_SC41_dT_Raw[bPlasTDCIDMain] = (bPlasTDC_TS_Raw[0] - bPlasTDC_TS_Raw[bPlasTDCIDMain]);
//            bPlas_SC41_dT_Calib[bPlasTDCIDMain] =  bPlas_SC41_dT_Raw[bPlasTDCIDMain] + fCal->BplasTDC_SC41dT[bPlasTDCIDMain];
//            // cout<<" bPlas_SC41_dT_Calib[bPlasTDCIDMain]  " <<  bPlas_SC41_dT_Calib[bPlasTDCIDMain] << endl;
//             pOutput-> pbPlas_SC41_dT[bPlasTDCIDMain] = bPlas_SC41_dT_Calib[bPlasTDCIDMain]; //Output SC41-bPlas dT
//             if(bPlasTDCIDMain!=0){
//                 hPLAS_TimeDiff_SC41_Raw[bPlasTDCIDMain] ->Fill(bPlas_SC41_dT_Raw[bPlasTDCIDMain]);
//                 hPLAS_TimeDiff_SC41_Calib[bPlasTDCIDMain] ->Fill(bPlas_SC41_dT_Calib[bPlasTDCIDMain]);
//                 hPLAS_TimeDiff_SC41_Sum ->Fill(bPlas_SC41_dT_Calib[bPlasTDCIDMain]);
//                 ///Fill for Multiplicites
// if (bPlas_TDC_Multiplicity[bPlasTDCIDMain] ==1) hPLAS_TimeDiff_SC41_Sum_M1 ->Fill(bPlas_SC41_dT_Calib[bPlasTDCIDMain]);
// if (bPlas_TDC_Multiplicity[bPlasTDCIDMain] ==2) hPLAS_TimeDiff_SC41_Sum_M2 ->Fill(bPlas_SC41_dT_Calib[bPlasTDCIDMain]);
// if (bPlas_TDC_Multiplicity[bPlasTDCIDMain] >2)  hPLAS_TimeDiff_SC41_Sum_M3P ->Fill(bPlas_SC41_dT_Calib[bPlasTDCIDMain]);
//
//             }
//           }
//   /**-----------------------bPlastic Energy gated Timing (gates defined in Correlations.dat)-----------------------------------------**/
//            /// NOTE: TDC ID = QDC ID +1 (TDC Ch.0 =SC41)
//
//
//                 ///Energy-Time matrices
//                 if(bPlasQDCGainMatch_i[bPlasTDCIDMain-1]>50){
//            hPLAS_CoincE_dTSiPM_Ch[bPlasTDCIDMain-1] -> Fill(bPlasQDCGainMatch_i[bPlasTDCIDMain-1],  pOutput ->  pbPlas_SiPM_dT_Calib[bPlasTDCIDMain]);
//           //   cout <<"event " <<event_number<<" dT "<<    pOutput ->  pbPlas_SiPM_dT_Calib[bPlasTDCIDMain] << " E " <<bPlasQDCGainMatch_i[bPlasTDCIDMain-1] <<endl;
//
//                 }
//
//                         for (int j=0; j<bPlasQDCFired; j++){
//                     ///Energy Gate
//                 if(bPlasQDCGainMatch_i[bPlasQDCID[j]] > fCorrel->GbPlas_Egate_low &&  bPlasQDCGainMatch_i[bPlasQDCID[j]] < fCorrel->GbPlas_Egate_high   ){
//                             ///Energy gated SiPM - SiPM x.
//                     if(bPlasTDCIDMain!=1){
//                     hPLAS_TimeDiffSiPM_Ch_Calib_Egated[bPlasTDCIDMain] -> Fill(bPlas_SiPM_dT_Calib[bPlasTDCIDMain]);
//                     hPLAS_TimeDiffSiPM_Ch_Sum_Egated -> Fill(bPlas_SiPM_dT_Calib[bPlasTDCIDMain]);
//                             }
//
//                     if(bPlasTDC_TS_Raw[0]>0 && bPlasTDC_TS_Raw[bPlasTDCIDMain]>0){
//                             ///Energy gated SC41 - SiPM Ch.x
//                         if(bPlasTDCIDMain!=0 && bPlasTDC_TS_Raw[bPlasTDCIDMain]> 2580 && bPlasTDC_TS_Raw[bPlasTDCIDMain]<2640){
//                     hPLAS_TimeDiff_SC41_Calib_Egated[bPlasTDCIDMain] ->Fill(bPlas_SC41_dT_Calib[bPlasTDCIDMain]);
//                     hPLAS_TimeDiff_SC41_Sum_Egated ->Fill(bPlas_SC41_dT_Calib[bPlasTDCIDMain]);
//               }
//             }
//           }
//         }
//       }
//   }

 /**----------------------------------------------------------------------------------------------**/   
     /**--------------------------------------  bPlastic and FATIMA Combined TAMEX ----------------------------------------------   
     //**----------------------------------------------------------------------------------------------**/  
        
    void EventAnlProc::Make_Plastic_Tamex_Histos(){ 
        Text_t chis[256];   
        Text_t chead[256];  
        
         for (int i =1; i<3; i++)   
      { 
          for(int j=0; j<16; j++){  
           sprintf(chis,"bPlastic/ToT/ToT Det.%2d Ch.%2d", i,j);    
           sprintf(chead,"TOT Detector %2d Ch. %2d", i,j);  
           hbPlas_ToT_det[i][j] = MakeTH1('I', chis,chead, 25000, 0., 150000.); 
            
           hbPlas_Energy_Calib[i][j] = MakeTH1('D', Form("bPlastic/Energy_Calib/Energy Calib Plas Det. %2d Ch. %2d",  i,j), Form("Energy Calib Det. %2d Ch. %2d", i,j),150000, 0., 150000.);    
          } 
      } 
        
            
        for(int i =0; i<bPLAS_TAMEX_NUM; i++){      
        for (int j =0; j<32; j++)   
      { 
          if(i==bPLAS_TAMEX_ID){    
    //        sprintf(chis,"bPlastic/ToT/ToT Det. %2d Ch.%2d", i,j);    
    //        sprintf(chead,"TOT Detector %2d Ch. %2d", i,j);   
    //        hbPlas_ToT[i][j] = MakeTH1('I', chis,chead, 25000, 0., 150000.);  
            
    //        sprintf(chis,"bPlastic/ToT_coinc_2chans/ToT_coinc_2Chans. %2d Ch.%2d", i,j);  
    //        sprintf(chead,"TOT coinc (min 2 channels fired) Detector %2d Ch. %2d", i,j);  
    //        hbPlas_ToT_coinc[i][j] = MakeTH1('I', chis,chead, 2500, 0., 15000.);  
            
    //         sprintf(chis,"bPlastic/ToT_coinc_2dets/ToT_coinc 2Dets. %2d Ch.%2d", i,j);   
    //        sprintf(chead,"TOT coinc (min 2 channels and 2 detectors fired) Detector %2d Ch. %2d", i,j);  
    //        hbPlas_ToT_coinc_dets[i][j] = MakeTH1('I', chis,chead, 2500, 0., 15000.); 
            
            
    //        sprintf(chis,"bPlastic/ToT_coinc_2dets_perchan/ToT_coinc 2Dets_perchan. %2d Ch.%2d", i,j);    
    //        sprintf(chead,"TOT coinc (coincident channels and 2 detectors fired) Detector %2d Ch. %2d", i,j); 
    //        hbPlas_ToT_coinc_dets_perchan[i][j] = MakeTH1('I', chis,chead, 2500, 0., 15000.); 
       // hbPlas_ToTCoin[i][j] = MakeTH1('I',Form("bPlastic/TOT_COIN/ToTCoin %2d Ch.%2d",  i,j), 2500, 0., 10000.); 
            
        hbPlas_lead_lead[i][j] = MakeTH1('D', Form("bPlastic/Lead-Lead/Lead-Lead Plas Det. %2d Ch.%2d",  i,j), Form("Lead - Lead Det %2d Ch. %2d", i,j),2500, -50000., 50000);  
        
        hbPlas_lead_lead_ref[i][j] = MakeTH1('D', Form("bPlastic/Lead-Lead_Ref/Lead-Lead Plas Det. %2d RefCh. %2d", i,j), Form("Lead Ref Ch.0 - Lead Det.%2d Ch. %2d", i,j),2500, -50000., 50000.); 
            
        hbPlas_lead_lead_gated[i][j] = MakeTH1('D', Form("bPlastic/Lead-Lead_Egated/Lead-Lead Egated Plas Det. %2d Ch. %2d",  i,j), Form("Lead - Lead Energy gated Det. %2d Ch.  %2d", i,j),2500, -50000., 50000.); 
            
    //     hbPlas_trail_trail[i][j] = MakeTH1('D', Form("bPlastic/Trail-Trail/trail-trail Plas Det. %2d Ch. %2d",  i,j), Form("Trail - Trail Det. %2d Ch. %2d", i,j), 2500, -50000., 50000.);   
    //      
        
            
       // hbPlas_SC41L_lead[i][j] = MakeTH1('D', Form("bPlastic/SC41-Lead_Plas/SC41_Lead Plas Det. %2d Ch.%02d", i,j), Form("SC41 Lead - bPlas Lead Det. %2d Lead Ch. %2d ", i,j), 4002, -100000., 100000.);    
         hbPlas_SC41L_Anal_lead[i][j] = MakeTH1('D', Form("bPlastic/SC41L_Anal-Lead_bPlas/SC41L_Anal_Lead bPlas Det. %2d Ch.%02d", i,j), Form("SC41L Analogue Lead - bPlas Lead Det. %2d Lead Ch. %2d ", i,j), 4002, -100000., 100000.);    
            
         hbPlas_SC41R_Anal_lead[i][j] = MakeTH1('D', Form("bPlastic/SC41R_Anal-Lead_bPlas/SC41R_Anal_Lead bPlas Det. %2d Ch.%02d", i,j), Form("SC41R Analogue Lead - bPlas Lead Det. %2d Lead Ch. %2d ", i,j), 4002, -100000., 100000.);    
            
         hbPlas_SC41L_Digi_lead[i][j] = MakeTH1('D', Form("bPlastic/SC41L_Digi-Lead_bPlas/SC41L_Digi_Lead bPlas Det. %2d Ch.%02d", i,j), Form("SC41L Digital Lead - bPlas Lead Det. %2d Lead Ch. %2d ", i,j), 4002, -100000., 100000.); 
            
         hbPlas_SC41R_Digi_lead[i][j] = MakeTH1('D', Form("bPlastic/SC41R_Digi-Lead_bPlas/SC41R_Digi_Lead bPlas Det. %2d Ch.%02d", i,j), Form("SC41R Digital Lead - bPlas Lead Det. %2d Lead Ch. %2d ", i,j), 4002, -100000., 100000.);             
        
        }   
                        
        if (i==FAT_TAMEX_ID){   
           sprintf(chis,"FATIMA_TAMEX/ToT/ToT Fat Det. %2d Ch.%2d", i,j);   
           sprintf(chead,"TOT Detector %2d Ch. %2d", i,j);  
           hFATTEMP_ToT[i][j] = MakeTH1('I', chis,chead, 25000, 0., 150000.);   
                
    //        sprintf(chis,"FATIMA_TAMEX/ToT_coinc_2dets_perchan/ToT_coinc 2Dets_perchan. %2d Ch.%2d", i,j);    
    //        sprintf(chead,"TOT coinc (coincident channels and 2 detectors fired) Detector %2d Ch. %2d", i,j); 
    //        hFATTEMP_ToT_coinc_dets_perchan[i][j] = MakeTH1('I', chis,chead, 2500, 0., 15000.);   
       // hbPlas_ToTCoin[i][j] = MakeTH1('I',Form("bPlastic/TOT_COIN/ToTCoin %2d Ch.%2d",  i,j), 2500, 0., 10000.); 
            
        hFATTEMP_lead_lead[i][j] = MakeTH1('D', Form("FATIMA_TAMEX/Lead-Lead/Lead-Lead Fat Det. %2d Ch.%2d",  i,j), Form("Lead - Lead Det %2d Ch. %2d", i,j),2500, -50000., 50000); 
        
            
        
        hFATTEMP_lead_lead_gated[i][j] = MakeTH1('D', Form("FATIMA_TAMEX/Lead-Lead_Egated/Lead-Lead Fat Egated Det. %2d Ch. %2d",  i,j), Form("Lead - Lead Energy gated Det. %2d Ch.  %2d", i,j),2500, -50000., 50000.);    
            
    //     hFATTEMP_trail_trail[i][j] = MakeTH1('D', Form("FATIMA_TAMEX/Trail-Trail/trail-trail Fat Det. %2d Ch. %2d",  i,j), Form("Trail - Trail Det. %2d Ch. %2d", i,j), 2500, -50000., 50000.);  
    //      
        hFATTEMP_Energy_Calib[i][j] = MakeTH1('D', Form("FATIMA_TAMEX/Energy_Calib/Energy_Calib Fat Det. %2d Ch. %2d",  i,j), Form("Energy Calib Det. %2d Ch. %2d", i,j),2500, -50000., 50000.);        
            
        hFATTEMP_SC41L_Anal_lead[i][j] = MakeTH1('D', Form("FATIMA_TAMEX/SC41L_Anal-Lead_Fat/SC41L_Anal_Lead Fat Det. %2d Ch.%02d", i,j), Form("SC41L Analogue Lead - Fat Lead Det. %2d Lead Ch. %2d ", i,j), 4002, -100000., 100000.); 
            
        hFATTEMP_SC41R_Anal_lead[i][j] = MakeTH1('D', Form("FATIMA_TAMEX/SC41R_Anal-Lead_Fat/SC41R_Anal_Lead Fat Det. %2d Ch.%02d", i,j), Form("SC41R Analogue Lead - Fat Lead Det. %2d Lead Ch. %2d ", i,j), 4002, -100000., 100000.); 
            
        hFATTEMP_SC41L_Digi_lead[i][j] = MakeTH1('D', Form("FATIMA_TAMEX/SC41L_Digi-Lead_Fat/SC41L_Digi_Lead Fat Det. %2d Ch.%02d", i,j), Form("SC41L Digital Lead - Fat Lead Det. %2d Lead Ch. %2d ", i,j), 4002, -100000., 100000.);  
            
        hFATTEMP_SC41R_Digi_lead[i][j] = MakeTH1('D', Form("FATIMA_TAMEX/SC41R_Digi-Lead_Fat/SC41R_Digi_Lead Fat Det. %2d Ch.%02d", i,j), Form("SC41R Digital Lead - Fat Lead Det. %2d Lead Ch. %2d ", i,j), 4002, -100000., 100000.);  
                    }   
        
            }   
        }   
            
        for(int j=0; j<16; j++){    
        hFATTEMP_lead_lead_ref[j] = MakeTH1('D', Form("FATIMA_TAMEX/Lead-Lead_Ref/Lead-Lead Fat RefCh. %2d",j), Form("Lead Ref Ch.0 - Lead Ch. %2d", j),2500, -50000., 50000.); 
        
        }   
        hSC41_Analogue_Tamex = MakeTH1('D',"bPlastic/SC41/Analogue L-R","SC41 Analogue L - R",4002, -100000., 100000.); 
        hSC41_Digital_Tamex = MakeTH1('D',"bPlastic/SC41/Digital L-R","SC41 Analogue L - R",4002, -100000., 100000.);   
            
        hbPlas_ToT_Sum_Det1 = MakeTH1('I',"bPlastic/ToT_Sum_Det1","bPlastic Sum Gainmatched ToT Detector 1 ",250000, 0., 1500000.); 
            
        hbPlas_ToT_Sum_Det2 = MakeTH1('I',"bPlastic/ToT_Sum_Det2","bPlastic Sum Gainmatched ToT Detector 2 ",250000, 0., 1500000.); 
            
        hbPlas_hit_pattern = MakeTH1('D',"bPlastic/Stats/HitPat","Scalar Hit pattern",32,0,32); 
        hbPlas_num_fired_chans = MakeTH1('D',"bPlastic/Stats/Hits_Chan","bPlastic Multiplicity",32,0,32);   
        
    }   
    /////////////////////////////////////////////////// 
    void EventAnlProc::Do_Plastic_Tamex_Histos(EventUnpackStore* pInput, EventAnlStore* pOutput){   
           
        bool fired_det1=false, fired_det2=false;    
        int bPlas_tot_hits; 
       
        for (int i = 0; i < 48; i++)    
            {   
                for(int j=0; j<100;j++){    
                 lead_lead_bplas[i][j]=0;   
                 lead_lead_fat[i][j]=0; 
                 lead_lead_bplas_Ref0[i][j]=0;  
                 lead_lead_fat_Ref0[i][j]=0;    
                 trail_trail_bplas[i][j]=0; 
                 ToT_bplas[i][j] = 0;   
                 trail_bplas1[i][j] = 0;    
                 trail_bplas2[i][j] = 0;    
                 lead_bplas1[i][j]=0;   
                 lead_bplas2[i][j] =0;  
    //              SC41_lead_bplas[i][j] = 0;  
    //              SC41_lead_fat[i][j] = 0;    
                 SC41L_ANA_lead_bPlas[i][j] = 0;    
                 SC41R_ANA_lead_bPlas[i][j] = 0;    
                 SC41L_DIG_lead_bPlas[i][j] = 0;    
                 SC41R_DIG_lead_bPlas[i][j] = 0;    
                 SC41L_ANA_lead_fat[i][j] = 0;  
                 SC41R_ANA_lead_fat[i][j] = 0;  
                 SC41L_DIG_lead_fat[i][j] = 0;  
                 SC41R_DIG_lead_fat[i][j] = 0;  
                }   
            }
          
            
        
     ///**---------------------------------------------LEAD -------------------------------------------------**/        
           ///Loop on channels First    
              
                                  
              for (int i = 0; i < 48; i++)  
            {   
                for(int j=0; j< pInput->fbPlas_PMT_Lead_N[i]; j++){ 
              Fat_RefCh[j] = pInput->fbPlas_Lead_PMT[1][j]; 
              bPlas_RefCh[j] = pInput->fbPlas_Lead_PMT[16][j];  
                }   
            }   
            ////////////////////////////    
              ///Loop over channels 
          for (int i = 0; i < 48; i++)  
            {   
                 int bplas_detnum =-1;  
                       if(i>15 && i<32) bplas_detnum=1; 
                       if(i>31)  bplas_detnum=2;    
                        
                pOutput->pbPlas_PMT_Lead_N[i] = pInput->fbPlas_PMT_Lead_N[i];   
                //Lead T    
                for(int j=0; j< pInput->fbPlas_PMT_Lead_N[i]; j++){ ///Hits 
                 //cout<<"j " << j << " pInput->fbPlas_Lead_PMT[1][j] " << pInput->fbPlas_Lead_PMT[1][j]<<" Fat_RefCh[j] " <<Fat_RefCh[j] <<  " i " << i <<endl;    
                        
                    
                    lead_bplas1[i][j] = pInput->fbPlas_Lead_PMT[i][j];      
                   // cout <<"lead_bplas1[i][j] " <<lead_bplas1[i][j] << " i " << i << " j " << j <<endl;   
                    hits_bplas_lead++;  
                    pOutput->pbPlas_LeadT[i][j] = lead_bplas1[i][j];    
                    pOutput->pbPlas_LeadHits = hits_bplas_lead; 
                    pOutput->pbPlas_LeadT_Avg = lead_bplas1[i][j]/hits_bplas_lead;  
        
        
    ///**---------------------------------------------Ref channel -------------------------------------------------**/          
            //Plastic               
                    if(i>15 && pInput->fbPlas_Lead_PMT[16][j]>0 && pInput->fbPlas_Lead_PMT[i][j]>0) {   
                     lead_lead_bplas_Ref0[i][j] = (bPlas_RefCh[j] -  pInput->fbPlas_Lead_PMT[i][j]);    
                     hbPlas_lead_lead_ref[bPLAS_TAMEX_ID][i-16] ->Fill(lead_lead_bplas_Ref0[i][j]);  
//                      if(j>0){
//                      cout<<"lead_lead_bplas_Ref0[i][j] " << lead_lead_bplas_Ref0[i][j] << " bPlas_RefCh[j] " <<bPlas_RefCh[j] << " pInput->fbPlas_Lead_PMT[i][j] " << pInput->fbPlas_Lead_PMT[i][j] << " i " << i << " j " << j << " bPLAS_TAMEX_ID " <<bPLAS_TAMEX_ID << endl;    }
                    }   
    //                  
                    
                   // if(i>15 && pInput->fbPlas_Lead_PMT[31][j]>0 && pInput->fbPlas_Lead_PMT[15][j]>0) {    
                
                
                    // lead_lead_bplas_Ref0[i][j] = (pInput->fbPlas_Lead_PMT[15][j] -  pInput->fbPlas_Lead_PMT[31][j]);   
            
        
                    ///Fatima   
                        
                    
                
                    if(i<8 && Fat_RefCh[j]>0 && pInput->fbPlas_Lead_PMT[i][j]>0){   
                
                  lead_lead_fat_Ref0[i][j] = (Fat_RefCh[j] -  pInput->fbPlas_Lead_PMT[i][j]);   
                    
                  if(lead_lead_fat_Ref0[i][j]!=0){  
                     //cout<<"lead_lead_fat_Ref0[i][j] " << lead_lead_fat_Ref0[i][j] << " Fat_RefCh[j] " <<Fat_RefCh[j]<< " pInput->fbPlas_Lead_PMT[i][j] " << pInput->fbPlas_Lead_PMT[i][j] <<" pInput->fbPlas_Lead_PMT[1][j]  " << pInput->fbPlas_Lead_PMT[1][j]  <<" i " << i << " j " << j <<   endl;   
                        
                     hFATTEMP_lead_lead_ref[i] ->Fill(lead_lead_fat_Ref0[i][j]);    
                        }   
                    }   
                        
     ///**--------------------------------------------SCI41 Comparisons -------------------------------------------------**/                     
                           hSC41_Analogue_Tamex->Fill(bPlas_TAM_SC41L_ANA[j] - bPlas_TAM_SC41R_ANA[j]); 
                           hSC41_Digital_Tamex->Fill(bPlas_TAM_SC41L_DIG[j] - bPlas_TAM_SC41R_DIG[j]);  
                            
                      if (i<8){ 
                            ///For Fatima SC41L Analogue    
                                
               if (bPlas_TAM_SC41L_ANA[j] >0 && lead_bplas1[i][j] >0) SC41L_ANA_lead_fat[i][j] = (bPlas_TAM_SC41L_ANA[j] - lead_bplas1[i][j]);  
                    
                if(SC41L_ANA_lead_fat[i][j]>0)  hFATTEMP_SC41L_Anal_lead[FAT_TAMEX_ID][i]->Fill(SC41L_ANA_lead_fat[i][j]);  
                            
                                
                            ///For Fatima SC41R Analogue    
                            if (bPlas_TAM_SC41R_ANA[j] >0 && lead_bplas1[i][j] >0)SC41R_ANA_lead_fat[i][j] = (bPlas_TAM_SC41R_ANA[j] - lead_bplas1[i][j]);  
                             hFATTEMP_SC41R_Anal_lead[FAT_TAMEX_ID][i]->Fill(SC41R_ANA_lead_fat[i][j]); 
                                
                             ///For Fatima SC41L Digital    
                            if (bPlas_TAM_SC41L_DIG[j] >0)SC41L_DIG_lead_fat[i][j] = (bPlas_TAM_SC41L_DIG[j] - lead_bplas1[i][j]);  
                           //  hFATTEMP_SC41L_Anal_lead[FAT_TAMEX_ID][i]->Fill(SC41L_DIG_lead_fat[i][j]);   
                                
                              ///For Fatima SC41R Digital   
                            if (bPlas_TAM_SC41R_DIG[j] >0)SC41R_DIG_lead_fat[i][j] = (bPlas_TAM_SC41R_DIG[j] - lead_bplas1[i][j]);  
                          //   hFATTEMP_SC41R_Anal_lead[FAT_TAMEX_ID][i]->Fill(SC41R_DIG_lead_fat[i][j]);   
                             }  
                                
                          if(i>15){ 
                         ///For bPlas SC41L Analogue    
                           if (bPlas_TAM_SC41L_ANA[j] >0)SC41L_ANA_lead_bPlas[i][j] = (bPlas_TAM_SC41L_ANA[j] - lead_bplas1[i][j]);     
                           // hbPlas_SC41L_Anal_lead[bPLAS_TAMEX_ID][i]->Fill(SC41L_ANA_lead_bPlas[i][j]);  
                                
                         ///For bPlas SC41R Analogue    
                           if (bPlas_TAM_SC41R_ANA[j] >0)SC41R_ANA_lead_bPlas[i][j] = (bPlas_TAM_SC41R_ANA[j] - lead_bplas1[i][j]);     
                            // hbPlas_SC41R_Anal_lead[bPLAS_TAMEX_ID][i]->Fill(SC41R_ANA_lead_bPlas[i][j]); 
                                
                          ///For bPlas SC41L Digital    
                           if (bPlas_TAM_SC41L_DIG[j] >0)SC41L_DIG_lead_bPlas[i][j] = (bPlas_TAM_SC41L_DIG[j] - lead_bplas1[i][j]);     
                           //  hbPlas_SC41L_Digi_lead[bPLAS_TAMEX_ID][i]->Fill(SC41L_DIG_lead_bPlas[i][j]); 
                                
                          ///For bPlas SC41R Digital    
                           if (bPlas_TAM_SC41R_DIG[j] >0)SC41R_DIG_lead_bPlas[i][j] = (bPlas_TAM_SC41R_DIG[j] - lead_bplas1[i][j]);     
                          //   hbPlas_SC41R_Digi_lead[bPLAS_TAMEX_ID][i]->Fill(SC41R_DIG_lead_bPlas[i][j]);     
                                
                          } 
                            
                        
                    ///all Channels compare per detector system 
                    for (int k=0; k< 48; k++){  
                        if(i != k && hits_bplas_lead>1){    
                                
                             lead_bplas2[k][j] = pInput->fbPlas_Lead_PMT[k][j];     
        
        
        
        
                             if (i>15 && k>15 &&  lead_bplas1[i][j] >0 && lead_bplas2[k][j]>0) lead_lead_bplas[i][j] = (lead_bplas1[i][j] -lead_bplas2[k][j]); ///For PLASTIC   
                            // cout<<"event_number " << pInput->fevent_number << " lead_bplas1[i][j] " << lead_bplas1[i][j] << " lead_bplas1[k][j] " << lead_bplas1[k][j] << " i " << i << " k " << k <<  " lead_lead_bplas[i][j] " << lead_lead_bplas[i][j] << endl;}  
                                
                             if (i<8 && k<8 ) lead_lead_fat[i][j] = (lead_bplas1[i][j] -lead_bplas2[k][j]); ///For FATIMA   
        
                             lead_lead_bplas[i][j]  = CYCLE_TIME*lead_lead_bplas[i][j]; 
        
                             if(lead_lead_bplas[i][j]!=0 && lead_bplas1[i][j]>0 && lead_bplas2[k][j]>0)  {  
                     if(i>15) hbPlas_lead_lead[bPLAS_TAMEX_ID][i-16] -> Fill(lead_lead_bplas[i][j]);    
                  } 
                                if(lead_lead_fat[i][j]!=0 && lead_bplas1[i][j]>0 && lead_bplas2[k][j]>0 && i<8)  {  
                              hFATTEMP_lead_lead[FAT_TAMEX_ID][i] -> Fill(lead_lead_fat[i][j]);         
                                }                               
                            }   
                          } 
                       }    
        
                //Trail T   
                for(int j=0; j< pInput->fbPlas_PMT_Trail_N[i]; j++){    
                    trail_bplas1[i][j] = pInput->fbPlas_Trail_PMT[i][j];    
                    hits_bplas_trail++; 
                    pOutput->pbPlas_TrailT[i][j] = trail_bplas1[i][j];  
                    pOutput->pbPlas_TrailHits = hits_bplas_trail;   
                    for (int k=0; k< 32; k++){  
                        if(i != k && hits_bplas_trail>1){   
                             trail_bplas2[i][j] = pInput->fbPlas_Trail_PMT[k][j];   
                             if(trail_bplas2[i][j]>0 && j>0){   
        
                trail_trail_bplas[i][j] = (trail_bplas1[i][j] - trail_bplas2[k][j]);    
    //              if (ABS(trail_trail_bplas[i][j]) > (double)(COARSE_CT_RANGE>>1))        // overflow 
    //                             {    
    //                             trail_trail_bplas[i][j] = CYCLE_TIME*(trail_trail_bplas[i][j] + COARSE_CT_RANGE) ;   
    //                         }    
    //                         else {   
                                 trail_trail_bplas[i][j]  = CYCLE_TIME*trail_trail_bplas[i][j] ;    
    //                         }    
                    
              /*  if(i>15) hbPlas_trail_trail[bPLAS_TAMEX_ID][i-16] -> Fill(trail_trail_bplas[i][j]);       
                else  hbPlas_trail_trail[FAT_TAMEX_ID][i] -> Fill(trail_trail_bplas[i][j]);    */   
                            }   
                        }   
                    }   
                }   
                //ToT (~Energy) 
              for(int j=0; j< pInput->fbPlas_PMT_Lead_N[i]; j++){   
                  if(pInput->fbPlas_Trail_PMT[i][j] >0 && pInput->fbPlas_Lead_PMT[i][j]>0){ 
                
                      ToT_bplas[i][j] = (pInput->fbPlas_Trail_PMT[i][j] - pInput->fbPlas_Lead_PMT[i][j]);   
                    
                ///Correction for overflows 
                if(ABS( ToT_bplas[i][j]) >(double)(COARSE_CT_RANGE>>1)) {   
                    
                       ToT_bplas[i][j] = CYCLE_TIME*(ToT_bplas[i][j] + COARSE_CT_RANGE);    
                      } 
                 else{  
                           ToT_bplas[i][j]= CYCLE_TIME*ToT_bplas[i][j];                         
                       }    
                                
                      if(ToT_bplas[i][j]>0) {   
                          hbPlas_hit_pattern->Fill(i);  
                          bPlas_tot_hits++; 
                
                        //   if(i>31)      cout<<"i " << i <<" ToT_bplas[i][j] " << ToT_bplas[i][j]<<  endl;    
                    
                        
                        
             if(bplas_detnum==1)   hbPlas_ToT_det[bplas_detnum][i-16] ->Fill(ToT_bplas[i][j]);  
             if(bplas_detnum==2) hbPlas_ToT_det[bplas_detnum][i-32] ->Fill(ToT_bplas[i][j]);    
                        
                   // else {    
                   //     if(ToT_bplas[i][j]>0){    
                            
                      if(i<16)  hFATTEMP_ToT[FAT_TAMEX_ID][i] ->Fill(ToT_bplas[i][j]); //raw Energy 
                      //  if(i==10 || i ==11){  
                      //    cout<<"ANALYSIS SCI " <<pOutput->pEvent_Number <<  " ToT_bplas[i][j] " << ToT_bplas[i][j] << " i " << i << endl;}   
                        }   
                  //  } 
                        
                  ///Gain matching  
                  pOutput-> pbPlas_ToTCalib[i] = fCal->Abplas_TAMEX[i]* ToT_bplas[i][j] + fCal->Bbplas_TAMEX[i];    
            //cout<<" ToT_bplas[i][j] " << ToT_bplas[i][j] << " pbPlas_ToTCalib[i] " <<pOutput->pbPlas_ToTCalib[i] << " fCal->Abplas_TAMEX[i] " <<  fCal->Abplas_TAMEX[i] << " fCal->Bbplas_TAMEX[i] " <<fCal->Bbplas_TAMEX[i] <<" i " << i <<  endl; 
            
                  if(i>15) {    
                        
               if(bplas_detnum==1)  {   
                   hbPlas_Energy_Calib[bplas_detnum][i-16]->Fill(pOutput-> pbPlas_ToTCalib[i]); 
               hbPlas_ToT_Sum_Det1->Fill(pOutput-> pbPlas_ToTCalib[i]); 
               }    
                
               if(bplas_detnum==2) {    
                   hbPlas_Energy_Calib[bplas_detnum][i-32]->Fill(pOutput-> pbPlas_ToTCalib[i]);                 
                   hbPlas_ToT_Sum_Det2->Fill(pOutput-> pbPlas_ToTCalib[i]); 
               }    
                    //cout<<"1 event " << pInput->fevent_number<<" pOutput-> pbPlas_ToTCalib[i]/5 " << pOutput-> pbPlas_ToTCalib[i]/5 << " i " << i << endl;    
                            
                  } 
                  else {    
                        
                    hFATTEMP_Energy_Calib[FAT_TAMEX_ID][i]->Fill(pOutput-> pbPlas_ToTCalib[i]); 
        
                  } 
                  ///Lead-Lead Energy gating    
               for (int k=0; k< 32; k++){   
                   if(i != k && pOutput->pbPlas_LeadHits>1 ){   
                       if(pOutput->pbPlas_LeadT[k][j]>0){   
                             ///Temp E gate 
                             if(ToT_bplas[0][j]>27806 && ToT_bplas[0][j]<28537){    
                            
                  if(pOutput->pbPlas_LeadT[i][j] >0 && pOutput->pbPlas_LeadT[k][j] >0){ 
                lead_lead_bplas_gated[i][j] = (pOutput->pbPlas_LeadT[i][j]  -  pOutput->pbPlas_LeadT[k][j])*5000; ///into ps    
                    
               if(i>15) hbPlas_lead_lead_gated[bPLAS_TAMEX_ID][i-16] -> Fill(lead_lead_bplas_gated[i][j]);  
               else hFATTEMP_lead_lead_gated[FAT_TAMEX_ID][i] -> Fill(lead_lead_bplas_gated[i][j]); 
                     }  
                   }    
                 }  
               }        
             }  
           }    
         }      
       }    
       /*for(int i=0; i<32; i++){   
        if(i>15 &&bPlas_tot_hits>1) fired_det1=true;    
        if(i<16 &&bPlas_tot_hits>1) fired_det2=true;    
    //     cout<<"event " << pInput->fevent_number << " pOutput-> pbPlas_ToTCalib[i] "<<pOutput-> pbPlas_ToTCalib[i] << " i " << i <<" fired_det1 " << fired_det1 << " fired_det2 " << fired_det2 <<" bPlas_tot_hits " << bPlas_tot_hits<<endl; 
    //     if(fired_det1==true && pOutput-> pbPlas_ToTCalib[i]>0)   
    //         hbPlas_ToT_coinc[bPLAS_TAMEX_ID][i] ->Fill(pOutput-> pbPlas_ToTCalib[i]/5); //coincident Energy  
    //      
    //     if(fired_det2==true && pOutput-> pbPlas_ToTCalib[i]>0)   
    //         hbPlas_ToT_coinc[FAT_TAMEX_ID][i-16] ->Fill(pOutput-> pbPlas_ToTCalib[i]/5); //coincident Energy 
    //      
    //     if(fired_det1==true &&fired_det2==true ){    
    //       if(i<16)  hbPlas_ToT_coinc_dets[bPLAS_TAMEX_ID][i] ->Fill(pOutput-> pbPlas_ToTCalib[i]/5); 
    //       else  hbPlas_ToT_coinc_dets[FAT_TAMEX_ID][i-16] ->Fill(pOutput-> pbPlas_ToTCalib[i]/5);    
    //          
    //     }    
    //     if(pOutput-> pbPlas_ToTCalib[i]/5>0 && pOutput-> pbPlas_ToTCalib[i+16]/5>0){ 
    //        if(i<16) hbPlas_ToT_coinc_dets_perchan[bPLAS_TAMEX_ID][i] ->Fill(pOutput-> pbPlas_ToTCalib[i]/5); 
    //        else hbPlas_ToT_coinc_dets_perchan[FAT_TAMEX_ID][i-16] ->Fill(pOutput-> pbPlas_ToTCalib[i]/5);    
    //          
    //        // cout<<"3 event " << pInput->fevent_number<< " i " << i <<" ToTCalib[i] " << pOutput-> pbPlas_ToTCalib[i]/5 << " ToTCalib[i+16] " << pOutput-> pbPlas_ToTCalib[i+16]/5 <<  endl;    
    //     }    
       }*/  
        
      hbPlas_num_fired_chans->Fill(bPlas_tot_hits); 
        
            
       hits_bplas_lead=0;   
       //bPlas_tot_hits=0;  
         
     }
  /**----------------------------------------------------------------------------------------------**/
 /**--------------------------------------  FATIMA VME (QDC) AND TAMEX -----------------------------**/
/**----------------------------------------------------------------------------------------------**/
// void EventAnlProc::Make_Fatima_VME_Tamex_Histos(){
//
//
//     for (int i =0; i<50; i++)
//   {
//     hFAT_ToT[i] = MakeTH1('D', Form("FATIMA/TAMEX/ToT/ToT_Ch.%02d", i), Form("TOT Detector %2d", i), 16250, 0., 65000.);
//     hFAT_lead_lead[i] = MakeTH1('D', Form("FATIMA/TAMEX/Lead-Lead/Lead-LeadCh.%02d", i), Form("Lead - Lead Ch. %2d", i),2500, -50000., 50000);
//     hFAT_lead_lead_ref[i] = MakeTH1('D', Form("FATIMA/TAMEX/Lead-Lead_Ref/Lead-LeadRefCh.%02d", i), Form("Lead Ref Ch.0 - Lead Ch. %2d", i),2500, -50000., 50000.);
//
//
//     hFAT_lead_lead_QDC_Gate[i] = MakeTH1('D', Form("FATIMA/Combined/Lead-Lead_QDCGate/Lead-Lead_QDCGateCh.%02d", i), Form("Lead - Lead QDC energy Gated Ch. %2d", i),2500, -50000., 50000);
// //     hFAT_lead_lead_gated[i] = MakeTH1('D', Form("FATIMA/TAMEX/Lead-Lead_Egated/Lead-Lead_Egated_Ch.%02d", i), Form("Lead - Lead Energy gated Ch.  %2d", i),2500, -50000., 50000.);
// //     hFAT_lead_lead_gated1[i] = MakeTH1('D', Form("FATIMA/TAMEX/Lead-Lead_Egated1/Lead-Lead_Egated1_Ch.%02d", i), Form("Lead - Lead Energy gated Ch.  %2d", i),2500, -50000., 50000.);
// //     hFAT_lead_lead_gated2[i] = MakeTH1('D', Form("FATIMA/TAMEX/Lead-Lead_Egated2/Lead-Lead_Egated2_Ch.%02d", i), Form("Lead - Lead Energy gated Ch.  %2d", i),2500, -50000., 50000.);
//
//     hFAT_lead_lead_energy[i] = MakeTH2('D',Form("FATIMA/TAMEX/Lead-Lead_energy/Energy_vs._dT/Energy_vs._dT_Ch.%02d", i),Form("Fatima Energy vs Lead-Lead.%02d", i),500,0,2000, 160, -40000., 40000);
//
//     hFAT_trail_trail[i] = MakeTH1('D', Form("FATIMA/TAMEX/Trail-Trail/trail-trailCh.%02d", i), Form("Trail - Trail Ch. %2d", i), 2500, -50000., 50000.);
//
//    // hFAT_Sc41lead_leadmaxtot[i] = MakeTH1('D', Form("FATIMA/Sc41-LeadMaxToT/SC41Lead_LeadCh.%02d", i), Form("SC41 Lead - (max ToT chan) Lead %2d ", i), 4002, -100000., 100000.);
//
//     ///VME
//
//     hFAT_QDCCalib1[i] =  MakeTH1('D', Form("FATIMA/VME/Energy/EnergyCalib/LaBr_ECalib_Ch.%2d",i), Form("QDC Calib Ch. %2d",i), 4000,0,4000);
//     }
//
//     hFAT_QDCCalib1Sum = MakeTH1('D', "FATIMA/VME/Energy/Fat_VME_EnergySum", "LaBr Energy (all detectors)",40000,0,40000);
//     hFAT_hits_QDC       = MakeTH1('D', "FATIMA/VME/Stats/QDC_FAThits", "bPlastic hit pattern QDC1",50,0,50);
//
//
//     hFAT_gamma_gamma = MakeTH2('D', "FATIMA/Gamma-Gamma/Sum", "FATIMA Gamma-Gamma (all detectors)",6500,0,65000, 6500,0,65000);
//     }
 /**----------------------------------------------------------------------------------------------**/

// void EventAnlProc::Do_Fatima_VME_Tamex_Histos(EventUnpackStore* pInput, EventAnlStore* pOutput){
//     double lead_lead_fat[50], trail_trail_fat[50];
//     double ToT_fat[50];
//     double lead_fat1[50], lead_fat2[50] ;
//     double trail_fat1[50], trail_fat2[50];
//     double  lead_lead_fat_Ref0[50];
//     double lead_lead_fat_gated[50], lead_lead_fat_gated1[50], lead_lead_fat_gated2[50];
//     int hits_fat_lead = 0, hits_fat_trail=0;
//     double Fat_QDC_GainMatch[50];
//
//
//     int Fat_QDC_IDMain_i = 0;
//     for (int i = 0; i < 50; i++)
//         {
//              Fat_QDC_GainMatch[i] = 0;
//              lead_lead_fat[i]=0;
//              lead_lead_fat_Ref0[i]=0;
//              trail_trail_fat[i]=0;
//              ToT_fat[i] = 0;
//              trail_fat1[i] = 0;
//              trail_fat2[i] = 0;
//              lead_fat1[i] =0;
//              lead_fat2[i] =0;
//
//              lead_lead_fat_gated[i] = 0;
//              lead_lead_fat_gated1[i] = 0;
//              lead_lead_fat_gated2[i] = 0;
//
//         }
//         pOutput->pFat_firedQDC_Comb = pInput->fFat_firedQDC;
//         ///VME QDC
//         for (int i=0; i< pInput->fFat_firedQDC; i++){
//
//         Fat_QDC_IDMain_i = pInput->fFat_QDC_ID[i];
//         pOutput->pFat_QDC_ID_Comb[i] = pInput->fFat_QDC_ID[i];
//        // cout<<"pInput->fFat_QDC_E[i] " << pInput->fFat_QDC_E[i] << " pInput->fFat_firedQDC " << pInput->fFat_firedQDC <<" Fat_QDC_IDMain_i " << Fat_QDC_IDMain_i <<endl;
//
//         ///Gainmatched Energy
//         Fat_QDC_GainMatch[Fat_QDC_IDMain_i] = fCal->Afat[Fat_QDC_IDMain_i]* pow( pInput->fFat_QDC_E[i],3) + fCal->Bfat[Fat_QDC_IDMain_i]* pow( pInput->fFat_QDC_E[i],2) + fCal->Cfat[Fat_QDC_IDMain_i]* pInput->fFat_QDC_E[i] + fCal->Dfat[Fat_QDC_IDMain_i];
//
//          hFAT_QDCCalib1[i]->Fill(Fat_QDC_GainMatch[Fat_QDC_IDMain_i]);
//          hFAT_QDCCalib1Sum->Fill(Fat_QDC_GainMatch[Fat_QDC_IDMain_i]);
//          hFAT_hits_QDC->Fill(Fat_QDC_IDMain_i);
//         pOutput-> pFat_QDC_E_Comb[Fat_QDC_IDMain_i] = Fat_QDC_GainMatch[Fat_QDC_IDMain_i];
//         }
//         ///Loop on channels
//       for (int i = 0; i < 50; i++)
//         {
//             ///Lead T
//             for(int j=0; j< pInput->fFat_PMT_Lead_N[i]; j++){
//                 lead_fat1[i] = pInput->fFat_Lead_PMT[i][j];
//                 hits_fat_lead++;
//                 pOutput->pFat_LeadT[i][j] = lead_fat1[i];
//                 pOutput->pFat_LeadHits = hits_fat_lead;
//
//                 ///Ref channel 0
//                 if(pInput->fFat_Lead_PMT[0][j]>0 && pInput->fFat_Lead_PMT[i][j]>0&& hits_fat_lead>1){
//                 lead_lead_fat_Ref0[i] = (pInput->fFat_Lead_PMT[0][j] -  pInput->fFat_Lead_PMT[i][j])*5000;
//
//                 hFAT_lead_lead_ref[i] ->Fill(lead_lead_fat_Ref0[i]);
//                 }
//
//                 ///all Dets compare
//                 for (int k=0; k< 50; k++){
//                     if(i != k && hits_fat_lead>1){
//                          lead_fat2[k] = pInput->fFat_Lead_PMT[k][j];
//
//                     lead_lead_fat[i] = (lead_fat1[i] - lead_fat2[k]);
//                         if (ABS(lead_lead_fat[i]) > (double)(COARSE_CT_RANGE>>1))        // overflow
//                             {
//                             lead_lead_fat[i] = CYCLE_TIME*(lead_lead_fat[i]  + COARSE_CT_RANGE);
//                         }
//                          else {
//                              lead_lead_fat[i]  = CYCLE_TIME*lead_lead_fat[i];
//                          }
//
//                          if(lead_lead_fat[i]!=0 && lead_fat1[i]>1 && lead_fat2[k]>1 ) {
//                              ///Fill Lead-Lead
//                              hFAT_lead_lead[i] -> Fill(lead_lead_fat[i]);
//
//                        ///      QDC Energy Gates
//                     if( pOutput-> pFat_QDC_E_Comb[0]>200 && pOutput-> pFat_QDC_E_Comb[0]<220 && pOutput-> pFat_QDC_E_Comb[1]>235 && pOutput-> pFat_QDC_E_Comb[1]<255){
//                                  hFAT_lead_lead_QDC_Gate[i] -> Fill(lead_lead_fat[i]);
//
//                          }
//          //cout<<"event " << event_number << " lead_lead_fat[i] " <<lead_lead_fat[i]<<" lead_fat1[i] " << lead_fat1[i] << " lead_fat2[k] "<<lead_fat2[k]  << " i " << i <<endl;
//                       }
//                     }
//                 }
//             }
//             ///Trail T
//             for(int j=0; j< pInput->fFat_PMT_Trail_N[i]; j++){
//                 trail_fat1[i] = pInput->fFat_Trail_PMT[i][j];
//                 hits_fat_trail++;
//                 pOutput->pFat_TrailT[i][j] = trail_fat1[i];
//                 pOutput->pFat_TrailHits = hits_fat_trail;
//                 for (int k=0; k< 50; k++){
//                     if(i != k && hits_fat_trail>1){
//                          trail_fat2[i] = pInput->fFat_Trail_PMT[k][j];
//                          if(trail_fat2>0 && j>0){
//
//             trail_trail_fat[i] = (trail_fat1[i] - trail_fat2[i]);
//             if (ABS(trail_trail_fat[i]) > (double)(COARSE_CT_RANGE>>1))        // overflow
//                             {
//                             trail_trail_fat[i] = CYCLE_TIME*(trail_trail_fat[i] + COARSE_CT_RANGE);
//                         }
//              else {
//                             trail_trail_fat[i]  = CYCLE_TIME*trail_trail_fat[i];
//                         }
//
//             hFAT_trail_trail[i] -> Fill(trail_trail_fat[i]);
//                         }
//                     }
//                 }
//             }
//             ///ToT (~Energy)
//           for(int j=0; j< pInput->fFat_PMT_Lead_N[i]; j++){
//            //   if(pInput->fFat_Trail_PMT[i][j] >0 && pInput->fFat_Lead_PMT[i][j]>0){
//                 ToT_fat[i] = (pInput->fFat_Trail_PMT[i][j] - pInput->fFat_Lead_PMT[i][j]);
//
//                 ///Correction for overflows
//             if(ABS(ToT_fat[i]) >(double)(COARSE_CT_RANGE>>1)) {
//                    ToT_fat[i] = CYCLE_TIME*(ToT_fat[i] + COARSE_CT_RANGE);
//             }
//              else{
//                         ToT_fat[i]= CYCLE_TIME*ToT_fat[i];
//                    }
// //
//  //cout<<"pInput->fFat_Trail_PMT[i][j]  " <<pInput->fFat_Trail_PMT[i][j] << " pInput->fFat_Lead_PMT[i][j] " << pInput->fFat_Lead_PMT[i][j] <<" ToT_fat[i] " << ToT_fat[i] << endl;
//                     hFAT_ToT[i] ->Fill(ToT_fat[i]); //raw Energy
//
//       }
//    }
// }


 /**----------------------------------------------------------------------------------------------**/
 /**--------------------------------------  FATIMA VME ----------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
void EventAnlProc::FAT_det_pos_setup(){

    FAT_positions   = new double*[36];
    FAT_neighbour_check = new bool*[36];
    FAT_angle_diffs = new double*[36];

    for(int i = 0; i < 36; ++i){
    FAT_positions[i] = new double[3];
    FAT_angle_diffs[i] = new double[36];
    FAT_neighbour_check[i] = new bool[36];
    for (int j = 0; j < 3; ++j) FAT_positions[i][j] = -1;
    for (int k = 0; k < 36; ++k){

        FAT_neighbour_check[i][k] = true;

        FAT_angle_diffs[i][k] = -1;

            }
    }

    const char* format = "%d %lf %lf %lf";

    ifstream file("Configuration_Files/FATIMA_Detector_Positions.txt");

    if(file.fail()){
        cerr << "Could not find FATIMA Detector Positions File!" << endl;
        exit(0);
    }

    string line;
    int pos_num;
    double r, theta, phi;

    while(file.good()){
        getline(file,line,'\n');
        if(line[0] == '#') continue;
        sscanf(line.c_str(),format, &pos_num, &r, &theta, &phi);

        FAT_positions[pos_num][0] = r;
        FAT_positions[pos_num][1] = theta;
        FAT_positions[pos_num][2] = phi;

    }

    if(FAT_nearest_neighbour_exclusion){

    for(int i = 0; i < 36; ++i){

        if(i%12 == 11) FAT_neighbour_check[i][(i-11)] = false; // Same Ring Rignt
        else FAT_neighbour_check[i][(i+1)] = false; // Same Ring Right
        if(i%12 == 0) FAT_neighbour_check[i][i+11] = false; // Same Ring Left
        else FAT_neighbour_check[i][(i-1)] = false; // Same Ring Left

        if(!same_ring_exclusion){

        if(i < 12){

            FAT_neighbour_check[i][(i+12)] = false; // Middle Ring Left

            if(i == 12) FAT_neighbour_check[i][(i+1)] = false; // Middle Ring Left for 11

            else FAT_neighbour_check[i][(i+13)] = false; // Middle Ring Right
        }
        if(i > 11 && i < 24){

            FAT_neighbour_check[i][(i+12)] = false; // Upper Outer Ring
            FAT_neighbour_check[i][(i-12)] = false; // Lower Outer Ring

            if(i == 12){
             FAT_neighbour_check[i][(i+23)] = false; // Upper Outer Ring
             FAT_neighbour_check[i][(i-1)] = false; // Lower Outer Ring
            }
            else{
            FAT_neighbour_check[i][(i+11)] = false; // Upper Outer Ring
            FAT_neighbour_check[i][(i-13)] = false; // Lower Outer Ring
            }
        }
        if(i > 23){

            FAT_neighbour_check[i][(i-12)] = false; // Middle Ring Left

            if(i == 35) FAT_neighbour_check[i][(i-23)] = false; // Middle Ring Left for 35

            else FAT_neighbour_check[i][(i-11)] = false; // Middle Ring Right
                }
            }
        }
    }

    ofstream output_position_matrix_file;
    output_position_matrix_file.open ("Configuration_Files/FATIMA_Exclusion_Matrix.txt");
    cout<<endl;
    cout << "============================================================" << endl;
    cout << "A Matrix of excluded detector pairings can be found in" << endl;
    cout << "'Configuration_Files/FATIMA_Exclusion_Matrix.txt'"<<endl;
    cout << "============================================================" << endl;
    cout<<endl;

    if (output_position_matrix) output_position_matrix_file <<"        "<<"0 "<<"1 "<<"2 "<<"3 "<<"4 "<<"5 "<<"6 "<<"7 "<<"8 "
            <<"9 "<<"10 "<<"11 "<<"12 "<<"13 "<<"14 "<<"15 "<<"16 "<<"17 "
            <<"18 "<<"19 "<<"20 "<<"21 "<<"22 "<<"23 "<<"24 "<<"25 "<<"26 "
            <<"27 "<<"28 "<<"29 "<<"30 "<<"31 "<<"32 "<<"33 "<<"34 "<<"35 "<<endl;

    for(int i = 0; i < 36; ++i){

    if (i >= 10 && output_position_matrix) output_position_matrix_file <<"Det "<<i<<": ";
    if (i < 10  && output_position_matrix) output_position_matrix_file <<"Det "<<i<<" : ";

    for (int k = 0; k < 36; ++k){

        if(k > 9 && output_position_matrix) output_position_matrix_file<<" ";

        double dist = distance_between_detectors( FAT_positions[i][0],  FAT_positions[i][1],  FAT_positions[i][2],
                              FAT_positions[k][0],  FAT_positions[k][1],  FAT_positions[k][2]);

        double angle = angle_between_detectors(FAT_positions[i][0], FAT_positions[k][0], dist);

        FAT_angle_diffs[i][k] = angle;

        if((dist < FAT_exclusion_dist && (((i < 12 && k < 12) ||
                        (i < 24 && i > 11 && k < 24 && k > 11) ||
                        (i > 23 && k > 23)) || !same_ring_exclusion )) || i == k ){


         FAT_neighbour_check[i][k] = false;

        }


        if (output_position_matrix && !FAT_neighbour_check[i][k]) output_position_matrix_file<<"X ";

        else if(output_position_matrix && FAT_neighbour_check[i][k]) output_position_matrix_file<<"0 ";

    }

    if (output_position_matrix) output_position_matrix_file<<endl;

    }

    output_position_matrix_file.close();

}
//-----------------------------------------------------------------------------------------------------------------------------//

double EventAnlProc::distance_between_detectors(double _r, double _theta, double _phi, double r_, double theta_, double phi_){

    _theta = _theta * M_PI/180.0;
    theta_ = theta_ * M_PI/180.0;

    _phi = _phi * M_PI/180.0;
    phi_ = phi_ * M_PI/180.0;

    double dist = sqrt(_r*_r + r_*r_ - 2.0*_r*r_*(sin(_theta)*sin(theta_)*cos(_phi - phi_) + cos(_theta)*cos(theta_)));

    return dist;


}
//-----------------------------------------------------------------------------------------------------------------------------//
double EventAnlProc::angle_between_detectors(double _r, double r_, double dist_){


    double angle_diff = acos((_r*_r + r_*r_ - dist_*dist_)/(2.0*_r*r_));

    angle_diff = angle_diff * 180.0/M_PI;

    return angle_diff;

}
//-----------------------------------------------------------------------------------------------------------------------------//

void EventAnlProc::Make_Fatima_Histos(){


  for (int i=0; i<50; i++){
    hFAT_QDCCalib1[i] =  MakeTH1('D', Form("FATIMA_VME/Energy/EnergyCalib/LaBr_ECalib_Ch.%2d",i), Form("QDC Calib Ch. %2d",i), 4000,0,4000);
   // hFAT_QDCdt[i]   = MakeTH1('D', Form("FATIMA_VME/Timing/QDCdt/QDCdt%2d",i), Form("QDCdT Ch.%2d",i), 3201,-40,40);
    //hFAT_TDCCalib1[i] =  MakeTH1('D', Form("FATIMA/Timing/TDCCalib/LaBr_Tcalib%2d",i), Form("TDC channel Calib %2d",i), 1E5,0,2E5);
    hFAT_TDCdt_refSC41[i] = MakeTH1('D', Form("FATIMA_VME/Timing/TDCdt_SC41-FatTDC/TDCdT_SC41_LaBr%02d", i), Form("TDC dtSC41 All Multip SC41- LaBr%02d", i),4E4,0,2E6);
//     hFAT_TDCdt_refSC41_M1[i] = MakeTH1('D', Form("FATIMA/Timing/TDCdt_SC41-FatTDC_M1/TDCdT_SC41_M1_LaBr%02d", i), Form("TDC dtSC41 Multip 1 SC41- LaBr%02d", i),25000,-50000,50000);
//     hFAT_TDCdt_refSC41_M2[i] = MakeTH1('D', Form("FATIMA/Timing/TDCdt_SC41-FatTDC_M2/TDCdT_SC41_M2_LaBr%02d", i), Form("TDC dtSC41 Multip 2 SC41- LaBr%02d", i),25000,-50000,50000);
//     hFAT_TDCdt_refSC41_M3[i] = MakeTH1('D', Form("FATIMA/Timing/TDCdt_SC41-FatTDC_M3/TDCdT_SC41_M2+_LaBr%02d", i), Form("TDC dtSC41 Multip 2+ SC41- LaBr%02d", i),25000,-50000,50000);
//
    hFAT_TDCdt_refSC41_gated[i] = MakeTH1('D', Form("FATIMA_VME/Timing/TDCdt_SC41-FatTDC_EGated/TDCdT_Egated_SC41_LaBr%02d", i), Form("TDC Gamma gated dtSC41 SC41- LaBr%02d", i),4E4,0,2E6);
/*    hFAT_TDCdt_refSC41_M1_gated[i] = MakeTH1('D', Form("FATIMA/Timing/EGated/TDCdt_SC41-FatTDC/TDCdT_Egated_SC41_M1_LaBr%02d", i), Form("TDC Gamma gated dtSC41 SC41- LaBr%02d", i),4000,-1000,1000);
    hFAT_TDCdt_refSC41_M2_gated[i] = MakeTH1('D', Form("FATIMA/Timing/EGated/TDCdt_SC41-FatTDC/TDCdT_Egated_SC41_M2_LaBr%02d", i), Form("TDC Gamma gated dtSC41 SC41- LaBr%02d", i),4000,-1000,1000);
    hFAT_TDCdt_refSC41_M3_gated[i] = MakeTH1('D', Form("FATIMA/Timing/EGated/TDCdt_SC41-FatTDC/TDCdT_Egated_SC41_M2+_LaBr%02d", i), Form("TDC Gamma gated dtSC41 SC41- LaBr%02d", i),4000,-1000,1000);
   */

    hFAT_TDCdt_refCha[i] = MakeTH1('D', Form("FATIMA_VME/Timing/TDCdT_TDC0-TDC/TDCdT_Cha_LaBr%02d_LaBr%02d", 0, i), Form("TDC dt Channel All Multip LaBr%02d - LaBr%02d",0 , i),250,-2E4,2E4);
    hFAT_TDCdt_refCha_gated[i] = MakeTH1('D', Form("FATIMA_VME/Timing/TDCdT_TDC0-TDC_EGated/TDCdT_Cha_EGated_LaBr%02d_LaBr%02d", 0, i), Form("TDC dt Channel All Multip LaBr%02d - LaBr%02d",0 , i),250,-2E4,2E4);
    //     hFAT_TDCdt_refCha_M1[i] = MakeTH1('D', Form("FATIMA/Timing/TDCdT_TDC0-TDC_M1/TDCdT_Cha_M1_LaBr%02d_LaBr%02d", 0, i), Form("TDC dt Channel Multip 1 LaBr%02d - LaBr%02d",0 , i),4000,-1000,1000);
//     hFAT_TDCdt_refCha_M2[i] = MakeTH1('D', Form("FATIMA/Timing/TDCdT_TDC0-TDC_M2/TDCdT_Cha_M2_LaBr%02d_LaBr%02d", 0, i), Form("TDC dt Channel Multip 2 LaBr%02d - LaBr%02d",0 , i),4000,-1000,1000);
//     hFAT_TDCdt_refCha_M3[i] = MakeTH1('D', Form("FATIMA/Timing/TDCdT_TDC0-TDC_M2+/TDCdT_Cha_M2+_LaBr%02d_LaBr%02d", 0, i), Form("TDC dt Channel Multip 2+ LaBr%02d - LaBr%02d",0 , i),4000,-1000,1000);
//
    hFAT_TDC_Multipl_ch[i] = MakeTH1('D', Form("FATIMA_VME/Stats/TDC_MultiplCh/TDCM_Ch_LaBr%2d",i), Form("TDC channel Multi Fatima %2d",i), 50, 0, 50);

  }
//      hFAT_QDC_vs_TDC_PMT_dT_Ch[7] = MakeTH2('D',Form("FATIMA_VME/Timing/Energy_vs._Time_SiPMdT_Ch/Energy_vs._Time_SiPMdT_Ch.%02d", 7),Form("Fatima Energy vs SiPMCh.0-SiPMCh.%02d", 7),4000,0,4000, 250,-2E4,2E4);
//      hFAT_QDC_vs_TDC_SC41dT_Ch[7] = MakeTH2('D',Form("FATIMA_VME/Timing/Energy_vs._Time_SC41_Ch/Energy_vs._Time_SC41_Ch.%02d", 7),Form("Fatima Energy vs SC41-SiPMCh.%02d", 7),4000,0,4000, 4E2,0,2E6);

    hFAT_QDCCalib1Sum = MakeTH1('D', "FATIMA_VME/Energy/Fat_VME__EnergySum", "LaBr Energy (all detectors)",40000,0,40000);
    hFAT_hits_QDC       = MakeTH1('D', "FATIMA_VME/Stats/QDC_FAThits", "bPlastic hit pattern QDC1",50,0,50);
   // hFAT_E_Mat_Sum = MakeTH2('D', "FATIMA_VME/Energy/Gam-GamSum", "FATIMA Gamma-Gamma (all detectors)",4000,0,4000, 4000,0,4000);
    hFAT_hits_TDC       = MakeTH1('D', "FATIMA_VME/Stats/TDC_FAThits", "FATIMA TDC statistics",50,0,50);
    hFAT_TDC_Multipl_PerChan       = MakeTH1('D', "FATIMA_VME/Stats/TDC_FAT_Multiplicity_perCh", "FATIMA TDC Multiplicity (hits per channel)",50,0,50);
    hFAT_TDC_Multipl       = MakeTH1('D', "FATIMA_VME/Stats/TDC_FAT_Multiplicity", "FATIMA TDC Multiplicity",50,0,50);

    hFAT_TDCdt_refSC41_Sum       = MakeTH1('D', "FATIMA_VME/Timing/TDCdt_refSC41_Sum", "TDC dT Ref SC41(all detectors)",4E4,0,2E6);
    hFAT_TDCdt_refSC41_Sum_gated       = MakeTH1('D', "FATIMA_VME/Timing/TDCdt_refSC41_Sum_EGated", "TDC dT (all detectors) Energy gated", 4E4,0,2E6);
    hFAT_TDCdt_refCha_Sum       = MakeTH1('D', "FATIMA_VME/Timing/TDCdt_ref0_AllM_Sum", "TDC dT LaBr0 - LaBr Multip All (all detectors)", 250,-2E4,2E4);
//     hFAT_TDCdt_refCha_Sum_M1       = MakeTH1('D', "FATIMA_VME/Timing/TDCdt_ref0_M1_Sum", "TDC dT LaBr0 - LaBr Multip 1 (all detectors)", 250,-2E4,2E4);
//     hFAT_TDCdt_refCha_Sum_M2       = MakeTH1('D', "FATIMA_VME/Timing/TDCdt_ref0_M2_Sum", "TDC dT LaBr0 - LaBr Multip 2 (all detectors)", 250,-2E4,2E4);
//     hFAT_TDCdt_refCha_Sum_M3       = MakeTH1('D', "FATIMA_VME/Timing/TDCdt_ref0_M2+_Sum", "TDC dT LaBr0 - LaBr Multip 2+ (all detectors)", 250,-2E4,2E4);

    hFAT_TDCdt_refCha_Sum_gated     = MakeTH1('D', "FATIMA_VME/Timing/TDCdt_ref0_Sum_EGated","TDC dT LaBr0 Gamma gated (all detectors)",250,-2E4,2E4);
   // hFAT_QDC_vs_TDC_PMT_dT = MakeTH2('D',"FATIMA_VME/Energy_vs._Time_PMdT","Energy_vs.Time_PMT_dT",4000,0,4000, 250,-2E4,2E4);
   // hFAT_QDC_vs_TDC_SC41dT = MakeTH2('D',"FATIMA_VME/Energy_vs._Time_S41dT","Energy_vs._Time_S41dT",4000,0,4000,4E2,0,2E6);

    hFAT_SC41_check      = MakeTH1('D', "FATIMA_VME/Timing/dt_SC41_L_R", "SC41 diff L vs R", 4E4,-2E6,2E6);

}
///-----------------------------------------------------------------------------------------------------------------------------------------------------------------------///
void EventAnlProc::Do_Fatima_Histos(EventAnlStore* pOutput){
    double Fat_QDC_i[50];
    double Fat_QDC_GainMatch[50], Fat_QDCGainMatch_j[50];
    double FATgate1_low, FATgate1_high;
    double Fat_TDC_T_Main[50], Fat_SC41_dT_Raw[50], Fat_SC41_dT_Calib[50],  Fat_Ch_dT[50], Fat_Ch_dT_Calib[50];
    int Fat_QDC_IDMain_i, Fat_QDC_IDMain_j, Fat_TDC_IDMain;
    double Fat_QDC_dt, Fat_QDCtime1, Fat_QDCtime2;
    int Fat_TDC_Incr;
    int  Fat_TDC_Multipl_perCh[50] ={0};


    Fat_QDC_IDMain_i = -1;
    Fat_QDC_IDMain_j = -1;
    Fat_TDC_IDMain = -1;
    Fat_QDC_dt = 0;
    Fat_QDCtime1 = 0;
    Fat_QDCtime2 = 0;


    for(int i=0; i<50; i++){
        Fat_QDC_i[i] = -1;
        //Fat_QDC_j[i] = -1;
        Fat_QDC_GainMatch[i] = 0;
        Fat_QDCGainMatch_j[i] = 0;
        Fat_Ch_dT_Calib[i] = 0;


        Fat_TDC_T_Main[i] = 0;
        Fat_SC41_dT_Raw[i] = 0;
        Fat_SC41_dT_Calib[i] = 0;
        Fat_Ch_dT[i] = 0;
  }
      //Fatima Energy gates
    FATgate1_low  = fCorrel->GFat_Egate_low;
    FATgate1_high = fCorrel->GFat_Egate_high;
   // Fat_E_gate1 = FATgate1_low + (FATgate1_high - FATgate1_low)/2.;

    /**------------------------------FATIMA Energy -----------------------------------------**/
       
      
        pOutput->pFAT_WR = Fat_WR;

    for (int i=0; i<FatQDCFired; i++){
  
        Fat_QDC_IDMain_i = FatQDCID[i]; //Channel ID
        pOutput->pFat_QDCFired = FatQDCID[i];
       // if(Fat_QDC_IDMain_i<40){
        pOutput->pFat_QDCID[i] = FatQDCID[i];
        hFAT_hits_QDC->Fill(Fat_QDC_IDMain_i);

        Fat_QDC_i[Fat_QDC_IDMain_i] = FatQDC[i];  //Calibrated energy
        Fat_QDCtime1 = FatQDC_T[i];
if(Fat_QDC_IDMain_i<40){
          ///FATIMA Calibrated Energy Singles
        //Fat_QDC_GainMatch[Fat_QDC_IDMain_i] = fCal->Afat[Fat_QDC_IDMain_i]* pow(Fat_QDC_i[Fat_QDC_IDMain_i],3) + fCal->Bfat[Fat_QDC_IDMain_i]* pow(FatQDC[i],2) + fCal->Cfat[Fat_QDC_IDMain_i]*FatQDC[i] + fCal->Dfat[Fat_QDC_IDMain_i];
        Fat_QDC_GainMatch[Fat_QDC_IDMain_i] = Fat_QDC_i[Fat_QDC_IDMain_i];
        hFAT_QDCCalib1[Fat_QDC_IDMain_i]->Fill(Fat_QDC_GainMatch[Fat_QDC_IDMain_i]);
        pOutput->pFat_QDCGainMatch[Fat_QDC_IDMain_i] = Fat_QDC_GainMatch[Fat_QDC_IDMain_i];
        hFAT_QDCCalib1Sum->Fill(Fat_QDC_GainMatch[Fat_QDC_IDMain_i]);
}
               ///Gamma-Gamma Fatima
          for (int j=0; j<FatQDCFired; j++){
            Fat_QDC_IDMain_j= FatQDCID[j];
            Fat_QDCGainMatch_j[Fat_QDC_IDMain_j] = FatQDC[j];
            Fat_QDCtime2 = FatQDC_T[j];

            ///Dont loop on the first hit again:
           if(Fat_QDC_IDMain_i < Fat_QDC_IDMain_j){
                Fat_QDC_dt = Fat_QDCtime1 - Fat_QDCtime2;
               // hFAT_QDCdt[Fat_QDC_IDMain_i] ->Fill(Fat_QDC_dt);
                //Fill Energy-Energy matrix (NOTE: turned off making 2D matrices for each channel for now to speed things up)
               // hFAT_Chan_E_Mat[Fat_QDC_IDMain_i]->Fill(Fat_QDC_GainMatch[Fat_QDC_IDMain_i], Fat_QDCGainMatch_j[Fat_QDC_IDMain_j]);

              //  hFAT_E_Mat_Sum->Fill(Fat_QDC_GainMatch[Fat_QDC_IDMain_i], Fat_QDCGainMatch_j[Fat_QDC_IDMain_j]);
                }
            }
       // }
  }

 /**---------------------------------FATIMA TIMING -----------------------------------------**/
          hFAT_TDC_Multipl->Fill(FatTDCFired);
          pOutput -> pFat_TDCFired = FatTDCFired;
         for (int i=0; i<FatTDCFired; i++){

           ///FAT TDC ID and Raw TDC data
          Fat_TDC_IDMain = FatTDCID[i];
          pOutput -> pFat_TDCID[i] = FatTDCID[i];
          Fat_TDC_T_Main[Fat_TDC_IDMain] = FatTDC_TS[i][Fat_TDC_IDMain]; //In ps
          pOutput ->  pFat_TDC_T[Fat_TDC_IDMain] =  Fat_TDC_T_Main[Fat_TDC_IDMain];

          ///FAT Multiplicity
           Fat_TDC_Multipl_perCh[Fat_TDC_IDMain]++;
           pOutput -> pFat_TDC_Multipl_perCh[Fat_TDC_IDMain] =  Fat_TDC_Multipl_perCh[Fat_TDC_IDMain];
           hFAT_TDC_Multipl_ch[i] -> Fill(Fat_TDC_Multipl_perCh[Fat_TDC_IDMain]);
           hFAT_TDC_Multipl_PerChan -> Fill(Fat_TDC_Multipl_perCh[Fat_TDC_IDMain]);

           hFAT_SC41_check->Fill(Fat_TDC_T_Main[35]-Fat_TDC_T_Main[34]);
        // cout <<"event " << event_number << " Fat_TDC_T_Main[41] " << Fat_TDC_T_Main[40] << " Fat_TDC_T_Main[41] " << Fat_TDC_T_Main[41]<< endl;
            ///Hit pattern
           hFAT_hits_TDC->Fill(Fat_TDC_IDMain);
            ///SC41 (TDC Ch.40) - FAT SiPM Ch.x
          // if( Fat_TDC_T_Main[Fat_TDC_IDMain]>fCorrel->GFat_TRawgate_low && Fat_TDC_T_Main[Fat_TDC_IDMain]<fCorrel->GFat_TRawgate_high){

            Fat_SC41_dT_Raw[Fat_TDC_IDMain] = (SC41 -  Fat_TDC_T_Main[Fat_TDC_IDMain]); //ps
            Fat_SC41_dT_Calib[Fat_TDC_IDMain]  = Fat_SC41_dT_Raw[Fat_TDC_IDMain] + fCal-> TFatTDC_SC41dT[Fat_TDC_IDMain];
           // pOutput ->  pFat_SC41_dT_Calib[Fat_TDC_IDMain] =  Fat_SC41_dT_Calib[Fat_TDC_IDMain];
        //   cout<<"1  Fat_SC41_dT_Raw[Fat_TDC_IDMain] "<< Fat_SC41_dT_Raw[Fat_TDC_IDMain]  <<endl;

            if( FatQDCID[i] == Fat_TDC_IDMain){
            double Mindiff = Fat_SC41_dT_Calib[0];
               if(Mindiff>Fat_SC41_dT_Calib[i]) Mindiff =Fat_SC41_dT_Calib[i];

        // cout<<"event " << event_number <<" Fat_SC41_dT_Calib[Fat_TDC_IDMain]  " << Fat_SC41_dT_Calib[Fat_TDC_IDMain]  << " SC41 " <<SC41 <<  "  Fat_TDC_T_Main[Fat_TDC_IDMain] " <<  Fat_TDC_T_Main[Fat_TDC_IDMain] <<" Fat_TDC_IDMain " << Fat_TDC_IDMain << " FatQDCID[i] " << FatQDCID[i] << " Mindiff " << Mindiff<<endl;

            }
         ///SC41 - Fatima TDC Only take the first hit per channel (Sultan)
           Fat_TDC_Incr = 0;
            for(int j=0; j<=i; j++){
             if (Fat_TDC_IDMain != FatTDCID[j] ) Fat_TDC_Incr++;



            }
            //if(Fat_TDC_Incr == i  && Fat_TDC_IDMain < 40 ){
             for (int j = 0; j< FatQDCFired; j++){
            if( FatQDCID[j] == Fat_TDC_IDMain){
         //   cout<<"Fat_SC41_dT_Calib[Fat_TDC_IDMain] " <<Fat_SC41_dT_Calib[Fat_TDC_IDMain] << endl;
              hFAT_TDCdt_refSC41[Fat_TDC_IDMain] -> Fill(Fat_SC41_dT_Calib[Fat_TDC_IDMain]);
              hFAT_TDCdt_refSC41_Sum -> Fill(Fat_SC41_dT_Calib[Fat_TDC_IDMain]);
              //
              //pOutput->pFat_SC41_dT_Calib[Fat_TDC_IDMain] = Fat_SC41_dT_Calib[Fat_TDC_IDMain];
            }

              /// Fatima Time SiPM 0 - SiPM Ch.x (Ch. 0 used as the reference)
              if(Fat_TDC_IDMain < 40 && Fat_CHA_0_TDC>0  && Fat_TDC_T_Main[Fat_TDC_IDMain] > 0&& FatQDCID[j] == Fat_TDC_IDMain){
                    Fat_Ch_dT[Fat_TDC_IDMain] =  (Fat_CHA_0_TDC - Fat_TDC_T_Main[Fat_TDC_IDMain]);
                  //  Fat_Ch_dT_Calib[Fat_TDC_IDMain] =  Fat_Ch_dT[Fat_TDC_IDMain] + fCal-> TFatTDC_Chref_dT[Fat_TDC_IDMain];
                    //pOutput->pFat_Ch_dT[Fat_TDC_IDMain] =  Fat_Ch_dT_Calib[Fat_TDC_IDMain] ;

                    if(Fat_TDC_IDMain!=0){
                        hFAT_TDCdt_refCha[Fat_TDC_IDMain]->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);
                        hFAT_TDCdt_refCha_Sum  ->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);

                    /// Histogram TDC cha 0 - TDC Ch.x
                    if(Fat_TDC_IDMain!=0 && Fat_TDC_Multipl_perCh[Fat_TDC_IDMain]==1 ){
                    //hFAT_TDCdt_refCha_M1[Fat_TDC_IDMain]->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);
                   // hFAT_TDCdt_refCha_Sum_M1  ->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);
                        }
                    if(Fat_TDC_IDMain!=0 && Fat_TDC_Multipl_perCh[Fat_TDC_IDMain]==2 ){
                   // hFAT_TDCdt_refCha_M2[Fat_TDC_IDMain]->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);
                   // hFAT_TDCdt_refCha_Sum_M2  ->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);
                        }
                    if(Fat_TDC_IDMain!=0 && Fat_TDC_Multipl_perCh[Fat_TDC_IDMain]>2 ){
                   // hFAT_TDCdt_refCha_M3[Fat_TDC_IDMain]->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);
                   // hFAT_TDCdt_refCha_Sum_M3  ->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);
                            }
                        }
                   // }
              }
                        ///Energy Time matrix (just one channel (ch. 7) for now to speed things up)

                       if(FatQDCID[j] == Fat_TDC_IDMain && Fat_QDC_GainMatch[FatQDCID[i]]>0  ){
                  //  hFAT_QDC_vs_TDC_PMT_dT_Ch[7] ->Fill( Fat_QDC_GainMatch[7],Fat_Ch_dT[Fat_TDC_IDMain]);
                  //  hFAT_QDC_vs_TDC_SC41dT_Ch[7] ->Fill( Fat_QDC_GainMatch[7],Fat_SC41_dT_Calib[Fat_TDC_IDMain]);

                  //  hFAT_QDC_vs_TDC_PMT_dT ->Fill( Fat_QDC_GainMatch[FatQDCID[j]],Fat_Ch_dT[Fat_TDC_IDMain]);
                  //  hFAT_QDC_vs_TDC_SC41dT ->Fill( Fat_QDC_GainMatch[FatQDCID[j]],Fat_SC41_dT_Calib[Fat_TDC_IDMain]);
                            }


            ///Gamma energy gates
             if(Fat_QDC_GainMatch[FatQDCID[j]] > FATgate1_low && Fat_QDC_GainMatch[FatQDCID[j]] < FATgate1_high){
                 /// Fatima Time SiPM 0 - SiPM Ch.x Energy gated
                if(Fat_TDC_IDMain!=0 && Fat_Ch_dT[Fat_TDC_IDMain]!=0) {
                   hFAT_TDCdt_refCha_gated[Fat_TDC_IDMain] ->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);
                   hFAT_TDCdt_refCha_Sum_gated ->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);
            ///SC41 - Fatima TDC Energy Gated
            hFAT_TDCdt_refSC41_gated[Fat_TDC_IDMain] -> Fill(Fat_SC41_dT_Calib[Fat_TDC_IDMain]);
         //  if (Fat_TDC_Multipl_perCh[Fat_TDC_IDMain]==1) hFAT_TDCdt_refSC41_M1_gated[i] -> Fill(Fat_SC41_dT_Calib[Fat_TDC_IDMain]);
         //    if (Fat_TDC_Multipl_perCh[Fat_TDC_IDMain]==2) hFAT_TDCdt_refSC41_M2_gated[i] -> Fill(Fat_SC41_dT_Calib[Fat_TDC_IDMain]);
          //   if (Fat_TDC_Multipl_perCh[Fat_TDC_IDMain]==3) hFAT_TDCdt_refSC41_M3_gated[i] -> Fill(Fat_SC41_dT_Calib[Fat_TDC_IDMain]);
            hFAT_TDCdt_refSC41_Sum_gated -> Fill(Fat_SC41_dT_Calib[Fat_TDC_IDMain]);

                   }
                }
             }
          }
        }
    


/**----------------------------------------------------------------------------------------------**/
/**--------------------------------------  FATIMA TAMEX ----------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
// void EventAnlProc::Make_Fatima_Tamex_Histos(){
//
//     for (int i =0; i<50; i++)
//   {
//     hFAT_ToT[i] = MakeTH1('D', Form("FATIMA/ToT/ToT_Ch.%02d", i), Form("TOT Detector %2d", i), 16250, 0., 65000.);
//     hFAT_lead_lead[i] = MakeTH1('D', Form("FATIMA/Lead-Lead/Lead-LeadCh.%02d", i), Form("Lead - Lead Ch. %2d", i),2500, -50000., 50000);
//     hFAT_lead_lead_ref[i] = MakeTH1('D', Form("FATIMA/Lead-Lead_Ref/Lead-LeadRefCh.%02d", i), Form("Lead Ref Ch.0 - Lead Ch. %2d", i),2500, -50000., 50000.);
//
//     hFAT_lead_lead_gated[i] = MakeTH1('D', Form("FATIMA/Lead-Lead_Egated/Lead-Lead_Egated_Ch.%02d", i), Form("Lead - Lead Energy gated Ch.  %2d", i),2500, -50000., 50000.);
//     hFAT_lead_lead_gated1[i] = MakeTH1('D', Form("FATIMA/Lead-Lead_Egated1/Lead-Lead_Egated1_Ch.%02d", i), Form("Lead - Lead Energy gated Ch.  %2d", i),2500, -50000., 50000.);
//     hFAT_lead_lead_gated2[i] = MakeTH1('D', Form("FATIMA/Lead-Lead_Egated2/Lead-Lead_Egated2_Ch.%02d", i), Form("Lead - Lead Energy gated Ch.  %2d", i),2500, -50000., 50000.);
//
//     hFAT_lead_lead_energy[i] = MakeTH2('D',Form("FATIMA/Lead-Lead_energy/Energy_vs._dT/Energy_vs._dT_Ch.%02d", i),Form("Fatima Energy vs Lead-Lead.%02d", i),500,0,2000, 160, -40000., 40000);
//
//     hFAT_trail_trail[i] = MakeTH1('D', Form("FATIMA/Trail-Trail/trail-trailCh.%02d", i), Form("Trail - Trail Ch. %2d", i), 2500, -50000., 50000.);
//
//    // hFAT_Sc41lead_leadmaxtot[i] = MakeTH1('D', Form("FATIMA/Sc41-LeadMaxToT/SC41Lead_LeadCh.%02d", i), Form("SC41 Lead - (max ToT chan) Lead %2d ", i), 4002, -100000., 100000.);
//
//     }
//
//     hFAT_gamma_gamma = MakeTH2('D', "FATIMA/Gamma-Gamma/Sum", "FATIMA Gamma-Gamma (all detectors)",6500,0,65000, 6500,0,65000);
// }
//-----------------------------------------------------------------------------------------------//
// void EventAnlProc::Do_Fatima_Tamex_Histos(EventUnpackStore* pInput, EventAnlStore* pOutput){
//     double lead_lead_fat[50], trail_trail_fat[50];
//     double ToT_fat[50];
//     double lead_fat1[50], lead_fat2[50] ;
//     double trail_fat1[50], trail_fat2[50];
//     double  lead_lead_fat_Ref0[50];
//     double lead_lead_fat_gated[50], lead_lead_fat_gated1[50], lead_lead_fat_gated2[50];
//     int hits_fat_lead = 0, hits_fat_trail=0;
//     for (int i = 0; i < 50; i++)
//         {
//              lead_lead_fat[i]=0;
//              lead_lead_fat_Ref0[i]=0;
//              trail_trail_fat[i]=0;
//              ToT_fat[i] = 0;
//              trail_fat1[i] = 0;
//              trail_fat2[i] = 0;
//              lead_fat1[i] =0;
//              lead_fat2[i] =0;
//
//              lead_lead_fat_gated[i] = 0;
//              lead_lead_fat_gated1[i] = 0;
//              lead_lead_fat_gated2[i] = 0;
//
//         }
//
//         ///Loop on channels
//       for (int i = 0; i < 50; i++)
//         {
//             ///Lead T
//             for(int j=0; j< pInput->fFat_PMT_Lead_N[i]; j++){
//                 lead_fat1[i] = pInput->fFat_Lead_PMT[i][j];
//                 hits_fat_lead++;
//                 pOutput->pFat_LeadT[i][j] = lead_fat1[i];
//                 pOutput->pFat_LeadHits = hits_fat_lead;
//
//                 ///Ref channel 0
//                 if(pInput->fFat_Lead_PMT[0][j]>0 && pInput->fFat_Lead_PMT[i][j]>0&& hits_fat_lead>1){
//                 lead_lead_fat_Ref0[i] = (pInput->fFat_Lead_PMT[0][j] -  pInput->fFat_Lead_PMT[i][j])*5000;
//                // cout<<"event 1 " << event_number << " lead N " << pInput->fFat_PMT_Lead_N[i] << " fFat_Lead_PMT[0][j] " << pInput->fFat_Lead_PMT[0][j]*5000 << " Fat_Lead_PMT[i][j] " << pInput->fFat_Lead_PMT[i][j]*5000 << " lead_lead_fat_Ref0[i] " <<lead_lead_fat_Ref0[i]<< " i " << i << " j " << j <<  endl;
//                 hFAT_lead_lead_ref[i] ->Fill(lead_lead_fat_Ref0[i]);
//                 }
//
//                 ///all Dets compare
//                 for (int k=0; k< 50; k++){
//                     if(i != k && hits_fat_lead>1){
//                          lead_fat2[k] = pInput->fFat_Lead_PMT[k][j];
//
//                     lead_lead_fat[i] = (lead_fat1[i] - lead_fat2[k]);
//                         if (ABS(lead_lead_fat[i]) > (double)(COARSE_CT_RANGE>>1))        // overflow
//                             {
//                             lead_lead_fat[i] = CYCLE_TIME*(lead_lead_fat[i]  + COARSE_CT_RANGE);
//                         }
//                          else {
//                              lead_lead_fat[i]  = CYCLE_TIME*lead_lead_fat[i];
//                          }
//                          if(lead_lead_fat[i]!=0 && lead_fat1[i]>1 && lead_fat2[k]>1 ) {
//                              hFAT_lead_lead[i] -> Fill(lead_lead_fat[i]);
//          //cout<<"event " << event_number << " lead_lead_fat[i] " <<lead_lead_fat[i]<<" lead_fat1[i] " << lead_fat1[i] << " lead_fat2[k] "<<lead_fat2[k]  << " i " << i <<endl;
//                       }
//                     }
//                 }
//             }
//             ///Trail T
//             for(int j=0; j< pInput->fFat_PMT_Trail_N[i]; j++){
//                 trail_fat1[i] = pInput->fFat_Trail_PMT[i][j];
//                 hits_fat_trail++;
//                 pOutput->pFat_TrailT[i][j] = trail_fat1[i];
//                 pOutput->pFat_TrailHits = hits_fat_trail;
//                 for (int k=0; k< 50; k++){
//                     if(i != k && hits_fat_trail>1){
//                          trail_fat2[i] = pInput->fFat_Trail_PMT[k][j];
//                          if(trail_fat2>0 && j>0){
//
//             trail_trail_fat[i] = (trail_fat1[i] - trail_fat2[i]);
//             if (ABS(trail_trail_fat[i]) > (double)(COARSE_CT_RANGE>>1))        // overflow
//                             {
//                             trail_trail_fat[i] = CYCLE_TIME*(trail_trail_fat[i] + COARSE_CT_RANGE);
//                         }
//              else {
//                             trail_trail_fat[i]  = CYCLE_TIME*trail_trail_fat[i];
//                         }
//
//             hFAT_trail_trail[i] -> Fill(trail_trail_fat[i]);
//                         }
//                     }
//                 }
//             }
//             ///ToT (~Energy)
//           for(int j=0; j< pInput->fFat_PMT_Lead_N[i]; j++){
//            //   if(pInput->fFat_Trail_PMT[i][j] >0 && pInput->fFat_Lead_PMT[i][j]>0){
//                 ToT_fat[i] = (pInput->fFat_Trail_PMT[i][j] - pInput->fFat_Lead_PMT[i][j]);
//
//                 ///Correction for overflows
//             if(ABS(ToT_fat[i]) >(double)(COARSE_CT_RANGE>>1)) {
//                    ToT_fat[i] = CYCLE_TIME*(ToT_fat[i] + COARSE_CT_RANGE);
//             }
//              else{
//                         ToT_fat[i]= CYCLE_TIME*ToT_fat[i];
//                    }
// //
//  //cout<<"pInput->fFat_Trail_PMT[i][j]  " <<pInput->fFat_Trail_PMT[i][j] << " pInput->fFat_Lead_PMT[i][j] " << pInput->fFat_Lead_PMT[i][j] <<" ToT_fat[i] " << ToT_fat[i] << endl;
//                     hFAT_ToT[i] ->Fill(ToT_fat[i]); //raw Energy
//
//
//               ///Gain matching
//               pOutput-> pFat_ToTCalib[i] = fCal->Afat_TAMEX[i]* pow(ToT_fat[i],3) + fCal->Bfat_TAMEX[i]* pow(ToT_fat[i],2) + fCal->Cfat_TAMEX[i]*ToT_fat[i] + fCal->Dfat_TAMEX[i];
//
//               ///Lead-Lead Energy gating
//            for (int k=0; k< 50; k++){
//                if(i != k ){
//                 if(ToT_fat[i]>0 && ToT_fat[k]>0 ) hFAT_gamma_gamma->Fill(ToT_fat[i],ToT_fat[k]);
//
//                  if(pOutput->pFat_LeadT[k][j]>0 && pOutput->pFat_LeadHits>1 ){
//                     lead_lead_fat_gated[i] = (pOutput->pFat_LeadT[i][j]  -  pOutput->pFat_LeadT[k][j])*5000; ///into ps
//
//            if( ToT_fat[1]>32700&&ToT_fat[1]<34000){
//             hFAT_lead_lead_energy[i]->Fill(ToT_fat[0]/27,lead_lead_fat_gated[i]);
//
//             }
//                          ///Temp E gate
//                          if(ToT_fat[1]<31400){
//
//              hFAT_lead_lead_gated[i] -> Fill(lead_lead_fat_gated[i]);
//                  }
//                        if(ToT_fat[1]>27000 && ToT_fat[1]<27716){
//
//             hFAT_lead_lead_gated1[i] -> Fill(lead_lead_fat_gated[i]);
//                  }
//                   if(ToT_fat[1]>27806 && ToT_fat[1]<28537){
//             hFAT_lead_lead_gated2[i] -> Fill(lead_lead_fat_gated[i]);
//                  }
//                }
//              }
//           }
//        // }
//      }
//    }
// }


/**----------------------------------------------------------------------------------------------**/
/**--------------------------------------  GALILEO  ---------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
void EventAnlProc::Make_Galileo_Histos()
{
  hGAL_ESum = MakeTH1('I',"GALILEO/Sum/GALILEO_ESum","GALILEO Energy Sum",20000,0,20000);
  //hGAL_ESum_largerange_all = MakeTH1('I',"GALILEO/Sum/GALILEO_ESum_largerange","GALILEO Energy Sum",20000,0,20000);
  hGAL_ESum_largerange_OF = MakeTH1('I',"GALILEO/Sum/hGAL_ESum_largerange_OF","GALILEO Energy Sum (Overflow)",20000,0,20000);
  hGAL_ESum_largerange_PU = MakeTH1('I',"GALILEO/Sum/hGAL_ESum_largerange_PU","GALILEO Energy Sum (Pileup)",20000,0,20000);
  hGAL_Hit_Pat = MakeTH1('I',"GALILEO/Stats/GALILEO_Hit_Pat","GALILEO Hit Pattern",36,0,36);
  //hGAL_Multi_1 = MakeTH1('I',"GALILEO/Stats/GALILEO_Multiplicity_1","GALILEO Multiplicity 1",50,0,50);
  //hGAL_Multi_2 = MakeTH1('I',"GALILEO/Stats/GALILEO_Multiplicity_2","GALILEO Multiplicity 2",50,0,50);
  //hGAL_Multi_3 = MakeTH1('I',"GALILEO/Stats/GALILEO_Multiplicity_3","GALILEO Multiplicity 3",50,0,50);
  //hGAL_Multi_4 = MakeTH1('I',"GALILEO/Stats/GALILEO_Multiplicity_4","GALILEO Multiplicity 4",50,0,50);
    hGAL_Chan_E_Mat = MakeTH2('D',"GALILEO/GALILEO_E_Mat","GALILEO Energy-Energy Matrix",2500,0,10000,2500,0,10000);
  //hGAL_Chan_E_M1= MakeTH1('I',"GALILEO/Stats/GALILEO_multiplicity_1","GALILEO Channel Energy",5000,0,5000);
  //hGAL_Chan_E_M2= MakeTH1('I',"GALILEO/Stats/GALILEO_multiplicity_2","GALILEO Channel Energy",5000,0,5000);
  hGAL_AddbackSum = MakeTH1('I',"GALILEO/Sum/GALILEO_Addback","GALILEO Addback Energy Sum",20000,0,20000);

  for (int i=0; i<GALILEO_MAX_DETS; i++)
  {
    for (int j = 0; j < GALILEO_CRYSTALS; j++)
    {
      hGAL_Chan_E[i][j] = MakeTH1('D',Form("GALILEO/Energy_Ch./GALILEO_E_Det_%2d_%1d",i, j), Form("GALILEO Channel Energy Detector %2d Crystal %1d",i, j),5000,0,5000);
    }
  //  hGAL_FatdT[j] = MakeTH1('I',Form("Correlations/Fatima_Galilieo/Fat_GAldT%2d",j),Form("GALILEO Fatima dT Ch. %2d",j),2000,-1000,1000);
   // hGAL_Chan_E2[j] = MakeTH1('D',Form("GALILEO/GALILEO_Energy2/GALILEO_E2%2d",j), Form("GALILEO Channel Energy Channel %2d",j),5000,0,5000);
    //hGAL_Chan_Egate[j] = MakeTH1('D',Form("GALILEO/gated energy/GALILEO_Egate%2d",j), Form("GALILEO Channel Energy Channel %2d",j),5000,0,5000);
    }
  //for (int k=0; k<32; k++){

    //hGAL_Chan_Time_Diff[k] = MakeTH1('D',Form("GALILEO/Time_diff/GALILEO_Chan_Time_Diff%2d",k), Form("GALILEO Channel Time Difference for %2d",k),2000,-1000,1000);
    //hGAL_Chan_Timedifference_new[k] = MakeTH1('D',Form("GALILEO/Timediff_new/GALILEO_Chan_T_Diff%2d",k), Form("GALILEO Channel T Difference for %2d",k),2000,-1000,1000);
   // hGAL_Time_Diff_vs_Energy[k] = MakeTH2('D',Form("GALILEO/GALILEO_dT_vs_Energy_Spectra/GALILEO_dT_vs_E%2d",k), Form("GALILEO Time Difference Vs Channel Energy Channel %2d",k),5000,0,5000,100,-1000,1000);
  //}
    }
///-----------------------------------------------------------------------------------------------------------------------------------------------------------------------///
void EventAnlProc::Do_Galileo_Histos(EventAnlStore* pOutput)
{
    // Process hits once
    for (int i = 0; i < GalFired; i++)
    {
       // Skip pileup/overflow events
       if (GalPileUp[i])
       {
          hGAL_ESum_largerange_PU->Fill(GalE_Cal[i]);
          continue;
        }

       if (GalOverFlow[i])
       {
         hGAL_ESum_largerange_OF->Fill(GalE_Cal[i]);
         continue;
        }

       int det = GalDet[i];
       int crys = GalCrys[i];
       pOutput->pGal_T[det][crys] = GalT[i];
       pOutput->pGal_E[det][crys] = GalE_Cal[i];
  
       hGAL_Hit_Pat->Fill(det * GALILEO_CRYSTALS + crys);
       hGAL_ESum->Fill(GalE_Cal[i]);
       hGAL_Chan_E[det][crys]->Fill(GalE_Cal[i]);

       // 2D Matrix generation
       for (int j = 0; j < GalFired; j++)
       {
          if (i == j) continue;
          hGAL_Chan_E_Mat->Fill(GalE_Cal[i], GalE_Cal[j]);
       }
    }
            
    static const long dT_addback = 50;
            
    // Detector addback
    for (int i = 0; i < GALILEO_MAX_DETS; i++)
    {
        double E[GALILEO_CRYSTALS] = { 0 };
        long T[GALILEO_CRYSTALS] = { 0 };
        int n[GALILEO_CRYSTALS] = { 0 };
        int v = 0;
        for (int j = 0; j < GALILEO_CRYSTALS; j++)
            {
          if (pOutput->pGal_E[i][j] == 0) continue;
          bool added = false;
          // Try to addback to an existing hit
          for (int k = 0; k < v; k++)
          {
            if (T[k] - pOutput->pGal_T[i][j] < dT_addback)
            {
              E[k] += pOutput->pGal_E[i][j];
              T[k] = (T[k] + pOutput->pGal_T[i][j]) / (n[k] + 1);
              n[k] += 1;
              added = true;
            }
            }
           
          // Add to a new hit
          if (!added)
          {
            T[v] = pOutput->pGal_T[i][j];
            E[v] = pOutput->pGal_E[i][j];
            n[v] = 1;
            v++;
          }
        }
           
        // Fill and write to Tree the addback energies
        for (int j = 0; j < v; j++)
        {
            pOutput->pGal_EAddback[i][j] = E[j];
            hGAL_AddbackSum->Fill(E[j]);
       }
    }
} //end of Do_Galileo_Histos()

/**----------------------------------------------------------------------------------------------**/
/**--------------------------------------  FINGER  ----------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
void EventAnlProc::Make_Finger_Histos(){
for (int i =0; i<52; i++)
  {

    hFING_ToT_Strip[i] = MakeTH1('D', Form("FINGER/ToT_Strip/ToT_Strip%02d", i), Form("TOT Strip %2d", i), 4001, -100000., 100000.);
    hFING_MaxToT_Strip[i] = MakeTH1('D', Form("FINGER/MaxToT_Strip/MaxToT_Strip%02d", i), Form("Max TOT Strip %2d", i), 4001,-100000., 100000.);
    hFING_ToT_PMT[i] = MakeTH1('D', Form("FINGER/ToT_PMT/ToT_PMT%02d", i), Form("TOT PMT %2d", i), 4001,-100000., 100000.);
    hFING_lead_lead[i] = MakeTH1('D', Form("FINGER/Lead-Lead/Lead-LeadCh.%02d", i), Form("Lead - Lead %2d", i), 4001, -100000., 100000.);
    hFING_trail_trail[i] = MakeTH1('D', Form("FINGER/trail-trail/trail-trailCh.%02d", i), Form("trail - trail %2d", i), 4001, -100000., 100000.);

    hFING_Sc41lead_leadmaxtot[i] = MakeTH1('D', Form("FINGER/Sc41-LeadMaxToT/SC41Lead_LeadCh.%02d", i), Form("SC41 Lead - (max ToT chan) Lead %2d ", i), 4001, -100000., 100000.);

}

hFING_ToT_v_Strip = MakeTH2('I',"FINGER/TOT_vs_Strip","ToT vs Strip number", 52, 0, 52, 2001, 0, 100000., "Strip", "ToT");
hFING_MaxToT_v_Strip = MakeTH2('I',"FINGER/MaxTOT_vs_Strip","MaxToT vs Strip number", 52, 0, 52, 2001, 0., 100000., "Strip", "Max ToT");
hFING_ToT_v_PMT = MakeTH2('I',"FINGER/TOT_vs_PMT","ToT vs PMT number", 52, 0, 52, 2001, 0., 100000., "PMT", "Max ToT");

hFing_ToTRatio_v_Strip = MakeTH2('D', "FINGER/ToT_Ratio_v_Strip", "TOT Ratio vs Strip", 52, 0, 52, 2001, 0, 1., "Strip", "ToT Ratio");
hFing_MaxToTRatio_v_Strip = MakeTH2('D', "FINGER/MaxToT_Ratio_v_Strip", "Max TOT Ratio vs Strip", 52, 0, 52, 2001, 0, 1., "Strip", "Max ToT Ratio");

hFING_SC41_lead_lead = MakeTH1('I',"FINGER/SC41_leadlead","SC41 lead-lead", 500, -500., 500.);
hFING_SC41_trail_trail = MakeTH1('I',"FINGER/SC41_trailtrail","SC41 trail-trail", 500, -500., 500.);

hFING_SC41_tot = MakeTH1('I',"FINGER/SC41_tot","SC41 trail-trail", 4002, 0., 100000.);
hFING_Multiplicity =  MakeTH1('I',"FINGER/Multiplicity","Finger Multiplicity", 500, 0., 500.);

/*
  hFING_Hit_Pat = MakeTH1('I',"FINGER/Stats/FINGER_Hit_Pat","FINGER Hit Pattern",52,0,52);
  hFING_ToT_StripID = MakeTH2('I',"FINGER/TOT_vs_PMT","ToT vs Strip number", 2001, -100000., 100000., 52, 0, 52);
  hFING_MaxToT_StripID = MakeTH2('I',"FINGER/MaxTOT_vs_PMT","MaxToT vs Strip number", 2001, -100000., 100000., 52, 0, 52);
  hFING_Pos = MakeTH2('D',"FINGER/position","Time ratio vs Strips",51,1,51, 1000, -10., 10.);
  hFING_Pos_ToT = MakeTH2('D',"FINGER/positionToT","ToT ratio vs Strips",51,1,51, 5000, -1., 1.);
  hFING_Pos_ToT_Max = MakeTH2('D',"FINGER/positionToTMax","ToT ratio vs Strips Max",51,1,51, 5000, -1., 1.);

  hFING_ToT_StripID_Exp= MakeTH2('I',"FINGER/TOT_vs_PMT_Exp","ToT exponential vs Strip number", 1000, 0., 10000., 52, 0, 52);
  hFING_MaxToTExp_StripID = MakeTH2('I',"FINGER/MaxTOTExp_vs_PMT","MaxToT exponential vs Strip number", 1000, 0., 100000., 52, 0, 52);

  hFING_ToT_StripID_UpDown = MakeTH2('I',"FINGER/TOT_vs_PMT_sumpmt","ToT vs Strip number sum PMT", 2001, -100000., 100000., 52, 0, 52);

  hFING_ToT_StripID_UpDown_max = MakeTH2('I',"FINGER/TOT_vs_PMT_sumpmt_MAX","ToT vs Strip number sum PMT MAX", 2001, -100000., 100000., 52, 0, 52);

  hFING_LeadLead_StripID = MakeTH2('I', "FINGER/LeadLead_vs_Strip", "Lead-Lead vs Strip ID", 1000, -50000, 50000, 52, 0, 52);
  */
}
///-----------------------------------------------------------------------------------------------------------------------------------------------------------------------///

    void EventAnlProc::Do_Finger_Histos(EventUnpackStore* pInput, EventAnlStore* pOutput){
  // Verify Finger Data Sanity
  // Need to investigate if lead/trail mismatches are fine or not
  int maxtot = 0;
  int maxtotup = 0;
  int maxtotdown = 0;
  int maxtotchan = 0;
  double lead_lead = 0;
  double trail_trail = 0;
  double SC41_lead_lead = 0;
  double SC41_trail_trail = 0;
  double SC41_tot=0;
  double maxtotSc41lead=0;

 // for (int i=0; i<2; i++){
    for(int l=0; l<50; l++){
  SC41_lead_lead = (pInput->fFing_SC41_lead[0][l] - pInput->fFing_SC41_lead[1][l]);
  SC41_trail_trail = (pInput->fFing_SC41_trail[0][l] - pInput->fFing_SC41_trail[1][l]);
  SC41_tot =   (pInput->fFing_SC41_trail[0][l]-pInput->fFing_SC41_lead[0][l]);

  if(SC41_lead_lead!=0) hFING_SC41_lead_lead->Fill(SC41_lead_lead);
  if(SC41_trail_trail!=0)  hFING_SC41_trail_trail->Fill(SC41_trail_trail);
  if(SC41_tot!=0) hFING_SC41_tot ->Fill(SC41_tot);
    }

//  cout << "event " << event_number <<" lead chan 0 " << pInput->fFing_SC41_lead[0][l] <<" lead chan 1 "<<  pInput->fFing_SC41_lead[1][l] << " SC41_lead_lead " <<SC41_lead_lead <<"  SC41_trail_trail " << SC41_trail_trail << endl;
   // }
  //}



  // PMTs
  for (int i = 0; i < 52; i++)
  {
    if (pInput->fFing_PMT_Lead_N[i] != pInput->fFing_PMT_Trail_N[i])
    {
      continue;
    }


    for (int j = 0; j < pInput->fFing_PMT_Lead_N[i]; j++)
    {
        hFING_Multiplicity->Fill(pInput->fFing_PMT_Lead_N[i]);
      int tot =-( pInput->fFing_Trail_PMT[i][j] - pInput->fFing_Lead_PMT[i][j]);
      // cout<<"pInput->fFing_PMT_Lead_N[i] " << pInput->fFing_PMT_Lead_N[i] <<" pInput->fFing_Trail_PMT[i][j] " << pInput->fFing_Trail_PMT[i][j] <<  " i " << i <<endl;
//cout<<"event " << event_number <<" tot "<< tot<<" pInput->fFing_Trail_PMT[i][j] " << pInput->fFing_Trail_PMT[i][j] << " pInput->fFing_Lead_PMT[i][j] " << pInput->fFing_Lead_PMT[i][j] <<" i " << i  << " j " << j <<" fFing_Trail_Up[i][j] " <<pInput->fFing_Trail_Up[i][j] <<" fFing_Trail_Down[i][j] "<<  pInput->fFing_Trail_Down[i][j] << endl;
      // Fill tot
      hFING_ToT_PMT[i]->Fill(tot);
      hFING_ToT_v_PMT->Fill(i, tot);
    }
  }

  for (int i = 0; i < 52; i++)
  {
    if (pInput->fFing_Strip_N_LU[i] != pInput->fFing_Strip_N_TU[i])
    {
      //std::cout << event_number << "!!!!!!! Strip : " << i << " Up PMT Leading and Trailing Error " << pInput->fFing_Strip_N_LU[i] << "," << pInput->fFing_Strip_N_TU[i] << std::endl;
      continue;
    }
    if (pInput->fFing_Strip_N_LD[i] != pInput->fFing_Strip_N_TD[i])
    {
      //std::cout << event_number << "!!!!!!! Strip : " << i << " Down PMT Leading and Trailing Error " << pInput->fFing_Strip_N_LD[i] << "," << pInput->fFing_Strip_N_TD[i] << std::endl;
      continue;
    }
    if (pInput->fFing_Strip_N_LD[i] != pInput->fFing_Strip_N_LU[i])
    {
      //std::cout << "Strip : " << i << " Down PMT and Up PMT Mismatch (intended)?" << std::endl;
      continue;
    }

    for (int j = 0; j < pInput->fFing_Strip_N_LD[i]; j++)
    {
      //int pmtUp = i - (i % 2 == 1 ? 0 : 1);
      //int pmtDown = i + (i % 2 == 1 ? 1 : 0);

      //cout << "Strip ID is " << i << ", hit is " << j << endl;
      //cout << "Top PMT is: " << pmtUp << " Times: " << pInput->fFing_Lead_Up[i][j] << " , " << pInput->fFing_Trail_Up[i][j] << std::endl;
      //cout << "Dwn PMT is: " << pmtDown << " Times: " << pInput->fFing_Lead_Down[i][j] << " , " << pInput->fFing_Trail_Down[i][j] << std::endl;
     pOutput->pFing_Lead_Up[i][j] = pInput->fFing_Lead_Up[i][j];
     pOutput->pFing_Lead_Down[i][j] = pInput->fFing_Lead_Down[i][j];
     pOutput->pFing_Trail_Up[i][j] = pInput->fFing_Trail_Up[i][j];
     pOutput->pFing_Trail_Down[i][j] = pInput->fFing_Trail_Down[i][j];
        
     lead_lead = (pInput->fFing_Lead_Up[i][j] - pInput->fFing_Lead_Down[i][j]);
      
      hFING_lead_lead[i] -> Fill(lead_lead);
      pOutput->pFing_lead_lead[i] = lead_lead;
      trail_trail = (pInput->fFing_Trail_Up[i][j]-pInput->fFing_Trail_Down[i][j]);
      hFING_trail_trail[i]->Fill(trail_trail);

      int tot_up = pInput->fFing_Trail_Up[i][j] - pInput->fFing_Lead_Up[i][j];
      int tot_down = pInput->fFing_Trail_Down[i][j] - pInput->fFing_Lead_Down[i][j];
      int tot = tot_up + tot_down;

   //   cout<<"event " << event_number << " tot_up " << tot_up << " tot_down " <<  tot_down << " i " <<i << " j " << j <<endl;

      pOutput->pFing_tot = tot;
      pOutput->pFing_stripID = i;
      hFING_ToT_Strip[i]->Fill(tot);
      hFING_ToT_v_Strip->Fill(i, tot);
      hFing_ToTRatio_v_Strip->Fill(i, (double)tot_up / (tot_up + tot_down));
      if (i > 17 && tot > maxtot)
      {
        maxtot = tot;
        maxtotup = tot_up;
        maxtotdown = tot_down;
        maxtotchan = i;

      }
       maxtotSc41lead = SC41_lead_lead - (lead_lead);

       if(maxtotSc41lead!=0)  hFING_Sc41lead_leadmaxtot[maxtotchan]->Fill(maxtotSc41lead);

    }
  }

  if (maxtotchan > 0)
  {


    hFING_MaxToT_Strip[maxtotchan]->Fill(maxtot);
    hFING_MaxToT_v_Strip->Fill(maxtotchan, maxtot);
    hFing_MaxToTRatio_v_Strip->Fill(maxtotchan, (double)maxtotup / (maxtotup + maxtotdown));
    pOutput->pFing_maxtotchan = maxtotchan;
    pOutput->pFing_maxtot= maxtot;
        }
  return;
    }
//--------------------------------------------------------------------------------------------------------------------//
TH1I* EventAnlProc::MakeH1I(const char* fname,
                            const char* hname,
                            Int_t nbinsx,
                            Float_t xmin, Float_t xmax,
                            const char* xtitle,
                            Color_t linecolor,
                            Color_t fillcolor,
                            const char* ytitle) {
//    TNamed* res = TestObject((getfunc)&TGo4EventProcessor::GetHistogram, fname, hname);
//    if (res!=0) return dynamic_cast<TH1I*>(res);

   TH1I* histo = new TH1I(hname, hname, nbinsx, xmin, xmax);
   histo->SetXTitle(xtitle);
   if (ytitle) histo->SetYTitle(ytitle);
   histo->SetLineColor(linecolor);
   histo->SetFillColor(fillcolor);
   AddHistogram(histo, fname);
   return histo;
}
//-----------------------------------------------------------------------------------------------------------------------------//

TH2I* EventAnlProc::MakeH2I(const char* fname,
                             const char* hname,
                             Int_t nbinsx, Float_t xmin, Float_t xmax,
                             Int_t nbinsy, Float_t ymin, Float_t ymax,
                             const char* xtitle, const char* ytitle,
                             Color_t markercolor) {
//    TNamed* res = TestObject((getfunc)&TGo4EventProcessor::GetHistogram, fname, hname);
//    if (res!=0) return dynamic_cast<TH2I*>(res);

   TH2I* histo = new TH2I(hname, hname, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
   histo->SetMarkerColor(markercolor);
   histo->SetXTitle(xtitle);
   histo->SetYTitle(ytitle);
   AddHistogram(histo, fname);
   return histo;
}
//-----------------------------------------------------------------------------------------------------------------------------//

TGo4WinCond* EventAnlProc::MakeWindowCond(const char* fname,
                                           const char* cname,
                                           float left,
                                           float right,
                                           const char* HistoName) {
  // TNamed* res = TestObject((getfunc)&TGo4EventProcessor::GetAnalysisCondition, fname, cname);
   //if (res!=0) return dynamic_cast<TGo4WinCond*>(res);

   TGo4WinCond* cond = new TGo4WinCond((Text_t*)cname);
   cond->SetValues(left, right);
   cond->Enable();
   if (HistoName!=0)
     cond->SetHistogram(HistoName);
   AddAnalysisCondition(cond, fname);
   return cond;
}
//-----------------------------------------------------------------------------------------------------------------------------//

TGo4PolyCond* EventAnlProc::MakePolyCond(const char* fname,
                                          const char* cname,
                                          Int_t size,
                                          Float_t (*points)[2],
                                          const char* HistoName) {
   //TNamed* res = TestObject((getfunc)&TGo4EventProcessor::GetAnalysisCondition, fname, cname);
   //if (res!=0) return dynamic_cast<TGo4PolyCond*>(res);

   Float_t fullx[size+1], fully[size+1];
   int numpoints = size;

   for (int i=0;i<numpoints;i++) {
     fullx[i] = points[i][0];
     fully[i] = points[i][1];
   }

   // connect first and last points
   if ((fullx[0]!=fullx[numpoints-1]) || (fully[0]!=fully[numpoints-1])) {
      fullx[numpoints] = fullx[0];
      fully[numpoints] = fully[0];
      numpoints++;
   }

   TCutG mycat("initialcut", numpoints, fullx, fully);
   TGo4PolyCond* cond = new TGo4PolyCond((Text_t*)cname);
   cond->SetValues(&mycat);
   cond->Enable();
   if (HistoName!=0)
     cond->SetHistogram(HistoName);
   AddAnalysisCondition(cond, fname);
   return cond;
}

void EventAnlProc::checkTAMEXorVME(){

  std::ifstream PL_FILE("Configuration_Files/TAMEX_or_VME.txt");

  std::string line;

  if(PL_FILE.fail()){
    std::cerr << "Could not find Configuration_Files/TAMEX_or_VME.txt file" << std::endl;
    exit(1);
  }
  bool T_or_V_bPlas = false;
  bool T_or_V_Fatima = false;
  bool T_and_V_Fatima = false;
  while(std::getline(PL_FILE,line)){
    if(line[0] == '#') continue;

    if(line == "VME_bPlas") T_or_V_bPlas = true;
    if(line == "TAMEX_bPlas") T_or_V_bPlas = false;

    if(line == "VME_Fatima") T_or_V_Fatima = true;
    if(line == "TAMEX_Fatima") T_or_V_Fatima = false;

    if(line == "VME_AND_TAMEX_Fatima") T_or_V_Fatima = false;
    if(line == "VME_AND_TAMEX_Fatima") T_and_V_Fatima = true;

    if(line == "VME_Fatima") T_and_V_Fatima = false;
    if(line == "TAMEX_Fatima") T_and_V_Fatima = false;


//     if(line != "VME_bPlas" && line != "TAMEX_bPlas"){
//       std::cerr << line << " module of PLASTIC not known!" <<std::endl;
//       exit(1);
//     }
  }

  VME_TAMEX_bPlas = T_or_V_bPlas;
  VME_TAMEX_Fatima = T_or_V_Fatima;
  VME_AND_TAMEX_Fatima = T_and_V_Fatima;

}
//-----------------------------------------------------------------------------------------------------------------------------//
//                                                            END                                                              //
//-----------------------------------------------------------------------------------------------------------------------------//
