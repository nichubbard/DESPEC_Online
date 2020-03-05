#include "EventUnpackProc.h"
#include "EventUnpackStore.h"
#include "Riostream.h"

// Root Includes //
#include "TROOT.h"
// #include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TCutG.h"
#include "TArc.h"
#include "TTree.h"

#include <time.h>
#include <math.h>
#include <iomanip>

// Go4 Includes //
#include "TGo4UserException.h"
#include "TGo4Picture.h"
#include "Go4StatusBase/TGo4Picture.h"
#include "TGo4MbsEvent.h"

#include "TGo4MbsSubEvent.h"

// General Includes //
#include <fstream>
#include <vector>

#include "Detector_System.cxx"
#include "AIDA_Detector_System.h"
#include "FATIMA_Detector_System.h"
#include "FATIMA_TAMEX_Detector_System.h"
#include "PLASTIC_TAMEX_Detector_System.h"
#include "PLASTIC_VME_Detector_System.h"
#include "FINGER_Detector_System.h"
#include "GALILEO_Detector_System_TEST.h"
#include "FRS_Detector_System.h"

#include "TAidaConfiguration.h"

#include "CalibParameter.h"
#include "CorrelParameter.h"

#include "Data_Stream.cxx"
#include "White_Rabbit.h"

#include <string>


using namespace std;
// std::unordered_map<int, std::string> names =
// {
//   { 0x100, "FRS" },
//   { 0x700, "AIDA" },
//   { 0x500, "bPlas" },
//   { 0x1500, "FATIMA" },
//   { 0x400, "GALILEO" },
// };

//***********************************************************
EventUnpackProc::EventUnpackProc() :TGo4EventProcessor("Proc")
{
  cout << "**** EventUnpackProc: Create instance " << endl;


}
//***********************************************************
// standard factory
EventUnpackProc::EventUnpackProc(const char* name) : TGo4EventProcessor(name)
{


  cout << "**** EventUnpackProc: Create" << endl;

  //  input_data_path_old = "old";
   // WR_out.open ("WR_diff_store_270289.txt");
  WR_used = false;

  //used_systems
  get_used_Systems();
  get_WR_Config();
    //read_setup_parameters();

  //  FAT_det_pos_setup();
    checkTAMEXorVME();
    checkPADI_or_PADIWA();
  //create White Rabbit obj
  WR = new White_Rabbit();

  fCal = new CalibParameter("CalibPar");
  AddParameter(fCal);
  if (fCal) fCal->PrintParameter(0,0);
  else cout << "**** ERRR - CalibPar doesn't exist - program will crash.\n";

  fCorrel = new CorrelParameter("CorrelPar");
  AddParameter(fCorrel);
  if (fCorrel) fCorrel->PrintParameter(0,0);
  else cout << "**** ERRR - CorrelPar doesn't exist - program will crash.\n";

  //create Detector Systems
  Detector_Systems = new Detector_System*[7];

  // all non used systems intialized as NULL
  //-> calling uninitialized system will cause an error !

  Detector_Systems[0] = !Used_Systems[0] ? nullptr : new FRS_Detector_System();
  Detector_Systems[1] = !Used_Systems[1] ? nullptr : new AIDA_Detector_System();
 // Detector_Systems[2] = !Used_Systems[2] ? nullptr : new PLASTIC_VME_Detector_System();
  Detector_Systems[5] = !Used_Systems[5] ? nullptr : new GALILEO_Detector_System();
  Detector_Systems[6] = !Used_Systems[6] ? nullptr : new FINGER_Detector_System();

  if(VME_TAMEX_bPlas==true) Detector_Systems[2] = !Used_Systems[2] ? nullptr : new PLASTIC_VME_Detector_System();
  if(VME_TAMEX_bPlas==false) Detector_Systems[2] = !Used_Systems[2] ? nullptr : new PLASTIC_TAMEX_Detector_System();

  if(VME_TAMEX_Fatima==true || VME_AND_TAMEX_Fatima==true) Detector_Systems[3] =  new FATIMA_Detector_System();


  if(VME_TAMEX_Fatima==false|| VME_AND_TAMEX_Fatima==true) Detector_Systems[4] = new FATIMA_TAMEX_Detector_System();
//   cout << "VME_TAMEX_Fatima "<<VME_TAMEX_Fatima <<  " Used_Systems[3] " << Used_Systems[3] << " Used_Systems[4] " << Used_Systems[4] << " VME_AND_TAMEX_Fatima " << VME_AND_TAMEX_Fatima << endl;
//   if(VME_AND_TAMEX_Fatima==true){
//    Detector_Systems[4] = !Used_Systems[4] ? nullptr : new FATIMA_Detector_System();
//    Detector_Systems[3] = !Used_Systems[3] ? nullptr : new FATIMA_TAMEX_Detector_System();
//
// }

  //          else Detector_Systems[5] = nullptr;
  // Detector_Systems[5] = (!Used_Systems[2] && VME_TAMEX==false) ? nullptr : new PLASTIC_VME_Detector_System();

  //for(int i = 0;i < 7;++i) if(!Used_Systems[i]) Detector_Systems[i] = nullptr;

    PLASTIC_CALIBRATION = Used_Systems[2] ? Check_Cal_Plastic() : false; //Not needed anymore

  //Only create histograms if system is used
  if(Used_Systems[0]) Make_FRS_Histos();

  if(Used_Systems[1]) Make_AIDA_Histos();


  if(Used_Systems[4] && VME_TAMEX_Fatima==false && VME_AND_TAMEX_Fatima==false) Make_FATIMA_TAMEX_Histos();

  if(Used_Systems[3] && VME_TAMEX_Fatima==true && VME_AND_TAMEX_Fatima==false)  Make_FATIMA_Histos();


  if((Used_Systems[3] || Used_Systems[4]) && VME_TAMEX_Fatima==false && VME_AND_TAMEX_Fatima==true) Make_FATIMA_VME_TAMEX_Histos();

  if(Used_Systems[5]) Make_GALILEO_Histos();

  RAW = new Raw_Event();

  load_PrcID_File();

  load_FingerID_File();
  load_FatTamex_Allocationfile();
  read_setup_parameters();

    WR_count = 0;
    count = 0;
    array_count = 0;
    iterator = 0;
    val_it = 0;


  Cout_counter = 0;

  ///Clear for AIDA
  lastTime = 0;
  ID = 0;
  totalEvents = 0;
  startTime = 0;
  stopTime = 0;
  fAida.ImplantEvents.clear();
  fAida.DecayEvents.clear();
  fAida.Implants.clear();
  fAida.Decays.clear();
  /// Setup AIDA arrays
  if(Used_Systems[1])
  {
    TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
    adcLastTimestamp.resize(conf->FEEs());
    adcCounts.resize(conf->FEEs());
  }

}

void EventUnpackProc::UserPostLoop()
{
  if (Used_Systems[1])
  {
    ((AIDA_Detector_System*)Detector_Systems[1])->PrintStatistics();
  }
}


//----------------------------------------------------------


EventUnpackProc::~EventUnpackProc()
{

  delete[] Detector_Systems;
  delete RAW;
  delete WR;
  cout << "**** EventUnpackProc: Delete instance" << endl;
}

//----------------------------------------------------------
Bool_t EventUnpackProc::BuildEvent(TGo4EventElement* dest)
{
  bool skip;
  static TString oldfname = "";
  TGo4MbsEvent *fMbsEvent = dynamic_cast<TGo4MbsEvent*>    (GetInputEvent("Unpack"));
// s_filhe* fileheader=fMbsEvent->GetMbsSourceHeader();
  //s_bufhe* head = GetMbsBufferHeader();

  EventUnpackStore* fOutput = (EventUnpackStore*) dest;
  TGo4MbsEvent* fInput = (TGo4MbsEvent*) GetInputEvent();

  if(fOutput==0){
    cout << "UnpackProc: no unpack output event !";
    return false;
  }
 // input_data_path = fileheader->filhe_file;

  count++;

  if (count % 100000 == 0){

    cout << "\r";
    cout << "Event " << count << " Reached!!!"<<"    Data File Number : "<<data_file_number;
    cout <<"\t\t\t\t";
    cout.flush();
  }

  Bool_t isValid=kFALSE; // validity of output event //

  if (fInput==0) // Ensures that there is data in the event //
  {
    cout << "EventUnpackProc: no input event !"<< endl;
    fOutput->SetValid(isValid);
    return isValid;
  }
  isValid=kTRUE;
  event_number=fInput->GetCount();
  fOutput-> fevent_number = event_number;

  fInput->ResetIterator();
  TGo4MbsSubEvent* psubevt(0);


  // ------------------------------------------------------ //
  // |                                                    | //
  // |               START OF EVENT ANALYSIS              | //
  // |                                                    | //
  // ------------------------------------------------------ //


  if (event_number>0){
  
      int subevent_iter = 0;
      Int_t PrcID_Conv = 0;

      Int_t* pdata = nullptr;
      Int_t lwords = 0;
      Int_t PrcID = 0;
      Int_t sub_evt_length = 0;
      WR_tmp = 0;
      WR_d=0;
      AIDA_Loop = 0;
      WR_main=0;

      while ((psubevt = fInput->NextSubEvent()) != 0) // subevent loop //
      {
        subevent_iter++;
        pdata = psubevt->GetDataField();
        lwords = psubevt->GetIntLen();
        PrcID = psubevt->GetProcid();
        PrcID_Conv = get_Conversion(PrcID);
        fOutput -> fProcID[PrcID_Conv] = PrcID_Conv;

        for (int i =0; i<7;i++){
          fOutput->fUsed_Systems[i] = Used_Systems[i];
     }

        sub_evt_length  = (psubevt->GetDlen() - 2) / 2;

    ///------------------------------WHITE RABBIT --------------------------------------////
        if(WHITE_RABBIT_USED){
          //sub_evt_length = sub_evt_length - 5;

            //Pulls it straight from White_Rabbit class
            WR_tmp = WR->get_White_Rabbit(pdata);
            WR_d = WR->get_Detector_id();
            pdata = WR->get_pdata();

            ///Temp Wr detector fix
            if(PrcID ==10|| PrcID == 30 || PrcID==20 || PrcID == 25)WR_d=0;
          //cout<<"WR_tmp " << WR_tmp << " WR_d " << WR_d<<" PrcID " << PrcID<< endl;

            //FRS
           // cout<<"WR_d " << WR_d << " WR_tmp "<< WR_tmp << endl;
           if(WR_d==0) fOutput->fFRS_WR = WR_tmp; //FRS
           if(WR_d==1) fOutput->fAIDA_WR = WR_tmp; //AIDA
           if(WR_d==2) fOutput->fbPlas_WR = WR_tmp; //bPlas (TAMEX)
           if(WR_d==3) fOutput->fFat_WR = WR_tmp; //Fatima (VME)
           if(WR_d==4) fOutput->fGal_WR = WR_tmp; //Galileo
           if(WR_d==5) fOutput->fFinger_WR = WR_tmp; //FINGER
           ///if(WR_d==xxx) fOutput->fMon_WR = WR_tmp; ///Monster WR

            WR_main = WR_tmp;
            //if(fOutput->fGal_WR>0 && fOutput->fbPlas_WR>0){
          //cout<<" event " << event_number << " WR_d " <<WR_d<<" WR_main " <<  WR_main <<" fOutput->fGal_WR " << fOutput->fGal_WR << " fOutput->fbPlas_WR " << fOutput->fbPlas_WR << "dT " << fOutput->fbPlas_WR-fOutput->fGal_WR <<   endl;}
           ///NOTE: IMPLEMENT WR Sync Checks
           /// if (wr1 != 0x03e1 || wr2 != 0x04e1 || wr3 != 0x05e1 || wr4 != 0x06e1) continue;
            
            
           // Fill_WRSyncs_Histos();        
            
            
            
        }

///-----------------------------------------------------------------------------------------------------------///
        //if necessary, directly print MBS for wanted Detector_System
        if(PrcID_Conv == AIDA && false) print_MBS(pdata,lwords);
        if(PrcID_Conv == FATIMA && false) print_MBS(pdata,lwords);
        if(PrcID_Conv == PLASTIC && false) print_MBS(pdata,lwords);
        if(PrcID_Conv == GALILEO && false) print_MBS(pdata,lwords);
        // if(PrcID_Conv == FINGER && false) print_MBS(pdata,lwords);

        //=================================================================
        //UNPACKING
        ///send subevent to respective unpacker

        if(Detector_Systems[PrcID_Conv] !=0){
        Detector_Systems[PrcID_Conv]->Process_MBS(psubevt);
        Detector_Systems[PrcID_Conv]->Process_MBS(pdata);
       // cout<<"Detector_Systems[PrcID_Conv] " <<Detector_Systems[PrcID_Conv] <<" PrcID_Conv " <<PrcID_Conv<<endl;
        ///get mbs stream data from unpacker (pointer copy solution)
        pdata = Detector_Systems[PrcID_Conv]->get_pdata();

        ///get data from subevent and send to RAW
        Detector_Systems[PrcID_Conv]->get_Event_data(RAW);
        }


        //=================================================================
        //HISTOGRAM FILLING (only singles)
        FILL_HISTOGRAMS(PrcID_Conv);
        //=================================================================

        pdata = nullptr;

        ///--------------------------------------------------------------------------------------------///
                                /** Unpack Tree for each detector subsystem**/
        ///--------------------------------------------------------------------------------------------///
                                                /** Output FRS **/
        ///--------------------------------------------------------------------------------------------///

   if (Used_Systems[0] && PrcID_Conv==0){

        ///MUSIC
//            for(int i =0; i<2; ++i){
//             //fOutput->fFRS_Music_dE[i] = RAW->get_FRS_MusicdE(i);
//           //  fOutput->fFRS_Music_dE_corr[i] = RAW->get_FRS_MusicdE_corr(i);
//            
//            }
        ///SCI
//        if(RAW->get_FRS_sci_l(2)>0){
//        fOutput->fFRS_sci_l2=RAW->get_FRS_sci_l(2);
//      fOutput->fFRS_sci_l[2] = RAW->get_FRS_sci_l(2);
//      cout<<"event " << event_number <<"  fOutput->fFRS_sci_l[2] " << fOutput->fFRS_sci_l[2] << " RAW->get_FRS_sci_l(2) " << RAW->get_FRS_sci_l(2) << endl;
//        }
           for(int l=0;l<12;++l){
   if(RAW->get_FRS_sci_e(l)>0) fOutput->fFRS_sci_e[l] = RAW->get_FRS_sci_e(l);
   if(RAW->get_FRS_sci_l(l)>0) fOutput->fFRS_sci_l[l] = RAW->get_FRS_sci_l(l);
   if(RAW->get_FRS_sci_r(l)>0) fOutput->fFRS_sci_r[l] = RAW->get_FRS_sci_r(l);
            
   if(RAW->get_FRS_sci_tx(l)>0) fOutput->fFRS_sci_tx[l] = RAW->get_FRS_sci_tx(l);
   if(RAW->get_FRS_sci_x(l)>0) fOutput->fFRS_sci_x[l] = RAW->get_FRS_sci_x(l);
             }
           
            ///SCI TOF
//         fOutput->fFRS_sci_tofll2 = RAW->get_FRS_tofll2();
//         fOutput->fFRS_sci_tofll3 = RAW->get_FRS_tofll3();
//         fOutput->fFRS_sci_tof2 = RAW->get_FRS_tof2();
//         fOutput->fFRS_sci_tofrr2 = RAW->get_FRS_tofrr2();
//         fOutput->fFRS_sci_tofrr3 = RAW->get_FRS_tofrr3();
//         fOutput->fFRS_sci_tof3 = RAW->get_FRS_tof3();
        ///ID 2 4
        fOutput->fFRS_ID_x2 = RAW->get_FRS_x2();
        
        //fOutput->fFRS_ID_y2 = RAW->get_FRS_y2();
//         fOutput->fFRS_ID_a2 = RAW->get_FRS_a2();
//         fOutput->fFRS_ID_b2 = RAW->get_FRS_b2();
        
        fOutput->fFRS_ID_x4 = RAW->get_FRS_x4();
//         fOutput->fFRS_ID_y4 = RAW->get_FRS_y4();
//         fOutput->fFRS_ID_a4 = RAW->get_FRS_a4();
//         fOutput->fFRS_ID_b4 = RAW->get_FRS_b4();
            ///SCI dT
//         fOutput->fFRS_sci_dt_21l_21r = RAW->get_FRS_dt_21l_21r();
//         fOutput->fFRS_sci_dt_41l_41r = RAW->get_FRS_dt_41l_41r();
//         fOutput->fFRS_sci_dt_42l_42r = RAW->get_FRS_dt_42l_42r();
//         fOutput->fFRS_sci_dt_43l_43r = RAW->get_FRS_dt_43l_43r();
// 
//         fOutput->fFRS_sci_dt_21l_41l = RAW->get_FRS_dt_21l_41l();
//         fOutput->fFRS_sci_dt_21r_41r = RAW->get_FRS_dt_21r_41r();
// 
//         fOutput->fFRS_sci_dt_21l_42l = RAW->get_FRS_dt_21l_42l();
//         fOutput->fFRS_sci_dt_21r_42r = RAW->get_FRS_dt_21r_42r();
            ///ID Beta Rho
//         for(int i =0; i<2; ++i){
//             fOutput->fFRS_ID_brho[i] = RAW->get_FRS_brho(i);
//             fOutput->fFRS_ID_rho[i] = RAW->get_FRS_rho(i);
//         }
//         fOutput->fFRS_beta = RAW->get_FRS_beta();
//         fOutput->fFRS_beta3 = RAW->get_FRS_beta3();
//         fOutput->fFRS_gamma = RAW->get_FRS_gamma();
            ///ID Z AoQ
       // cout<<"1111event_ " << event_number << " RAW->get_FRS_AoQ_corr() " <<RAW->get_FRS_AoQ_corr() <<" Unpack PRCID " << PrcID << endl;
       // cout<<" " << endl;
        if(PrcID==25){
        fOutput->fFRS_AoQ = RAW->get_FRS_AoQ();
        fOutput->fFRS_AoQ_corr = RAW->get_FRS_AoQ_corr();
        
        //cout <<"2222event_ " << event_number << " RAW->get_FRS_AoQ_corr() " <<RAW->get_FRS_AoQ_corr() << endl;
        //    cout<<" " << endl;
        }
     
        fOutput->fFRS_z = RAW->get_FRS_z();
        fOutput->fFRS_z2 = RAW->get_FRS_z2();
        //fOutput->fFRS_z3 = RAW->get_FRS_z3();
            ///ID Timestamp
//         fOutput->fFRS_timestamp = RAW->get_FRS_timestamp();
//         fOutput->fFRS_ts = RAW->get_FRS_ts();
//         fOutput->fFRS_ts2 = RAW->get_FRS_ts2();
         }
         ///--------------------------------------------------------------------------------------------///
                                            /** Output AIDA **/
        ///--------------------------------------------------------------------------------------------///

        if (Used_Systems[AIDA] && PrcID_Conv==AIDA){
          TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
          AIDA_Hits = RAW->get_AIDA_HITS();

          AidaEvent evt;
          fOutput->fAIDAHits += AIDA_Hits;

          for(int i = 0; i<AIDA_Hits; i++){

            AIDA_Energy[i] = RAW->get_AIDA_Energy(i);

            AIDA_FEE[i] = RAW-> get_AIDA_FEE_ID(i);
            AIDA_ChID[i] = RAW-> get_AIDA_CHA_ID(i);
            AIDA_Time[i] = RAW-> get_AIDA_WR(i);
            AIDA_HighE_veto[i] = RAW-> get_AIDA_HighE_VETO(i);
            AIDA_Side[i] = RAW-> get_AIDA_SIDE(i);
            AIDA_Strip[i] = RAW-> get_AIDA_STRIP(i);
            AIDA_evtID[i] = RAW-> get_AIDA_EVTID(i);

            evt.Channel = AIDA_ChID[i];
            evt.Module = AIDA_FEE[i];
            evt.Time = AIDA_Time[i];
            evt.HighEnergy =  AIDA_HighE_veto[i];
            evt.DSSD = conf->FEE(AIDA_FEE[i]).DSSD;
            evt.Side = AIDA_Side[i];
            evt.Strip = AIDA_Strip[i];
            evt.ID = AIDA_evtID[i];
            evt.Energy = AIDA_Energy[i];
            evt.FastTime = RAW->get_AIDA_FastTime(i);
            //if(evt.HighEnergy>0){
       // cout<<"Event AIDA " << event_number << " evt.Energy " << evt.Energy <<" evt.HighEnergy " << evt.HighEnergy <<endl;}
            if (!startTime) startTime = evt.Time;
            stopTime = evt.Time;
            /// Build events from everything until there's a gap of 2000 �s (event window)

            /// If lastTime is 0 it's the first event
            /// New event for timewarps
            if (lastTime > 0 && (evt.Time - lastTime > conf->EventWindow()))
            {
              // if event happened too late, redo the event again with a new out_event
              lastTime = 0;
              ResetMultiplexer();

              totalEvents++;

              fOutput->Aida.push_back(fAida);
              fAida.ImplantEvents.clear();
              fAida.DecayEvents.clear();
              fAida.AIDATime = 0;
            }

            lastTime = evt.Time;
            CorrectTimeForMultiplexer(evt);


            if (evt.HighEnergy==1)
            {
              fAida.ImplantEvents.push_back(evt);
            }
            else
            {
              fAida.DecayEvents.push_back(evt);
            }

            if (fAida.AIDATime == 0)
            {
              fAida.AIDATime = evt.Time;
            }
          }

          if (conf->ucesb())
          {
            lastTime = 0;
            ResetMultiplexer();
            totalEvents++;

            fOutput->Aida.push_back(fAida);
            fAida.ImplantEvents.clear();
            fAida.DecayEvents.clear();
            fAida.AIDATime = 0;
          }
        }
        ///--------------------------------------------------------------------------------------------///
                                        /** Output bPlastic + SCALAR **/
        ///--------------------------------------------------------------------------------------------///

        //  int Plas_VME_QDC_ID;
        int Plas_VME_TDC_ID_sing, Plas_VME_TDC_ID[50], bPlasTDCmulti[50];
        Plas_VME_TDC_ID_sing=-1;

        for (int i =0; i<50; i++){
          Plas_VME_TDC_ID[i] = -1;
          bPlasTDCmulti[i] = 0;
        }

        if (Used_Systems[2] && PrcID_Conv==2){
          fOutput->fbPlas_VME_firedQDC = RAW->get_plastic_VME_QDC_fired();

          for (int i = 0; i < RAW->get_plastic_VME_QDC_fired();  i++){
            // Plas_VME_QDC_ID = RAW->get_plastic_VME_QDC_cha(i);
            // bPlasQDCmulti[Plas_VME_QDC_ID]++;

            fOutput-> fbPlas_VME_QDC_ID[i] = RAW->get_plastic_VME_QDC_cha(i);
            fOutput-> fbPlas_VME_QDC_E[i] = RAW->get_plastic_VME_QDC_dat1(i);
            // fOutput-> fbPlas_VME_QDC_Multiplicity = bPlasQDCmulti[Plas_VME_QDC_ID];

            }
          fOutput->fbPlas_VME_firedTDC = RAW->get_plastic_VME_TDC_fired();

          for (int j=0; j < RAW->get_plastic_VME_TDC_fired(); j++){
            Plas_VME_TDC_ID[j] = RAW->get_plastic_VME_TDC_cha(j); //Plastic TDC ID
            Plas_VME_TDC_ID_sing = RAW->get_plastic_VME_TDC_cha(j);

            bPlasTDCmulti[Plas_VME_TDC_ID_sing]++;
            fOutput->fbPlas_VME_TDC_ID[j] = Plas_VME_TDC_ID[j];
            fOutput->fbPlas_VME_TDC_TS[j][Plas_VME_TDC_ID[j]] = RAW->get_plastic_VME_TDC_dat(Plas_VME_TDC_ID[j]);//25ps //Plastic TDC Data

            fOutput->fbPlas_VME_TDC_Multiplicity[Plas_VME_TDC_ID_sing] =  bPlasTDCmulti[Plas_VME_TDC_ID_sing];

          }

//           fOutput->fScalar_fired = RAW->get_scalar_iterator();
//           for (int g=0; g < RAW->get_scalar_iterator(); g++){
// //            fOutput-> fScalar_ID = RAW->get_scalar_chan(g);
// 
//           }
        }
          ///--------------------------------------------------------------------------------------------///
                                                /**Output bPLASTIC TAMEX and FATIMA **/
        ///--------------------------------------------------------------------------------------------///
        int Fatfired[4];
        int  bPlasfired[2];
        int Phys_Channel_Lead_Fat[4][256];
        int Phys_Channel_Trail_Fat[4][256];
        int Phys_Channel_Lead_bPlas[2][256];
        int Phys_Channel_Trail_bPlas[2][256];
     if (Used_Systems[2]&& PrcID_Conv==2 && VME_TAMEX_bPlas == false){
          
         for(int f=FATIMA_TAMEX_MODULES;f<bPLASTIC_TAMEX_MODULES;f++){       
           bPlasfired[f]  = 0;
           for(int g=0; g<bPLASTIC_TAMEX_HITS; g++){
           Phys_Channel_Lead_bPlas[f][g] = 0;
           Phys_Channel_Trail_bPlas[f][g] = 0;
           }
         }
           for(int p=0;p<FATIMA_TAMEX_MODULES;p++){
              Fatfired[p] = 0;         
             for(int q=0; q<FATIMA_TAMEX_HITS; q++){
                 Phys_Channel_Lead_Fat[p][q] = 0;
                 Phys_Channel_Trail_Fat[p][q] = 0;
                 
             }
         }

          for (int i=0; i<RAW->get_PLASTIC_tamex_hits(); i++){///Loop over tamex ID's
         ///--------------------------------------------------------------------------------------------///
                                                /**Output FATIMA TAMEX **/
        ///--------------------------------------------------------------------------------------------///
               
            if(RAW->get_PLASTIC_TAMEX_ID(i) < FATIMA_TAMEX_MODULES){
               
                Fatfired[i] = RAW->get_PLASTIC_am_Fired(i);

            for(int j = 0;j < Fatfired[i];j++){

              if(RAW->get_PLASTIC_CH_ID(i,j) % 2 == 1){ //Lead odd j
                Phys_Channel_Lead_Fat[i][j] =TAMEX_bPlasFat_ID[i][RAW->get_PLASTIC_physical_channel(i, j)]; //From allocation file in future
                int chan_fat = Phys_Channel_Lead_Fat[i][j];
            
                int N1 = fOutput->fFat_PMT_Lead_N[chan_fat]++;
                fOutput->fFat_Lead_PMT[chan_fat][N1] = RAW->get_PLASTIC_lead_T(i,j);
               
              }
              else{ //Trail even j
                Phys_Channel_Trail_Fat[i][j] = TAMEX_bPlasFat_ID[i][RAW->get_PLASTIC_physical_channel(i, j)];

                int chan_fat = Phys_Channel_Trail_Fat[i][j];

                // PMT allocation succeeded
                int N1 = fOutput->fFat_PMT_Trail_N[chan_fat]++;
                fOutput->fFat_Trail_PMT[chan_fat][N1] = RAW->get_PLASTIC_trail_T(i,j);
              // cout<<"fOutput->fFat_Trail_PMT[chan][N1] " << fOutput->fFat_Trail_PMT[chan][N1] << " chan " << chan << " N1 "<< N1 <<endl;
                        }
                    }
                }
        
            
        
        ///--------------------------------------------------------------------------------------------///
                                                /**Output bPLASTIC TAMEX  **/
        ///--------------------------------------------------------------------------------------------///
        else{
            int chan=-1;
            
           //fOutput->fbPlas_TAMEX_ID = i;
            bPlasfired[i] = RAW->get_PLASTIC_am_Fired(i); ///Iterator
            int bPlasdetnum=i;
            for(int j = 0;j < bPlasfired[i];j++){

              if(RAW->get_PLASTIC_CH_ID(i,j) % 2 == 1){ //Lead odd j
                Phys_Channel_Lead_bPlas[i][j] = RAW->get_PLASTIC_physical_channel(i, j); //From allocation file in future
                
                chan = (Phys_Channel_Lead_bPlas[i][j])/bPlasdetnum-16;
                fOutput->fbPlaschan=  chan; 
                // PMT allocation succeeded
                int N1 = fOutput->fbPlas_PMT_Lead_N[bPlasdetnum][chan]++;
                fOutput->fbPlas_Lead_PMT[bPlasdetnum][chan][N1] = RAW->get_PLASTIC_lead_T(i,j);
            
            //cout<<"event " << event_number << " chan " << chan <<" bPlasdetnum " << bPlasdetnum<<" Phys_Channel_Lead_bPlas[i][j] " <<Phys_Channel_Lead_bPlas[i][j]<<" fOutput->fbPlas_Lead_PMT[bPlasdetnum][chan][N1]  " << fOutput->fbPlas_Lead_PMT[bPlasdetnum][chan][N1]  << " N1 " << N1<<  " i " << i << " j " << j << endl;
                
              }
              else{ //Trail even j
                  
                Phys_Channel_Trail_bPlas[i][j] = RAW->get_PLASTIC_physical_channel(i,j);

                 chan = (Phys_Channel_Trail_bPlas[i][j])-(bPlasdetnum*16);
                 
                // PMT allocation succeeded
                
                   
                int N1 = fOutput->fbPlas_PMT_Trail_N[bPlasdetnum][chan]++;
                
               
               
              fOutput->fbPlas_Trail_PMT[bPlasdetnum][chan][N1] = RAW->get_PLASTIC_trail_T(i,j);
              // cout<<"event " << event_number << " chan " << chan <<" bPlasdetnum " << bPlasdetnum<<" Phys_Channel_Trail_bPlas[i][j] " <<Phys_Channel_Trail_bPlas[i][j]<<" fOutput->fbPlas_Trail_PMT[bPlasdetnum][chan][N1]  " << fOutput->fbPlas_Trail_PMT[bPlasdetnum][chan][N1]   << " N1 " << N1<<  " i " << i << " j " << j << endl;
                
              }
            }
          }
        }
       }
     
        ///--------------------------------------------------------------------------------------------///
                                                /**Output FATIMA VME **/
        ///--------------------------------------------------------------------------------------------///

        int Fat_QDC_ID;
        int Fat_TDC_ID_sing;
        int Fat_TDC_ID[48];
        int Fat_TDC_multi[48];

        for (int i = 0; i<48; i++){
          Fat_TDC_multi[i] = 0;

        }
        if (Used_Systems[3]&& PrcID_Conv==3){
          fOutput->fFat_firedQDC = RAW->get_FAT_QDCs_fired();

          for (int i=0; i<RAW->get_FAT_QDCs_fired(); i++){

            fOutput->fFat_QDC_ID[i] =  RAW->get_FAT_QDC_id(i);
            fOutput->fFat_QDC_E[i] = RAW->get_FAT_QLong(i);
            fOutput->fFat_QDC_T[i] = RAW->get_FAT_t_qdc(i);
          }
          fOutput->fFat_firedTDC = RAW->get_FAT_TDCs_fired();
          for (int j=0; j<RAW->get_FAT_TDCs_fired(); j++){
            Fat_TDC_ID[j] =  RAW->get_FAT_TDC_id(j);
            Fat_TDC_ID_sing = RAW->get_FAT_TDC_id(j);

            Fat_TDC_multi[Fat_TDC_ID_sing]++;

            fOutput->fFat_TDC_ID[j] =  RAW->get_FAT_TDC_id(j);
            fOutput->fFat_TDC_TS[j][Fat_TDC_ID[j]] = RAW->get_FAT_TDC_timestamp(j); //In 25ps
            fOutput->fFat_TDC_Multiplicity[Fat_TDC_ID_sing] =  Fat_TDC_multi[Fat_TDC_ID_sing];
            
            /** SC41 Trigger **/
            if(Fat_TDC_ID[j]==34){
            fOutput->fSC41[j] =  RAW->get_FAT_TDC_timestamp(j);
            }
          }
        }

       
//         cout<<"VME_TAMEX_Fatima "<< VME_TAMEX_Fatima <<endl;
      if (Used_Systems[2]&& PrcID_Conv==2 && VME_TAMEX_Fatima==false){

         
      }
        ///--------------------------------------------------------------------------------------------///
                                            /**Output GALILEO **/
        ///--------------------------------------------------------------------------------------------///
        if (Used_Systems[5]&& PrcID_Conv==5){
         for (int i=fOutput->fGal_fired; i<RAW->get_GALILEO_am_Fired() && i < GALILEO_MAX_HITS; i++){
                fOutput->fGal_Detector[i] =  RAW->get_GALILEO_Det_id(i);
                fOutput->fGal_Crystal[i] =  RAW->get_GALILEO_Crystal_id(i);
                fOutput->fGal_E[i] = RAW->get_GALILEO_Chan_E(i);
                fOutput->fGal_T[i] = RAW->get_GALILEO_Chan_T(i);
                fOutput->fGal_Pileup[i] = RAW->get_GALILEO_Pileup(i);
                fOutput->fGal_Overflow[i] = RAW->get_GALILEO_Overflow(i);
                fOutput->fGal_fired++;
            
            
          }
        }
        ///--------------------------------------------------------------------------------------------///
                                        /** Output FINGER **/
  ///--------------------------------------------------------------------------------------------///
      if (Used_Systems[6]&& PrcID_Conv==6){

          int Phys_Channel_Lead[FINGER_TAMEX_MODULES][FINGER_TAMEX_HITS] = {0,0};
          int Phys_Channel_Trail[FINGER_TAMEX_MODULES][FINGER_TAMEX_HITS] = {0,0};

          int fingfired[FINGER_TAMEX_MODULES] = {0};

          for (int i=0; i<RAW->get_FINGER_tamex_hits(); i++){
            fingfired[i] = RAW->get_FINGER_am_Fired(i);
     
            for(int j = 0;j < fingfired[i];j++){
        
              if(RAW->get_FINGER_CH_ID(i,j) % 2 == 1){ //Lead odd j
                Phys_Channel_Lead[i][j] = fingID[i][RAW->get_FINGER_physical_channel(i, j)]; //From allocation file
                int chan = Phys_Channel_Lead[i][j];

                if (chan < 0)
                  continue;

                // PMT allocation succeeded
                int N1 = fOutput->fFing_PMT_Lead_N[chan]++;
                fOutput->fFing_Lead_PMT[chan][N1] = RAW->get_FINGER_lead_T(i,j);
  
                // PMT "0" is the trigger
                if (chan == 0 || chan == 1){
                    fOutput->fFing_SC41_lead[chan][N1] = RAW->get_FINGER_lead_T(i,j);
                  continue;
                }
                // chan = "PMT" number
                // this maps to two strips to fill in
                if (chan % 2 == 0) // even PMT = up pmts
                {
                  int strip1 = chan;
                  int strip2 = chan + 1;
                  int N1 = fOutput->fFing_Strip_N_LU[strip1]++;
                  int N2 = fOutput->fFing_Strip_N_LU[strip2]++;
                  fOutput->fFing_Lead_Up[strip1][N1] = RAW->get_FINGER_lead_T(i,j);
                  fOutput->fFing_Lead_Up[strip2][N2] = RAW->get_FINGER_lead_T(i,j);
                  fOutput->fFing_Strip_N[strip1]++;
                  fOutput->fFing_Strip_N[strip2]++;
                                }
                else // odd = lower PMT
                {
                  int strip1 = chan + 1;
                  int strip2 = chan;
                  int N1 = fOutput->fFing_Strip_N_LD[strip1]++;
                  int N2 = fOutput->fFing_Strip_N_LD[strip2]++;
                  fOutput->fFing_Lead_Down[strip1][N1] = RAW->get_FINGER_lead_T(i,j);
                  fOutput->fFing_Lead_Down[strip2][N2] = RAW->get_FINGER_lead_T(i,j);
                      }
              }
              else{ //Trail even j
                Phys_Channel_Trail[i][j] = fingID[i][RAW->get_FINGER_physical_channel(i,j)];
             
                int chan = Phys_Channel_Trail[i][j];
                if (chan < 0)
                  continue;

                // PMT allocation succeeded
                int N1 = fOutput->fFing_PMT_Trail_N[chan]++;
                fOutput->fFing_Trail_PMT[chan][N1] = RAW->get_FINGER_trail_T(i,j);
                 // PMT "0" is the trigger
                if (chan == 0 || chan == 1){
                    fOutput->fFing_SC41_trail[chan][N1] = RAW->get_FINGER_trail_T(i,j);

                  continue;
                }
                if (chan % 2 == 0) // even PMT = up pmts
                {
                  int strip1 = chan + 1;
                  int strip2 = chan;
                  int N1 = fOutput->fFing_Strip_N_TU[strip1]++;
                  int N2 = fOutput->fFing_Strip_N_TU[strip2]++;
                  fOutput->fFing_Trail_Up[strip1][N1] = RAW->get_FINGER_trail_T(i,j);
                  fOutput->fFing_Trail_Up[strip2][N2] = RAW->get_FINGER_trail_T(i,j);
                 }
                else // odd = lower PMT
                {
                  int strip1 = chan + 1;
                  int strip2 = chan;
                  int N1 = fOutput->fFing_Strip_N_TD[strip1]++;
                  int N2 = fOutput->fFing_Strip_N_TD[strip2]++;
                  fOutput->fFing_Trail_Down[strip1][N1] = RAW->get_FINGER_trail_T(i,j);
                  fOutput->fFing_Trail_Down[strip2][N2] = RAW->get_FINGER_trail_T(i,j);
                }
              }
            }
          }
        }      
        ///--------------------------------------------------------------------------------------------///

      } //End of subevent loop


      fOutput->SetValid(isValid);

      pdata = nullptr;
    //} //End of Skip
  }
  //
  return isValid;

}

void EventUnpackProc::FILL_HISTOGRAMS(int PrcID_Conv){
    
 // switch(PrcID_Conv){
  //  case 0:
    ///WARNING USED SYSTEMS!!!
  if(PrcID_Conv==0 && Used_Systems[0])  Fill_FRS_Histos();
  //else break;
   // break;
   // case 1:
    //   Process_AIDA_Event(fOutput);
   if(PrcID_Conv==1)  Fill_AIDA_Histos();
  //  break;
  //  case 2:
   // if(!PLASTIC_CALIBRATION && VME_TAMEX_bPlas==false && PrcID_Conv==2) Fill_Plastic_Histos();
    if(!PLASTIC_CALIBRATION && VME_TAMEX_bPlas==true && PrcID_Conv==2 ) Fill_Plastic_VME_Histos();
//     break;
//     case 3:
  //  Fill_FATIMA_Histos();
    if(PrcID_Conv==3) Fill_FATIMA_Histos();

    if(VME_TAMEX_Fatima==false && VME_AND_TAMEX_Fatima==false && PrcID_Conv==4) Fill_FATIMA_TAMEX_Histos();
    //case 4:

    if(VME_AND_TAMEX_Fatima==true && PrcID_Conv==4 || PrcID_Conv==3){ Fill_FATIMA_VME_TAMEX_Histos();

    }
    //cout <<"VME_AND_TAMEX_Fatima " <<VME_AND_TAMEX_Fatima << " PrcID_Conv " << PrcID_Conv<<endl;
//     break;
//     case 5:
    if(PrcID_Conv==5) Fill_GALILEO_Histos();
//     break;
//     case 6:
    //Fill_Finger_Histos();
//     break;
//     default:
  //  cerr << "PrcID_Conv " << PrcID_Conv << " not known" << endl;
    //exit(0);
 // }
}


//-----------------------------------------------------------------------------------------------------------------------------//
void EventUnpackProc::ResetMultiplexer()
{
  for (int i = 0; i < 12; ++i)
  {
    for (int j = 0; j < 4; ++j)
    {
      adcLastTimestamp[i][j] = 0;
      adcCounts[i][j] = 0;
    }
  }
}



void EventUnpackProc::CorrectTimeForMultiplexer(AidaEvent &evt)
{
  int fee = evt.Module;
  int adc = evt.Channel / 16;
  int64_t time = evt.Time;

  if ((time - adcLastTimestamp[fee][adc] > 2500) && adcLastTimestamp[fee][adc] != 0)
  adcCounts[fee][adc] = 0;

  adcLastTimestamp[fee][adc] = time;

  evt.Time = time - (2000 * adcCounts[fee][adc]++);
  if (evt.HighEnergy) evt.FastTime = evt.Time;
}

//-----------------------------------------------------------------------------------------------------------------------------//


void EventUnpackProc::load_PrcID_File(){
  ifstream data("Configuration_Files/PrcID_to_Det_Sys.txt");
  if(data.fail()){
    cerr << "Could not find PrcID config file!" << endl;
    exit(0);
  }
  int id[7] = {0,0,0,0,0,0,0};
  int i = 0;
  string line;
  char s_tmp[100];
  while(data.good()){
    getline(data,line,'\n');
    if(line[0] == '#') continue;
    sscanf(line.c_str(),"%s %d %d %d %d %d",s_tmp,&id[0],&id[1],&id[2],&id[3],&id[4],&id[5],&id[6]);
    for(int j = 0; j < 6; ++j) PrcID_Array[i][j] = id[j];
    i++;
  }
}
//---------------------------------------------------------------------------------------------------
void EventUnpackProc::load_FatTamex_Allocationfile(){

  const char* format = "%d %d %d";
  ifstream data("Configuration_Files/Fatima_TAMEX_allocation.txt");
  if(data.fail()){
    cerr << "Could not find Fatima_TAMEX_allocation config file!" << endl;
    exit(0);
  }
  //     int id[5] = {0,0,0,0,0};
  //int i = 0;
  int TamID = 0;
  int TamCh = 0;
  int Sys_ch =0;
  string line;
  //char s_tmp[100];
  while(data.good()){

    getline(data,line,'\n');
    if(line[0] == '#') continue;
    sscanf(line.c_str(),format,&TamID,&TamCh,&Sys_ch);
    TAMEX_bPlasFat_ID[TamID][TamCh] = Sys_ch;
  }
}
//---------------------------------------------------------------------------------------------------
void EventUnpackProc::load_FingerID_File(){

  const char* format = "%d %d %d";
  ifstream data("Configuration_Files/Finger_allocation.txt");
  if(data.fail()){
    cerr << "Could not find Finger_allocation config file!" << endl;
    exit(0);
  }
  //     int id[5] = {0,0,0,0,0};
  //int i = 0;
  int tamid = 0;
  int tamch = 0;
  int fingid = 0;
  string line;
  //char s_tmp[100];
  while(data.good()){

    getline(data,line,'\n');
    if(line[0] == '#') continue;
    sscanf(line.c_str(),format,&tamid,&tamch,&fingid);
    fingID[tamid][tamch] = fingid;
  }
}

//-----------------------------------------------------------------------------------------------------------------------------//

void EventUnpackProc::read_setup_parameters(){

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
  file >> WHITE_RABBIT_USED;//dummy_var;


  cout<<endl;
  cout<<endl;
  cout<<"////////////////////////////////////////////////////////////////////////"<<endl;
    cout<<"Setup Parameters List Unpack Proc: "<<endl;
  if(WHITE_RABBIT_USED) cout<<"White Rabbit: Enabled"<<endl;
  else if(!WHITE_RABBIT_USED) cout<<"White Rabbit: Disabled"<<endl;
    cout<<"////////////////////////////////////////////////////////////////////////"<<endl;
  cout<<endl;
  cout<<endl;
 
}

//-----------------------------------------------------------------------------------------------------------------------------//
Int_t EventUnpackProc::get_Conversion(Int_t PrcID){

  for(int i = 0;i < 7;++i){
    for(int j = 0;j < 5;++j){
          
      if(PrcID == PrcID_Array[i][j]) return i;
    }
  }
  cerr << "ProcID " << PrcID << " not known!" << endl;
  exit(0);
}

void EventUnpackProc::get_used_Systems(){
    for(int i = 0;i < 6;i++) Used_Systems[i] = false;

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

    cout << "\n=====================================================" << endl;
    cout << "USED SYSTEMS" << endl;
    cout << "-----------------------------------------------------" << endl;
    for(int j = 0;j < 6;++j){
        if(Used_Systems[j]) cout << DET_NAME[j] << endl;
    }
    cout << "=====================================================" << endl;


}

//-----------------------------------------------------------------------------------------------------------------------------//

void EventUnpackProc::get_WR_Config(){
  ifstream data("Configuration_Files/White_Rabbit.txt");
  if(data.fail()){
    cerr << "Could not find White_Rabbit config file!" << endl;
    exit(0);
  }

  int id = 0;
  string line;
  char s_tmp[100];
  while(data.good()){
    getline(data,line,'\n');
    if(line[0] == '#') continue;
    sscanf(line.c_str(),"%s %d",s_tmp,&id);
    WR_used = (id == 1);
  }
}

  //-----------------------------------------------------------------------------------------------------------------------------//
  // ################################################################## //
  // ################################################################## //
  // ################# Raw Histogram Filling Section ################## //
  // ################################################################## //
  // ################################################################## //
  
  
  
  /**----------------------------------------------------------------------------------------------**/
  /**---------------------------------------------  FRS  ------------------------------------------**/
  /**----------------------------------------------------------------------------------------------**/

  void EventUnpackProc::Make_FRS_Histos(){
 char fname[50], name[50], title[60];//, title2[60];

  const char *count_title1[12]={"(0:1)", "(1:1)", "(2:1)",
                "(2:2)", "(3:1)", "(4:1)",
                "(4:2)", "(4:3)", "(6:1)",
                "(6:2)", "(8:1)"};
  const char *fext1[12]={"0", "1", "2", "2", "3", "4", "4", "4", "6", "6", "8", "8"};
  const char *fext2[12]={"01", "11", "21", "22","31", "41",
             "42", "43", "61",
             "62", "81", "82"};
             
  ///FRS Scalars           
  bool scaler_enable_hist[64];
  char scaler_name[64][256];
  scaler_ch_1kHz=39; //ch7 of 2nd scaler
  scaler_ch_spillstart=8; //ch8 of 1st scaler 
  scaler_check_first_event=1;      
  for(int ii=0; ii<64; ii++){
    sprintf(scaler_name[ii],"scaler_ch%d",ii);//default name
    scaler_enable_hist[ii]=false;
  }
  sprintf(scaler_name[0],"IC01curr-old"); 
  sprintf(scaler_name[1],"SEETRAM-old");
  sprintf(scaler_name[2],"SEETRAM-new");
  sprintf(scaler_name[3],"IC01curr-new");
  sprintf(scaler_name[4],"IC01 count");
  sprintf(scaler_name[5],"SCI00");
  sprintf(scaler_name[6],"SCI01");
  sprintf(scaler_name[7],"SCI02");
  sprintf(scaler_name[8],"Start Extr");
  sprintf(scaler_name[9],"Stop Extr");
  sprintf(scaler_name[10],"Beam Transformer");
  
  sprintf(scaler_name[32],"Free Trigger");
  sprintf(scaler_name[33],"Accept Trigger");
  sprintf(scaler_name[34],"Spill Counter");
  sprintf(scaler_name[35],"1 Hz clock");
  sprintf(scaler_name[36],"10 Hz clock");
  sprintf(scaler_name[37],"100 kHz X veto dead-time");
  sprintf(scaler_name[38],"100 kHz clock");
  sprintf(scaler_name[39],"1 kHz clock");
  
  sprintf(scaler_name[48],"SCI21L");
  sprintf(scaler_name[49],"SCI41L");
  sprintf(scaler_name[50],"SCI42L");
  sprintf(scaler_name[51],"SCI43L");
  sprintf(scaler_name[52],"SCI81L");
  sprintf(scaler_name[53],"SCI21R");
  sprintf(scaler_name[54],"SCI41R");
  sprintf(scaler_name[55],"SCI42R");
  sprintf(scaler_name[56],"SCI43R");
  sprintf(scaler_name[57],"SCI81R");
  sprintf(scaler_name[58],"SCI31L");
  sprintf(scaler_name[59],"SCI31R");
  sprintf(scaler_name[60],"SCI11");
  sprintf(scaler_name[61],"SCI51");
  

  for(int ii=0; ii<64; ii++){
    hScaler_per_s[ii]     = MakeH1I("FRS/Scaler/Scaler_per_1s",Form("%s_per_1s",scaler_name[ii]),1000,0,1000,"Time (s)", 2,5, "Count per second");
    hScaler_per_100ms[ii] = MakeH1I("FRS/Scaler/Scaler_per_0.1s",Form("%s_per_0.1s",scaler_name[ii]),4000,0,400,"Time (s)", 2,5, "Count per 0.1 second");
    hScaler_per_spill[ii] = MakeH1I("FRS/Scaler/Scaler_per_spill",Form("%s_per_spill",scaler_name[ii]),1000,0,1000,"Spill", 2,5, "Count per spill");
  }

  for (int cnt = 0; cnt<7; cnt++) //changed from 3 to 6 04.07.2018
    {
      int index = 0;
      switch(cnt)
    {
        case 0: index = 2; break;
        case 1: index = 3; break;
        case 2: index = 5; break;
        case 3: index = 6; break;
        case 4: index = 7; break;
        case 5: index = 10; break;
        case 6: index = 8; break;
    }
      sprintf(fname,"FRS/SCI/SCI%s/SCI%s",fext1[index],fext2[index]);
      sprintf(name, "SCI%s_L", count_title1[index]);
      sprintf(title, "Sc%s L dE [ch]", count_title1[index]);
      hSCI_L[index] = MakeH1I(fname,name,4096,0,4096,title,2,3);

      sprintf(name, "SCI%s_R", count_title1[index]);
      sprintf(title, "Sc%s R dE [ch]", count_title1[index]);
      hSCI_R[index] = MakeH1I(fname,name,4096,0,4096,title,2,3);

      sprintf(name, "SCI%s_E", count_title1[index]);
      sprintf(title, "Sc%s Energy [ch]", count_title1[index]);
      hSCI_E[index] = MakeH1I(fname,name,4096,0,4096,title,2,3);

      sprintf(name, "SCI%s_Tx", count_title1[index]);
      sprintf(title, "Sc%s t_lr [ch] TAC", count_title1[index]);
      hSCI_Tx[index] = MakeH1I(fname,name,4096,0,4096,title,2,3);

      sprintf(name, "SCI%s_X", count_title1[index]);
      sprintf(title, "Sc%s x-pos [mm]", count_title1[index]);
      hSCI_X[index] = MakeH1I(fname,name,240,-120,120,title,2,3);

    }
      hSCI_dE24 = MakeH2I("FRS/SCI","SCI_dE21-41", 100,10,4000,100,10,4000,"SC21 dE","SC41 dE",2);
     // ToF SC21-SC41
        sprintf(fname,"FRS/SCI/TOF/TOF(%d)",2);
        sprintf(name,"SCI_21_41_TofLL");
        hSCI_TofLL2 = MakeH1I(fname,name,1500,0,62000,"TAC SC41L-SC21L [ps]",2,3);

        sprintf(name,"SCI_21_41_TofRR");
        hSCI_TofRR2 = MakeH1I(fname,name,1500,0,62000,"TAC SC41R-SC21R [ps]",2,3);

        hSCIdE41_TPC42X= MakeH2I("FRS/SCI_TPC/","SCIdE41_TPC42X", 1024,0,4096, 400,-100.,100, "SC41 dE", "TPC42 X[mm]", 2);
        hSCIdE41L_TPC42X= MakeH2I("FRS/SCI_TPC/","SCIdE41L_TPC42X", 1024,0,4096, 400,-100.,100, "SC41L dE", "TPC42 X[mm]", 2);
        hSCIdE41L_TPC41X= MakeH2I("FRS/SCI_TPC/","SCIdE41L_TPC41X", 1024,0,4096, 400,-100.,100, "SC41L dE", "TPC41 X[mm]", 2);
        hSCIdE41R_TPC42X= MakeH2I("FRS/SCI_TPC/","SCIdE41R_TPC42X", 1024,0,4096, 400,-100.,100, "SC41R dE", "TPC42 X[mm]", 2);
        hSCIdE41R_TPC41X= MakeH2I("FRS/SCI_TPC/","SCIdE41R_TPC41X", 1024,0,4096, 400,-100.,100, "SC41R dE", "TPC41 X[mm]", 2);

        hSCIdE21_TPC42X= MakeH2I("FRS/SCI_TPC/","SCIdE21_TPC42X", 1024,0,4096, 400,-100.,100, "SC41 dE", "TPC42 X[mm]", 2);
        hSCIdE21L_TPC42X= MakeH2I("FRS/SCI_TPC/","SCIdE21L_TPC42X", 1024,0,4096, 400,-100.,100, "SC21L dE", "TPC42 X[mm]", 2);
        hSCIdE21L_TPC41X= MakeH2I("FRS/SCI_TPC/","SCIdE21L_TPC41X", 1024,0,4096, 400,-100.,100, "SC21L dE", "TPC41 X[mm]", 2);
        hSCIdE21R_TPC42X= MakeH2I("FRS/SCI_TPC/","SCIdE21R_TPC42X", 1024,0,4096, 400,-100.,100, "SC21R dE", "TPC42 X[mm]", 2);
        hSCIdE21R_TPC41X= MakeH2I("FRS/SCI_TPC/","SCIdE21R_TPC41X", 1024,0,4096, 400,-100.,100, "SC21R dE", "TPC41 X[mm]", 2);

        // ToF SC21-SC42 changed on 03.07.2018 SB
        sprintf(fname,"FRS/SCI/TOF/TOF(%d)",3);
        sprintf(name,"SCI_21_42_TofLL");
        hSCI_TofLL3 = MakeH1I(fname,name,1500,0,62000,"TAC SC42L-SC21L [ps]",2,3);

        sprintf(name,"SCI_21_42_TofRR");
        hSCI_TofRR3 = MakeH1I(fname,name,1500,0,62000,"TAC SC42R-SC21R [ps]",2,3);

        sprintf(name,"SCI_21_42_Tof3");
        hSCI_Tof3 = MakeH1I(fname,name,1000,0,62000,"TAC SC42-SC21 [ps] (pos.corr.)",2,3);

    hSCI_dT_21l_41l = MakeTH1('D',"FRS/SCI/dT/SCI_dt_21l_41l","hSCI_dT_21l_41l",5001,0,5000); //from Multihit TDCS
    hSCI_dT_21r_41r = MakeTH1('D',"FRS/SCI/dT/SCI_dt_21r_41r","hSCI_dT_21r_41r",5001,0,5000);

    hSCI_dT_21l_42l = MakeTH1('D',"FRS/SCI/dT/SCI_dt_21l_42l","hSCI_dT_21l_42l",5001,0,5000);
    hSCI_dT_21r_42r = MakeTH1('D',"FRS/SCI/dT/SCI_dt_21r_42r","hSCI_dT_21r_42r",5001,0,5000);

    //ID
    hID_AoQ = MakeH1I("FRS/ID","ID_AoQ",2000,1.4,5.0,"A/Q S2-S4",2,6);
    hID_AoQ_corr = MakeH1I("FRS/ID","ID_AoQ_corr",2000,1.4,3.0,"A/Q S2-S4",2,6);
  //   hID_Z = MakeH1I("ID",Form("ID_Z, gain=%f",music->e1_gain[0]),1000,10,93,"Z s2-s4",2,6);
    hID_Z = MakeH1I("FRS/ID","ID_Z",1000,0,93,"Z s2-s4",2,6);
    hID_Z2 = MakeH1I("FRS/ID","ID_Z2",1000,0,93,"Z2 s2-s4",2,6);
   // hID_Z3 = MakeH1I("FRS/ID","ID_Z3",1000,10,93,"Z3 s2-s4",2,6);
    ////////////////////////////////////////////////////////////

    hID_Z_dE2 = MakeH2I("FRS/ID","ID_Z_dE2", 250,1,30, 250,0.,4000.,
              "Z", "MUSIC2_dE", 2);

    hID_Z_Sc21E = MakeH2I("FRS/ID","ID_Z_Sc21E", 300,0,25.,400,0,4000.,
            "Z s2-s4", "sqrt(Sc21_L*sC21_R)", 2);
    hID_x2z = MakeH2I("FRS/ID","ID_x2z", 300,1.,30., 200,-100.,100., "Z s2-s4", "X at S2 [mm]", 2);
    hID_x4z = MakeH2I("FRS/ID","ID_x4z", 300,1.,30., 200,-100.,100., "Z s2-s4", "X at S4 [mm]", 2);
    hID_E_Xs4 = MakeH2I("FRS/ID","ID_E_Xs4", 200,-100.,100., 400,0.,4000., "X s4 [mm]", "Delta E", 2);
    hID_E_Xs2 = MakeH2I("FRS/ID","ID_E_Xs2", 200,-100.,100., 400,0.,4000., "X s2 [mm]", "Delta E", 2);
    hID_x2a2 = MakeH2I("FRS/ID", "ID_x2_a2", 200, -100., 100., 200, -100., 100., "X s2 [mm]", "AngleX s2 [mrad]", 2);
    hID_y2b2 = MakeH2I("FRS/ID", "ID_y2_b2", 200, -100., 100., 200, -100., 100., "Y s2 [mm]", "AngleY s2 [mrad]", 2);
    hID_x4a4 = MakeH2I("FRS/ID", "ID_x4_a4", 200, -100., 100., 200, -100., 100., "X s4 [mm]", "AngleX s4 [mrad]", 2);
    hID_y4b4 = MakeH2I("FRS/ID", "ID_y4_b4", 200, -100., 100., 200, -100., 100., "Y s4 [mm]", "AngleY s4 [mrad]", 2);
    hID_x2x4 = MakeH2I("FRS/ID","ID_x2_x4",200,-100,100,200,-100,100,"x2 mm","x4 mm",2);
    hID_SC41dE_AoQ = MakeH2I("FRS/ID","ID_SC41dE_AoQ", 300,1.2,3.0, 800,0.,4000.,"A/Q s2-s4", "SC41 dE", 2);
  ///////////////////////////////////////////////////////////////////////////////////////////////
for(int i=0;i<7;i++)
    {
      char fname[100];
      char name[100];
      sprintf(fname,"FRS/TPC/%s/",tpc_folder_ext1[i]);

      hTPC_X[i]=MakeH1I_TPC(fname,"X",i,800,-100.,100.,
                "x[mm]",2,3);
      hTPC_Y[i]=MakeH1I_TPC(fname,"Y",i,800,-100.,100.,
                "y[mm]",2,3);


      sprintf(name,"%s%s",tpc_name_ext1[i],"XY");
      hcTPC_XY[i]=MakeH2I(fname,name, 120,-120.,120., 120,-120.,120.,
              "X [mm] ","Y [mm] ", 2);

      sprintf(name,"%s%s",tpc_name_ext1[i],"LTRT");
      hTPC_LTRT[i]=MakeH2I(fname,name, 2048,0,4095, 2048,0,4095,
               "LT [ch]","RT[ch] ", 2);
      hTPC_DELTAX[i]=MakeH1I_TPC(fname,"x0-x1",i,100,-10.,10.,
                 "x0-x1[mm]",2,3);

    }
    hID_x2 = MakeTH1('D',"FRS/TPC/S2_X","ID_x2",3000,0,5000);
    hID_y2 = MakeTH1('D',"FRS/TPC/S2_Y","ID_y2",3000,0,5000);
    hID_a2 = MakeTH1('D',"FRS/TPC/S2_angA","ID_a2",3000,0,5000);
    hID_b2 = MakeTH1('D',"FRS/TPC/S2_angB","ID_b2",3000,0,5000);

    hID_x4 = MakeTH1('D',"FRS/TPC/S4_X","ID_x4",800,-100,100);
    hID_y4 = MakeTH1('D',"FRS/TPC/S4_Y","ID_y4",3000,0,5000);
    hID_a4 = MakeTH1('D',"FRS/TPC/S4_angA","ID_a4",3000,0,5000);
    hID_b4 = MakeTH1('D',"FRS/TPC/S4_angA","ID_b4",3000,0,5000);

    htpc_X2 = MakeTH1('D',"FRS/TPC/S2_TPCX","tpc_x S21",800,-100,100);
    htpc_Y2 = MakeTH1('D',"FRS/TPC/S2_TPCY","tpc_y S21",800,-100,100);
    htpc_X4 = MakeTH1('D',"FRS/TPC/S4_TPCX","tpc_x S41",800,-100,100);
    htpc_Y4 = MakeTH1('D',"FRS/TPC/S4_TPCY","tpc_y S41",800,-100,100);

    hID_beta = MakeH1I("FRS/ID","ID_beta",2000,0,2000,"id.beta(2)*1000",2,6);

    hID_dEToF = MakeH2I("FRS/ID","ID_dEToF", 2000, 00000.,70000.,400,0,4000, "tof S2-S4 Sci.Tof(2)", "Music_dE(1)", 2);
    hID_BRho[0] = MakeH1I("FRS/ID","ID_BRho0",5000,2.5,14.5,"BRho of 1. Stage [Tm]",2,6);
    hID_BRho[1] = MakeH1I("FRS/ID","ID_BRho1",5000,2.5,14.5,"BRho of 2. Stage [Tm]",2,6);

  // char name[80], xtitle[80];
   for(int i=0;i<8;i++)
     {
       hMUSIC1_E[i] = MakeTH1('D', Form("FRS/MUSIC/MUSIC(1)/Energy/EnergyM1%2d",i), Form("Music 1 E%2d",i), 4096,0,4096);
       hMUSIC2_E[i] = MakeTH1('D', Form("FRS/MUSIC/MUSIC(2)/Energy/EnergyM2%2d",i), Form("Music 2 E%2d",i), 4096,0,4096);
       hMUSIC1_T[i] = MakeTH1('D', Form("FRS/MUSIC/MUSIC(1)/Time/TimeM1%2d",i), Form("Music 1 T%2d",i), 4096,0,4096);
       hMUSIC2_T[i] = MakeTH1('D', Form("FRS/MUSIC/MUSIC(2)/Time/TimeM2%2d",i), Form("Music 2 T%2d",i), 4096,0,4096);
     }
 
  //  hMUSIC1_dE1dE2 = MakeTH2('D',"FRS/MUSIC/MUSIC(1)/E1E2","E1_E2", 1024,0,4096,1024,0,4096);
   
  //  hMUSIC1_MUSIC2 = MakeTH2('D',"FRS/MUSIC/MUSIC1_MUSIC2","dE1_dE2", 1024,0,4096,1024,0,4096);  
   
    htimestamp = MakeTH1('D',"FRS/timestamp","timestamp",30,0.,300.);
    hts = MakeTH1('D',"FRS/ts","ts",30,0.,300.);
    hts2 = MakeTH1('D',"FRS/ts2","ts2",30,0.,300.);

  }
  //-----------------------------------------------------------------------------------------------------------------------------//
  void EventUnpackProc::Fill_FRS_Histos(){

     time_in_ms = 0;
     spill_count = 0;
     ibin_for_s = 0;
     ibin_for_100ms = 0;
     ibin_for_spill = 0;


    for(int i =0; i<8; i++){
        Music_E1[i] = 0;
        Music_E2[i] = 0;
        Music_T1[i] = 0;
        Music_T2[i] = 0;
    }

    for(int i =0; i<3; i++){
        Music_dE[i] = RAW->get_FRS_MusicdE(i);
        Music_dE_corr[i] = RAW->get_FRS_MusicdE_corr(i);
    }
    for(int i=0; i<8; i++){
        Music_E1[i] = RAW->get_FRS_MusicE1(i);
        Music_E2[i] = RAW->get_FRS_MusicE2(i);
        Music_T1[i] = RAW->get_FRS_MusicT1(i);
        Music_T2[i] = RAW->get_FRS_MusicT2(i);

    }

    for(int l=0;l<12;++l){
        sci_l[l] = RAW->get_FRS_sci_l(l);
        sci_r[l] = RAW->get_FRS_sci_r(l);
        sci_e[l] = RAW->get_FRS_sci_e(l);
        sci_tx[l] = RAW->get_FRS_sci_tx(l);
        sci_x[l] = RAW->get_FRS_sci_x(l);
    }
    sci_tofll2 = RAW->get_FRS_tofll2();
    sci_tofll3 = RAW->get_FRS_tofll3();
    sci_tof2 = RAW->get_FRS_tof2();
    sci_tofrr2 = RAW->get_FRS_tofrr2();
    sci_tofrr3 = RAW->get_FRS_tofrr3();
    sci_tof3 = RAW->get_FRS_tof3();

    ID_x2 = RAW->get_FRS_x2();
    ID_y2 = RAW->get_FRS_y2();
    ID_a2 = RAW->get_FRS_a2();
    ID_b2 = RAW->get_FRS_b2();

    ID_x4 = RAW->get_FRS_x4();
    ID_y4 = RAW->get_FRS_y4();
    ID_a4 = RAW->get_FRS_a4();
    ID_b4 = RAW->get_FRS_b4();

    for(int i =0; i<7; i++){
    TPC_X[i] = RAW-> get_FRS_tpcX(i);
    TPC_Y[i] = RAW-> get_FRS_tpcY(i);
    for(int j=0; j<2; j++){
    TPC_LT[i][j] = RAW->get_FRS_tpclt(i,j);
    TPC_RT[i][j] = RAW->get_FRS_tpcrt(i,j);
        }

    }
    TPC_X0 = RAW->get_FRS_tpcx0();
    TPC_X1 = RAW->get_FRS_tpcx0();

    sci_dt_21l_21r = RAW->get_FRS_dt_21l_21r();
    sci_dt_41l_41r = RAW->get_FRS_dt_41l_41r();
    sci_dt_42l_42r = RAW->get_FRS_dt_42l_42r();
    sci_dt_43l_43r = RAW->get_FRS_dt_43l_43r();

    sci_dt_21l_41l = RAW->get_FRS_dt_21l_41l();
    sci_dt_21r_41r = RAW->get_FRS_dt_21r_41r();

    sci_dt_21l_42l = RAW->get_FRS_dt_21l_42l();
    sci_dt_21r_42r = RAW->get_FRS_dt_21r_42r();

    for(int k =0; k<2; ++k){
        ID_brho[k] = RAW->get_FRS_brho(k);
        ID_rho = RAW->get_FRS_rho(k);
    }
    beta = RAW->get_FRS_beta();
    beta3 = RAW->get_FRS_beta3();
    gamma = RAW->get_FRS_gamma();

    AoQ = RAW->get_FRS_AoQ();
    AoQ_corr = RAW->get_FRS_AoQ_corr();

    ID_z = RAW->get_FRS_z();
    ID_z2 = RAW->get_FRS_z2();
    ID_z3 = RAW->get_FRS_z3();
   

    timestamp = RAW->get_FRS_timestamp();
    ts = RAW->get_FRS_ts(); //Spill time structrue
    ts2 = RAW->get_FRS_ts2();
    
     /// --------FRS SCALARS-------------------------------- //

    
    time_in_ms           = RAW->get_FRS_time_in_ms();
    spill_count          = RAW->get_FRS_spill_count();
    ibin_for_s           = RAW->get_FRS_ibin_for_s();
    ibin_for_100ms       = RAW->get_FRS_ibin_for_100ms();
    ibin_for_spill       = RAW->get_FRS_ibin_for_spill();
    ibin_clean_for_s     = RAW->get_FRS_ibin_clean_for_s();
    ibin_clean_for_100ms = RAW->get_FRS_ibin_clean_for_100ms();
    ibin_clean_for_spill = RAW->get_FRS_ibin_clean_for_spill();
    increase_scaler_temp = RAW->get_FRS_increase_scaler_temp();
    

    /// ------------MUSIC---------------------------- //

    
     //MUSIC 1 is TUM MUSIC (8 anodes). MUSIC 3 not required
    for(int i=0; i<3; i++){
   //  hMUSIC1_MUSIC2->Fill(Music_dE[0],Music_dE[1]);
   //  hMUSIC1_dE1dE2->Fill(Music_E1[0],Music_E1[1]);
    }
    for(int i=0; i<8; i++){
   if(Music_E1[i]!=0) hMUSIC1_E[i]->Fill(Music_E1[i]);
   if(Music_E2[i]!=0) hMUSIC2_E[i]->Fill(Music_E2[i]);
   if(Music_T1[i]!=0) hMUSIC1_T[i]->Fill(Music_T1[i]);
   if(Music_T2[i]!=0) hMUSIC2_T[i]->Fill(Music_T2[i]);
 }
    //SCI
    for (int cnt=0;cnt<6;cnt++) //
     {
       int idx = 0 ;
       //int mw_idx = 0;
       //Float_t mwx = 0;
       switch(cnt)
     {
     case 0:        /* SC21 */
       idx = 2;
       //mw_idx = 2;
       //mwx = clb.sc21_x;
       break;
     case 1:        /* SC21 delayed */
       idx = 3;
       //mw_idx = 2;
       //mwx = clb.sc21_x;
       break;
     case 2:        /* SC41 */
       idx = 5;
       //mw_idx = 5;
       //mwx = clb.tpc_sc41_x;
       break;
     case 3:        /* SC42 */
           idx = 6;
       break;
     case 4:
       idx = 7;     /* SC43 */
       break;
     case 5:
       idx = 10;    /* SC81 */
       break;
     default: idx = 2;
     }
    /*-------------------------------------------------------------------------*/
    /* focus index: detector number                  tof index  tof path       */
    /*       0:     Sc01                                0:     TA - S1         */
    /*       1:     Sc11                                1:     S1 - S2         */
    /*       2:     Sc21                                2:     S2 - S41        */
    /*       3:     Sc21                                3:     S2 - S42        */
    /*       4:     Sc31                                4:     S2 - 81         */
    /*       5:     Sc41                                5:     S2 - E1         */
    /*                                                                         */
    /*       6:     Sc42                              tof index not used up to */
    /*       7:     Sc43 (previously Sc51)             now, only separate      */
    /*       8:     Sc61                              variables for S2-S41 and */
    /*       9:     ScE1 (ESR)                                S2-S42           */
    /*      10:     Sc81                                                       */
    /*      11:     Sc82                                                       */
    /*-------------------------------------------------------------------------*/
      if(sci_l[idx]!=0)   hSCI_L[idx]->Fill(sci_l[idx]);
      if(sci_r[idx]!=0)   hSCI_R[idx]->Fill(sci_r[idx]);
      if(sci_e[idx]!=0)   hSCI_E[idx]->Fill(sci_e[idx]);
      if(sci_tx[idx]!=0)   hSCI_Tx[idx]->Fill(sci_tx[idx]);
      if(sci_x[idx]!=0)   hSCI_X[idx]->Fill(sci_x[idx]);
       
     if(sci_e[2]!=0 && sci_e[5]!=0)    hSCI_dE24->Fill(sci_e[2],sci_e[5]);
     }
     if(sci_e[5]!=0 && TPC_X[5]!=0) hSCIdE41_TPC42X->Fill(sci_e[5],TPC_X[5]); //dE_SCI_41 vs TPC_42X
     if(sci_l[5]!=0 &&TPC_X[5]!=0) hSCIdE41L_TPC42X->Fill(sci_l[5],TPC_X[5]); //dE_SCI_41L vs TPC_42X
     if(sci_l[5]!=0 &&TPC_X[4]!=0) hSCIdE41L_TPC41X->Fill(sci_l[5],TPC_X[4]); //dE_SCI_41L vs TPC_41X
     if(sci_r[5]!=0 &&TPC_X[5]!=0) hSCIdE41R_TPC42X->Fill(sci_r[5],TPC_X[5]); //dE_SCI_41R vs TPC_42X
     if(sci_r[5]!=0 &&TPC_X[4]!=0) hSCIdE41R_TPC41X->Fill(sci_r[5],TPC_X[4]); //dE_SCI_41R vs TPC_41X
      //cout<<"sci_e[0] " << sci_e[0] <<endl;
     if(sci_e[0]!=0 &&TPC_X[5]!=0) hSCIdE21_TPC42X->Fill(sci_e[0],TPC_X[5]); //dE_SCI_21 vs TPC_42X
     if(sci_l[0]!=0 &&TPC_X[5]!=0) hSCIdE21L_TPC42X->Fill(sci_l[0],TPC_X[5]); //dE_SCI_21L vs TPC_42X
     if(sci_l[0]!=0 &&TPC_X[4]!=0) hSCIdE21L_TPC41X->Fill(sci_l[0],TPC_X[4]); //dE_SCI_21L vs TPC_41X
     if(sci_l[0]!=0 &&TPC_X[4]!=0) hSCIdE21R_TPC42X->Fill(sci_l[0],TPC_X[4]); //dE_SCI_21R vs TPC_42X
     if(sci_r[0]!=0 &&TPC_X[4]!=0) hSCIdE21R_TPC41X->Fill(sci_r[0],TPC_X[4]); //dE_SCI_21R vs TPC_41X


     if(sci_tofll2!=0)  hSCI_TofLL2->Fill(sci_tofll2);
     if(sci_tofll3!=0) hSCI_TofLL3->Fill(sci_tofll3);
  //   hSCI_Tof2->Fill(sci_tof2);
     if(sci_tofrr2!=0)   hSCI_TofRR2->Fill(sci_tofrr2);
     if(sci_tofrr3!=0)  hSCI_TofRR3->Fill(sci_tofrr3);
     if(sci_tof3!=0)  hSCI_Tof3->Fill(sci_tof3);


    if(ID_x2!=0) hID_x2->Fill(ID_x2);
    if(ID_y2!=0) hID_y2->Fill(ID_y2);
    if(ID_a2!=0) hID_a2->Fill(ID_a2);
    if(ID_b2!=0) hID_b2->Fill(ID_b2);

    if(ID_x4!=0) hID_x4->Fill(ID_x4);
    if(ID_y4!=0) hID_y4->Fill(ID_y4);
    if(ID_a4!=0) hID_a4->Fill(ID_a4);
    if(ID_b4!=0) hID_b4->Fill(ID_b4);

    for(int i=0; i<7; i++){
       if(TPC_X[i]!=0)   hTPC_X[i]->Fill(TPC_X[i]);
       if(TPC_Y[i]!=0)    hTPC_Y[i]->Fill(TPC_Y[i]);
       if(TPC_X[i]!=0 && TPC_Y[i]!=0)    hcTPC_XY[i]->Fill(TPC_X[i],TPC_Y[i]);
       if(TPC_LT[i][0]!=0 && TPC_RT[i][1]!=0)    hTPC_LTRT[i]->Fill(TPC_LT[i][0],TPC_RT[i][1]);
       if(TPC_X0-TPC_X1!=0)    hTPC_DELTAX[i]->Fill(TPC_X0-TPC_X1);
         // cout<<"TPC_LT[i][0] " << TPC_LT[i][0] << " i " << i << endl;
    }

//     htpc_X2->Fill(TPC_X[2]);
//     htpc_Y2->Fill(TPC_Y[2]);
//     htpc_X4->Fill(TPC_X[4]);
//     htpc_Y4->Fill(TPC_Y[4]);
    for(int i=0;i<2;i++){
      if(ID_brho[i]!=0)hID_BRho[i]->Fill(ID_brho[i]);
    }

    //SCI tx
    //if(sci_dt_21l_21r) hSCI_dT_21l_21r->Fill(sci_dt_21l_21r);
//     if(sci_dt_41l_41r) hSCI_dT_41l_41r->Fill(sci_dt_41l_41r);
//     if(sci_dt_42l_42r) hSCI_dT_42l_42r->Fill(sci_dt_42l_42r);

    if(sci_dt_21l_41l!=0) hSCI_dT_21l_41l->Fill(sci_dt_21l_41l);
    if(sci_dt_21r_41r!=0) hSCI_dT_21r_41r->Fill(sci_dt_21r_41r);

    if(sci_dt_21l_42l!=0) hSCI_dT_21l_42l->Fill(sci_dt_21l_42l);
    if(sci_dt_21r_42r!=0) hSCI_dT_21r_42r->Fill(sci_dt_21r_42r);



    if(beta!=0) hID_beta->Fill(beta);
   // if(beta3) hbeta3->Fill(beta3);

    //  if(gamma) hgamma->Fill(gamma);
    //cout<<"count" <<count <<endl;
    //if(count>10000 && count<40000)cout<<"event "<< event_number  <<" AoQ sanity check!! " << AoQ << endl;

    if(AoQ!=0) hID_AoQ->Fill(AoQ);
    if(AoQ_corr!=0) hID_AoQ_corr->Fill(AoQ_corr);

    /****  S4  (MUSIC 1)   */
     if(ID_z!=0)hID_Z->Fill(ID_z);
     /****  S4  (MUSIC 2)   */
     if(ID_z2!=0) hID_Z2->Fill(ID_z2);
     /****  S4  (MUSIC OLD)   */
     //hID_Z3->Fill(ID_z3);

//      hID_Z_Z2->Fill(ID_z,ID_z2);
     if(ID_z!=0 && Music_dE[1]!=0)hID_Z_dE2->Fill(ID_z,Music_dE[1]);
    // hID_Z_Z3->Fill(ID_z,ID_z3);
     if(ID_z!=0 && sci_l[2]!=0 && sci_r[2]!=0)hID_Z_Sc21E->Fill(ID_z, sqrt(sci_l[2]*sci_r[2]));

     if(ID_x2!=0&&ID_x4!=0 ) hID_x2x4->Fill(ID_x2, ID_x4);
     if(AoQ!=0 && sci_e[5]!=0) hID_SC41dE_AoQ->Fill(AoQ, sci_e[5]);

     if(sci_tof2!=0 && Music_dE[0]!=0) hID_dEToF->Fill(sci_tof2, Music_dE[0]);

     if(ID_z!=0 && ID_x2!=0) hID_x2z->Fill(ID_z, ID_x2);// MUSIC1
     if(ID_z!=0 && ID_x4!=0) hID_x4z->Fill(ID_z, ID_x4);// MUSIC1

     if(ID_x4!=0 && Music_dE[0]!=0) hID_E_Xs4->Fill(ID_x4,Music_dE[0]);
     if(ID_x4!=0 && Music_dE[0]!=0)hID_E_Xs2->Fill(ID_x2,Music_dE[0]);

     if(ID_x2!=0 && ID_a2!=0)hID_x2a2->Fill(ID_x2,ID_a2);
     if(ID_y2!=0 && ID_b2!=0)hID_y2b2->Fill(ID_y2,ID_b2);
     if(ID_x4!=0 && ID_a4!=0)hID_x4a4->Fill(ID_x4,ID_a4);
     if(ID_x4!=0 && ID_b4!=0) hID_y4b4->Fill(ID_y4,ID_b4);

    if(timestamp) htimestamp->Fill(timestamp);
    if(ts) hts->Fill(ts);
    if(ts2) hts2->Fill(ts2);
    
    for(int ii=0; ii<64; ii++){
    
    //    printf("ch %d: this event = %lld, increase =%lld\n",ii,src.sc_long[ii],increase_scaler_temp);
    hScaler_per_s[ii]->AddBinContent(ibin_for_s, increase_scaler_temp);
    hScaler_per_100ms[ii]->AddBinContent(ibin_for_100ms, increase_scaler_temp);
    hScaler_per_spill[ii]->AddBinContent(ibin_for_spill, increase_scaler_temp);    
   // if(ii=50)cout<<"ibin_clean_for_s " << ibin_clean_for_s << " increase_scaler_temp " << increase_scaler_temp<< endl;
  }
  
   for(int ii=0; ii<64; ii++){
    hScaler_per_s[ii]->SetBinContent(ibin_clean_for_s, 0);
    hScaler_per_100ms[ii]->SetBinContent(ibin_clean_for_100ms, 0);
    hScaler_per_spill[ii]->SetBinContent(ibin_clean_for_spill, 0);
  }

  }

  /**----------------------------------------------------------------------------------------------**/
  /**-------------------------------------------  AIDA   ------------------------------------------**/
  /**----------------------------------------------------------------------------------------------**/

  void EventUnpackProc::Make_AIDA_Histos(){

    TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
    hAIDA_ADC.resize(conf->FEEs());

    for (int i = 0; i < conf->FEEs(); i++)
    {
      for (int j = 0; j < 64; j++)
      {
        hAIDA_ADC[i][j][0] = MakeTH1('I',
          Form("AIDA/Unpacker/FEE%d/Fee%d_L_Channel%02d", i+1, i+1, j+1),
          Form("FEE %d Channel %2d (Low Energy)", i+1, j+1),
          2000, -32768, 32767
        );
      }
    }

    for (int i = 0; i < conf->FEEs(); i++)
    {
      for (int j = 0; j < 64; j++)
      {
        hAIDA_ADC[i][j][1] = MakeTH1('I',
          Form("AIDA/Unpacker/FEE%d/Fee%d_H_Channel%02d", i+1, i+1, j+1),
          Form("FEE %d Channel %2d (High Energy)", i+1, j+1),
          2000, -32768, 32767
        );
      }
    }
}

void EventUnpackProc::Fill_AIDA_Histos() {
  AIDA_Hits = RAW->get_AIDA_HITS();

  for(int i = 0; i<AIDA_Hits; i++) {
    int fee = RAW-> get_AIDA_FEE_ID(i);
    int chan = RAW-> get_AIDA_CHA_ID(i);
    int adc = RAW->get_AIDA_ADC(i);
    int veto = RAW->get_AIDA_HighE_VETO(i) ? 1 : 0;

    hAIDA_ADC[fee][chan][veto]->Fill(adc - 32767);

    //cout<<"chan " << chan << endl;
  }
}


  /**----------------------------------------------------------------------------------------------**/
  /**---------------------------------------  bPLASTIC TAMEX  -------------------------------------**/
  /**----------------------------------------------------------------------------------------------**/
//   void EventUnpackProc::Make_Plastic_Histos(){
//
//     TOT_TOT = new TH1***[100];
//     TOT_Single = new TH1**[100];
//     TRAIL_TRAIL = new TH1***[100];
//     LEAD_LEAD = new TH1***[100];
//
//     for(int i = 0;i < 100;++i){
//       TOT_Single[i] = new TH1*[100];
//       TOT_TOT[i] = new TH1**[100];
//       TRAIL_TRAIL[i] = new TH1**[100];
//       LEAD_LEAD[i] = new TH1**[100];
//
//       for(int j = 0;j < 100;++j){
//         TOT_TOT[i][j] = new TH1*[100];
//         TRAIL_TRAIL[i][j] = new TH1*[100];
//         LEAD_LEAD[i][j] = new TH1*[100];
//
//         for(int k = 0;k < 100;++k){
//           TOT_TOT[i][j][k] = nullptr;
//           TRAIL_TRAIL[i][j][k] = nullptr;
//           LEAD_LEAD[i][j][k]   = nullptr;
//         }
//
//         TOT_Single[i][j] = nullptr;
//
//       }
//     }
//   }
  //-----------------------------------------------------------------------------------------------------------------------------//

  

  /**----------------------------------------------------------------------------------------------**/
  /**--------------------------------------  bPLASTIC VME (+Scalar)  ----------------------------------------**/
  /**----------------------------------------------------------------------------------------------**/
//   void EventUnpackProc::Make_Plastic_VME_Histos(){
//     for (int i = 0; i<32; i++){
// 
//   hPLAS_QDCRaw1[i] =  MakeTH1('D', Form("bPlastic/Energy/Raw/QDC1Raw/QDC1_Ch.%2d",i), Form("QDC1 Ch. %2d",i), 2500, 0., 5000.);
//   hPLAS_QDCRaw2[i] =  MakeTH1('D', Form("bPlastic/Energy/Raw/QDC2Raw/QDC2_Ch.%2d",i), Form("QDC2 Ch. %2d",i), 2500, 0., 5000.);
//   hPLAS_TDCRaw[i] =  MakeTH1('D', Form("bPlastic/Timing/Raw/TDCRaw/TDC_Ch.%2d",i), Form("TDC Ch. %2d",i), 1000, 0, 10000);
//     }
//    
// }


  //-----------------------------------------------------------------------------------------------------------------------------//
  void EventUnpackProc::Fill_Plastic_VME_Histos(){

    double QDC1[32], QDC2[32];
    double bplasTDC[32];
    int bPlas_QCDID;
    int bplasTDCID, bPlasIT;

    int bPlasmulti[32];
    
   // double Scalar_Data[32];

    //Reset arrays and variables://
//     Scalar_Chan = 0;
    bplasTDCID=0;
    bPlasIT = 0;

    for (int i =0; i<32; i++){
      QDC1[i] = 0;
      QDC2[i] = 0;
      bplasTDC[i] = 0;
    //  Scalar_Data[i] = 0;

    }

    /**------------------PLASTIC Energy -----------------------------------------**/

    for (int jj = 0; jj < RAW->get_plastic_VME_QDC_fired();  jj++){
      bPlas_QCDID = RAW->get_plastic_VME_QDC_cha(jj);
      QDC1[bPlas_QCDID] = RAW->get_plastic_VME_QDC_dat1(jj);
      hPLAS_QDCRaw1[bPlas_QCDID]->Fill(QDC1[bPlas_QCDID]);
      hPLAS_QDCRaw2[bPlas_QCDID]->Fill(QDC2[bPlas_QCDID]);


    }
    /**------------------PLASTIC Time -----------------------------------------**/

    bPlasIT = RAW->get_plastic_VME_TDC_fired(); //Fired TDC channels

    for (int i=0; i<bPlasIT; i++){

      bplasTDCID = RAW->get_plastic_VME_TDC_cha(i); //Plastic TDC ID
      bplasTDC[bplasTDCID] = RAW->get_plastic_VME_TDC_dat(bplasTDCID)*0.025;//1ns //Plastic TDC Data

      bPlasmulti[bplasTDCID]++;


    if(bplasTDC[bplasTDCID]>0){

            hPLAS_TDCRaw[bplasTDCID]->Fill(bplasTDC[bplasTDCID]);  //1ns

             }
            }
              //Scalar
//         Scalar_iterator = RAW->get_scalar_iterator();
//         for (int g=0; g<Scalar_iterator; g++){
//             Scalar_Chan = RAW->get_scalar_chan(g);
//             hScalar_hit_pattern->Fill(Scalar_Chan);
//             }
         }

/**----------------------------------------------------------------------------------------------**/
/**-----------------------------------------  FATIMA TAMEX  ------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
 void EventUnpackProc::Make_FATIMA_TAMEX_Histos(){


        for(int i=0;i<FAT_MAX_DET;i++){
             hFATlead_Coarse[i]= MakeTH1('D', Form("FATIMA/Lead_Coarse/Lead-CoarseCh.%02d", i), Form("Lead Coarse %2d", i), 50000, -100000., 100000.);
             hFATlead_Fine[i]= MakeTH1('D', Form("FATIMA/Lead_Fine/Lead-FineCh.%02d", i), Form("Lead Fine %2d", i), 601, -1., 60.);
             hFATtrail_Coarse[i]= MakeTH1('D', Form("FATIMA/Trail_Coarse/Trail-CoarseCh.%02d", i), Form("Trail Coarse %2d", i), 50000, -100000., 100000.);
             hFATtrail_Fine[i]= MakeTH1('D', Form("FATIMA/Trail_Fine/Trail-FineCh.%02d", i), Form("Trail Fine %2d", i), 601, -1., 60.);


            }
        }
  //-----------------------------------------------------------------------------------------------------------------------------//

  void EventUnpackProc::Fill_FATIMA_TAMEX_Histos(){

      ///TAMEX
    //get amount of fired Tamex modules
    int TamexHits_Fatima = RAW->get_FATIMA_tamex_hits();

    int Lead_Fatima_Coarse[4][32];
    double Lead_Fatima_Fine[4][32];
    int Trail_Fatima_Coarse[4][32];
    double Trail_Fatima_Fine[4][32];
    int Phys_Channel_Fatima[32];
    int leadHits_Fatima = 0,leadHitsCh_Fatima = 0;
    int trailHits_Fatima = 0,trailHitsCh_Fatima = 0;


    for(int i=0; i<32; i++){
     Phys_Channel_Fatima[i] = 0;
    }

    for(int i =0; i<4; i++){
        for(int j=0; j<32;j++){
     Lead_Fatima_Coarse[i][j] = 0;
     Lead_Fatima_Fine[i][j] = 0;
     Trail_Fatima_Coarse[i][j] = 0;
     Trail_Fatima_Fine[i][j] = 0;
        }

    for(int i = 0;i < TamexHits_Fatima;++i){

      leadHits_Fatima = RAW->get_FATIMA_lead_hits(i);
      trailHits_Fatima = RAW->get_FATIMA_trail_hits(i);


        //Box diagrams for leading and trailing
      for(int j = 0;j < RAW->get_FATIMA_am_Fired(i);j++){
         if(RAW->get_FATIMA_CH_ID(i,j) % 2 == 1){
            Phys_Channel_Fatima[j] = RAW->get_FATIMA_physical_channel(i,j);
          //  Lead_Fatima[i][j] = RAW->get_FATIMA_lead_T(i,j);
            Lead_Fatima_Coarse[i][j] = RAW->get_FATIMA_coarse_lead(i,j);
            Lead_Fatima_Fine[i][j] = RAW->get_FATIMA_fine_lead(i,j);
            hFATlead_Coarse[Phys_Channel_Fatima[j]]->Fill(Lead_Fatima_Coarse[i][j]);
            hFATlead_Fine[Phys_Channel_Fatima[j]]->Fill(Lead_Fatima_Fine[i][j]);

        }
         if(RAW->get_FATIMA_CH_ID(i,j) % 2 == 0){
            Phys_Channel_Fatima[j] = RAW->get_FATIMA_physical_channel(i,j);
          //  Trail_Fatima[Phys_Channel_Fatima[j]] = RAW->get_FATIMA_trail_T(i,j);
            Trail_Fatima_Coarse[i][j] = RAW->get_FATIMA_coarse_trail(i,j);
            Trail_Fatima_Fine[i][j] = RAW->get_FATIMA_fine_trail(i,j);
            hFATtrail_Coarse[Phys_Channel_Fatima[j]]->Fill(Trail_Fatima_Coarse[i][j]);
            hFATtrail_Fine[Phys_Channel_Fatima[j]]->Fill(Trail_Fatima_Fine[i][j]);
                }
            }
        }
    }

  }
  /**----------------------------------------------------------------------------------------------**/
/**-----------------------------------------  FATIMA VME  ------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
void EventUnpackProc::Make_FATIMA_Histos(){


  for (int i=0; i<FAT_MAX_DET; i++){
    hFAT_Eraw_VME[i] = MakeTH1('D', Form("FATIMA_VME/Unpacker/Energy/Raw/E_Raw_LaBrCh. %02d", i),
    Form("LaBr%02d energy (raw)", i),2000,0,40000);

    hFAT_Traw[i] = MakeTH1('D', Form("FATIMA_VME/Unpacker/Timing/Raw/Traw_LaBrCh. %02d", i),
                    Form("LaBr%02d energy", i),5000,0,4E6);

            }
    hScalar_hit_pattern = MakeTH1('D',"Scalar/HitPat","Scalar Hit pattern",32,0,32);
        }
//-----------------------------------------------------------------------------------------------------------------------------//

void EventUnpackProc::Fill_FATIMA_Histos(){
  double FAT_E[50],FAT_T[50];
  double En_i;
  int detQDC, detTDC;
  detTDC = 0;
  detQDC = 0;
  for (int k=0; k<50; k++){
    FAT_T[k] = 0;
    FAT_E[k] = 0;
  }
  int Scalar_iterator = 0;
  int Scalar_Chan = 0;


  /**------------------FATIMA Energy -----------------------------------------**/
  for (int i=0; i<RAW->get_FAT_QDCs_fired(); i++){ /** Loops over only channels in the QDC **/

    detQDC = RAW->get_FAT_QDC_id(i); /**FAT ID QDC*/
    En_i = RAW->get_FAT_QLong_Raw(i); /**Raw FAT Energy*/
  if (detQDC<40){
   hFAT_Eraw_VME[detQDC]->Fill(En_i);
    //cout << " EN " << En_i<<endl;
  }

  }
  /**------------------FATIMA TIMING -----------------------------------------**/
  for (int i=0; i<RAW->get_FAT_TDCs_fired(); i++){ /** Loops over only channels in the TDC 1-4 **/

   
    detTDC = (RAW->get_FAT_TDC_id(i));
     if(detTDC<50){
        FAT_T[detTDC] = (RAW->get_FAT_TDC_timestamp(i));
        hFAT_Traw[detTDC]->Fill(FAT_T[detTDC]*25); //in ps
    
                }
         }
  /**------------------Scaler TIMING -----------------------------------------**/  
         Scalar_iterator = RAW->get_scalar_iterator();
         
          for (int g=0; g<Scalar_iterator; g++){
              if(RAW->get_scalar_data(g)>0){
                hScalar_hit_pattern->Fill(g);
           }
        }
    }


/**----------------------------------------------------------------------------------------------**/
/**-----------------------------------------  FATIMA VME AND TAMEX  ------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
 void EventUnpackProc::Make_FATIMA_VME_TAMEX_Histos(){


        for(int i=0;i<40;i++){
             hFATlead_Coarse[i]= MakeTH1('D', Form("FATIMA/Unpacker/Lead_Coarse/Lead-CoarseCh.%02d", i), Form("Lead Coarse %2d", i), 50000, -100000., 100000.);
             hFATlead_Fine[i]= MakeTH1('D', Form("FATIMA/Unpacker/Lead_Fine/Lead-FineCh.%02d", i), Form("Lead Fine %2d", i), 601, -1., 60.);
             hFATtrail_Coarse[i]= MakeTH1('D', Form("FATIMA/Unpacker/Trail_Coarse/Trail-CoarseCh.%02d", i), Form("Trail Coarse %2d", i), 50000, -100000., 100000.);
             hFATtrail_Fine[i]= MakeTH1('D', Form("FATIMA/Unpacker/Trail_Fine/Trail-FineCh.%02d", i), Form("Trail Fine %2d", i), 601, -1., 60.);


            hFAT_Eraw[i] = MakeTH1('D', Form("FATIMA/Unpacker/Energy/Raw/E_Raw_LaBr%02d", i),
            Form("LaBr%02d energy (raw)", i),2000,0,40000);



            }
         
        }
  //-----------------------------------------------------------------------------------------------------------------------------//

  void EventUnpackProc::Fill_FATIMA_VME_TAMEX_Histos(){

    //get amount of fired Tamex modules
    int TamexHits_Fatima = RAW->get_FATIMA_tamex_hits();

    int Lead_Fatima_Coarse[4][32];
    double Lead_Fatima_Fine[4][32];
    int Trail_Fatima_Coarse[4][32];
    double Trail_Fatima_Fine[4][32];
    int Phys_Channel_Fatima[32];
    int leadHits_Fatima = 0,leadHitsCh_Fatima = 0;
    int trailHits_Fatima = 0,trailHitsCh_Fatima = 0;

    for(int i=0; i<32; i++){
     Phys_Channel_Fatima[i] = 0;
    }

    for(int i =0; i<4; i++){
        for(int j=0; j<32;j++){
     Lead_Fatima_Coarse[i][j] = 0;
     Lead_Fatima_Fine[i][j] = 0;
     Trail_Fatima_Coarse[i][j] = 0;
     Trail_Fatima_Fine[i][j] = 0;
        }
    }
        ///QDC
    double FAT_E[50];
    double En_i;
    int detQDC;
    detQDC = 0;
    for (int k=0; k<50; k++){

        FAT_E[k] = 0;
  }
 //cout<<"TamexHits_Fatima "<<TamexHits_Fatima<<endl;
    for(int i = 0;i < TamexHits_Fatima;++i){

      leadHits_Fatima = RAW->get_FATIMA_lead_hits(i);
      trailHits_Fatima = RAW->get_FATIMA_trail_hits(i);


        //Box diagrams for leading and trailing
      for(int j = 0;j < RAW->get_FATIMA_am_Fired(i);j++){
         if(RAW->get_FATIMA_CH_ID(i,j) % 2 == 1){
            Phys_Channel_Fatima[j] = RAW->get_FATIMA_physical_channel(i,j);
          //  Lead_Fatima[i][j] = RAW->get_FATIMA_lead_T(i,j);
            Lead_Fatima_Coarse[i][j] = RAW->get_FATIMA_coarse_lead(i,j);
            Lead_Fatima_Fine[i][j] = RAW->get_FATIMA_fine_lead(i,j);
            hFATlead_Coarse[Phys_Channel_Fatima[j]]->Fill(Lead_Fatima_Coarse[i][j]);
            hFATlead_Fine[Phys_Channel_Fatima[j]]->Fill(Lead_Fatima_Fine[i][j]);

        }
         if(RAW->get_FATIMA_CH_ID(i,j) % 2 == 0){
            Phys_Channel_Fatima[j] = RAW->get_FATIMA_physical_channel(i,j);
          //  Trail_Fatima[Phys_Channel_Fatima[j]] = RAW->get_FATIMA_trail_T(i,j);
            Trail_Fatima_Coarse[i][j] = RAW->get_FATIMA_coarse_trail(i,j);
            Trail_Fatima_Fine[i][j] = RAW->get_FATIMA_fine_trail(i,j);
            hFATtrail_Coarse[Phys_Channel_Fatima[j]]->Fill(Trail_Fatima_Coarse[i][j]);
            hFATtrail_Fine[Phys_Channel_Fatima[j]]->Fill(Trail_Fatima_Fine[i][j]);
                }
            }
        }

    /**------------------FATIMA QDC Energy -----------------------------------------**/

  for (int i=0; i<RAW->get_FAT_QDCs_fired(); i++){ /** Loops over only channels in the QDC **/


   detQDC = RAW->get_FAT_QDC_id(i); /**FAT ID QDC*/
   En_i = RAW->get_FAT_QLong_Raw(i); /**Raw FAT Energy*/
   if(detQDC<40){

  // cout<<"detQDC " << detQDC << " En_i " << En_i << " i " << i << endl;

  // hFAT_Eraw[detQDC]->Fill(En_i);

      }
    }
  }
/**----------------------------------------------------------------------------------------------**/
/**----------------------------------------   GALILEO   -----------------------------------------**/
/**----------------------------------------------------------------------------------------------**/


void EventUnpackProc::Make_GALILEO_Histos(){
  for (int j; j<GALILEO_MAX_HITS; j++){
        hGAL_Raw_E[j] = MakeTH1('D',Form("GALILEO/Raw/GALILEO_Energy_Spectra/GALILEO_Raw_E%2d",j),
                            Form("GALILEO Channel Energy Channel Raw %2d",j),20000,0,20000);

                    }
                }
//-----------------------------------------------------------------------------------------------------------------------------//
void EventUnpackProc::Fill_GALILEO_Histos(){

    double tmpGAL[32];
    int  GALILEO_hits, GalID;

     /**------------------GALILEO Raw Energy -----------------------------------------**/
      GALILEO_hits = RAW->get_GALILEO_am_Fired();
         for(int i=0; i<GALILEO_hits; i++){
         GalID = RAW->get_GALILEO_Det_id(i) * 3 + RAW->get_GALILEO_Crystal_id(i);
        tmpGAL[GalID] = RAW->get_GALILEO_Chan_E(i);
        hGAL_Raw_E[GalID]->Fill(tmpGAL[GalID]);
         }
   }





//-----------------------------------------------------------------------------------------------------------------------------//

void EventUnpackProc::checkPADI_or_PADIWA(){

  std::ifstream PADIFILE("Configuration_Files/PADI_or_PADIWA.txt");

  std::string line;

  if(PADIFILE.fail()){
    std::cerr << "Could not find Configuration_Files/PADI_or_PADIWA.txt file" << std::endl;
    exit(1);
  }
  bool P_or_PW = false;
  while(std::getline(PADIFILE,line)){
    if(line[0] == '#') continue;

    if(line == "PADI") P_or_PW = true;
    if(line == "PADIWA") P_or_PW = false;

    if(line != "PADIWA" && line != "PADI"){
      std::cerr << line << " module of PLASTIC not known!" <<std::endl;
      exit(1);
    }
  }

  PADI_OR_PADIWA = P_or_PW;

}

//-----------------------------------------------------------------------------------------------------------------------------//

void EventUnpackProc::checkTAMEXorVME(){

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

bool EventUnpackProc::Check_Cal_Plastic(){
  ifstream data("Configuration_Files/PLASTIC_CALIB_FILE.txt");
  if(data.fail()){
    cerr << "Could not find Calibration type file for PLASTIC" << endl;
    exit(0);
  }
  string line;
  const char* format = "%s %d";
  char s[100];
  int val;
  bool CALIBRATE = false;

  while(data.good()){
    getline(data,line,'\n');
    if(line[0] == '#') continue;
    sscanf(line.c_str(),format,&s,&val);
    if(string(s) == string("ONLINE")) CALIBRATE = (val == 1);
  }

  return CALIBRATE;

}

//-----------------------------------------------------------------------------------------------------------------------------//

void EventUnpackProc::print_MBS(int* pdata,int lwords){
  cout << "---------------------\n";
  for(int i = 0;i < lwords;++i){
    cout << hex << *(pdata + i) << " ";
    if(i % 5 == 0 && i > 0) cout << endl;
  }
  cout << "\n---------------------\n";
}
//-----------------------------------------------------------------------------------------------------------------------------//

TH1I* EventUnpackProc::MakeH1I(const char* fname,
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

TH2I* EventUnpackProc::MakeH2I(const char* fname,
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
TH1I* EventUnpackProc::MakeH1I_TPC(const char* foldername, const char* name, int nameindex,
                  Int_t nbinsx, Float_t xmin, Float_t xmax,
                  const char* xtitle, Color_t linecolor, Color_t fillcolor)
{
  char fullname[100];
  if(nameindex>=0)
    sprintf(fullname,"%s%s",tpc_name_ext1[nameindex],name);
  else
    strcpy(fullname, name);
  return MakeH1I(foldername, fullname, nbinsx, xmin, xmax, xtitle,
         linecolor, fillcolor);
}
const  char* EventUnpackProc::tpc_name_ext1[7]={"TPC21_","TPC22_","TPC23_","TPC24_","TPC41_","TPC42_", "TPC31_"};
const  char* EventUnpackProc::tpc_folder_ext1[7]={"TPC21","TPC22","TPC23","TPC24","TPC41","TPC42","TPC31"};


//-----------------------------------------------------------------------------------------------------------------------------//
//                                                            END                                                              //
//-----------------------------------------------------------------------------------------------------------------------------//
