#ifndef PLASTIC_VME_DETECTOR_SYSTEM_H
#define PLASTIC_VME_DETECTOR_SYSTEM_H


#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <TFile.h>
#include "TH1.h"
#include "TObjArray.h"

#include "TGo4MbsEvent.h"
#include "Data_Stream.cxx"
#include "EventBuilder.cxx"
//#include "TSCNUnpackEvent.h"
#include "EventUnpackStore.h"

// #include "QDC_751.h"
// #include "TDC_1290.h"
//#include "QDC.h"
#include "Detector_System.cxx"
//#include "Raw_Event.h"
#include "TGo4EventProcessor.h"

typedef unsigned long long ULong64_t;
//typedef unsigned long ULong;

class PLASTIC_VME_Detector_System : public Detector_System {


private:

    // ================= THE SUDIPTA VME VARIABLES ================== //
    
    
    int fthr; //100
    
    int *pl_data,diff,first_e_sample, sub_evt_len;

    char exit_flg;
    
    int* pdata, *data_field;
    
    int lwords;
    
    bool next_sub_evt;

    int sub_evt_length; 
    
    long raw_data, energy_data, first_lw;

    long card, chan, rchan, samples, energy, trg_cnt, ft_cnt, trailer;

    long hit_pattern, nof_hits, i;

    long raw_data_old;
    
    unsigned short *pl_data16, first_sample;

    char index_pat;
  
    long base_sum, base2_sum, top_sum;

    char pi_flg, rt_flg;
  
    int value, value_t, chNo;
    
    long d_cnt;
    
    unsigned long type, geo, caen_header, data;

    unsigned long chan1;


    const unsigned int  OV_MASK = 0x00001000;  /* overflow bit */
    const unsigned int  UN_MASK = 0x00002000;  /* underthreshold bit */
    const unsigned int   V_MASK = 0x00004000;  /* data valid bit */
    const unsigned int FCH_MASK = 0x00003F00;  /* number of fired channels mask */
    const unsigned int  CH_MASK = 0x003F0000;  /* channel mask */
    const unsigned int  DA_MASK = 0x00000FFF;  /* data mask */
    const unsigned int  TY_MASK = 0x07000000;  /* type mask */
    const unsigned int GEO_MASK = 0xF8000000;  /* geo address mask */

    const unsigned int CH_MASK2 = 0x03E00000;  // channel mask for V1290
    const unsigned int DA_MASK2 = 0x001FFFFF;  // data mask for V1290*/
    
    // ===== THE ONE'S I ADDED ======== //
    
//     std::vector<double> VME_QDC_Data;
//     std::vector<double> VME_QDC_Channels;
//     std::vector<double> VME_TDC_Data;
//     std::vector<double> VME_TDC_Channels;

     int QDC_iterator;
     int TDC_iterator;   
     double* VME_QDC_Data1;
     double* VME_QDC_Data2;
     int* VME_QDC_Channels;
     double* VME_TDC_Data;
     int* VME_TDC_Channels;
     
//      int Scalar_iterator;
//      double* Scalar_Data;
//      int* Scalar_Channels;
     
     int unknown_header_counter;
     
     void reset_fired_channels();
     
    
     int num_TDC_modules = 0;  //Set in the constructor
          
    
public:
    PLASTIC_VME_Detector_System();
    ~PLASTIC_VME_Detector_System();

    void Process_MBS(TGo4MbsSubEvent* psubevt);
    void Process_MBS(TGo4MbsEvent* test);
   
    //functions from abstract class Detector_System
    void Process_MBS(int*);

    void get_Event_data(Raw_Event*);

    int* get_pdata();
    
    bool calibration_done(){return false;}

    void write(){return;};
    void set_Gain_Match_Filename(std::string){return;};
    void make_histos (Int_t);
    
  //  virtual Bool_t BuildEvent(TGo4EventElement* dest);
    

};



#endif

