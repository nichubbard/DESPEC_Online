#ifndef GALILEO_DETECTOR_SYSTEM_H
#define GALILEO_DETECTOR_SYSTEM_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <map>


#include "GALILEO_Energy_Calibration.h"
#include "GALILEO_Time_Calibration.h"

#include "FEBEX.h"

#include "Detector_System.cxx"

typedef unsigned long long ULong64_t;

class GALILEO_Detector_System : public Detector_System{


private:
    int max_am_dets;
    
    int board_id;
    int num_channels;
    
    int* pdata;

    ULong64_t tmp_Sum_Time;
    int tmp_Pileup;
    int tmp_Hit_Pattern;
    
    int* det_ids;
    int* crystal_ids;
    ULong64_t* Sum_Time;
    int* Hit_Pattern;
    ULong64_t* Chan_Time;
    double* Chan_Energy;
  bool* Overflow;
  bool* Pileup;
    

    std::map<std::pair<int,int>, std::pair<int,int>> GALILEO_map;
    
    int num_channels_fired = 0;

    int pileup_flags[16];
    int Ge_channels[16];
    int fired_FEBEX_amount;
    void load_board_channel_file();
    void reset_fired_channels();
    void Calibrate_FEBEX(int);

    GALILEO_Time_Calibration* GALILEO_T_CALIB;
    GALILEO_Energy_Calibration* GALILEO_E_CALIB;


public:
    GALILEO_Detector_System();
    ~GALILEO_Detector_System();
    
    void Process_MBS(TGo4MbsSubEvent* psubevt){};

    void Process_MBS(int*);
    void get_Event_data(Raw_Event*);
    int* get_pdata();


    // Useless but needed for Detector_System
    bool calibration_done(){return false;}
    void write(){return;};
    void set_Gain_Match_Filename(std::string){return;};


};



#endif
