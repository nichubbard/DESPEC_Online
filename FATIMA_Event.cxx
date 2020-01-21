#include "FATIMA_Event.h"

using namespace std;

//---------------------------------------------------------------

FATIMA_Event::FATIMA_Event(int* positions_tmp,int length,Raw_Event* RAW) : Events(positions_tmp,length){
//     set_DATA(RAW);
  
    for(int i = 0;i < 100;++i) DATA.Energy[i] = 0;
    DATA.amountHits = 0;

    //set_DATA(RAW);
}


//---------------------------------------------------------------

FATIMA_Event::~FATIMA_Event(){}

//---------------------------------------------------------------

inline void FATIMA_Event::set_DATA(Raw_Event* RAW){
    int am_fired = RAW->get_FAT_det_fired();
    double energy_tot = 0;
    //for(int i = 0;i < am_fired;++i) energy_tot += RAW->get_FAT_E(i);
    
    DATA.amountHits = am_fired;
//     FATIMA_DataStruct* P_Tmp = RAW->PassFATIMA();
//     //set event data vectors
//     DATA.SetDATA_Directly(P_Tmp->FAT_DET_FIRED,P_Tmp->FAT_id,P_Tmp->FAT_E,P_Tmp->FAT_ratio,P_Tmp->FAT_t,
//                  P_Tmp->FAT_t_qdc,P_Tmp->FAT_QDCs_FIRED,P_Tmp->FAT_QDC_id,P_Tmp->FAT_QLong,
//                  P_Tmp->FAT_QLong_Raw,P_Tmp->FAT_QShort_Raw,P_Tmp->FAT_QDC_t_coarse,
//                  P_Tmp->FAT_QDC_t_fine,P_Tmp->FAT_TDCs_FIRED,P_Tmp->FAT_TDC_id,P_Tmp->FAT_TDC_timestamp);
//     P_Tmp = nullptr;

}

//---------------------------------------------------------------
//vector 
//double FATIMA_Event::get_energy(){return 0;}
double FATIMA_Event::get_energy(){return (am_fired > 0) ? energy_tot : -1;}

//---------------------------------------------------------------

void FATIMA_Event::Write_Event(Tree_Creator* Tree){
	//User defined what variables should be written
}

//---------------------------------------------------------------

// FATIMA_DataStruct* FATIMA_Event::GET_FATIMA(){
//     return &DATA;
// }

//---------------------------------------------------------------