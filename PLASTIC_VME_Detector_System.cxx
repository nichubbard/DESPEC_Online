#include "PLASTIC_VME_Detector_System.h"
#include "TGo4MbsEvent.h"
#include <fstream>
#include "TGo4MbsSubEvent.h"

#include "QDC.h"
#include <map>
#include <cstdlib>

using namespace std;

//---------------------------------------------------------------

PLASTIC_VME_Detector_System::PLASTIC_VME_Detector_System(){

    // ================= THE SUDIPTA VME VARIABLES ================== //
    
    
    fthr = 50; //100
    
    exit_flg=0;

    raw_data=0;
    energy_data=0;

    card=0;
    energy=0;

    nof_hits=0;

    raw_data_old=0;

    value = 0;  
    value_t = 0;
    
    chNo = 0;
    
//     lwords = 0;
//     next_sub_evt = 0;
//     data_field=0;
//     
//     pdata=0;
    
    
    
    // ===== THE ONE'S I ADDED ======== //
    
//     VME_QDC_Data = std::vector<double>(100,0);
//     VME_QDC_Channels = std::vector<double>(100,0);
//     VME_TDC_Data = std::vector<double>(100,0);
//     VME_TDC_Channels = std::vector<double>(100,0);
     //TDC_iterator = new int;
     VME_QDC_Data1 = new double[32];
     VME_QDC_Data2 = new double[32];
     VME_QDC_Channels = new int[32];
     VME_TDC_Data = new double[50];
     VME_TDC_Channels = new int[50];
         
//     Scalar_Data = new double[50];
  //   Scalar_Channels = new int[50];
     
     TDC_iterator = 0;
     for(int i; i<50; i++){
     VME_TDC_Data[i] = 0;
     VME_TDC_Channels[i] = 0;
//     Scalar_Data[i] = 0;
  //   Scalar_Channels[i] = 0;
     }
     
     for(int i; i<32; i++){
     VME_QDC_Data1[i] = 0;
     VME_QDC_Data2[i] = 0;
     VME_QDC_Channels[i] = 0;
     
    
     }
     

}

//--------------------destructor-------------------------------------------

PLASTIC_VME_Detector_System::~PLASTIC_VME_Detector_System(){
    
  
    delete[]VME_QDC_Data1;
    delete[]VME_QDC_Data2;
    delete[]VME_QDC_Channels;
    delete[]VME_TDC_Data;
    delete[]VME_TDC_Channels;
//    delete[]Scalar_Data;
   // delete[]Scalar_Channels;

}

//---------------------------------------------------------------

void PLASTIC_VME_Detector_System::get_Event_data(Raw_Event* RAW){
    
    RAW->set_DATA_PLASTIC_VME(QDC_iterator, TDC_iterator,VME_QDC_Data1, VME_QDC_Data2, VME_QDC_Channels,VME_TDC_Data,VME_TDC_Channels);
   // RAW->set_DATA_SCALAR(Scalar_iterator,Scalar_Data, Scalar_M);

}

//---------------------------------------------------------------
 void PLASTIC_VME_Detector_System::Process_MBS(TGo4MbsSubEvent* psubevt){
    
   lwords = psubevt->GetIntLen();
   pdata = psubevt->GetDataField();

 }
 
 void PLASTIC_VME_Detector_System::Process_MBS(TGo4MbsEvent* test){
    next_sub_evt = test-> NextSubEvent();
   
 }
 
   
 void PLASTIC_VME_Detector_System::Process_MBS(int* hah){

     reset_fired_channels();
     pl_data=pdata;

   
       
    //***********************************
    //***  CAEN V785 or V775 or V792  ***
 // printf("0x%08x\n", (unsigned int) caen_header);       
    // module loop
     int position = 0;
  //   while(true){
       diff =  pdata-pl_data;
      
         caen_header = *pdata++;
         type = (caen_header & TY_MASK) >> 24;    
              
            position++;
                
         if (type == 2){        
            // get geographical address 
            geo  = (caen_header & GEO_MASK) >> 27;      
            //QDC_iterator = d_cnt;    
                   
            // get number of fired channels
            d_cnt= (caen_header & FCH_MASK) >> 8;   
            for(int jj = 0;jj < d_cnt;++jj){
           
              data=*pdata++;
                
             //  printf("data QDc 0x%08x\n", (unsigned int) data);
                position++;
                
                chan = (data & CH_MASK) >> 16;     // [0..31] 
                data = data & DA_MASK;
               
                if(geo==6){
       
                 
                 if (chan >=0 && chan <=31){
                VME_QDC_Channels[chan] = chan;                 
                VME_QDC_Data1[chan] = data;
                          
                                            }
                       QDC_iterator++;           
                   
                                    }                                
                                }                          
                          }
      
            pdata++;   // skip EOB word 
            position++;
                     
       

    //************************
    //***  CAEN V1290 TDC  ***
    //************************
    pdata++; //skip empty word
   // cout << "Position " << position<<endl;
       // read header word 
        
     caen_header=*pdata++;
     type = (caen_header >> 27) & 0x1F;   // get global header ID 
 
     geo  = caen_header & 0x1F;            // get geographical address 
      
  // printf("caen_header1 cut0x%08x\n", (unsigned int) caen_header & 0xFFFF0000); 
  
    while(true){
        
         data=*pdata++; 
         //data = *pdata++;  
        // trigger time tag
        //if ((data & 0xF8000000) == 0x88000000){ //data = *pdata++;
        
       // printf("data head cut 0x%08x\n", (unsigned int) data ); 
       if ((data & 0xF8000000) == 0x88000000){  // trigger time tag
                  data=*pdata++;
       }
      //  printf("data 0x%08x\n", (unsigned int) data); 
       
        // global trailer
        if ((data & 0xF8000000 ) == 0x80000000){
          //  if (data & 0x07000000) printf("V1290 Global Trailer contains Status Info: 0x%08x\n",(unsigned int)data);
            break;  // exit loop
          //  cout << "trailer TDC " << endl;           
        }
      
        chan1 = (data & CH_MASK2) >> 21;     // [0..31] 
        data = data & DA_MASK2;
      
       
      // cout << "TDC_iterator "<< TDC_iterator << endl;
         if (chan1 >=0 && chan1 <=31)
      {
          
            VME_TDC_Channels[TDC_iterator] = chan1;      
            VME_TDC_Data[chan1] = data;
       //   cout <<"chan1 " <<chan1 <<" VME_TDC_Channels[TDC_iterator] "<< VME_TDC_Channels[TDC_iterator] <<  " TDC_iterator " << TDC_iterator <<endl;
     
      //    cout <<" TDC_iterator " << TDC_iterator <<" VME_TDC_Channels[TDC_iterator] " << VME_TDC_Channels[TDC_iterator]<< " chan1 " << chan1 << " VME_TDC_Data[TDC_iterator] " <<VME_TDC_Data[TDC_iterator] << endl;
     //   cout <<"TDC_iterator "<<TDC_iterator << " chan1 " << chan1 <<" VME_TDC_Channels[TDC_iterator] " <<VME_TDC_Channels[TDC_iterator] <<" VME_TDC_Data[chan1] "<< VME_TDC_Data[chan1] << " data " << data<< endl;
              
            }
         TDC_iterator++;   
    }
    
     
   //************************
    //***  CAEN V1290 TDC (SCALAR) ***
    //************************
    
       caen_header = *pdata++;    
    //  printf("caen_header scalar 0x%08x\n", (unsigned int) caen_header& 0xFF000000);
     while(true){    
         data = *pdata++;  

        // global trailer
        if ((data & 0xFFFF0000) == 0x4c000000){
            break;  // exit loop
        //cout << "scalar break "<<endl;
//            
        }
//       //   printf("0x%08x\n", (unsigned int) data);
        int chan = (data & CH_MASK2) >> 21;     // [0..31] 
        data = data & DA_MASK2;
        type = (data & TY_MASK) >> 24; 
        geo  = caen_header & 0x1F;       
      //  cout << "type " << type << " geo " << geo<< " chan "<< chan << " data "<< data<< endl;
//        Scalar_Data[Scalar_iterator] = data;
   //     Scalar_Channels[Scalar_iterator] = chan;
   //     Scalar_iterator++;
     }
 }
 
   

void PLASTIC_VME_Detector_System::reset_fired_channels(){
    TDC_iterator = 0;
    QDC_iterator = 0;
//    Scalar_iterator = 0;
     for (int i=0; i<QDC_iterator; i++){
        VME_QDC_Data1[i] = 0 ;
        VME_QDC_Data2[i] = 0 ;
        VME_QDC_Channels[i] =0;
    }
    for (int i=0; i<TDC_iterator; i++){
        VME_TDC_Data[i] = 0 ;
        VME_TDC_Channels[i] =0;
    }
   // for (int i=0; i<Scalar_iterator; i++){
   //     Scalar_Data[i] = 0;
   //     Scalar_Channels[i] = 0;
   // }
       
}

int* PLASTIC_VME_Detector_System::get_pdata(){return pdata;}

//---------------------------------------------------------------
