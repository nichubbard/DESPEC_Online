#include "GALILEO_Detector_System_TEST.h"

#include <cstdlib>
#include <map>

#include "FEBEX.h"
#include "DESPEC_Array_Sizes.h"



using namespace std;

//---------------------------------------------------------------

GALILEO_Detector_System::GALILEO_Detector_System(){

    //set amount of Germanium Detectors

  max_am_dets = GALILEO_MAX_HITS;

    fired_FEBEX_amount = 0;
   
    Sum_Time = new ULong64_t[max_am_dets];
    Hit_Pattern = new int[max_am_dets];

  det_ids = new int[max_am_dets];
  crystal_ids = new int[max_am_dets];
  Chan_Time = new ULong64_t[max_am_dets];
  Chan_Energy = new double[max_am_dets];
  Pileup = new bool[max_am_dets];
  Overflow = new bool[max_am_dets];

    GALILEO_T_CALIB = new GALILEO_Time_Calibration();
    GALILEO_E_CALIB = new GALILEO_Energy_Calibration();

    load_board_channel_file();

    
}

//---------------------------------------------------------------

GALILEO_Detector_System::~GALILEO_Detector_System(){

  //GALILEO_map.clear(); - not necessary
    
  delete[] det_ids;
  delete[] crystal_ids;
    delete[] Sum_Time;
    delete[] Pileup;
    delete[] Hit_Pattern;
    delete[] Chan_Time;
    delete[] Chan_Energy;
    delete[] Pileup;
    delete[] Overflow;
    
    delete GALILEO_T_CALIB;
    delete GALILEO_E_CALIB;
}

//---------------------------------------------------------------

void GALILEO_Detector_System::load_board_channel_file()
{
  ifstream file("Configuration_Files/GALILEO_Detector_Map.txt");

  if(file.fail())
  {
    cerr << "Could not find FEBEX GALILEO Mapping File!" << endl;
        exit(0);
    }

  constexpr auto ignore = std::numeric_limits<std::streamsize>::max();

  while(file.good())
  {
    if (file.peek() == '#')
    {
      file.ignore(ignore, '\n');
      continue;
    }
    //#Board ID   Ch ID    Det Num    Crysal Num
    int board_id, channel_num, det_num, crystal_num;
    file >> board_id >> channel_num >> det_num >> crystal_num;
    file.ignore(ignore, '\n');
    // Map board,channel to detector,crystal
    GALILEO_map[std::make_pair(board_id,channel_num)] = std::make_pair(det_num, crystal_num);
    }
}

//---------------------------------------------------------------

void GALILEO_Detector_System::get_Event_data(Raw_Event* RAW){
    
  RAW->set_DATA_GALILEO(fired_FEBEX_amount,Sum_Time,Ge_channels,Chan_Time,Chan_Energy,det_ids,crystal_ids,Pileup,Overflow);
    
}

//---------------------------------------------------------------

void GALILEO_Detector_System::Process_MBS(int* pdata){

    reset_fired_channels();
    
  std::pair<int,int> current_det;
    
    this->pdata = pdata;
    
    bool FEBEX_data_loop = true;
    
  int num_modules = GALILEO_FEBEX_MODULES; // Number of modules in the data //
    
    fired_FEBEX_amount = 0;
    
    FEBEX_Add* FEBEX_add  = (FEBEX_Add*) this->pdata;

  while (FEBEX_add->add == 0xadd){
        this->pdata++;
        FEBEX_add = (FEBEX_Add*) this->pdata;
    }

    FEBEX_Header* FEBEXhead  = (FEBEX_Header*) this->pdata;

    while(FEBEX_data_loop){

    if (FEBEXhead->ff == 0xFF){ // FEBEX module idicator //

            // FEBEXhead->ff;
            // FEBEXhead->chan_head;
            // FEBEXhead->three_four;
    
            board_id = FEBEXhead->chan_head;
            
            this->pdata++; // Moves to Channel Size //
            
            FEBEX_Chan_Size *fbx_size=(FEBEX_Chan_Size*) this->pdata;    
            
            num_channels = ((fbx_size->chan_size)/4) - 1;
    
            if (num_channels == 0) num_modules--;
        
            this->pdata++; // Moves to Event Timestamp //

            FEBEX_Half_Time *fbx_hT=(FEBEX_Half_Time*) this->pdata;
        
            this->pdata++; // Moves to rest of Event Timestamp //
            
            FEBEX_Evt_Time *fbx_time=(FEBEX_Evt_Time*) this->pdata;

            ULong64_t tmp_ext_time = ((fbx_hT->ext_time));
            
            tmp_Sum_Time = (fbx_time->evt_time)+ (tmp_ext_time<<32);//((fbx_hT->ext_time)<<32);
            
            
            this->pdata++; // Moves to Pileup & Hit Pattern //
                        
            FEBEX_Flag_Hits *fbx_flag=(FEBEX_Flag_Hits*) this->pdata;
        
            tmp_Pileup = fbx_flag->pile_flags;
        
            tmp_Hit_Pattern = fbx_flag->hit_pattern;
        
      for(int j = 15; j >= 0; j--)
      {
        if(tmp_Pileup & (1 << j)) {
                    pileup_flags[j] = 1;
                }
        if(tmp_Hit_Pattern & (1 << j)) {
                    Ge_channels[j] = j;
                    num_channels_fired++;
                }
            }
            this->pdata++; // Moves to DEADBEEF //
        }
    else if (FEBEXhead->ff == 0xF0){ // FEBEX channel indicator //
            this->pdata--; // Moves back to DEADBEEF so channel loop functions properly //

            for(int i=0; i<num_channels; ++i){
                this->pdata++; // Moves to channel header //

                FEBEX_Chan_Header *fbx_Ch=(FEBEX_Chan_Header*) this->pdata;
        
                int tmp_Ch_ID = fbx_Ch->Ch_ID;
        
        //if(false && pileup_flags[tmp_Ch_ID] == 1) this->pdata += 3;
            
        //else
        {
                        
          auto idx = std::make_pair(board_id, tmp_Ch_ID);
            
          if(GALILEO_map.find(idx) != GALILEO_map.end())
          {
            
            //det_ids[i] = current_det;
            Sum_Time[fired_FEBEX_amount] = tmp_Sum_Time;
            this->pdata++; // Moves to rest of channel timestamp //
        
            FEBEX_TS *fbx_Ch_TS=(FEBEX_TS*) this->pdata; 
            ULong64_t tmp_ext_chan_ts = (fbx_Ch->ext_chan_ts);
    
            Chan_Time[fired_FEBEX_amount] = ((fbx_Ch_TS->chan_ts)+(tmp_ext_chan_ts<<32))*10; // in nanoseconds
            this->pdata++; // Moves to Channel Energy //
    
            FEBEX_En *fbx_Ch_En=(FEBEX_En*) this->pdata; 
            
            Chan_Energy[fired_FEBEX_amount] = fbx_Ch_En->chan_en;
            Pileup[fired_FEBEX_amount] = fbx_Ch_En->pileup != 0;
            Overflow[fired_FEBEX_amount] = fbx_Ch_En->overflow != 0;
            det_ids[fired_FEBEX_amount] = GALILEO_map[idx].first;
            crystal_ids[fired_FEBEX_amount] = GALILEO_map[idx].second;
            //cout<<" det_ids[fired_FEBEX_amount] " <<  det_ids[fired_FEBEX_amount] << " fired_FEBEX_amount " << fired_FEBEX_amount <<endl;
            //cout <<  "1) Chan_Energy[current_det] "<<  Chan_Energy[fired_FEBEX_amount] << " current_det " <<  det_ids[fired_FEBEX_amount] <<endl;
            //Calibrate_FEBEX(current_det);
            //
            //std::cout << "pdata = " << std::hex << *this->pdata << std::dec << std::endl;
            //if ((*this->pdata & 0x80000000) != 0) cout<<" Bit 31 set "<<endl;
            //if ((*this->pdata & 0x40000000) != 0) cout<<" Bit 30 set "<<endl;
            //std::cout << "Pileup = " << fbx_Ch_En->pileup << ", OF = " << fbx_Ch_En->overflow << std::endl;
    
            this->pdata++; // Moves to Future Use //
    
            fired_FEBEX_amount++;
            
            // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ //
            // @@@@ Traces Would Go Here @@@@@ //
            // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ //
            
                    }   
          else
          {
           // std::cout << "Unknown detector with channel: " << board_id << ", " << tmp_Ch_ID << std::endl;
            this->pdata += 3; // Moves to Future Use //
          }
                }
                        
            }
            num_modules--;
        }
    
        if (num_modules != 0){
            this->pdata++; 
            FEBEXhead  = (FEBEX_Header*) this->pdata;
        }
        else FEBEX_data_loop = false; // Exits FEBEX Loop //
    }
}

//---------------------------------------------------------------

void GALILEO_Detector_System::reset_fired_channels(){
    
    fired_FEBEX_amount = 0;
    num_channels_fired = 0;
    
    for(int i = 0;i < max_am_dets;++i){
        Sum_Time[i] = -1;
        pileup_flags[i] = -1;
        Ge_channels[i] = 0;
        Pileup[i] = -1;
        Hit_Pattern[i] = 0;
        Chan_Time[i] = 0;
        Chan_Energy[i] = 0;
    }
}

//---------------------------------------------------------------

void GALILEO_Detector_System::Calibrate_FEBEX(int id){

    Sum_Time[id] = GALILEO_T_CALIB->Calibrate_FEBEX_Sum_T(Sum_Time[id],id);

    Chan_Time[id] = GALILEO_T_CALIB->Calibrate_FEBEX_Chan_T(Chan_Time[id],id);
   // Chan_Energy[id] = GALILEO_E_CALIB->Calibrate_FEBEX_E(Chan_Energy[id],id);
   // cout << "2) Chan_Energy[id] " <<Chan_Energy[id] << " id " << id << endl;
}

//---------------------------------------------------------------

int* GALILEO_Detector_System::get_pdata(){return pdata;}

//---------------------------------------------------------------
