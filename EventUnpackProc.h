// $Id: TSCNUnpackProc.h 755 2011-05-20 08:04:11Z linev $
//-----------------------------------------------------------------------
//       The GSI Online Offline Object Oriented (Go4) Project
//         Experiment Data Processing at EE department, GSI
//-----------------------------------------------------------------------
// Copyright (C) 2000- GSI Helmholtzzentrum für Schwerionenforschung GmbH
//                     Planckstr. 1, 64291 Darmstadt, Germany
// Contact:            http://go4.gsi.de
//-----------------------------------------------------------------------
// This software can be used under the license agreements as stated
// in Go4License.txt file which is part of the distribution.
//-----------------------------------------------------------------------

#ifndef EVENTUNPACKPROC_H
#define EVENTUNPACKPROC_H

#include "Riostream.h"

// Root Includes //
#include "TROOT.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TCutG.h"
#include "TArc.h"
#include "TTree.h"


// Go4 Includes //
#include "TGo4UserException.h"
#include "TGo4Picture.h"
#include "TGo4MbsEvent.h"


#include "EventUnpackStore.h"
#include "CalibParameter.h"
#include "CorrelParameter.h"

#include "Detector_System.cxx"
#include "AIDA_Detector_System.h"
#include "FATIMA_Detector_System.h"
#include "PLASTIC_VME_Detector_System.h"
#include "PLASTIC_TAMEX_Detector_System.h"
#include "GALILEO_Detector_System_TEST.h"
#include "DESPECAnalysis.h"
#include "FRS_Detector_System.h"

#include "Data_Stream.cxx"
#include "White_Rabbit.h"
#include "AIDA_Event.h"
#include "AIDA_Decay_Event_Store.h"
#include "PLASTIC_Data_Stream.h"

#include "Raw_Event.h"

#include "EventBuilder.cxx"
#include "Time_EventBuilder.h"
#include "AidaMover.h"
#include <string>

#include <array>
#include <vector>

#include "AIDA_Headers.h"
#include "AIDA_Event.h"
#include "AIDA_Data_Types.h"
#include "TGo4MbsEvent.h"
#include "Data_Stream.cxx"
#include "EventBuilder.cxx"
//#include "TSCNUnpackEvent.h"
#include "EventUnpackStore.h"
//#include "AIDA_Processor.h"
//#include "FRS_Basic.h"
#include "Detector_System.cxx"

///////////////////////////////

using namespace std;

#include "TGo4EventProcessor.h"

	class EventUnpackStore;


	class EventUnpackProc : public TGo4EventProcessor
	{
		public:

            EventUnpackProc();
			EventUnpackProc(const char* name);
            virtual ~EventUnpackProc();
						virtual void UserPostLoop();
			//FATIMA variables
			double FATgate1_low, FATgate1_high;
			double FATgate2_low, FATgate2_high;

            Bool_t BuildEvent(TGo4EventElement* dest);
             void Process_AIDA_Event(EventUnpackStore* event);

             ofstream WR_out;

//            std::vector<AidaCluster> EventsToClusters(std::vector<AidaEvent> const&);
//             AidaHit ClusterPairToHit(std::pair<AidaCluster, AidaCluster> const&);

            //parameters
            CalibParameter *fCal;
            CorrelParameter *fCorrel;
            EventUnpackStore *fOutput;
            AidaUnpackData    fAida;

		protected:
//                 TGo4MbsEvent *fInput;

            //FRS histogram settings
      TH1I* MakeH1I(const char* foldername,
        const char* histoname,
        Int_t nbins,
        Float_t xmin, Float_t xmax,
        const char* xtitle = "channels",
        Color_t linecolor = 2,
        Color_t fillcolor = 6,
        const char* ytitle = 0);

      TH2I* MakeH2I(const char* foldername,
        const char* histoname,
        Int_t nbinsx, Float_t xmin, Float_t xmax,
        Int_t nbinsy, Float_t ymin, Float_t ymax,
        const char* xtitle, const char* ytitle,
        Color_t marker);

      TH1I* MakeH1I_TPC(const char* foldername, const char* name, int nameindex,
            Int_t nbins, Float_t xmin, Float_t xmax,
            const char* xtitle = "channels", Color_t linecolor = 2, Color_t fillcolor = 6);

        const static char* tpc_name_ext1[7];
        const static char* tpc_folder_ext1[7];
        
        //--------Scaler Graphs----------
  TH1I          *hScaler_per_s[64];
  TH1I          *hScaler_per_100ms[64];
  TH1I          *hScaler_per_spill[64];  
  int            scaler_ch_1kHz=0;
  int            scaler_ch_spillstart=0;
  UInt_t         scaler_initial[64];
  UInt_t         scaler_previous[64];
  int            scaler_check_first_event = 1;
  
   Float_t Music_dE[3], Music_dE_corr[3];

    Int_t Music_E1[8], Music_E2[8], Music_T1[8], Music_T2[8];

    Float_t sci_l[12], sci_r[12], sci_e[12], sci_tx[12], sci_x[12];

    Float_t sci_tofll2, sci_tofll3, sci_tof2, sci_tofrr2, sci_tofrr3, sci_tof3;

    Float_t ID_x2, ID_y2, ID_a2, ID_b2;

    Float_t ID_x4, ID_y4, ID_a4, ID_b4;

    Float_t TPC_X[7], TPC_Y[7];
    Int_t TPC_LT[7][2], TPC_RT[7][2];
    Float_t TPC_X0, TPC_X1;


    Int_t sci_dt_21l_21r, sci_dt_41l_41r, sci_dt_42l_42r, sci_dt_43l_43r;

    Int_t sci_dt_21l_41l, sci_dt_21r_41r, sci_dt_21l_42l, sci_dt_21r_42r;

    Float_t ID_brho[2], ID_rho;

    Float_t beta, beta3, gamma;

    Float_t AoQ, AoQ_corr;

    Float_t ID_z, ID_z2, ID_z3;

    Float_t timestamp, ts, ts2;
    
    Int_t time_in_ms;
    Int_t spill_count;
    Int_t ibin_for_s;
    Int_t ibin_for_100ms;
    Int_t ibin_for_spill;
    Int_t ibin_clean_for_s;
    Int_t ibin_clean_for_100ms;
    Int_t ibin_clean_for_spill;
    UInt_t increase_scaler_temp;

			TH1* hID_x2;
			TH1* hID_y2;
			TH1* hID_a2;
			TH1* hID_b2;

			TH1* hID_x4;
			TH1* hID_y4;
			TH1* hID_a4;
			TH1* hID_b4;


            /// positions and control sum

            TH1I *hTPC_X[7];
            TH1I *hTPC_Y[7];
            TH2I *hcTPC_XY[7];
            TH2I *hTPC_LTRT[7];
            TH1I *hTPC_DELTAX[7];

            // CSUM[index][anode_no]
            TH1I *hTPC_CSUM[7][4];

            TH1* htpc_X2;
            TH1* htpc_Y2;
            TH1* htpc_X4;
            TH1* htpc_Y4;

            TH1I *hSCI_L[12];
            TH1I *hSCI_R[12];
            TH1I *hSCI_E[12];
            TH1I *hSCI_Tx[12];
            TH1I *hSCI_X[12];
            TH2I *hSCI_dE24;

            TH2I *hSCIdE41_TPC42X;
            TH2I *hSCIdE41L_TPC41X;
            TH2I *hSCIdE41L_TPC42X;
            TH2I *hSCIdE41R_TPC42X;
            TH2I *hSCIdE41R_TPC41X;
            TH2I *hSCIdE21_TPC42X;
            TH2I *hSCIdE21L_TPC42X;
            TH2I *hSCIdE21L_TPC41X;
            TH2I *hSCIdE21R_TPC42X;
            TH2I *hSCIdE21R_TPC41X;

            TH1I *hSCI_TofLL2;
            TH1I *hSCI_TofLL3;
            TH1I *hSCI_Tof2;
            TH1I *hSCI_TofRR2;
            TH1I *hSCI_TofRR3;
            TH1I *hSCI_Tof3;

			TH1* hSCI_dT_21l_41l;
			TH1* hSCI_dT_21r_41r;

			TH1* hSCI_dT_21l_42l;
			TH1* hSCI_dT_21r_42r;

			TH1* hSCI_dT_21l_81l;
			TH1* hSCI_dT_21r_81r;

			TH1I *hID_AoQ;
            TH1I *hID_AoQ_corr;
//             TH2I *hID_x2AoQ;
//             TH2I *hID_Z_AoQ;
//             TH2I *hID_Z_AoQ_zsame;
//             TH2I *hID_Z_AoQ_corr;
			TH1I *hID_Z;
            TH1I *hID_Z2;
            TH1I *hID_Z3;

            TH2I *hID_x2a2;
            TH2I *hID_y2b2;
            TH2I *hID_x4a4;
            TH2I *hID_y4b4;

          //  TH2I *hID_Z_Z2;
           // TH2I *hID_Z_Z3;

            TH2I *hID_Z_dE2;

            TH2I *hID_Z_Sc21E;
            TH2I *hID_SC41dE_AoQ;
            TH2I *hID_x2z;
            TH2I *hID_x4z;
            TH2I *hID_x4AoQ;
            TH2I *hID_E_Xs4;
            TH2I *hID_E_Xs2;

            TH2I *hID_x2x4;
            TH2I *hID_dEToF;
            TH1I *hID_beta;
            TH1I *hID_BRho[2];

			TH1* htimestamp;
			TH1* hts;
			TH1* hts2;

            TH1 *hMUSIC1_E[8];
            TH1 *hMUSIC1_T[8];
            TH1 *hMUSIC1_dE;

            TH1 *hMUSIC2_E[8];
            TH1 *hMUSIC2_T[8];
            TH1 *hMUSIC2_dE;

            TH2 *hMUSIC1_dE1dE2;

            TH2 *hMUSIC1_MUSIC2;

            //AIDA Histograms
//             TH1* hAIDA_EnergyImpSum;
//             TH1* hAIDA_EnergyDecSum;
//             TH1* hAIDA_EnergyDecHitsX;
//             TH1* hAIDA_EnergyImpHitsX;
//             TH1* hAIDA_EnergyDecHitsY;
//             TH1* hAIDA_EnergyImpHitsY;

			//Fatima histograms
			//-general
			TH1* hFAT_Esum;
			TH2* hFAT_gg;
            TH1* hFAT_TDCdt_ref_sum;
			TH1* hFAT_TDCdt_refcalib_sum;
            TH1* hFAT_TDCdt_ref0_sum;
			TH1* hFAT_QDCdtsum;
			TH1* hFAT_TDCdtsum_ref_gated;
            TH1* hFAT_TDCdtsum_ref0_gated;
			TH1* hFAT_QDCdtsum_ref_gated;   //for now...
			TH1* hFAT_Angular_Diff_ref_gated; // Histogram of Gated Angular Differences
			//-statistics
			TH1* hFAT_hits;		     //number of hits per detector id
			TH1* hFAT_hits_QDC;
			TH1* hFAT_hits_TDC;
            TH1* hFAT_TDC_multi;

			TH2* hFAT_QDC_TDC_hitmap; //hits of qdc and tdc in same event
			TH2* hFAT_correlations;   //det-det coincidence map
			//-energy
			TH1* hFAT_E[40];

			TH1*  hFAT_Traw[40];
			TH1*  hFAT_Eraw_VME[40];

            TH1*  hFAT_Eraw[40];

			TH1* hFAT_TDC_multich[48];

			TH2** hFAT_E_ratio;
			TH2** hFAT_gg_ref;
			//-timing
			TH1* hFAT_TDCdt_ref[48];
			TH1* hFAT_TDCdt_refcalib[48];
			TH1* hFAT_TDCdt_ref0[48];
			TH1** hFAT_QDCdt_ref;
			TH2** hFAT_TDC_QDC_dt;
			TH1*  hFAT_TDCdt_ref_gated[48];
			TH1*  hFAT_TDCdt_ref0_gated[48];
			TH2** hFAT_E_TDCdt_ref_gated;

			//Other histograms
			TH1* WR_HIST;
			TH1* WR_HIST2;
			TH1* C_t;
			TH1** tamex_Mult_lead;
			TH1** tamex_Mult_trail;

			TH1*** mat;
			TH1* all;
			TH1* all2;

			TH1* WR_F;

            TH1 *hPLAS_QDCRaw1[32];
            TH1 *hPLAS_QDCRaw2[32];
            TH1 *hPLAS_TDCRaw[32];
            TH1 *hPLAS_TimeDiff[32]; // Time difference Histogram from TDC channels
            TH1 *hPLAS_TimeDiffSiPM[32];
            TH1 *hPLAS_TimeDiffCalib[32]; // Time difference Histogram from TDC channels Calibrated
            TH2 *hPLAS_CoincQ1A1; //Coincidences with germanium
            TH2 *hPLAS_CoincQ1Q2Add; //Coincidences Addback
            TH2 *hPLAS_CoincQ1Q2[32];
            TH2 *hPLAS_EETimeGated[32];
            TH1 *hPLAS_TimeD1;
            TH1 *hPLAS_TimeD2;
            TH1 *hPLAS_QDC1_hits;
            TH1 *hPLAS_TDC_hits;
            TH1 *hPLAS_TDC_multich[32];
            TH1 *hPLAS_QDC1_multi;
            TH1 *hPLAS_TDC_multi;
            TH2 *hPLAS_CoincQ1T1[32]; //Coincidences with germaniun
            TH2 *hPLAS_TempCoincQ1T1[32]; //Coincidences with germaniun
            TH1 *hPLAS_QDCGain[32];
            TH1 *hPLAS_TDCGain[32];
            TH1 *hPLAS_TDC_FATgated[32];
            TH1 *hPLAS_QDCRaw1_FATgated[32];
            TH1 *hPLAS_TimeDiff_FATgated[32];
            TH1 *hPLAS_QDCCalib1_FATgated[32];

            TH1 *hScalar_hit_pattern;

            TH2 *hCoinc_FAT_PLAS_TDetc;
            TH2 *hCoinc_FAT_PLAS_ESum;
            TH2 *hCoinc_FAT_PLAS_T;
            TH2 *hCoinc_FAT_PLAS_E[32];
            TH1 *hCoinc_FAT_PLAS_TDiff[32];
            TH1 *hCoinc_FAT_PLAS_TDiff_gg[32];
          //  TH1 *hCoinc_FATID_PLAS_TDiff[34];

            TH1 *hCoinc_FAT_PLAS_TDiffSum;
            TH1 *hCoinc_FAT_PLAS_TDiffSum_gg;

            TH1 *hPLAS_TDCSumCalib;
            //TH1 *hPLAS_QDCSum;
            TH1 *hPLAS_GatedHist;

			//Plastic histograms
			//TH1*** tamex_Mult_Ch_lead;
			//TH1*** tamex_Mult_Ch_trail;
			//TH2** tamex_mult_mat_lead;
			//TH2** tamex_mult_mat_trail;
			//TH1*** Trail_LEAD;
			TH1**** TRAIL_TRAIL;
			TH1**** LEAD_LEAD;

			//TH1*** Coarse;
			//TH1** DIFF_ARR;
			TH1**** TOT_TOT;
			TH1*** TOT_Single;
			//TH1*** LEAD_LEAD_Total;

            TH1 *hFATlead_Coarse[52];
            TH1 *hFATlead_Fine[52];
            TH1 *hFATtrail_Coarse[52];
            TH1 *hFATtrail_Fine[52];



			TH1 *hFin_ToT[33];
			TH1 *hlead_lead[48];

			TH1* FAT_TDC_Diff; // ****NEWLY ADDED****

			//for the SIS modules

			// GALILEO Histograms //
			TH1 *hGAL_Raw_E[32];

		private:
        // typedef std::unordered_map<int, TTSInfo> align_t;
       // align_t aligns_0;
       // align_t aligns_1;
       // align_t aligns_3;


//             long adcLastTimestamp[12][4];
//             int adcCounts[12][4];
//             void ResetMultiplexer();
//             void CorrectTimeForMultiplexer(AidaEvent& evt);

            int AIDA_Hits=0;
            double AIDA_Energy[10000] = {0};
            int AIDA_FEE[10000] = {0};
            int AIDA_ChID[10000] = {0};
            ULong64_t AIDA_Time[10000] = {0};
            bool AIDA_HighE_veto[10000] = {false};
            int AIDA_Side[10000] = {0};
            int AIDA_Strip[10000] = {0};
            int AIDA_evtID[10000] = {0};

			const int FATIMA_reference_det = 0;
			const int FAT_MAX_DET = 50;

			const int FRS = 0;
			const int AIDA = 1;
			const int PLASTIC = 2;
			const int FATIMA = 3;
			const int GALILEO = 4;
			const int FINGER = 5;

			int FAT_REF_DET;


			float E_gate1,E_gate2;

			Bool_t ffill;
			Int_t fshift;
			ULong64_t White_Rabbbit_old;

			Int_t PrcID_Array[10][5];
			bool Used_Systems[10];

//             int tamID[4][100];
//             int tamCH[4][100];
            int fingID[4][16];
            
//             int TamID[2][256];
//             int TamCh[2][256];
            int TAMEX_bPlasFat_ID[2][256];

			bool SKIP_EVT_BUILDING;

			bool PADI_OR_PADIWA,VME_TAMEX_bPlas,VME_TAMEX_Fatima, VME_AND_TAMEX_Fatima;

        ///For AIDA

            long lastTime;
            int ID;
            AidaEvent old;
            AidaEvent evt;

            long startTime;
            long stopTime;

      /* Multiplexer correction */
						std::vector<std::array<uint64_t, 4>> adcLastTimestamp;
						std::vector<std::array<int, 4>> adcCounts;
            void ResetMultiplexer();
            void CorrectTimeForMultiplexer(AidaEvent& evt);

            int totalEvents;
            int implantEvents;
            int decayEvents;
            int pulserEvents;
            int nonsenseEvents;

						// AIDA histograms
						std::vector<std::array<std::array<TH1*, 2>, 64>>	hAIDA_ADC;
						TH2* hAIDA_ADC_unaligned;
						TH2* hAIDA_ADC_aligned;

      ///End AIDA
			double vals[100000];
			int val_it;
            int event_number;
			string input_data_path;
			string input_data_path_old;
            string test;

			bool cals_done,WR_used;
			bool WHITE_RABBIT_USED; // Read from General Setup File
			bool FAT_make_raw_histograms;


			int file_pwd, file_end;
			std::string gain_match_filename;
			int data_file_number = 0;

			Detector_System** Detector_Systems;
			Data_Stream** data_stream;
			White_Rabbit* WR;
			Raw_Event* RAW;
           // AIDA_Event* Aida_inp;
			EventBuilder** EvtBuilder;
            AidaEvent CallTheAida();

			int amount_interest;
			int* length_interest;
			int** interest_array;
            int FAT_GamGatelow;
            int FAT_GamGatehigh;

			//Event_Builder** EvtBuilder;

			double fatima_E_save[4];
			int am_FATIMA_hits;
			int num_full_FAT_evts;


			Int_t get_Conversion(Int_t);
			void get_used_Systems();
			void get_WR_Config();

			void read_setup_parameters();

            void load_FingerID_File();
            void load_FatTamex_Allocationfile();
			void load_PrcID_File();
			void get_interest_arrays();

			void Make_FRS_Histos();
			void Fill_FRS_Histos();

            void Make_AIDA_Histos();
            void Fill_AIDA_Histos();

			void Make_Plastic_Histos();
			void Fill_Plastic_Histos();

            void Make_Plastic_VME_Histos();
            void Fill_Plastic_VME_Histos();

			void Make_FATIMA_Histos();
			void Fill_FATIMA_Histos();

            void Make_FATIMA_TAMEX_Histos();
            void Fill_FATIMA_TAMEX_Histos();

            void Make_FATIMA_VME_TAMEX_Histos();
            void Fill_FATIMA_VME_TAMEX_Histos();

			void Make_GALILEO_Histos();
			void Fill_GALILEO_Histos();

            void Make_Finger_Histos();
            void Fill_Finger_Histos();

			void FILL_HISTOGRAMS(int);


			bool Check_Cal_Plastic();

			void checkPADI_or_PADIWA();
			void checkTAMEXorVME();

			bool PLASTIC_CALIBRATION;



// 			void FAT_det_pos_setup();
// 			double distance_between_detectors(double, double, double, double, double, double);
// 			double angle_between_detectors(double, double, double);

			void print_MBS(int*,int);

			int count;
            int array_count;
			int called[2];
			int iterator;
			int Cout_counter;

			ULong64_t WR_tmp;
            ULong64_t WR_main;
            int WR_count;
            int WR_d;
            ULong64_t WR_AIDA[10000];
            Long64_t  WR_diff[10000];
            int AIDA_Loop;
			ClassDef(EventUnpackProc,1)
	};

#endif
