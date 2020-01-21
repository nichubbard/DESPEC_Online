#include "TFile.h"
//#include "TTree.h"
//#include "TTreeReader.h"
//#include "TTreeReaderValue.h"

// Basic root macro to display a divided TCanvas of spectra for online monitoring... 
// H.M. Albers December 2019


void SpectrumCheck(char* filename)
{
    
    TFile *f1 = TFile::Open(filename);
    if (f1 == 0) {
      printf("Error: cannot open root file\n");
      return;
    }    
    
    //Histograms to plot here...    
    TH2I *AidaImpStripXY1 = (TH2I*) f1->Get("Histograms/AIDA/Implants/DSSD1_implants_strip_XY");
    TH2I *AidaImpStripXY2 = (TH2I*) f1->Get("Histograms/AIDA/Implants/DSSD2_implants_strip_XY");
    TH2I *AidaImpStripXY3 = (TH2I*) f1->Get("Histograms/AIDA/Implants/DSSD3_implants_strip_XY");
    TH1F *AidaImpEnergy1 = (TH1F*) f1->Get("Histograms/AIDA/Implants/DSSD1_implants_energy");
    TH1F *AidaImpEnergy2 = (TH1F*) f1->Get("Histograms/AIDA/Implants/DSSD2_implants_energy");
    TH1F *AidaImpEnergy3 = (TH1F*) f1->Get("Histograms/AIDA/Implants/DSSD3_implants_energy");

    TH2I *AidaDecStripXY1 = (TH2I*) f1->Get("Histograms/AIDA/Decays/DSSD1_decays_strip_XY");
    TH2I *AidaDecStripXY2 = (TH2I*) f1->Get("Histograms/AIDA/Decays/DSSD2_decays_strip_XY");
    TH2I *AidaDecStripXY3 = (TH2I*) f1->Get("Histograms/AIDA/Decays/DSSD3_decays_strip_XY");
    TH1F *AidaDecEnergy1 = (TH1F*) f1->Get("Histograms/AIDA/Decays/DSSD1_decays_energy");
    TH1F *AidaDecEnergy2 = (TH1F*) f1->Get("Histograms/AIDA/Decays/DSSD2_decays_energy");
    TH1F *AidaDecEnergy3 = (TH1F*) f1->Get("Histograms/AIDA/Decays/DSSD3_decays_energy");
    
    TH1I *bPlasToT1 = (TH1I*) f1->Get("Histograms/bPlastic/ToT/ToT Det.  1 Ch. 0");
    TH1I *bPlasToT5 = (TH1I*) f1->Get("Histograms/bPlastic/ToT/ToT Det.  1 Ch. 4");
    TH1I *bPlasToT9 = (TH1I*) f1->Get("Histograms/bPlastic/ToT/ToT Det.  1 Ch. 8");    
    TH1I *bPlasToT13 = (TH1I*) f1->Get("Histograms/bPlastic/ToT/ToT Det.  1 Ch.12");
    
    TH1D *GalEn1 = (TH1D*) f1->Get("Histograms/GALILEO/Energy_Ch./GALILEO_E_Ch. 0");
    TH1D *GalEn2 = (TH1D*) f1->Get("Histograms/GALILEO/Energy_Ch./GALILEO_E_Ch. 1");
    TH1D *GalEn3 = (TH1D*) f1->Get("Histograms/GALILEO/Energy_Ch./GALILEO_E_Ch. 2");
    TH1D *GalEn4 = (TH1D*) f1->Get("Histograms/GALILEO/Energy_Ch./GALILEO_E_Ch. 3");
    TH1D *GalEn5 = (TH1D*) f1->Get("Histograms/GALILEO/Energy_Ch./GALILEO_E_Ch. 4");
    TH1D *GalEn6 = (TH1D*) f1->Get("Histograms/GALILEO/Energy_Ch./GALILEO_E_Ch. 5");
    
    TH1D *FatTAMEX1 = (TH1D*) f1->Get("Histograms/FATIMA_TAMEX/Energy_Calib/Energy_Calib Fat Det.  0 Ch.  0");
    TH1D *FatTAMEX2 = (TH1D*) f1->Get("Histograms/FATIMA_TAMEX/Energy_Calib/Energy_Calib Fat Det.  0 Ch.  1");
    TH1D *FatTAMEX3 = (TH1D*) f1->Get("Histograms/FATIMA_TAMEX/Energy_Calib/Energy_Calib Fat Det.  0 Ch.  2");
    TH1D *FatVME1 = (TH1D*) f1->Get("Histograms/FATIMA_VME/Energy/EnergyCalib/LaBr_ECalib_Ch. 0");
    TH1D *FatVME2 = (TH1D*) f1->Get("Histograms/FATIMA_VME/Energy/EnergyCalib/LaBr_ECalib_Ch. 1");
    TH1D *FatVME3 = (TH1D*) f1->Get("Histograms/FATIMA_VME/Energy/EnergyCalib/LaBr_ECalib_Ch. 2");
    
    TH2I *CorrAidaFRSgate1 = (TH2I*) f1->Get("Histograms/Correlations/AIDA-FRS/Implants/Z1Z2_Gate/DSSD1_implants_strip_XY_Z1Z2g");
    TH2I *CorrAidaFRSgate2 = (TH2I*) f1->Get("Histograms/Correlations/AIDA-FRS/Implants/Z1Z2_x2x4AoQ_Gate/DSSD1_implants_strip_XY_Z1Z2IDx2x4AoQg");
    TH1I *CorrAidaFRS_WR = (TH1I*) f1->Get("Histograms/Correlations/AIDA-FRS/WR_timediff");
    TH1I *CorrAidabPlas_WR = (TH1I*) f1->Get("Histograms/Correlations/AIDA-bPlast/AidaDecay_bPlas_WR_timediff");
    TH1I *CorrAidaImpDecdT = (TH1I*) f1->Get("Histograms/Correlations/AIDA-bPlast/Implant_Decay_dT");
    TH1I *CorrAidabPlas_CoinEnergy = (TH1I*) f1->Get("Histograms/Correlations/AIDA-bPlast/Aida_GatedE_bPlas_CoincEnergy");
    
    
    
    
    TCanvas *c = new TCanvas();
    c->Divide(6,6);
    // AIDA IMPLANTS
    c->cd(1); AidaImpStripXY1->Draw("col");    //gPad->SetLogz();
    c->cd(2); AidaImpStripXY2->Draw("col");    
    c->cd(3); AidaImpStripXY3->Draw("col");
    c->cd(4); AidaImpEnergy1->Draw(); AidaImpEnergy1->GetXaxis()->SetRangeUser(0,2000);
    c->cd(5); AidaImpEnergy2->Draw(); AidaImpEnergy2->GetXaxis()->SetRangeUser(0,2000);   
    c->cd(6); AidaImpEnergy3->Draw(); AidaImpEnergy3->GetXaxis()->SetRangeUser(0,2000);
    // AIDA DECAYS
    c->cd(7); AidaDecStripXY1->Draw("col"); gPad->SetLogz();
    c->cd(8); AidaDecStripXY2->Draw("col"); gPad->SetLogz();   
    c->cd(9); AidaDecStripXY3->Draw("col"); gPad->SetLogz();
    c->cd(10); AidaDecEnergy1->Draw(); AidaDecEnergy1->GetXaxis()->SetRangeUser(0,3000);
    c->cd(11); AidaDecEnergy2->Draw(); AidaDecEnergy2->GetXaxis()->SetRangeUser(0,3000);   
    c->cd(12); AidaDecEnergy3->Draw(); AidaDecEnergy3->GetXaxis()->SetRangeUser(0,3000);
    // bPlast ToTs 
    c->cd(13); bPlasToT1->Draw();
    c->cd(14); bPlasToT5->Draw();
    c->cd(15); bPlasToT9->Draw();
    c->cd(16); bPlasToT13->Draw();    
    //c->cd(17); //FINGER here?
    //c->cd(18); //FINGER here?
    // GALILEO
    c->cd(19); GalEn1->Draw();
    c->cd(20); GalEn2->Draw();
    c->cd(21); GalEn3->Draw();
    c->cd(22); GalEn4->Draw();
    c->cd(23); GalEn5->Draw();
    c->cd(24); GalEn6->Draw();
    //FATIMA
    c->cd(25); FatTAMEX1->Draw();
    c->cd(26); FatTAMEX2->Draw();
    c->cd(27); FatTAMEX3->Draw();
    c->cd(28); FatVME1->Draw();
    c->cd(29); FatVME2->Draw();
    c->cd(30); FatVME3->Draw();    
    //CORRELATIONS
    c->cd(31); CorrAidaFRSgate1->Draw("col"); 
    c->cd(32); CorrAidaFRSgate2->Draw("col");
    c->cd(33); CorrAidaFRS_WR->Draw();
    c->cd(34); CorrAidabPlas_WR->Draw();
    c->cd(35); CorrAidaImpDecdT->Draw();
    c->cd(36); CorrAidabPlas_CoinEnergy->Draw();
    
    
    

    
          
   
}
