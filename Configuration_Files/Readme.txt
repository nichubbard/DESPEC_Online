DESPEC-Analysis-Frameworkv1.3 Notes

The code was tested with the data: 252_Co60_fp_trigwind2us (FATIMA-bPlastic) and fpg_Co60 (FATIMA-bPlastic-Germanium)

Inputs:
(i) Detectors: FATIMA, bPlastic(VME), Galileo; change between different setups in /Configuration_Files/Used_Systems.txt (The code will crash if the wrong detector systems are input).

(ii) FATIMA gamma gate (lower and upper energy range) is given in /Configuration_Files/Correlations.dat

(iii) Energy and time calibrations: Gains files located in /Configuration_Files/
Plastic_VME_Energy_Calibration.dat
Plastic_VME_Time_Calibration.dat
FATIMA_Energy_Calibration.txt
FATIMA_Time_Calibration.txt

Raw and calibrated FATIMA/bPlastic time and energy spectra in respective FATIMA/ and PLASTIC/ folders
'Info' in the histogram browser gives details. 
For FATIMA, the gamma gated time is given in the FATIMA/Timing/ folder.

FATIMA-bPlastic coincidences are in  /Correlations/FATIMA_PLASTIC_General/
The FATIMA (gamma) gated bPlastic spectra are in /Correlations/FATIMA_PLASTIC_GamGate/

bPlastic reference channel is 0
FATIMA reference channel is 40 

Next version will include: 
-AIDA 
-Walk correction
-GALILEO correlations (check plastic VME time data stream...)
-root tree output
