Setting up the pipeline:

This pipeline is meant to preprocess standard 10-20 EEG data (can be collected from various data acquisition systems) using EEGLAB. A recent version of EEGLAB is included in this repository, so no need to download anything from the web. In order to use this pipeline, you must have MATLAB 2017a or more recent version installed. It may work on older versions of MATLAB, but this has not been tested. 

Before using this pipeline, you must add some files/folders from this repository to your MATLAB path. The folders you need to add are:
  1. ~/dependencies
  2. ~/eeglab14_1_2b/functions/sigprocfunc/FastICA_25

To do so, navigate (in MATLAB) to the folder that you want to add, right click it, and click "add this folder to path". Do not click "add this folder and all subfolders". 

Before running the script, the user must launch eeglab manually. To do so, navigate to "~/eeglab14_1_2b/" and type 'eeglab' in your command prompt. This will launch EEGLAB and add the necessary plugins to your path. It is important to launch EEGLAB this way, rather than adding it to your path with all subfolders, because EEGLAB calls functions that are named the same as some lower-level code that MATLAB calls, and adding both to your path confuses MATLAB.



What the pipeline does:

This pipeline will only do EEG data preprocessing. It will not segment the data, nor will it compute metrics to be analyzed. The preprocessing steps of this pipeline are as follows:
  1. Convert photodiodes from 5 digit codes to 2-3 digit codes (this makes the event markers easier to interpret and handle)
  2. Remove external electrodes 
  3. Downsample the data to 1024 Hz
  4. Remove indexed bad electrodes (if any) and apply a spherical spline interpolation from neighboring electrodes
  5. Lowpass filter at a user-specified frequency (recommend 40 Hz)
  6. Highpass filter the data above 1 Hz
  7. Do an automatic ocular correction using ICA (if user select '1' on the prompt. Else, user will have to manually identify the eye-blink ICs. 
  8. Rereference to the average signal
  
Note: automatic ocular correction is done via the function 'fDestroyBlinks'. It takes compares each independent component (IC) to the spatio-temporal properties of a typical eye blink, and the IC that most closely resembles a characteristic ocular artifact will be removed. As with any automatic data cleaning technique, it is not perfect. It works in about 80-90% of cases, depending on the data quality. It is good practice to visually inspect your data after running it through this pipeline to be sure that you are getting desirable results.  
