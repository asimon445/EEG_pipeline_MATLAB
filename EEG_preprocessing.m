
%========================== PREPROCESSING STEPS ===========================%
%%% This script will completely preprocess your EEG data suited to your
%%% needs. Before running this script, the user should have already
%%% visually inspected his/her data to identify all files that contain bad
%%% electrodes. Before any processing begins, the user must list all 
%%% subjects with noisy channels and index all of those channels (this is 
%%% done so that the program can remove them interpolate data from 
%%% neighboring electrodes). 
%%%
%%% The preprocessing steps are as follows:
%%%     1. Convert photodiodes from 5 digit codes to 2-3 digit codes
%%%     2. Remove external electrodes 
%%%     3. Downsample the data to 1024 Hz
%%%     4. Remove indexed bad electrodes (if any) and apply a spherical 
%%%        spline interpolation from
%%%        neighboring electrodes
%%%     5. Lowpass filter at a user-specified frequency (recommend 40 Hz)
%%%     6. De-trend the data using a 1 Hz highpass filter
%%%     7. Do an automatic ocular correction using ICA (if user select '1'
%%%        on the prompt. Else, user will have to manually identify the
%%%        eye-blink ICs. 
%%%     8. Rereference to the average signal
%%%
%%%
%=========================== AJ Simon ====================================%


clc;clear;

%% setup the program
%load eeglab
%eeglab;

% tell it where to look for data and where to save results to
origfilepath='/Users/ajsimon/Desktop/Sample_EEG_data/raw/'; %Folder where all the data is stored
replacefilepath='/Users/ajsimon/Desktop/Sample_EEG_data/preprocessed_test/';
filetype='.bdf';   %specify file extention (.set or .bdf)

% create filelist to loop through
fileList = getAllFiles(origfilepath);
Fileidx=strfind(fileList,filetype); Fileidx=find(~cellfun(@isempty,Fileidx));   %just find the files of 'filetype' within fileList
fileList=fileList(Fileidx); %Extract just these files
ix_failures=0;

% find where eeglab is stored
tofind = 'mydir';
esctofind = regexptranslate('escape', tofind);   %in case it has special characters
dirs = regexp(path, pathsep,'split');          %cell of all individual paths

founddir=0; d=0;
while founddir==0
    d=d+1; 
    locs=strfind(dirs{d}(end-12:end),'eeglab14_1_2b');
    if ~isempty(locs)
        locpath=dirs{d};
        founddir=1;
    end
end
        
% tell it what cutoff freq to use for lowpass filtering
high_cutoff=input('Indicate lowpass cutoff frequency (type number only): ');
low_cutoff=[]; %detrending done using a different filter

%run autoamtic or manual ocular correction?
autorej=input('Would you like this program to automatically correct ocular artifacts? (yes=1/no=0): ');

% index which subjects have bad electrodes
subbad_elecs=input('List subjects that have bad electrodes (surround by brackets, empty brackets if no subjects have busted electrodes): ');

% indicate which electrodes are bad in each subject with bad electrodes
for subLen=1:length(subbad_elecs)
    
    prompt=sprintf('How many bad electrodes does subject %d have? ',subbad_elecs(subLen));
    numbad_elecs=input(prompt);
    
    for be=1:numbad_elecs
        prompt=sprintf('Enter the name of bad electrode #%d (ALL CAPS): ',be);
        thisbadelec = input(prompt,'s');
        bad_elecs{1,be} = thisbadelec;
    end
    
    %MAKE SURE THESE ELECS STAY IN THIS ORDER!!!!!
    elecs = {'Fp1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3','FC1','C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1','P3','P5','P7','P9','PO7','PO3','O1','Iz','Oz','POz','Pz','CPz','Fpz','Fp2','AF8','AF4','AFz','Fz','F2','F4','F6','F8','FT8','FC6','FC4','FC2','FCz','Cz','C2','C4','C6','T8','TP8','CP6','CP4','CP2','P2','P4','P6','P8','P10','PO8','PO4','O2'};
    elecs = upper(elecs);
    
    for i = 1:length(bad_elecs)
        for ii = 1:length(elecs)
            if strcmp(bad_elecs{i},elecs{ii});
                remove_elecs(subLen,i) = ii;
            end
        end
    end
    
    clear bad_elecs
end

%% Process
for fn=1:length(fileList)
    try
        
        % filepath and filename
        pathsep=strfind(fileList{fn},filesep);pathsep=pathsep(end); %find where the directory ends and filename starts
        filename=fileList{fn}(pathsep+1:end-length(filetype));
        
        filepath=fileList{fn}(1:pathsep);
        
        newfilepath=strrep(filepath, origfilepath, replacefilepath);
        
        if ~exist([newfilepath filename '_preprocessed.set'],'file')
            if ~exist(newfilepath,'dir');
                mkdir(newfilepath);
            end
            
            isset=strfind(filetype,'set');
            if ~isempty(isset)
                EEG=pop_loadset([filepath filename filetype]);   %use pop_loadset to load in .set files
                EEG=pop_chanedit(EEG, 'lookup',[locpath '/functions/resources/Standard-10-5-Cap385_witheog.elp']);   %get elec locations
                
            else
                EEG=pop_fileio([filepath filename filetype]);
                EEG=pop_chanedit(EEG, 'lookup',[locpath '/functions/resources/Standard-10-5-Cap385_witheog.elp']);   %get elec locations
            end
            
            %index number of unconverted diode markers
           ix_evt=0;
           if length(EEG.event)>3    %if there are 3 or fewer events, then this dataset is likely empty. Skip it. 
                for e=1:length(EEG.event)
                    if isnumeric(EEG.event(5).type)
                        if EEG.event(e).type > 1000   %if events are numbers, count how many are greater than 1000 (indicates they're unconverted)
                            ix_evt=ix_evt+1;
                        end
                    elseif ischar(EEG.event(5).type)
                        if length(EEG.event(e).type) > 3   %if events are strings, count how many are greater than 3 characters in length (indicates they're unconverted)
                            ix_evt=ix_evt+1;
                        end
                    end
                end
                
                % convert diodes if not done already
                if ix_evt > 10
                    EEG=PD_to_numeric(EEG);
                end
                
                for e=1:length(EEG.event)
                    if ~ischar(EEG.event(e).type)
                        EEG.event(e).type=num2str(EEG.event(e).type);    %convert events to strings
                    end
                end
                
                % remove external electrodes
                EEG.chanlocs=EEG.chanlocs(1:64); EEG.data=EEG.data(1:64,:); EEG.nbchan=64;
                
                % down sample
                EEG=pop_resample(EEG,1024);
                
                %lowpass filter
                EEG=pop_eegfiltnew(EEG,low_cutoff,high_cutoff);
                
                % highpass filter
                EEG=pop_eegfilt(EEG,1, 0, [], [0], 0, 0, 'fir1', 0);
                
                for r=1:length(subbad_elecs)
                    
                    current_sub=num2str(subbad_elecs(r));
                    thissub=strfind(filename,current_sub);
                    
                    if thissub==1
                        for c=1:length(remove_elecs(r,:))
                            if remove_elecs(r,c)~=0     %make sure there is an electrode in this column, not a 0
                                EEG=eeg_interp(EEG,remove_elecs(r,c));   %spherical
                            end
                        end
                    end
                end
                
                % run ICA
                EEG = pop_runica(EEG,'fastica');
                
                % ocular correction
                [EEG,compnums]=fDestroyBlinks(EEG,autorej);
                
                % index rejected components
                Rejcomp{fn,1}=filename;
                for cr=1:length(compnums)
                    Rejcomp{fn,cr+1}=compnums(cr);
                end
                
                % rereference to the average of all signals
                EEG=pop_reref(EEG, []);
                
                EEG=eeg_checkset(EEG);
                
                % save eeg
                pop_saveset(EEG,'filename',[filename '_preprocessed'],'filepath',newfilepath);
            else
                ix_failures=ix_failures+1;
                failures{ix_failures,1}=filename;
            end
        end
    catch
        ix_failures=ix_failures+1;
        failures{ix_failures,1}=filename;
    end
end

%save ocular correction info
d=datetime('Today'); d=datestr(d);
outfile=sprintf('%s/Ocularcorrection_info_%s',replacefilepath,d);
save(outfile,'Rejcomp');
