
function [EEG, compnum] = fDestroyBlinks(EEG,autorem)
%UNTITLED Summary of this function goes here
%   Write some cool shit here
%
%
% 
%

%autoremoval set to 'on'
if nargin==1 || autorem==1
    
    temp=EEG;
    temp.data=EEG.icaact;
    temp=pop_eegfiltnew(temp,2,9);   %bandpass filter ICA activations to increase SNR in eye blinks
    
    blcomp_ix=zeros(size(temp.data,1),1);
    compnum=[];
    
    for pc=1:size(temp.data,1)    %loop through components
        
        bl_ix_p=0; bl_ix_n=0; bl_ix=0;
        ix_p=0; pospeak=[];
        ix_n=0; negpeak=[];
        
        winvar=std(temp.data(pc,:));        %find standard deviation of component
        winavg=mean(temp.data(pc,:));       %find average of component
        
        for t=3:length(temp.data)      %find local peaks, can't find a peak on first data point    
            if temp.data(pc,t)<temp.data(pc,t-1) && temp.data(pc,t-1)>temp.data(pc,t-2)  %p-1 is a positive peak         
                ix_p=ix_p+1;
                pospeak(1,ix_p)=temp.data(pc,t-1);   %index positive peak      
            elseif temp.data(pc,t)>temp.data(pc,t-1) && temp.data(pc,t-1)<temp.data(pc,t-2)   %p-1 is a negative peak              
                ix_n=ix_n+1;
                negpeak(1,ix_n)=temp.data(pc,t-1);   %index negative peak
            end
        end
        
        % find and index significantly large peaks
        for p=1:length(pospeak)
            if pospeak(p) >= (winavg+(3.5*winvar))     
                bl_ix_p=bl_ix_p+1;
                bl_ix=bl_ix+1;
            end
        end
        
        for n=1:length(negpeak)
            if negpeak(n) <= (winavg-(3.5*winvar))
                bl_ix_n=bl_ix_n+1;
                bl_ix=bl_ix+1;
            end
        end 
        blcomp_ix(pc,1)=bl_ix;
        blcomp_ix(pc,2)=pc;
        blcomp_ix(pc,3)=bl_ix_p;
        blcomp_ix(pc,4)=bl_ix_n;
    end
    
    [s,I]=sort(blcomp_ix,'descend');    %sort components from most to least number of large peaks (possible blinks)
    
    id=0; c=0;
    while id==0; %loop through the components ordered from most to least peaks in 2-9 hz range and find the first one that's activation in FPZ>FZ>CZ
        c=c+1;
        
        if c>size(temp.data,1)
            
            idd=0; d=0;
            while idd==0
                d=d+1;
                
                if d>size(temp.data,1)
                    fprintf('no component removed \n');
                    compnum=0;
                    id=1; idd=1;
                else
                    
                    posblix=blcomp_ix(I(d,1),3); negblix=blcomp_ix(I(d,1),4);
                    
                    %index strength of ICA activation in frontal electrodes
                    FPZ=EEG.icawinv(33,I(d,1)); FP1=EEG.icawinv(1,I(d,1)); FP2=EEG.icawinv(1,I(d,1));
                    FZ=EEG.icawinv(37,I(d,1)); F3=EEG.icawinv(5,I(d,1)); F4=EEG.icawinv(40,I(d,1));
                    CZ=EEG.icawinv(48,I(d,1)); C3=EEG.icawinv(13,I(d,1)); C4=EEG.icawinv(50,I(d,1));
                    
                    Fronts=(FPZ+FP1+FP2)/3;
                    Midfr=(FZ+F3+F4)/3;
                    Mid=(CZ+C3+C4)/3;
                    
                    if posblix>negblix   %if blink related component has a positive polarity
                        if (Fronts>Midfr) && (Midfr>Mid)  % if component strength is highest in far frontal elecs and decreases as it moves more posterior                            
                            compnum=I(d,1); id=1; idd=1;
                            EEG=pop_subcomp(EEG,compnum,0);   %delete component from data
                            fprintf('removing component %d \n',I(d,1));
                        end         
                    elseif negblix>posblix
                        if (Fronts<Midfr) && (Midfr<Mid)
                            compnum=I(d,1); id=1; idd=1;
                            EEG=pop_subcomp(EEG,compnum,0);   %delete component from data
                            fprintf('removing component %d \n',I(d,1));               
                        end
                    end                  
                end
            end
        else
            posblix=blcomp_ix(I(c,1),3); negblix=blcomp_ix(I(c,1),4);
            
            %index strength of ICA activation in frontal electrodes
            FPZ=EEG.icawinv(33,I(c,1)); FP1=EEG.icawinv(1,I(c,1)); FP2=EEG.icawinv(1,I(c,1));
            FZ=EEG.icawinv(37,I(c,1)); F3=EEG.icawinv(5,I(c,1)); F4=EEG.icawinv(40,I(c,1));
            CZ=EEG.icawinv(48,I(c,1)); C3=EEG.icawinv(13,I(c,1)); C4=EEG.icawinv(50,I(c,1));
            
            Fronts=(FPZ+FP1+FP2)/3;
            Midfr=(FZ+F3+F4)/3;
            Mid=(CZ+C3+C4)/3;
            FPdiff=abs(FP2-FP1); Fdiff=abs(F4-F3); Cdiff=abs(C4-C3);
            FrontLeftdiff=abs(F3-FP1); FrontRightdiff=abs(F4-FP2);
            MidLeftdiff=abs(C3-F3); MidRightdiff=abs(C4-F4);
            
            if posblix>negblix   %if blink related component has a positive polarity
                if (Fronts>Midfr) && (Midfr>Mid)  % if component strength is highest in far frontal elecs and decreases as it moves more posterior
                    
                    if FrontLeftdiff > FPdiff && FrontLeftdiff > Fdiff && FrontRightdiff > FPdiff && FrontRightdiff > Fdiff ...   %check that differences in component strength are greater in the caudal-rostral plane than the saggital plane
                            && MidLeftdiff > Fdiff && MidLeftdiff > Cdiff && MidRightdiff > Fdiff && MidRightdiff > Cdiff
                        compnum=I(c,1); id=1;
                        EEG=pop_subcomp(EEG,compnum,0);   %delete component from data
                        fprintf('removing component %d \n',I(c,1));
                    end
                end
            elseif negblix>posblix
                if (Fronts<Midfr) && (Midfr<Mid)
                    if FrontLeftdiff > FPdiff && FrontLeftdiff > Fdiff && FrontRightdiff > FPdiff && FrontRightdiff > Fdiff ...   %check that differences in component strength are greater in the caudal-rostral plane than the saggital plane
                            && MidLeftdiff > Fdiff && MidLeftdiff > Cdiff && MidRightdiff > Fdiff && MidRightdiff > Cdiff
                        compnum=I(c,1); id=1;
                        EEG=pop_subcomp(EEG,compnum,0);   %delete component from data
                        fprintf('removing component %d \n',I(c,1));
                    end
                end
            end
        end 
    end
    clear temp
    
% manual eye blink removal
else 
    w=0;
    while w==0
        prompt={'View topoplots of components (1) or time series of components (2)?'};
        name = 'Visualization Method';
        defaultans = {'1'};
        method = inputdlg(prompt,name,[1 60],defaultans);
        method = cell2mat(method);
        method = str2num(method);

        if method==1
            w=1;
            %visualize
            last=size(EEG.icaact,1);
            pop_topoplot(EEG,0, [1:last] ,' resampled',[8 8] ,0,'electrodes','on');

            prompt={'Which component most closely resembled an eyeblink?'};
            name = 'Pick an eye blink component';
            defaultans = {''};
            compnum = inputdlg(prompt,name,[1 64],defaultans);
            compnum = cell2mat(compnum);
            compnum = str2num(compnum);
            
            %remove component
            EEG=pop_subcomp(EEG,compnum,0);
            close;

        elseif method==2
            w=1;
            figure;
            for c = 1:size(EEG.icaact,1)
                plot(EEG.data(c,:));
                comp=sprintf('Component: %d',c);
                title(comp,'FontSize',18);
                pause;
            end

            prompt={'Which component most closely resembled an eyeblink?'};
            name = 'Pick an eye blink component';
            defaultans = {''};
            compnum = inputdlg(prompt,name,[1 64],defaultans);
            compnum = cell2mat(compnum);
            compnum = str2num(compnum);
            
            %remove component
            EEG=pop_subcomp(EEG,compnum,0);

            close;   %close figure
        else
            disp('incorrect choice');
        end
    end
end
end



