function [peaks lats] = fPickPeakWin(data,times,win,type,area)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [peaks lats] = fPickPeakWin(data,times,win,type,area)
% 
% Picks peaks! If no peak can be identified in your specified time window,
% the peak returned will be the maximum or minimum value in the window.
% 
% Inputs: 
%       data = M x N matrix (M = electrodes, N = timepoints)
%       times = time values of N
%       win = [startwindow stopwindow] to find peak relative to times
%       type = 'max' or 'min' i.e. positive or negative peak to indentify
%       area = number of timepoints to average over (before and after peak)
%
%  Outputs: 
%       peaks = M x 1 vector of peaks (each value corresponding to an electrode)
%       lats = M x 1 vector of latencies (each value corresponding to an electrode)
%
% T.Zanto
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

tindx = find(times >= win(1) & times <= win(2));


%% Get peaks


for chan = 1:size(data,1)
    if type == 'min'
        temppeak = 1000; %some large number
        checkpeak = 'data(chan,tindx(t)) < data(chan,tindx(t)-1) && data(chan,tindx(t)) < data(chan,tindx(t)+1)';    
        checkpeak2 = 'data(chan,tindx(t)) < temppeak';    
        checkpeak3 = 'temppeak > data(chan,tindx(latindx)+adj+1)  || temppeak > data(chan,tindx(latindx)+adj-1)';
    elseif type == 'max'
        temppeak = -1000; %some small number
        checkpeak = 'data(chan,tindx(t)) > data(chan,tindx(t)-1) && data(chan,tindx(t)) > data(chan,tindx(t)+1)';
        checkpeak2 = 'data(chan,tindx(t)) > temppeak';
        checkpeak3 = 'temppeak < data(chan,tindx(latindx)+adj+1)  || temppeak < data(chan,tindx(latindx)+adj-1)';
    else
        error('Variable ''type'' must be either ''max'' or ''min''...jerky')
    end
    for t = 2:length(tindx)-1
        if eval(checkpeak)
            if eval(checkpeak2)
                temppeak = data(chan,tindx(t));
                tempindx = tindx(t);
            end
        end
    end
    if ~exist('tempindx') %if no peak identified, take maximum or minimum value
        eval(sprintf('indx = find(data(chan,tindx) == %s(data(chan,tindx)));',type));
        tempindx = tindx(indx);
    end
    tempindx = tempindx(1); %use the first value (in case the peak has more than 1 value)
    lats(chan) = times(tempindx);
    peaks(chan) = mean(data(chan,tempindx-area:tempindx+area));
    clear tempindx temppeak
end


