function [ variance ] = fErpVar( data,times,win )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% variance = fPickPeakWin(data,times,win)
% 
% Finds the standard deviation of ERPs within a specified time window
% 
% Inputs: 
%       data = M x N matrix (M = electrodes, N = timepoints)
%       times = time values of N
%       win = [startwindow stopwindow] to find peak relative to times
%
%  Outputs: 
%       variance = M x 1 vector of standard deviation values (each value corresponding to an electrode)       
%
% AJ Simon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


tindx = find(times >= win(1) & times <= win(2));

for chan = 1:size(data,1)
    variance(chan,1)=mean(data(chan,tindx(1):tindx(end)));
end

