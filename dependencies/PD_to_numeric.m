function signal = PD_to_numeric(varargin)
%
% arg_define(varargin, ...
% arg_norep({'signal','Signal'}));

signal=varargin{1,1};

for evc=1:length(signal.event)
    try
        if ischar(signal.event(evc).type)
            event=dec2bin(str2double(signal.event(evc).type),16);
        elseif isnumeric(signal.event(evc).type)
            event=dec2bin(signal.event(evc).type,16);
        end
        %event=dec2bin(str2double(evs(evc)),32);
        
        event=event(end-15:end-8); %This eliminates higher values
        event(5:8)='0000'; %This eliminates lower values
        signal.event(evc).type=bin2dec(event);
    catch
    end
end
end

