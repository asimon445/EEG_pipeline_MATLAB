function signal = flt_PD_to_numeric(varargin)
%

if ~exp_beginfun('filter') return; end

% used as a tool to select channel subsets before these ops are applied
declare_properties('name',{'EventConversion','events'}, 'precedes',{'flt_laplace','flt_ica','flt_reref'}, 'independent_trials',true, 'independent_channels',true);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}));

for evc=1:length(signal.event)
    try
    event=dec2bin(str2double(signal.event(evc).type),16);
        %event=dec2bin(str2double(evs(evc)),32);

    event=event(end-15:end-8); %This eliminates higher values
    event(5:8)='0000'; %This eliminates lower values
    signal.event(evc).type=bin2dec(event);
    catch
    end
end

exp_endfun;
