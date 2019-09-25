function [out] = JY_getreroutetrialtime(index, excludetimes, data)
% gets the duration of individual trials
% if a barrier is on

warning('OFF','MATLAB:divideByZero');
intertrialt=[];
if ~isempty(data{index(1)}{index(2)}.Run)
    
    for i=1:size(data{index(1)}{index(2)}.Run,1)-1;
        intertrial=[];
    intertrial(1,1)=data{index(1)}{index(2)}.Run(i,4)+1;
    intertrial(1,2)=data{index(1)}{index(2)}.Run(i+1,3)+1;
    intertrialt=[intertrialt;intertrial];
    end
    
out.trialduration=data{index(1)}{index(2)}.Run(:,5);
out.trialstartend=data{index(1)}{index(2)}.Run(:,3:4);
out.intertrialstartend=intertrialt;
out.ntrial=size(data{index(1)}{index(2)}.Run(:,5),1);
out.barrier=data{index(1)}{index(2)}.Run(:,6);
else
out.trialduration=[];
out.trialstartend=[];
out.intertrialstartend=[];
out.ntrial=0;
out.barrier=[];
end

