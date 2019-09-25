function [out] = JY_checkintersectioncrossing(index, excludetimes, data, linpos)
% checks if any intersection is skipped during linearise position
if ~isempty(linpos{index(1)}{index(2)})
    
warning('OFF','MATLAB:divideByZero');



trajmatrix=linpos{1,index(1)}{1,index(2)}.trajmatrix;

testmatrix=trajmatrix(2:end,2)-trajmatrix(1:end-1,3);

errors=trajmatrix(find(testmatrix~=0)+1,:);

% check number of trials 
out.datatrialn=size(data{1,index(1)}{1,index(2)}.Run,1);
out.linpostrialn=size(linpos{1,index(1)}{1,index(2)}.trialsegments,2);
out.result=sum(testmatrix);
out.errors=errors;

else
out.result=[];
out.errors=[];
out.datatrialn=[];
out.linpostrialn=[];

warning('ON','MATLAB:divideByZero');
end
end