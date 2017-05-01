function [out] = dfa_createbehaveperform(index, excludeperiods, rewardinfo,task, varargin) %oppseqrewardinfo



for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'learningtype'
            learningtype = varargin{option+1}; %options are "S1Acquisition" and "Switch"
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end


reward = rewardinfo{index(1)}{index(2)};

out.outboundperf = reward(reward(:,4)==11,3);
out.inboundperf = reward(reward(:,4)==10,3);
out.combinedperf = reward(:,3);
out.index =index;
out.learningtype = learningtype;
if isfield(task{index(1)}{index(2)},'sequence')
out.sequence = task{index(1)}{index(2)}.sequence;
end

if strcmp(learningtype,'Switch')
    oppseqreward = oppseqrewardinfo{index(1)}{index(2)};
    out.OPPoutboundperf = oppseqreward(oppseqreward(:,4)==11,3);
out.OPPinboundperf = oppseqreward(oppseqreward(:,4)==10,3);
out.OPPcombinedperf = oppseqreward(:,3);
end  

end