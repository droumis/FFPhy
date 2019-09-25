function out = calcripstats_new(index, excludetimes, rip, varargin)
% function out = calcripstats(index, excludetimes, rip, varargin)
%
%  Loads the rip structure and computes the average CA1 spectrum, the
%  average CA3 spectrum, the correlation between gamma power and ripple
%  power, and the cross correlation between gamma power and ripple power.
%

% set options
appendindex = 1;
min_separation = 1;

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'min_separation'
                min_separation = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end
rip = rip{index(1)}{index(2)};
out.starttime = [];
out.endtime = [];
out.time = [];
out.ca1_ripple_power = [];
out.peak = [];
out.ca1_power = [];
out.ca3_power = [];
out.ca1_ca1_coherence = [];
out.ca1_ca3_coherence = [];
out.ca3_ca3_coherence = [];
out.ca1_power_corr = [];
out.ca3_power_corr = [];
out.ca1_ca1_coherence_corr = [];
out.ca1_ca3_coherence_corr = [];
out.ca3_ca3_coherence_corr = [];
%Apply excludetimes
if ~isempty(rip)
    if ~isempty(rip.starttime)
        included = ~isExcluded(rip.starttime, excludetimes);
        valid = [1000; rip.endtime(2:end) - rip.starttime(1:end-1)];
        valid = valid > min_separation;
        included = included & valid;
        out.starttime = rip.starttime(included);
        out.endtime = rip.endtime(included);
        out.time = rip.time;

        %Determine peak of CA1 ripple power
        out.ca1_ripple_power = rip.ca1_ripple_power(:,included);
        out.peak = max(out.ca1_ripple_power);
        peak_ind = zeros(size(out.peak));
        for i = 1:length(out.peak)
            peak_ind(i) = find(out.peak(i) == out.ca1_ripple_power(:,i));        
        end

        %Initialize CA1 and CA3 gamma power, CA1-CA1, CA1-CA3, & CA3-CA3 coherence
        out.ca1_power = rip.ca1_gamma_power(:,included);
        out.ca3_power = rip.ca3_gamma_power(:,included);
        out.ca1_ca1_coherence = rip.ca1_ca1_gamma_coherence(:,included);
        out.ca1_ca3_coherence = rip.ca1_ca3_gamma_coherence(:,included);
        out.ca3_ca3_coherence = rip.ca3_ca3_gamma_coherence(:,included);

        %Exclude invalid ripples: max ripple power > 300ms away from starttime
        valid = peak_ind>26 & peak_ind<66;
        out.starttime = out.starttime(valid);
        out.endtime = out.endtime(valid);
        out.ca1_ripple_power = out.ca1_ripple_power(:,valid);
        out.peak = out.peak(:,valid);
        peak_ind = peak_ind(valid);
        out.ca1_power = out.ca1_power(:,valid);
        out.ca3_power = out.ca3_power(:,valid);
        out.ca1_ca1_coherence = out.ca1_ca1_coherence(:,valid);
        out.ca1_ca3_coherence = out.ca1_ca3_coherence(:,valid);
        out.ca3_ca3_coherence = out.ca3_ca3_coherence(:,valid);

        %Compute correlation between gamma power and coherence & ripple power
        tmp1 = zeros(10,length(peak_ind)); tmp3 = zeros(size(tmp1));
        tmp11 = zeros(size(tmp1)); tmp13 = zeros(size(tmp1)); tmp33 = zeros(size(tmp1));   
        for i = 1:length(peak_ind)
            tmp1(:,i) = out.ca1_power([6 16 26 36 46 56 66 76 86 peak_ind(i)],i);
            tmp3(:,i) = out.ca3_power([6 16 26 36 46 56 66 76 86 peak_ind(i)],i);
            tmp11(:,i) = out.ca1_ca1_coherence([6 16 26 36 46 56 66 76 86 peak_ind(i)],i);
            tmp13(:,i) = out.ca1_ca3_coherence([6 16 26 36 46 56 66 76 86 peak_ind(i)],i);
            tmp33(:,i) = out.ca3_ca3_coherence([6 16 26 36 46 56 66 76 86 peak_ind(i)],i);
        end
        out.ca1_power_corr = nan(size(tmp1,1),4);
        out.ca3_power_corr = nan(size(tmp3,1),4);
        out.ca1_ca1_coherence_corr = nan(size(tmp11,1),4);
        out.ca1_ca3_coherence_corr = nan(size(tmp13,1),4);
        out.ca3_ca3_coherence_corr = nan(size(tmp33,1),4);

        if size(tmp1,2) >10
            for i = 1:size(tmp1,1)
                [r p rL rU] = corrcoef(out.peak,tmp1(i,:));
                out.ca1_power_corr(i,:) = [r(1,2) p(1,2) rL(1,2) rU(1,2)];

                [r p rL rU] = corrcoef(out.peak,tmp3(i,:));
                out.ca3_power_corr(i,:) = [r(1,2) p(1,2) rL(1,2) rU(1,2)];

                [r p rL rU] = corrcoef(out.peak,tmp11(i,:));
                out.ca1_ca1_coherence_corr(i,:) = [r(1,2) p(1,2) rL(1,2) rU(1,2)];

                [r p rL rU] = corrcoef(out.peak,tmp13(i,:));
                out.ca1_ca3_coherence_corr(i,:) = [r(1,2) p(1,2) rL(1,2) rU(1,2)];

                [r p rL rU] = corrcoef(out.peak,tmp33(i,:));
                out.ca3_ca3_coherence_corr(i,:) = [r(1,2) p(1,2) rL(1,2) rU(1,2)];
            end
        else
            out.ca1_power_corr = [];
            out.ca3_power_corr = [];
            out.ca1_ca1_coherence_corr = [];
            out.ca1_ca3_coherence_corr = [];
            out.ca3_ca3_coherence_corr = [];
        end
    end
end
if appendindex
    out.index = index;
end

end