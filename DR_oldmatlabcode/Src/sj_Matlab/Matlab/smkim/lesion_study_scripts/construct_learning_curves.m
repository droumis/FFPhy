load('behavperform_2');

% learning_curves is a struct array with fields
%   subject: string
%   group: string, either 'lesion' or 'control'
%   prob: Nx3 matrix, one row per each trial
%   type: 1xN cell vector, elements are strings either 'inbound' or 'outbound'
%   correct: Nx1 vector of 1 or 0
%   dayrun: Nx2 matrix, first column is day of the trial and second number is session 1/2 of the trial
%   ld_inbound: learning day for inbound, a scalar, possibly NaN if the rat did not learn
%   ld_outbound: learning day for outbound, a scalar, possibly NaN
%   lt_inbound: (cumulative) learning trial
%   lt_outbound

for s = 1:length(behavperform)
    group = behavperform(s).group;
    subject = behavperform(s).subject;
    learning_curves(s).subject = subject;
    learning_curves(s).group = group;
    learning_curves(s).conf = behavperform(s).conf;
    learning_curves(s).criterion = behavperform(s).criterion;
    load([subject '_Wtrack_journeys.mat']);
    % cumulative journey number
    i = 0;
    incount = 0;
    outcount = 0;
    for j = 1:10
        for k = 1:2
            for l = 1:numel(journeys{j}{k}.startzone)
                if (journeys{j}{k}.startzone(l) == 2) & ...
                    ~isnan(journeys{j}{k}.endzone(l)) & ... 
                    ~isnan(journeys{j}{k}.correct(l))
                    i = i+1;
                    learning_curves(s).dayrun(i,:) = [j k];
                    learning_curves(s).correct(i) = journeys{j}{k}.correct(l);
                    learning_curves(s).type{i} = 'outbound';
                    outcount = outcount+1;
                    learning_curves(s).prob(i,:) = behavperform(s).outprobcorrect(outcount,:);
                elseif ((journeys{j}{k}.startzone(l) == 1) | ...
                    (journeys{j}{k}.startzone(l) == 3)) & ...
                    ~isnan(journeys{j}{k}.endzone(l)) & ...
                    ~isnan(journeys{j}{k}.correct(l))
                    i = i+1;
                    learning_curves(s).dayrun(i,:) = [j k];
                    learning_curves(s).correct(i) = journeys{j}{k}.correct(l);
                    learning_curves(s).type{i} = 'inbound';
                    incount = incount+1;
                    learning_curves(s).prob(i,:) = behavperform(s).inprobcorrect(incount,:);
                end
            end
        end
    end
    inbound_idx = find(strcmp(learning_curves(s).type,'inbound'));
    outbound_idx = find(strcmp(learning_curves(s).type,'outbound'));
    try
        inbound_lt = inbound_idx(1+find( ...
            learning_curves(s).prob(inbound_idx,2) <= learning_curves(s).criterion,1,'last'));
    catch
        inbound_lt = NaN;
    end
    try
        outbound_lt = outbound_idx(1+find( ...
            learning_curves(s).prob(outbound_idx,2) <= learning_curves(s).criterion,1,'last'));
    catch
        outbound_lt = NaN;
    end
    if ~isnan(inbound_lt) && (inbound_lt < size(learning_curves(s).prob,1))
        learning_curves(s).lt_inbound = inbound_lt;
        learning_curves(s).ld_inbound = learning_curves(s).dayrun(inbound_lt,1);
    else
        learning_curves(s).lt_inbound = NaN;
        learning_curves(s).ld_inbound = NaN;
    end
    if ~isnan(outbound_lt) && (outbound_lt < size(learning_curves(s).prob,1))
        learning_curves(s).lt_outbound = outbound_lt;
        learning_curves(s).ld_outbound = learning_curves(s).dayrun(outbound_lt,1);
    else
        learning_curves(s).lt_outbound = NaN;
        learning_curves(s).ld_outbound = NaN;
    end

end


