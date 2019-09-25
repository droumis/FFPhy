function [out ps] = runglm(g, colindex, colexcl, colout)
%out = runglm(g, colindex, colexcl, colout)
% g = matrix that includes colindex, columns to include in X of GLM and
%   colexcl
% colindex: columns identifying the index of each row, to separate into
%    sessions
% colexcl: columns to exclude from data in X.  also automatically excludes
%   colindex
% colout: identifies column of Y

ind = unique(g(:,colindex), 'rows');
ps=[];
for i = 1:size(ind,1) %for each epoch
    
    %get data for session
    data = g(ismember(g(:,colindex), ind(i,:), 'rows'),:);
    data = data(ismember(data(:,colout), [0 1]),:); %only include binaryoutcomes
    
    %get spiking data
    nancol = find(sum(isnan(data))./size(data,1)==1);
    nancol = nancol(~ismember(nancol, colexcl));
    if ~isempty(nancol)
        tdata = data(:,1:min(nancol)-1);
    else
        tdata = data;
    end
    colincl = 1:size(tdata,2);
    colincl = colincl(~ismember(colincl, [colindex colexcl]));
    %get rid of cases with NaN in ripspikes
    [exclr junk] = find(isnan(tdata(:,colincl)));
    rows = 1:size(tdata,1);
    rows = setdiff(rows, exclr);
    tdata = tdata(rows,:);  %each row is a trial, each column is one cell, value = spikes
    
    %get trial outcome data
    outcome = data(rows,colout);
    
    %run regression
    [b,dev, stats]  = glmfit([tdata(:,colincl)],outcome,'binomial');
    ps = [ps; repmat(ind(i,:), size(stats.p)) stats.p];
    
%     figure(1)
%     plot(sum(logical(tdata(:,colincl)),2),outcome, '*')
%     hold on
%     lsline
%     pause
%     clf(1)
    
    out{i}.index = ind(i,:);
    out{i}.b = b;
    out{i}.dev = dev;
    out{i}.stats = stats;
    out{i}.numcell = length(colincl);
    out{i}.numtrials = size(tdata,1);
    
end

