function out = JY_filterspikesbytrial(index, excludetimes, spikes,data)
% classifies each spike as occurring during a trial or between trials
% (intertrial)

if ~isempty(data{index(1)}{index(2)}.Run)
    
    % get trial times and intertrial trial times
    trialstartend=data{index(1)}{index(2)}.Run(:,3:4);
    
    intertrialt=[];
    for i=1:size(data{index(1)}{index(2)}.Run,1)-1;
        intertrial=[];
        intertrial(1,1)=data{index(1)}{index(2)}.Run(i,4)+1;
        intertrial(1,2)=data{index(1)}{index(2)}.Run(i+1,3)-1;
        intertrialt=[intertrialt;intertrial];
    end
    
    epochind = index(:,1:2);
    epochind = unique(epochind,'rows');
    if (size(epochind,1) ~= 1)
        error('Indices must have only one unique day/epoch');
    end
    
    for i = 1:size(index,1)
        if ~isempty(spikes{index(i,1)}{index(i,2)}{index(i,3)}{index(i,4)})
            if ~isempty(spikes{index(i,1)}{index(i,2)}{index(i,3)}{index(i,4)}.data)
                out = spikes{index(i,1)}{index(i,2)}{index(i,3)}{index(i,4)};
                out.data = out.data(find(~isExcluded(out.data(:,1), ...
                    excludetimes)),:);
                
            end
        end
    end
    
    % calculate which trial or intertrial interval the spike belongs to
    
    % trial
    out.trialspike=[];
    out.intertrialspike=[];
    out.trialtime=trialstartend;
    out.intertrialtime=intertrialt;
    
    for i=1:size(out.data,1)
        trialn=0;
        intertrialn=0;
        trialspike=isExcluded(out.data(i,1)*10000,trialstartend);
        if ~isempty(find(trialspike==1))
            trialn=find(sum(trialstartend>=out.data(i,1)*10000,2)==1);
        end
        out.trialspike(i,1)=trialn;
        
        intertrialspike=isExcluded(out.data(i,1)*10000,intertrialt);
        if ~isempty(find(intertrialspike==1))
            
            intertrialn=find(sum(intertrialt>out.data(i,1)*10000,2)==1);
            
        end
        
        out.intertrialspike(i,1)=intertrialn;
    end
    
    out.uniquetrialspike=unique(out.trialspike);
    out.uniqueintertrialspike=unique(out.intertrialspike);
else
    
    out.data=[];
    out.trialspike=[];
    out.uniquetrialspike=[];
    out.intertrialspike=[];
    out.uniqueintertrialspike=[];
end