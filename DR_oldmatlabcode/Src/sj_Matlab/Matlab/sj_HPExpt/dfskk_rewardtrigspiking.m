
% this function can be used for consensus ripples

runscript = 1

if runscript == 1
    
% Animal Selection
animals = {'Egypt'};

% Epoch Filter
epochfilter = '(isequal($type, ''run''))';

% Time Filter
timefilter = { {'kk_getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'}};
        % {'kk_get2dstate', '$velocity < 4'} } ;
        
ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && isequal($type, ''interneuron''))';    %% ($meanrate < 3))';
ca2cellfilter = '(isequal($area, ''CA2'') && ($numspikes > 100) && isequal($type, ''interneuron''))';   
ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && isequal($type, ''interneuron''))';

%ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && ($meanrate < 3))';
%ca2cellfilter = '(isequal($area, ''CA2'') && ($numspikes > 100) && ($meanrate < 3))';
%ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && ($meanrate < 3))';

% Iterator
iterator = 'singlecellanal';

% Filter Creation
ca1f = createfilter('animal', animals, 'epochs', epochfilter, 'cells', ca1cellfilter, 'excludetime', timefilter, 'iterator', iterator);
ca2f = createfilter('animal', animals, 'epochs', epochfilter, 'cells', ca2cellfilter, 'excludetime', timefilter, 'iterator', iterator);
ca3f = createfilter('animal', animals, 'epochs', epochfilter, 'cells', ca3cellfilter, 'excludetime', timefilter, 'iterator', iterator);

% Set Analysis Function
ca1f = setfilterfunction(ca1f, 'dfakk_getrewardtrigspiking', {'spikes','rewardinfo'},'window',[2 2],'binsize',0.001);
ca2f = setfilterfunction(ca2f, 'dfakk_getrewardtrigspiking', {'spikes','rewardinfo'},'window',[2 2],'binsize',0.001);
ca3f = setfilterfunction(ca3f, 'dfakk_getrewardtrigspiking', {'spikes','rewardinfo'},'window',[2 2],'binsize',0.001);

% Run Analysis
ca1f = runfilter(ca1f);
ca2f = runfilter(ca2f);
ca3f = runfilter(ca3f);

end





% Collect across epochs and plot

for reg=[1:3]       % iterate over regions
    if reg==1
        region='CA1';
        clr='k';
        f=ca1f;
    elseif reg==2
        region='CA2';
        clr='g';
        f=ca2f;
    else
        region='CA3';
        clr='r';
        f=ca3f;
    end

% Collapse cells across epochs within a day

%collect all cell indices, ignoring specific epochs, into daytetcell
indices = [];
for i=1:length(f.output{1})  % iterate over epochs
    indices = [indices ; f.output{1}(i).index];
end
daytetcell = unique(indices(:,[1 3 4]),'rows');

f.celloutput = struct;


for ind=1:size(daytetcell,1)        % iterate over each unique cell

    % initialize entry
    f.celloutput(ind).daytetcell = daytetcell(ind,:);
    f.celloutput(ind).times = f.output{1}(1).times;
    f.celloutput(ind).trigmatrix = f.output{1}(1).trigmatrix;
    f.celloutput(ind).trigmatrixdescript = f.output{1}(1).trigmatrixdescript;
    f.celloutput(ind).output = cell(1,size(f.output{1}(1).trigmatrix,1));
    for m=1:length(f.celloutput(ind).output)
        f.celloutput(ind).output{m}=[];    % initialize individual output matrices
        f.celloutput(ind).rewardtimes{m}=[];    % initialize individual output matrices
    end

    %collect data across epochs
    for c=1:size(indices,1)
        if daytetcell(ind,:)==indices(c,[1 3 4])
                for m=1:length(f.celloutput(ind).output)
                    % concatenate
                    f.celloutput(ind).output{m}=[f.celloutput(ind).output{m} ; f.output{1}(c).output{m}];    
                    f.celloutput(ind).rewardtimes{m}=[f.celloutput(ind).rewardtimes{m} ; f.output{1}(c).rewardtimes{m}];    
                end
        end
    end
end
    
super(reg) = f;
    
end



% plot reward and error trials

if 1
    
    
    for reg=[2]
        
        if reg==1
            region='CA1'; clr=[0 0 0];
        elseif reg==2
            region='CA2'; clr=[0 0.7 0];
        else
            region='CA3'; clr=[1 0 0];
        end
        

        
        for c=1:length(super(reg).celloutput)
            
            figure
       
            for eventtype = [0 1]
        
        if eventtype == 1
            flag = 'reward';
            subplot(1,2,1)
                        title({[region ' (' num2str(daytetcell) ')  ' eventtype ' rasters']; ...
                               [flag]},'fontsize',18,'fontweight','bold');
        elseif eventtype == 0
            flag = 'error';
            subplot(1,2,2)
                         title([flag],'fontsize',18,'fontweight','bold');
        end
                
                
            data = super(reg).celloutput(c);
            daytetcell = data.daytetcell;

            
            % collect all reward trial spike rasters, disregarding wells, in
            % chronological order
            taggedrasters = [];
            for trig = find(data.trigmatrix(:,2)==eventtype)'   % error triggers
                % rewardtime and raster
                taggedrasters = [ taggedrasters ; data.rewardtimes{trig}    data.output{trig} ] ;
            end
            taggedrasters = sortrows(taggedrasters,1);
            rasters = taggedrasters(:,2:end);
            
            % plot
            plotraster2(data.times,rasters,1,'Color',clr,'linewidth',2)
            
            end
            
        end
        
        
        
    end
    
end

































