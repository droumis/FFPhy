
% Plots output-DIO triggered spiking. Either reward or error.

    % relies on rewardinfo struct created by dioparser_wtrack (shantanu)
    % relies on "trigmatrix" outputted by dfakk_getrewardtrigspiking
    
    % Note that if your trials have a delay instituted between input and
    % output triggers at a well, then you will want to create a pseudotimes
    % data structure to simulate analogous times for error trials.
    

runscript = 0
if runscript
    
% Animal Selection
animals = {'Chapati'};

% Day Filter
dayfilter = 3:12;

% Epoch Filter
epochfilter = '(isequal($type, ''run''))';

% Time Filter
timefilter = { {'kk_getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'}};
             % {'kk_get2dstate', '$velocity < 4'} } ;
        
ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && isequal($type, ''principal''))';    %% ($meanrate < 3))';
ca2cellfilter = '(isequal($area, ''CA2'') && ($numspikes > 100) && isequal($type, ''principal''))';   
ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && isequal($type, ''principal''))';

%ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && ($meanrate < 3))';
%ca2cellfilter = '(isequal($area, ''CA2'') && ($numspikes > 100) && ($meanrate < 3))';
%ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && ($meanrate < 3))';

% Iterator
iterator = 'singlecellanal';

% Filter Creation
ca1f = createfilter('animal', animals, 'days',dayfilter,'epochs', epochfilter, 'cells', ca1cellfilter, 'excludetime', timefilter, 'iterator', iterator);
ca2f = createfilter('animal', animals, 'days',dayfilter,'epochs', epochfilter, 'cells', ca2cellfilter, 'excludetime', timefilter, 'iterator', iterator);
ca3f = createfilter('animal', animals, 'days',dayfilter,'epochs', epochfilter, 'cells', ca3cellfilter, 'excludetime', timefilter, 'iterator', iterator);

% Set Analysis Function
ca1f = setfilterfunction(ca1f, 'dfakk_getrewardtrigspiking', {'spikes','rewardinfo','ripples','pos'},'posflag',1,'window',[5 5],'binsize',0.001,'minthresh',0,'pseudo',.52);
ca2f = setfilterfunction(ca2f, 'dfakk_getrewardtrigspiking', {'spikes','rewardinfo','ripples','pos'},'posflag',1,'window',[5 5],'binsize',0.001,'minthresh',0,'pseudo',.52);
ca3f = setfilterfunction(ca3f, 'dfakk_getrewardtrigspiking', {'spikes','rewardinfo','ripples','pos'},'posflag',1,'window',[5 5],'binsize',0.001,'minthresh',0,'pseudo',.52);

% Run Analysis
ca1f = runfilter(ca1f);
ca2f = runfilter(ca2f);
ca3f = runfilter(ca3f);

end





% Collect across epochs and plot

for reg=1:3       % iterate over regions
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
    f.celloutput(ind).ripples = cell(1,size(f.output{1}(1).trigmatrix,1));
    f.celloutput(ind).output = cell(1,size(f.output{1}(1).trigmatrix,1));
    f.celloutput(ind).rewardoutcomes = cell(1,size(f.output{1}(1).trigmatrix,1));
    f.celloutput(ind).velocity = cell(1,size(f.output{1}(1).trigmatrix,1));
    % initialize individual output matrices
    for m=1:length(f.celloutput(ind).output)
        f.celloutput(ind).output{m}=[];         
        f.celloutput(ind).ripples{m}=[];        
        f.celloutput(ind).rewardtimes{m}=[];    
        f.celloutput(ind).rewardoutcomes{m}=[];
        f.celloutput(ind).velocity{m}=[]; 
    end

    %collect data across epochs
    for c=1:size(indices,1)
        if daytetcell(ind,:)==indices(c,[1 3 4])
                for m=1:length(f.celloutput(ind).output)
                    % concatenate
                    f.celloutput(ind).output{m}=[f.celloutput(ind).output{m} ; f.output{1}(c).output{m}];    
                    f.celloutput(ind).rewardtimes{m}=[f.celloutput(ind).rewardtimes{m} ; f.output{1}(c).rewardtimes{m}];
                    f.celloutput(ind).rewardoutcomes{m}=[f.celloutput(ind).rewardoutcomes{m} ; f.output{1}(c).rewardoutcomes{m}];
                    f.celloutput(ind).ripples{m}=[f.celloutput(ind).ripples{m} ; f.output{1}(c).ripples{m}];
                    f.celloutput(ind).velocity{m}=[f.celloutput(ind).velocity{m} ; f.output{1}(c).velocity{m}];
                end
        end
    end
end
    
super(reg) = f;
    
end



% plot reward vs error trials

if 0
    
    for reg=1:3
        
        if reg==1
            region='CA1'; clr=[0 0 0];
        elseif reg==2
            region='CA2'; clr=[0 0.7 0];
        else
            region='CA3'; clr=[1 0 0];
        end
        
        
        % plot cell-by-cell
        for c=1:length(super(reg).celloutput)
            
            data = super(reg).celloutput(c);
            daytetcell = data.daytetcell;
            figure
            hold on
            
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

                % collect all reward trial spike rasters (+ ripple 'rasters'), disregarding wells, in
                % chronological order
                taggedspikerasters = [];
                taggedriprasters = [];
                taggedvelocity = [];
                for trig = find(data.trigmatrix(:,2)==eventtype)'   % selected by reward or error
                    % rewardtime and raster
                    taggedspikerasters = [ taggedspikerasters ; data.rewardtimes{trig} data.output{trig} ] ;
                    taggedriprasters = [ taggedriprasters ; data.rewardtimes{trig} data.ripples{trig} ] ;
                    taggedvelocity = [ taggedvelocity ; data.rewardtimes{trig} data.velocity{trig} ] ;
                end
                % sort into chronological order + untag
                taggedspikerasters = sortrows(taggedspikerasters,1);
                spikerasters = taggedspikerasters(:,2:end);
                taggedriprasters = sortrows(taggedriprasters,1);
                riprasters = taggedriprasters(:,2:end);
                taggedvelocity = sortrows(taggedvelocity,1);
                velocitytraces = taggedvelocity(:,2:end);
                
                hold on
                % plot (no ripples)
                plotraster2(data.times,spikerasters,1,'Color',clr,'linewidth',2)
                % plot (with ripples)
                plotraster3(data.times,spikerasters,riprasters,velocitytraces,1,'Color',clr,'linewidth',2)
            end
        end
    end
end



% plot by well -- (disregard reward vs. error)

if 1
    
    for reg=1:3
        
        if reg==1
            region='CA1'; clr=[0 0 0];
        elseif reg==2
            region='CA2'; clr=[0 0.7 0];
        else
            region='CA3'; clr=[1 0 0];
        end

        % plot cell-by-cell
        for c=1:length(super(reg).celloutput)
            
            data = super(reg).celloutput(c);
            daytetcell = data.daytetcell;
                
                % cherry picking
                if 1


                cherrypick = [11 1 1; 10 14 1; 10 1 1; 8 14 3; 7 14 2 ; 6 14 2; 11 3 8; 11 3 7 ; 11 3 5 ; 11 3 4; ...
                              10 17 2; 9 3 1; 8 17 3; 8 3 1; 5 17 4; 5 17 2; 4 17 6 ; 4 17 2; 4 3 2; 3 17 1 ; 10 16 1; ...
                              10 2 2 ; 9 16 5 ; 9 13 1 ; 8 16 2 ; 8 16 1; 8 2 1; 7 16 2; 7 2 1; 6 16 9 ; 6 15 1; 5 16 8 ; ...
                              5 16 6 ; 5 16 5 ; 5 16 3 ; 5 16 2 ; 5 2 5];
                          
                cherrypick = [8 4 2; 4 9 7; 9 2 2; 4 9 3; 6 9 1 ; 4 9 5 ; 9 13 1; 9 10 1; 5 6 2; 9 3 1; 8 4 1; 9 9 2]                          
                if not(rowfind(daytetcell,cherrypick))        %~rowfind(daytetcell,cherrypick)
                    continue
                end
                end
                
            H = figure
            hold on
            
            for well = [0 1 2]
                
                if well == 0
                    flag = 'well 0';
                    subplot(1,3,1)
                    title({[animals{1} ' ' region ' (' num2str(daytetcell) ')  ' 'reward outs + error pseudos']; ...
                           [flag]},'fontsize',18,'fontweight','bold');
                elseif well == 1
                    flag = 'well 1';
                    subplot(1,3,2)
                    title([flag],'fontsize',18,'fontweight','bold');
                elseif well == 2
                    flag = 'well 2';
                    subplot(1,3,3)
                    title([flag],'fontsize',18,'fontweight','bold');                    
                end
                
                % collect all reward trial spike rasters (+ ripple 'rasters'), disregarding wells, in
                % chronological order
                taggedspikerasters = [];
                taggedriprasters = [];
                taggedvelocitytraces = [];
                for trig = find(data.trigmatrix(:,1)==well)'   % triggers selected by well number
                    taggedspikerasters = [ taggedspikerasters ; data.rewardtimes{trig} data.rewardoutcomes{trig} data.output{trig} ] ;
                    taggedriprasters = [ taggedriprasters ; data.rewardtimes{trig} data.rewardoutcomes{trig} data.ripples{trig} ] ;
                    taggedvelocitytraces = [ taggedvelocitytraces ; data.rewardtimes{trig} data.rewardoutcomes{trig} data.velocity{trig} ] ;
                end
                % sort into chronological order + untag
                taggedspikerasters = sortrows(taggedspikerasters,1);
                    spikerasters = taggedspikerasters(:,3:end);
                taggedriprasters = sortrows(taggedriprasters,1);
                    riprasters = taggedriprasters(:,3:end);
                    outcomes = taggedriprasters(:,2);
                taggedvelocitytraces = sortrows(taggedvelocitytraces,1);
                    velocitytraces = taggedvelocitytraces(:,3:end);
                
                hold on
                % plot long bar indicating error trial
                plotraster4(data.times,outcomes,1)
                
                % plot (no ripples)
                plotraster2(data.times,spikerasters,1,'Color',clr,'linewidth',2)
                % plot (with ripples)
                plotraster3(data.times,spikerasters,riprasters,velocitytraces,1,'Color',clr,'linewidth',2)
                
            end
            
                % if cherrypicked, then save file
                if 1
                str = sprintf('%d%d%dwv',daytetcell(1),daytetcell(2),daytetcell(3));
                saveas(H,str,'fig')
                clear H
                close all
                end
        end

        
    end
end





























