

%DR load and plot ripple triggered LFP across regions

%use multiple tetfilters to define regions to grab from
%Animal selection-------------
animals = {'D10', 'D12'}

%Filter creation----------------
epochfilter{1} = ['isequal($environment, ''wtrack'')'];
epochfilter{2} = ['isequal($environment, ''sleepbox'')'];
epochfilter{3} = ['isequal($environment, ''openfield'')'];

%cell & tet filter------------
%cellfilter = '($
ca1tetfilter = {?(isequal($area,??ca1??)?};
    mectetfilter = {?(isequal($area,??mec??)?};
    portetfilter = {?(isequal($area,??por??)?};
    v2ltetfilter = {?(isequal($area,??v2l??)?};
    
    %time filter----------
    %timefilter = {{'getriptimes', '($nripples > 0)', 'getlinstate', '(($traj ~= -1) & (abs($velocity) <= 4))'}};
    timefilter = {{'kk_getconstimes','($cons==1)','ripplescons',1,'consensus_numtets',3,...
    'minthresh',3,'exclusion_dur',0}};
    
    %iterator ----------
    iterator = 'eeganal';
    
    %filter creation -------------
    ca1f = createfilter('animal',animal,'epochs',epochfilter,...
    'excludetime', timefilter, 'iterator', iterator, 'eegtetrodes', ca1tetfilter);
    mecf = createfilter('animal',animal,'epochs',epochfilter,...
    'excludetime', timefilter, 'iterator', iterator, 'eegtetrodes', mectetfilter);
    porf = createfilter('animal',animal,'epochs',epochfilter,...
    'excludetime', timefilter, 'iterator', iterator, 'eegtetrodes', portetfilter);
    v2lf = createfilter('animal',animal,'epochs',epochfilter,...
    'excludetime', timefilter, 'iterator', iterator, 'eegtetrodes', v2ltetfilter);
    
    %set analysis function
    ca1f = setfilterfunction(ca1f, 'calctotalmeanrate', {'eeg', 'pos', 'ripples'});
    mecf = setfilterfunction(ca1f, 'calctotalmeanrate', {'eeg', 'pos', 'ripples'});
    porf = setfilterfunction(ca1f, 'calctotalmeanrate', {'eeg', 'pos', 'ripples'});
    v2lf = setfilterfunction(ca1f, 'calctotalmeanrate', {'eeg', 'pos', 'ripples'});
    
    %run analysis
    ca1f = runfilter(ca1f);
    mecf = runfilter(mecf);
    porf = runfilter(porf);
    v2lf = runfilter(v2lf);
    
%% scrap 
    %load Ripples
    animalinfo = animaldef(fileprefix);
    ripout = kk_getconstimes(animalinfo{2},animalinfo{3},[day ep], 'ripplesdons',1,'consensus_numtets',3);
    consvec_rip2 = ripout{day}{ep}.cons;
    consvectimes_rip2 = ripout{day}{ep}.time;
    [periodtimes_rip2 periodripsinds] = dr_vec2list(consvec_rip2,consvectimes_rip2);
    
    
    
    
%% 
    %plot ripples... from mari's plotEEGripples_rewardtimes.m
    
    
    %% plots selected tetrodes' EEG side-by-side around reward times, looking at ripples
    
    function out = plotEEGripples_rewardtimes(animdirectory,animprefix,day,epoch,tetrodes,datatype,windowsize)
%% tetrodes is a vector of tetrodes
%% windowsize is in seconds: window to look before or after reward time
%% datatype is 'eeg', 'ripple', etc.
% animdirectory = '/opt/data40/mari/Cha/';

%example: plotEEGripples_rewardtimes('/opt/data40/mari/Gov/','gov',8,6,[2 11 12 14 19 21],'ripple',5)


eegstruct = loadeegstruct(animdirectory,animprefix,datatype,day,epoch,tetrodes);

%>>>> WTF DOES LOADEEGSTRUCT OUPUT CONTAIN? ;::
%one animal with all the days, epochs, ntrodes specified. 
%eegstruct{day}{epoch}{ntrode}
    %-- .data
    %-- .samprate
    %-- .starttime


if (day < 10)
    rewardfile=sprintf('%s%srewardinfo0%d.mat',animdirectory,animprefix,day);
    ripplefile=sprintf('%s%sripples0%d.mat',animdirectory,animprefix,day);
else
    rewardfile=sprintf('%s%srewardinfo%d.mat',animdirectory,animprefix,day);
    ripplefile=sprintf('%s%sripples%d.mat',animdirectory,animprefix,day);
end


if isempty(tetrodes)
    tetrodes = 1:50;
end

load(rewardfile);
rewardtimes=rewardinfo{day}{epoch}(:,2);
num_tetrodes=length(tetrodes);
epoch_starttime=eegstruct{day}{epoch}{tetrodes(1)}.starttime;
no_samp=length(eegstruct{day}{epoch}{tetrodes(1)}.data);
samprate=eegstruct{day}{epoch}{tetrodes(1)}.samprate;
epoch_endtime=(no_samp-1) * (1 / samprate);
windowsize=floor(windowsize*samprate);

%get ripple times for ripples that are present on >minrip tetrodes
load(ripplefile);
load(sprintf('%s%scellinfo.mat',animdirectory,animprefix));
%[riptimes] = kk_getriptimes('/opt/data40/mari/Gov/','gov',[day epoch],tetrodes,'minthresh',3);
[out, ripplestdout] = getripples([day epoch], ripples, cellinfo,'cellfilter','(isequal($area, ''CA1'') && ($numspikes > 20))', 'minstd',3,'minrip',2);


%>>> WTF DOES GETRIPPLES OUTPUT LOOKLIKE?? DOES IT TAKE THE RIPPLES.MAT
%ALREADY CREATED?? AND THEN COMBNINE ACROSS VALIDRIP TETRODES?? 
xx
start=rewardtimes/10000;
%exclude rewardtimes too close to beginning or end of eegdata
%     while start(1)<windowsize
%             start(1) = [];
%         end
%         while ~isempty(start) && start(end)> epoch_endtime-windowsize
%             start(end) = [];
%         end
for s=1:length(start)
    epoch_times=epoch_starttime:1/samprate:(epoch_starttime+(1/samprate)*no_samp);
    if ~isnan(start(s))
        rewtime=lookup(start(s),epoch_times);
        riptimes=lookup(out,epoch_times);
        
        %>>> WTF DOES LOOKUP DO? .. takes each reward time range (start:end) and all the
        %epoch timestamps... and probably returns all the epoch timestamps
        %that fall in the reward range... same for ripple time ranges
        
        %THIS IS KEY
        
        riptimes=sort(riptimes);
        find1=find(riptimes<(rewtime+windowsize)); %can probably make this faster by specifying the find star/end ranges together rather than doing an ismember afterwards
        find2=find(riptimes>(rewtime-windowsize));
        ripind = find(ismember(find1,find2));
        ripstoplot = riptimes(ripind);
        figure %('units','normalized','outerposition',[0 0 1 1])
        hold on
        for t=1:length(tetrodes)
            subplot(round(length(tetrodes)/2),2,t)
            
            %this is the key part where it's grabbing the raw eeg traces
            %based on the determined time for specific tetrodes.. looks
            %like it's overlap plotting each trace for each specified
            %ntrode.. i can basically use this code but without the reward
            %times at first.. and also split up the plotting to search for
            %an 'area' tag and seperate into subplots.. 
            %also can add the ripple band traces to the plots..     
            plot((rewtime-windowsize):(rewtime+windowsize),eegstruct{day}{epoch}{tetrodes(t)}.data((rewtime-windowsize):(rewtime+windowsize)), ...
                'k','LineWidth',1);
            if ~isempty(ripstoplot)
                plotraster(ripstoplot,-100,200,[],'Color','m')
            end
            line(rewtime,-200:200,'Color','g','LineWidth',2)
            title(num2str(tetrodes(t)));
            axis xy
        end
        hold off
    end
end

    end