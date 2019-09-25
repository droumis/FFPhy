n=1
ripVeqn = '<1'
runVeqn = '>=2'
ripwelldist = 20
runwelldist = 20
PFwelldist = 20
lessthan = 1
maxstage = [10]%[10]%[1 2 3] %10 to concatenate stages
minVPF = 2
minPeakPF = 0  %if == 0 then no exclusion
mintime = 0 %=0 to not exclude times, 0.2-1 for ripples
comp = 1 % calctaskstage
subsampletwell = 0
subsamplepass = 0 % 1 to subsample by passlength
subsamplenpass = 0 %2,3,4
outripatwell = 0 %1 for spikes outside ripples at the well
vel2d = 0
wtrack = 5
nexttraj = 1
area = 'CA3CA1' %'CA3CA1', 'CA3', 'CA1'
figdir = []
% minnumspikes = 1
% subsample = 0
% nrips = 1
perday = 2 %1 to separate bu track, 2 to not separate
trk = 'TrackA' %'TrackA' or 'TrackB'
exp = 'exposureday' %'exposureday' or 'exposure

%Animal selection
%-----------------------------------------------------
if wtrack == 1
    animals = {'Bond', 'Frank', 'Nine'} %for more familiar track A, recording started on day 7
elseif wtrack == 2
    animals = {'Miles', 'Ten', 'Conley', 'Dudley'} %Mattias's animal: 3 day separation
    %conley seems to have linear velocity strangeness: peak at 34cm/s
elseif wtrack == 3 %all animals
    animals = {'Miles','Conley','Bond','Frank','Nine','Ten', 'Five', 'Eight', 'Coriander'}; %
elseif wtrack == 4 %maggie's animals all 3 day sep
    animals = {'Five', 'Eight', 'Coriander'}; %
elseif wtrack == 5 %maggie's and mattias' animals, all 3 day sep
    animals = {'Miles', 'Ten', 'Conley', 'Dudley', 'Eight', 'Five', 'Coriander'} ; %
end

%Filter creation
%--------------------------------------------------------
%only include cells with placefields
if minPeakPF >0
    if wtrack ~= 0
        epochfPF = ['(isequal($type, ''run'')) ']; %just analyze days where switching between tasks
    end
    includecells{1} = calcincludecells(minVPF, minPeakPF, animals, epochfPF, lessthan, 'CA3CA1');
    pf = '';
end

%set velcoity filter function
velocfcn = 'getlinvelocity';
if vel2d==1
    velocfcn = 'get2dvelocity';
elseif vel2d==2
    velocfcn = 'get2dstate';
end

k=3
for j = 1:length(maxstage)
    i = maxstage(j);
    if wtrack ~= 0
        
            if perday == 0
                epochfilter{1} = ['(isequal($type, ''run'')) & ($dailyexposure == 1)  & ( (isequal($environment, ''TrackA'')) | (isequal($description, ''TrackA'')) )'];
                epochfilter{2} = ['(isequal($type, ''run'')) & ($dailyexposure == 1)  & ( (isequal($environment, ''TrackB'')) | (isequal($description, ''TrackB'')) )'];
            elseif perday == 1 | 2
                maxd = 4;
                if isequal(trk, 'TrackA')
                    maxd = 7;
                end
                if isequal(exp, 'exposure')
                    maxd = 12;
                end
                if perday == 1
                    trkeqn = [' & ( (isequal($environment, ''',trk,''')) | (isequal($description, ''',trk,''')) )' ];
                elseif perday == 2
                trkeqn = '';
                end
                for d = 1:maxd
                    if isequal(exp, 'exposureday')
                        epochfilter{d} = ['(isequal($type, ''run''))  & ($', exp,' ==', num2str(d), ') & ($dailyexposure == 1)' , trkeqn];
                    elseif isequal(exp, 'exposure')
                        epochfilter{d} = ['(isequal($type, ''run''))  & ($', exp,' ==', num2str(d), ') & ( (isequal($environment, ''',trk,''')) | (isequal($description, ''',trk,''')) )'];
                    end
                end
            end
    end
        
    if isequal(area, 'CA3CA1') | isequal(area, 'CA1CA3')
        cellfilter =  ' ( (isequal($area, ''CA3'') | isequal($area, ''CA1'') ) && ($meanrate < 4) )';
        areaeqn = ' (isequal($area, ''CA1'') ) '  %to detect ripples, as of 8/10 should only be on CA1 tetrodes
    elseif isequal(area, 'CA3')
        cellfilter =  '(isequal($area, ''CA3'') && ($meanrate < 4))';
        areaeqn = ' (isequal($area, ''CA1'' ) )'%to detect ripples, as of 8/10 should only be on CA1 tetrodes
    elseif isequal(area, 'CA1')
        cellfilter =  '(isequal($area, ''CA1'') && ($meanrate < 4))';
        areaeqn = ' (isequal($area, ''CA1'' )) ' %to detect ripples, as of 8/10 should only be on CA1 tetrodes
    end
    
    iterator = 'multicellanal';
    
    %set timefilters
    timefilter{k}{i} = { };
    %timefilter{k}{i} = { {velocfcn, ['((abs($velocity) ',Veqn,'))']} };
    
    %if need to filter by stage
    if ismember(i, [ 1 2 3])
        timefilter{k}{i}{end+1} = {  'getcalctaskstage', ['($includebehave ==',num2str(i),')'], comp };
    end
    
    f{k}{i} = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter{k}{i}, 'iterator', iterator);
    if mintime>=0
        f{k}{i} = testexcludetimes(f{k}{i}, mintime); %removes epochs from analysis if all epoch excluded by excludetimes, mintime = 30
    end
    
    %only include cells with placefields
    if minPeakPF >0
f{k}{i} = excludecellsf(f{k}{i}, includecells{1});
    end
end

%% run filter
k=3
for j = 1:length(maxstage)
    i = maxstage(j);
    f{k}{i}= setfilterfunction(f{k}{i}, 'getspikesinripinrun', {'ripples','linpos', 'pos','spikes','task','cellinfo', 'DIO'} , ...
        areaeqn, ripVeqn, runVeqn, ripwelldist, runwelldist, 'appendindex', 1, 'appendincludetime', 1, ...
        'appendpasslength',1, 'subsampletwell', subsampletwell,'subsamplenpass', 0,'outripatwell',outripatwell,...
        'maxcelltet', 0, 'velocfcn', velocfcn, 'correcttraj',  [' ''correctorder'', [2 1 3]'], 'nexttraj', 1);%
    f{k}{i} = runfilter(f{k}{i});
end

%iteratorPFloc = 'singlecellanal';
%timefilterPFloc = {{'getlinvelocity', ['((abs($velocity) >=', num2str(minVPF),'))'] } };
%PFlocf = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilterPFloc, 'iterator', iteratorPFloc);

%% group
k=3;
for j = 1:length(maxstage)
    i = maxstage(j);
    for t=1:length(epochfilter) %for each task
        g{k}{i}{t}.ripspikes=[];
        g{k}{i}{t}.runspikes=[];
        g{k}{i}{t}.nrips=[];
        g{k}{i}{t}.rewardin=[];
        g{k}{i}{t}.passlength=[];
        g{k}{i}{t}.numtrials = [];
        for a = 1: size(animals,2)
            if size(f{k}{i}(a).output,2)>=t
            for e = 1:size(f{k}{i}(a).output{t},2)
                g{k}{i}{t}.ripspikes = stack( g{k}{i}{t}.ripspikes, [a*ones(size(f{k}{i}(a).output{t}(e).ripspikes,1),1) f{k}{i}(a).output{t}(e).ripspikes]);
                %[an day epoch passstart passend passcorrect traj  passlength ripincludetime spikespercell];
                % 1   2     3     4         5           6       7   8               9           10-end
                g{k}{i}{t}.runspikes = stack( g{k}{i}{t}.runspikes, [a*ones(size(f{k}{i}(a).output{t}(e).runspikes,1),1) f{k}{i}(a).output{t}(e).runspikes]);
                g{k}{i}{t}.nrips = stack( g{k}{i}{t}.nrips, [a*ones(size(f{k}{i}(a).output{t}(e).nrips,1),1) f{k}{i}(a).output{t}(e).nrips]);
                epnumtrials = [sum(f{k}{i}(a).output{t}(e).ripspikes(:,5)==0 )  sum(f{k}{i}(a).output{t}(e).ripspikes(:,5)==1 )];
                g{k}{i}{t}.numtrials = stack( g{k}{i}{t}.numtrials, repmat(epnumtrials, size(f{k}{i}(a).output{t}(e).nrips,1), 1) );
                
                %[index numberripples timewithinripples]
                g{k}{i}{t}.rewardin = stack( g{k}{i}{t}.rewardin, [a*ones(size(f{k}{i}(a).output{t}(e).rewardin,1),1) f{k}{i}(a).output{t}(e).rewardin]);
                g{k}{i}{t}.passlength = stack( g{k}{i}{t}.passlength, [a*ones(size(f{k}{i}(a).output{t}(e).passlength,1),1) f{k}{i}(a).output{t}(e).passlength]);
            end
            end
        end
    end
end

%% determine PF peak location & xclude
% if minPeakPF > 0 && PFwelldist > 0
%     PFlocf = excludecellsf(PFlocf, includecells{1});
%     PFlocf= setfilterfunction(PFlocf, 'calcpeaklocationpertraj', {'spikes','linpos'} , minPeakPF,'appendindex', 1, 'fromstart', 0);%
%     PFlocf = runfilter(PFlocf);
%     PFlocg = numericgroupcombine(PFlocf, 1);  %if output numeric, this combines across animals, first col is animal
%
%     exclude cells with PF near well
%     for j = 1:length(maxstage)
%         i = maxstage(j);
%         [g{3}{i}{1}.ripspikes cellexcl1]= excludewellPFcells(g{3}{i}{1}.ripspikes, PFlocg{1}, PFwelldist, 7, 9, 'noPF', noPFtraj);
%         [g{3}{i}{2}.ripspikes cellexcl2]= excludewellPFcells(g{3}{i}{2}.ripspikes, PFlocg{2}, PFwelldist, 7, 9, 'noPF', noPFtraj);
%     end
%
% end

%% separate into correct and incorrect & subsample as needed

for j = 1:length(maxstage)
    i = maxstage(j);
    for t = 1:length(g{3}{i})
        gripspk{1}{i}{t} = g{3}{i}{t}.ripspikes(g{3}{i}{t}.ripspikes(:,6)==0 & g{3}{i}{t}.ripspikes(:,7)>0 & ismember(g{3}{i}{t}.ripspikes(:,7), [1 3]) ,:);
        gripspk{2}{i}{t} = g{3}{i}{t}.ripspikes(g{3}{i}{t}.ripspikes(:,6)==1 & g{3}{i}{t}.ripspikes(:,7)>0 & ismember(g{3}{i}{t}.ripspikes(:,7), [1 3]),:);
        
        inclrows0 = g{3}{i}{t}.ripspikes(:,6)==0 & g{3}{i}{t}.ripspikes(:,7)>0 & ismember(g{3}{i}{t}.ripspikes(:,7), [1 3]) ;
        inclrows1 = g{3}{i}{t}.ripspikes(:,6)==1 & g{3}{i}{t}.ripspikes(:,7)>0 & ismember(g{3}{i}{t}.ripspikes(:,7), [1 3]) ;
        gdur{1}{i}{t} = [g{3}{i}{t}.nrips(inclrows0,:) g{3}{i}{t}.numtrials(inclrows0,:)]; %ripple include time
        gdur{2}{i}{t} = [g{3}{i}{t}.nrips(inclrows1,:) g{3}{i}{t}.numtrials(inclrows1,:)];
        
        
        %if need to subsample for # tri1als
        if subsamplenpass == 2
            eps = unique([gripspk{2}{i}{t}(:,1:3); gripspk{1}{i}{t}(:,1:3)] ,'rows');
            subgripspk{2}{i}{t} = []; subgripspk{1}{i}{t} = []; subgdur{2}{i}{t} = []; subgdur{1}{i}{t} = [];
            
            for e = 1:size(eps,1) %each epoch
                %find number of incorrect and randomly select = number of
                %correct
                data1 = gripspk{1}{i}{t}(ismember(gripspk{1}{i}{t}(:,1:3), eps(e,:), 'rows'),:);
                data2 = gripspk{2}{i}{t}(ismember(gripspk{2}{i}{t}(:,1:3), eps(e,:), 'rows'),:);
                %correct
                data3 = gdur{1}{i}{t}(ismember(gdur{1}{i}{t}(:,1:3), eps(e,:), 'rows'),:);
                data4 = gdur{2}{i}{t}(ismember(gdur{2}{i}{t}(:,1:3), eps(e,:), 'rows'),:);
                
                
                if ~isempty(data2) & size(data2,1)>size(data1,1)
                    corrind = randperm(size(data2,1)); %randomly shuffle them
                    inclind = corrind(1:size(data1,1));
                    tempdata2 = data2(inclind,:);
                    tempdata4 = data4(inclind,:);
                    subgripspk{2}{i}{t} = [ subgripspk{2}{i}{t}; tempdata2];
                    % subgripspk{1}{i}{t} = [ subgripspk{1}{i}{t}; data1];
                    subgdur{2}{i}{t} = [ subgdur{2}{i}{t}; tempdata4];
                    % subgdur{1}{i}{t} = [ subgdur{1}{i}{t}; data1];
                    %                 elseif ~isempty(data1) & size(data2,1)<size(data1,1) %if
                    %                 more inC than Correct
                    %                     incorrind = randperm(size(data1,1)); %randomly shuffle them
                    %                     inclind = incorrind(1:size(data2,1));
                    %                     tempdata1 = data1(inclind,:);
                    %                     subgripspk{1}{i}{t} = [ subgripspk{1}{i}{t}; tempdata1];
                    %                     subgripspk{2}{i}{t} = [ subgripspk{2}{i}{t}; data2];
                end
            end
            origgripspk = gripspk;
            gripspk{2} = subgripspk{2};
            gdur{2} = subgdur{2};
            %gripspk{1} = subgripspk{1};
            
            %if need to subsample for = L & R or 1 & 3 traj
        elseif subsamplenpass == 3
            eps = unique([gripspk{2}{i}{t}(:,1:3); gripspk{1}{i}{t}(:,1:3)] ,'rows');
            subgripspk{1}{i}{t} = []; subgripspk{2}{i}{t} = [];
            for e = 1:size(eps,1) %each epoch
                for k = 1:2
                    %find number of incorrect and randomly select = number of
                    %correct
                    data = gripspk{k}{i}{t}(ismember(gripspk{k}{i}{t}(:,1:3), eps(e,:), 'rows'),:);
                    data1 = data(data(:,7)==1,:);
                    data2 = data(data(:,7)==3,:);
                    
                    if ~isempty(data1) & ~isempty(data2)
                        shufind1 = randperm(size(data1,1)); %randomly shuffle them
                        shufind2 = randperm(size(data2,1)); %randomly shuffle them
                        maxsize = min([size(data1,1) size(data2,1)]);
                        subdata1 = data1(shufind1(1:maxsize),:);
                        subdata2 = data2(shufind2(1:maxsize),:);
                        
                        subgripspk{k}{i}{t} = [ subgripspk{k}{i}{t};subdata1; subdata2];
                    end
                end
            end
            origgripspk = gripspk;
            gripspk = subgripspk;
            
            %if need to subsample for 1st/correct vs 2nd/incorrect half
        elseif subsamplenpass == 4
            eps = unique([g{3}{i}{t}.ripspikes(:,1:3)] ,'rows');
            subgripspk{1}{i}{t} = []; subgripspk{2}{i}{t} = [];
            for e = 1:size(eps,1) %each epoch
                data = g{3}{i}{t}.ripspikes(ismember(g{3}{i}{t}.ripspikes(:,1:3), eps(e,:),'rows'),:);
                numtrials = size(data,1);
                %select 1st/2nd half of session
                data1 = data(floor(numtrials/2):end,:); %second half of sessions
                data2 = data(1:ceil(numtrials/2),:); %1st half of sessions
                %select for correct/incorrect and outbound trials
                data1 = data1(data1(:,6) == 0 & ismember(data1(:,7), [1 3]),:); %incorrect
                data2 = data2(data2(:,6) == 1 & ismember(data2(:,7), [1 3]),:); %correct
                subgripspk{1}{i}{t} = [ subgripspk{1}{i}{t};data1];
                subgripspk{2}{i}{t} = [ subgripspk{2}{i}{t};data2];
            end
            %origgripspk = gripspk;
            gripspk = subgripspk;
            
        end
    end
end

%% subsample by passlength
if subsamplepass == 1
    k =2;
    for j = 1:length(maxstage)
        i = maxstage(j);
        for t = 1:length(gripspk{1}{i})
            [gripspk{2}{i}{t} mindist{2}{i}{t}] = selectsubsample( gripspk{2}{i}{t}, gripspk{1}{i}{t}, [1:3], [8], 'zscore', 0); %col8 is passlength
        end
    end
end

%% activation prob & other measures
for k = 1:2 %1 = incorrect
    for j = 1:length(maxstage)
        i = maxstage(j);
        [ripout{k}{i}] = calcsubsamplemeasures2(gripspk{k}{i}, [1 2], [3 4 5 6 7 8 9], gdur{k}{i}, 5, 'minnumtrials', 5, 'minnumtrialsboth', [5 6 7], 'numcoactiv',0 ,'totalmeanrate', 1, 'proportioncellsactive',0,'numcellsactiveperrip', 0,  'coactivprob', 0, 'activprob', 0); %concat epochs, only for total mean rate
        [tmpout{k}{i}] = calcsubsamplemeasures2(gripspk{k}{i}, [1 2 3], [4 5 6 7 8 9], gdur{k}{i}, 5, 'minnumtrials', 5, 'minnumtrialsboth', [5 6 7],'totalmeanrate',0, 'numspikesperrip',0);
        for t = 1:length(ripout{k}{i})
            %calculate ripple rate & number of rips per pass too
            riprate = gdur{k}{i}{t}(:,4)./gripspk{k}{i}{t}(:,9); %number of ripples / time at well
            numrip = gdur{k}{i}{t}(:,4);
            ripout{k}{i}{t}.riprate = riprate;
            ripout{k}{i}{t}.numrip = numrip;
            
            if ~isempty(ripout{k}{i}{t})
                ripout{k}{i}{t}.propcellsactive = tmpout{k}{i}{t}.propcellsactive;
                ripout{k}{i}{t}.ncellsactiveperrip = tmpout{k}{i}{t}.ncellsactiveperrip;
                % ripout{k}{i}{t}.meanrates = tmpout{k}{i}{t}.meanrates;
                ripout{k}{i}{t}.coactivprob = tmpout{k}{i}{t}.coactivprob;
                ripout{k}{i}{t}.activprob = tmpout{k}{i}{t}.activprob;
            end
        end
    end
end

%% plot measures correct vs incorrect
names = fieldnames(ripout{1}{i}{1});
for j = 1:length(maxstage)
    i = maxstage(j);
    rprip{i} = [];
    for m =  size(names,1)
        nm = names{m};
        eval(['allmin = min([ripout{1}{i}{1}.',nm,';ripout{2}{i}{1}.',nm,';ripout{1}{i}{2}.',nm,';ripout{2}{i}{2}.',nm,']);']);
        eval(['allmax = max([ripout{1}{i}{1}.',nm,';ripout{2}{i}{1}.',nm,';ripout{1}{i}{2}.',nm,';ripout{2}{i}{2}.',nm,']);']);
        edges = linspace(allmin, allmax, 40);
        
        tmprp = zeros(1,length(ripout{1}{i})+1);
        for t = 1:length(ripout{1}{i})+1
            if t ==length(ripout{1}{i})+1
                eval(['data1=[ripout{1}{i}{1}.',nm,'; ripout{1}{i}{2}.',nm,'];'])
                eval(['data2=[ripout{2}{i}{1}.',nm,'; ripout{2}{i}{2}.',nm,'];'])
            else
                eval(['data1=[ripout{1}{i}{t}.',nm,'];'])
                eval(['data2=[ripout{2}{i}{t}.',nm,'];'])
            end
            if ~isempty(data1) && ~isempty(data2)
                data1 = data1(~isnan(data1));
                data2 = data2(~isnan(data2));
                %  data1 = data1(data1>0);
                %  data2 = data2(data2>0);
                
                
                label = nm;
                
                stitle = ['In/Correct Comparison Task', num2str(t),animals{1:end},' stage', num2str(i)];
                figname = [nm,'_correctVin_', pf ,'_tsk', num2str(t)];
                
                col = 1;
                plotmeanhistcumsum(data1, data2, col,label, stitle, n, edges)
                if ~isempty(figdir)
                    eval(['cd ',figdir])
                    print('-djpeg',figname)
                end
                n=n+1;
                [tmprp(t) h] = ranksum(data1,data2);
            end
        end
        rprip{i} = [rprip{i}; tmprp];
    end
end

% plot measures correct vs incorrect
names = fieldnames(ripout{1}{i}{1});
for j = 1:length(maxstage)
    i = maxstage(j);
    for m = size(names,1)
        nm = names{m};
        for t = length(ripout{k}{i})+1
            eval(['data1=[ripout{1}{i}{1}.',nm,'; ripout{1}{i}{2}.',nm,'];'])
            eval(['data2=[ripout{2}{i}{1}.',nm,'; ripout{2}{i}{2}.',nm,'];'])
            data1 = data1(~isnan(data1));
            data2 = data2(~isnan(data2));
            figure
            hold on
            bar([1 2], [mean(data1) mean(data2)], 0.6)
            errorbar2([1 2], [mean(data1) mean(data2)],  [stderr(data1) stderr(data2)] , 0.6)
            xlim([0.5 2.5])
            set(gca, 'fontsize', 24)
            set(gca, 'xtick', [1 2], 'xticklabel', {'Unrewarded', 'Rewarded'})
            n=n+1;
        end
    end
end

%% within cell


%% GLM per session per cell
% allps = []
% for j = 1:length(maxstage)
%     i = maxstage(j);
%     for t = 1:length(epochfilter)
%         tempg = g{3}{i}{t}.ripspikes(ismember(g{3}{i}{t}.ripspikes(:,6), [0 1]) & ismember(g{3}{i}{t}.ripspikes(:,7), [1 3]), :);
%         [out{i}{t} ps] = runglm(tempg, [1 2], [3:9], 6);
%         allps = [allps; ps];
%         %out = runglm(g, colindex, colexcl, colout)
%     end
% end

%% GLM over all sessions, prop cells active
allpropact = [];
for k = 1:2 %1 = incorrect
    for j = 1:length(maxstage)
        i = maxstage(j);
        [epout{k}{i}] = calcsubsamplemeasures(gripspk{k}{i}, [1 2 3], [ 4 5 6 7 8 9], 5, 'minnumtrials', 5, 'appendepoch', 1, 'numcoactiv',0 ,'totalmeanrate', 0, 'proportioncellsactive',1,'numcellsactiveperrip', 1,  'coactivprob', 0, 'activprob', 0); %concat epochs, only for total mean rate
        
        for t = 2%1:length(epochfilter)
            %eps = unique(epout{k}{i}{t}.propcellsactive(:,1:end-1), 'rows');
            eps = unique(epout{k}{i}{t}.ncellsactiveperrip(:,1:end-1), 'rows');
            outcome = gripspk{k}{i}{t}(ismember(gripspk{k}{i}{t}(:,1:size(eps,2)), eps, 'rows'),6:7);
            %allpropact = [allpropact; epout{k}{i}{t}.propcellsactive outcome];
            allpropact = [allpropact; epout{k}{i}{t}.ncellsactiveperrip outcome];
            % allpropact = [allpropact; epout{k}{i}{t}.propcellsactive(epout{k}{i}{t}.propcellsactive(:,end)>0,:) outcome(epout{k}{i}{t}.propcellsactive(:,end)>0,:)];
            %[ an day ep propcellsactive correct traj];
        end
    end
end

[b dev stats] = glmfit(allpropact(:,4), allpropact(:,5), 'binomial');

%% GLM 90% / 10%

propcorr = zeros(1000,1); propcorrI = propcorr; propcorrC = propcorr;
for n = 1:1000
    rind = randperm(length(allpropact)); %randomize indexes
    nind = sort( rind( 1:round(0.5*length(allpropact)) ) ); %randomly select 90% and sort
    tind = sort( rind(round(0.5*length(allpropact))+1:end) );%randomly select 10% and sort
    nallpropact = allpropact(nind,:); %apply 90% index to allpropact, so get random 90% of original data
    tallpropact = allpropact(tind,:); %apply 10% index to allpropact, so get random 10% of original data
    
    [nb ndev nstats] = glmfit(nallpropact(:,4), nallpropact(:,5), 'binomial');
    yhat = glmval(nb, tallpropact(:,4), 'logit');
    estout = round(yhat); %GLMS estimated outcome
    propcorr(n) = sum(estout == tallpropact(:,5))./length(estout); %proportion correct comparing GLM and real outcome
    incind = tallpropact(:,5)==0;
    propcorrI(n) = sum(estout(incind) == tallpropact(incind,5))./sum(incind); %proportion correct comparing GLM and real outcome
    Cind = tallpropact(:,5)==1;
    propcorrC(n) = sum(estout(Cind) == tallpropact(Cind,5))./sum(Cind); %proportion correct comparing GLM and real outcome
    
end
figure
subplot(3,1,1)
hist(propcorr)
xlabel('proportion correct')
title('comparing GLM prediction and outcome on 50% of data')
subplot(3,1,2)
hist(propcorrC)
xlabel('proportion correct')
title('comparing GLM prediction and outcome on Correct trials')
subplot(3,1,3)
hist(propcorrI)
xlabel('proportion correct')
title('comparing GLM prediction and outcome on incorrect trials')



%% path length
%calc activ prob and append pass and index
%minpl = 130; maxpl = 200;
minpl = -1; maxpl = 800;

for j = 1:length(maxstage)
    i = maxstage(j);
    for t = 1:length(g{3}{i})
        gpath{1}{i}{t} = gripspk{1}{i}{t}(: ,[1 2 3 8]);
        gpath{2}{i}{t} = gripspk{2}{i}{t}(: ,[1 2 3 8]);
    end
end

for j = 1:length(maxstage)
    i = maxstage(j);
    for t = 1:length(g{3}{i})
        if t == length(g{3}{i})+1
            pathl1 = [gpath{1}{i}{1}(:,4); gpath{1}{i}{2}(:,4)];
            pathl2 = [gpath{2}{i}{1}(:,4); gpath{2}{i}{2}(:,4)];
        else
            pathl1 = gpath{1}{i}{t}(:,4);
            pathl2 = gpath{2}{i}{t}(:,4);
        end
        %plot
        pathl1 = pathl1(pathl1>minpl & pathl1<maxpl);
        pathl2 = pathl2(pathl2>minpl & pathl2<maxpl);
        bins = linspace(minpl,maxpl,40);
        h1 = hist(pathl1, bins);
        h2= hist(pathl2, bins);
        figure
        plot(bins, h1/sum(h1))
        hold on
        plot(bins, h2/sum(h2), 'r')
        xlabel('path length')
        title(['stg', num2str(i), ' tsk', num2str(t)])
        rppl(t) = ranksum(pathl1, pathl2);
    end
end

%%
for j = 1:length(maxstage)
    i = maxstage(j);
    for t = 1:length(g{3}{i})
        if t == length(g{3}{i})+1
            pathl1 = [gpath{1}{i}{1}(:,4); gpath{1}{i}{2}(:,4)];
            priprate1 = [gdur{1}{i}{1}(:,4)./gdur{1}{i}{1}(:,5) ; gdur{1}{i}{2}(:,4)./gdur{1}{i}{2}(:,5) ]; % #rips/durripincludetime
            propact1 = [ripout{1}{i}{1}.propcellsactive; ripout{1}{i}{2}.propcellsactive];
            pathl2 = [gpath{2}{i}{1}(:,4); gpath{2}{i}{2}(:,4)];
            priprate2 = [gdur{2}{i}{1}(:,4)./gdur{2}{i}{1}(:,5) ; gdur{2}{i}{2}(:,4)./gdur{2}{i}{2}(:,5) ]; % #rips/durripincludetime
            propact2 = [ripout{2}{i}{1}.propcellsactive; ripout{2}{i}{2}.propcellsactive];
        else
            pathl1 = gpath{1}{i}{t}(:,4);
            priprate1 = gdur{1}{i}{t}(:,4)./gdur{1}{i}{t}(:,5); % #rips/durripincludetime
            propact1 = ripout{1}{i}{t}.propcellsactive;
            pathl2 = gpath{2}{i}{t}(:,4);
            priprate2 = gdur{2}{i}{t}(:,4)./gdur{2}{i}{t}(:,5); % #rips/durripincludetime
            propact2 = ripout{2}{i}{t}.propcellsactive;
        end
        %plot
        figure
        subplot(2,1,1);
        plot(pathl1(priprate1>=0 & pathl1>minpl & pathl1<maxpl), priprate1(priprate1>=0 & pathl1>minpl & pathl1<maxpl), '.')
        hold on
        lsline
        title('Incorrect')
        subplot(2,1,2);
        plot(pathl2(priprate2>=0 & pathl2>minpl & pathl2<maxpl), priprate2(priprate2>=0 & pathl2>minpl & pathl2<maxpl), '.')
        hold on
        lsline
        title('Correct')
        ylabel('riprate per pass')
        xlabel('path length')
        % xlim([100 400])
        subtitle(['stg', num2str(i), ' tsk', num2str(t)])
        
        figure
        subplot(2,1,1);
        plot(pathl1(pathl1>minpl & pathl1<maxpl), propact1(pathl1>minpl & pathl1<maxpl), '.')
        hold on
        lsline
        title('Incorrect')
        subplot(2,1,2);
        plot(pathl2(pathl2>minpl & pathl2<maxpl), propact2(pathl2>minpl & pathl2<maxpl), '.')
        hold on
        lsline
        title('Correct')
        ylabel('proportion cells active per pass')
        xlabel('path length')
        %xlim([100 400])
        subtitle(['stg', num2str(i), ' tsk', num2str(t)])
        
        %stats
        [rho{2}(1,t) corrp{2}(1,t)] = corr(pathl1(pathl1>minpl & pathl1<maxpl), propact1(pathl1>minpl & pathl1<maxpl));
        [rho{2}(2,t) corrp{2}(2,t)] = corr(pathl2(pathl2>minpl & pathl2<maxpl), propact2(pathl2>minpl & pathl2<maxpl));
        [rho{1}(1,t) corrp{1}(1,t)] = corr(pathl1(priprate1>=0 & pathl1>minpl & pathl1<maxpl), priprate1(priprate1>=0 & pathl1>minpl & pathl1<maxpl));
        [rho{1}(2,t) corrp{1}(2,t)] = corr(pathl2(priprate2>=0 & pathl2>minpl & pathl2<maxpl), priprate2(priprate2>=0 & pathl2>minpl & pathl2<maxpl));
    end
end


%% other measures
%plot duration of included time
for j = 1:length(maxstage)
    i = maxstage(j);
    rpdur{i} = [];
    allmin = min([gripspk{1}{i}{1}(:,9);gripspk{2}{i}{1}(:,9);gripspk{1}{i}{2}(:,9) ;gripspk{2}{i}{2}(:,9) ]);
    allmax = max([gripspk{1}{i}{1}(:,9) ;gripspk{2}{i}{1}(:,9) ;gripspk{1}{i}{2}(:,9) ;gripspk{2}{i}{2}(:,9) ]);
    edges = linspace(allmin, allmax, 40);
    
    tmprp = zeros(1,length(gripspk{1}{i})+1);
    for t = 1:length(gripspk{1}{i})+1
        if t ==length(gripspk{1}{i})+1
            data1=[gripspk{1}{i}{1}(:,9) ; gripspk{1}{i}{2}(:,9) ];
            data2=[gripspk{2}{i}{1}(:,9) ; gripspk{2}{i}{2}(:,9) ];
        else
            data1=[gripspk{1}{i}{t}(:,9) ];
            data2=[gripspk{2}{i}{t}(:,9) ];
        end
        if ~isempty(data1) && ~isempty(data2)
            data1 = data1(~isnan(data1));
            data2 = data2(~isnan(data2));
            label = 'Duration Included Time at Well';
            
            stitle = ['In/Correct Comparison Task', num2str(t),animals{1:end},' stage', num2str(i)];
            %figname = [nm,'_correctVin_', pf ,'_tsk', num2str(t)];
            
            col = 1;
            plotmeanhistcumsum(data1, data2, col,label, stitle, n, edges)
            if ~isempty(figdir)
                eval(['cd ',figdir])
                print('-djpeg',figname)
            end
            n=n+1;
            [tmprp(t) h] = ranksum(data1,data2);
        end
    end
    rpdur{i} = [rpdur{i}; tmprp];
end
%
% %mean number of ripples/pass
% for j = 1:length(maxstage)
%     i = maxstage(j);
%     data1 = [g{3}{i}{1}.nrips(g{3}{i}{1}.ripspikes(:,6)==0,4); g{3}{i}{2}.nrips(g{3}{i}{2}.ripspikes(:,6)==0,4)];
%     data2 = [g{3}{i}{1}.nrips(g{3}{i}{1}.ripspikes(:,6)==1,4); g{3}{i}{2}.nrips(g{3}{i}{2}.ripspikes(:,6)==1,4)];
%     figure
%     hold on
%     bar([1 2], [mean(data1) mean(data2)], 0.6)
%     errorbar2([1 2], [mean(data1) mean(data2)],  [stderr(data1) stderr(data2)] , 0.6)
%     xlim([0.5 2.5])
%     set(gca, 'fontsize', 24)
%     set(gca, 'xtick', [1 2], 'xticklabel', {'Unrewarded', 'Rewarded'})
%     ylabel('Ripples Per Pass')
%     rpnrip = ranksum(data1,data2);
% end

%% traj bias
 %[an day epoch passstart passend passcorrect traj  passlength ripincludetime spikespercell];
 % 1   2     3     4         5           6       7   8               9           10-end
 
 an = unique(gripspk{2}{i}{1}(:,1));
 figure
 for a=an'
     a
     %trajs1 = [gripspk{1}{i}{1}(gripspk{1}{i}{1}(:,1)==a,7)];
     %trajs2 = [gripspk{2}{i}{1}(gripspk{2}{i}{1}(:,1)==a,7)];
      trajs1 = [g{3}{i}{1}.ripspikes(g{3}{i}{1}.ripspikes(:,1)==a & g{3}{i}{1}.ripspikes(:,6)==0, 7)];
      trajs2 = [g{3}{i}{1}.ripspikes(g{3}{i}{1}.ripspikes(:,1)==a & g{3}{i}{1}.ripspikes(:,6)==1, 7)];
     h1 = histc(trajs1, [1 3]);
     h2 = histc(trajs2, [1 3])
     subplot(1,max(an), a)
     hold on
     bar(0.8:1.8, h1/sum(h1),0.5, 'b')
     bar(1.1:2.1, h2/sum(h2),0.5,'r')
     xlim([0 3])
     ylim([0 1])
 end
 
 subplot(1,max(an), 1)
     ylabel('fraction')
     xlabel('traj')

