% autocorrelations -- modified from dfskk_ca1ca3xcorr


runscript = 1
window = .5   % specific window to plot
maxvalue = .002  % max value of firing probability plot (see below)
regions = [2]   % 1 ca1, 2 ca2, 3 ca3

if runscript
%Animal selection
%-----------------------------------------------------
animals = {'Egypt'}
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter = [];
epochfilter{1} = '(isequal($type, ''run''))';           %|| isequal($type,''sleep''))';

%cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && ($meanrate < 7))';
%ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && isequal($type, ''principal''))';
ca2cellfilter = '(isequal($area, ''CA2'') && ($numspikes > 100) && isequal($type, ''principal''))';  
%ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && isequal($type, ''principal''))';  

%ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && ($meanrate > 7))';
%ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && ($meanrate > 7))';  

% linstate, moving
%timefilter = {{'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'},{'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 5))', 6}};

% 2d state, moving
timefilter = {{'kk_getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'} , ...
              {'kk_get2dstate', '((abs($velocity) < 3))'}, ... 
              {'gethighthetatimes2','($ntheta == 0)',[],'tetfilter','(isequal($area, ''CA1'') && ($numcells >= 1))','powerratio1',4,'powerratio2',.3,'mindur',.5}};
% not moving              
%timefilter = {{'kk_get2dstate', '((abs($velocity) > 5))'}};   %, ... 
              
              
              
eegfilter = {'geteegtet', 'thetagnd', 'sametet',1}
iterator = 'singlecelleeganal';

%ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetimefilter', timefilter, 'eegtetrodes',eegfilter,'iterator', iterator);
ca2f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca2cellfilter,'excludetimefilter', timefilter, 'eegtetrodes',eegfilter,'iterator', iterator);
%ca3f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetimefilter', timefilter, 'eegtetrodes',eegfilter,'iterator', iterator);

%f = setfilterfunction(f, 'dfakk_calcxcorrmeasures', {'spikes', 'linpos','theta'}, 'edgespikes', 1, 'calclinfields', 1, 'calctrajxcorr', 1);
%ca1f = setfilterfunction(ca1f, 'dfakk_autocorr', {'spikes', 'linpos','thetagnd'},'tmax',1,'bin',0.001);
ca2f = setfilterfunction(ca2f, 'dfakk_autocorr', {'spikes', 'linpos','thetagnd'},'tmax',1,'bin',0.001);
%ca3f = setfilterfunction(ca3f, 'dfakk_autocorr', {'spikes', 'linpos','thetagnd'},'tmax',1,'bin',0.001);

%ca1f = runfilter(ca1f);
ca2f = runfilter(ca2f);
%ca3f = runfilter(ca3f);


end





for reg=regions       % iterate over regions
    
    
    if reg==1
        region='CA1';
        clr=[0 0 0];
        f=ca1f;
    elseif reg==2
        region='CA2';
        clr=[0 1 0];
        f=ca2f;
    else
        region='CA3';
        clr=[1 0 0];
        f=ca3f;
    end



if 1
% consolidate across epochs to obtain .celloutput  
dummyindex=[];
for i=1:length(f.output{1})
    dummyindex = [dummyindex;f.output{1}(i).index];           % collect all day-epoch indices
end
f.celloutput=struct;

for i=1:size(f.output{1},2)
    daytetcell=f.output{1}(i).index([1 3 4]);
    ind=[];
    while rowfind(daytetcell,dummyindex(:,[1 3 4]))~=0          % collect all rows (epochs)
        ind = [ind rowfind(daytetcell,dummyindex(:,[1 3 4]))];
        dummyindex(rowfind(daytetcell,dummyindex(:,[1 3 4])),:)=[0 0 0 0];
    end
     
    %collect data for all epochs pair participates in
    ac1=[];
    for r=ind
        ac1=[ac1 ; f.output{1}(r).ac1];
        if ~isempty(f.output{1}(r).time)
            time = f.output{1}(r).time;
        end
    end

    if ~isempty(ac1) 
        if ~isfield(f.celloutput,'index')      % if initial entry
            f.celloutput(1).index=daytetcell;
            f.celloutput(1).time=time;
            f.celloutput(1).ac1=sum(ac1,1);
            % now that all epochs collapsed, compute normalized ac
            f.celloutput(1).ac1norm = sum(ac1,1)/sum(sum(ac1,1));
        else                                   % subsequent entries
            A.index=daytetcell;
            A.time=f.output{1}(r).time;
            A.ac1=sum(ac1,1);
            % now that all epochs collapsed, compute normalized ac
            A.ac1norm=A.ac1/sum(A.ac1);
            f.celloutput=horzcat(f.celloutput,A);
        end
    end
    
end    
end

super{reg} = f;

end
    
    
if 1
    for reg = regions
        
    if reg==1
        region='CA1';
    elseif reg==2
        region='CA2';
    else
        region='CA3';
    end
        
        
    % collect, score, order, plot auto-correlograms of sig units
    

    gaussiansd = 3  % in binsize, here since bins are 2 ms, this is 10 ms
    binsize = 1;    % in ms
    
    allac = [];
    
    %collect
    windowindices = find(abs(super{reg}.celloutput(1).time) < window);
    timevec = 1000*super{reg}.celloutput(1).time(windowindices);
    for p=1:length(super{reg}.celloutput)
        cellac = super{reg}.celloutput(p);
        %  **check if xcorr has at least 100 spikes**
        if sum(cellac.ac1(windowindices))>=100
            % note that want to append day and tet here, to keep in order
            allac = [allac ; [cellac.index(1:2) cellac.ac1norm(windowindices) ]];
        end
    end
    
    smoothedallac = [];
    
    % smooth cross-correlogram w/ Gaussian
    for r = 1:size(allac,1)
        smoothedallac(r,:) = smoothvect(allac(r,3:end),gaussian(gaussiansd,length(allac(r,:))));
    end
    % order by day and tetrode (so can see progression over days)
    smoothedallac = sortrows(smoothedallac,[1 2]);
    
    %   plot 
    
    %full scale
    figure
    imagesc(timevec,1:size(smoothedallac,1),smoothedallac,[0 maxvalue])
    hold on
    % plot center line
    plot(zeros(1,9999),1:9999,'Color',[1 1 1],'LineWidth',2)
    title({[animals{1} ' ' region ' norm autocorrelograms'] , ca2cellfilter},'FontSize',14,'FontWeight','bold') 
    

    superac{reg} = smoothedallac;

    end
end


   
    % plot each region's sem
    
    figure
    hold on
    for reg = regions
        
        if reg==1
            clr=[0 0 0];
        elseif reg==2
            clr=[0 1 0];
        else
            clr=[1 0 0];
        end
        
        seplot(timevec,superac{reg},'Color',clr,'FaceAlpha',0.5)
   
    end
    title({[animals{1} ' ' region ' sem autocorrelograms'] , timefilter{1}{2},'FontSize',14,'FontWeight','bold') 


