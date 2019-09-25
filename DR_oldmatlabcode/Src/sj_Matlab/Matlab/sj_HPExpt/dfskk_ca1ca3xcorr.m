
runscript = 1

if runscript
%Animal selection
%-----------------------------------------------------
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Dudley','Miles','Conley','Bond','Frank','Ten'};
%animals = {'Dudley'};
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
animals = {'Corriander'}
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter = [];
epochfilter{1} = '(isequal($type, ''run'') || isequal($type,''sleep''))';

cellpairfilter = {'allcomb', '(isequal($area, ''CA1'') && ($numspikes > 100) && ($meanrate < 7))', ...
                             '(isequal($area, ''CA1'') && ($numspikes > 100) && ($meanrate < 7))'};

%timefilter = {{'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'},{'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 5))', 6}};
timefilter = {{'kk_getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'}, ...
                  {'kk_get2dstate', '((abs($velocity) > 5))'}, ... 
                  {'gethighthetatimes2','($ntheta >= 2)',[],'tetfilter','(isequal($area, ''CA1'') && ($numcells >= 1))','powerratio1',4,'powerratio2',.3,'mindur',.5}};
              
eegfilter = {'geteegtet', 'thetagnd', 'sametet',1}
iterator = 'kk_paircelleeganal';

f = createfilter('animal',animals,'epochs',epochfilter,'cellpairs',cellpairfilter,'excludetimefilter', timefilter, 'eegtetrodes',eegfilter,'iterator', iterator);

%f = setfilterfunction(f, 'dfakk_calcxcorrmeasures', {'spikes', 'linpos','theta'}, 'edgespikes', 1, 'calclinfields', 1, 'calctrajxcorr', 1);
f = setfilterfunction(f, 'dfakk_calcxcorrmeasures', {'spikes', 'linpos','thetagnd'}, 'thetaphase',1,'edgespikes', 0);
f = runfilter(f);

end


if 1
% consolidate across epochs  
dummyindex=[];
for i=1:length(f.output{1})
    dummyindex = [dummyindex;f.output{1}(i).index];           % collect all day-epoch indices
end
f.pairoutput=struct;

for i=1:size(f.output{1},2)
    daytetpair=f.output{1}(i).index([1 3 4 5 6]);
    ind=[];
    while rowfind(daytetpair,dummyindex(:,[1 3 4 5 6]))~=0          % collect all rows (epochs)
        ind = [ind rowfind(daytetpair,dummyindex(:,[1 3 4 5 6]))];
        dummyindex(rowfind(daytetpair,dummyindex(:,[1 3 4 5 6])),:)=[0 0 0 0 0 0];
    end
     
    %collect data for all epochs pair participates in
    ac1=[]; ac2=[]; c1vsc2=[]; sph1=[]; sph2=[];
    for r=ind
        ac1=[ac1 ; f.output{1}(r).ac1];
        ac2=[ac2 ; f.output{1}(r).ac2];
        c1vsc2=[c1vsc2 ; f.output{1}(r).c1vsc2];
        sph1=[sph1 ; double(f.output{1}(r).sph1)/10000];
        sph2=[sph2 ; double(f.output{1}(r).sph2)/10000];
        if ~isempty(f.output{1}(r).time)
            time = f.output{1}(r).time;
        end
    end

    if ~isempty(c1vsc2) || ~isempty(ac1) || ~isempty(ac2)
        if ~isfield(f.pairoutput,'index')      % if initial entry
            f.pairoutput(1).index=daytetpair;
            f.pairoutput(1).time=time;
            f.pairoutput(1).ac1=ac1;
            f.pairoutput(1).ac2=ac2;
            f.pairoutput(1).c1vsc2=c1vsc2; 
            f.pairoutput(1).sph1=sph1; 
            f.pairoutput(1).sph2=sph2;
            if ~isempty(sph1)
                %sph1(sph1==3.1416)=-3.1415;  %bumps offending phase
                %sph1(sph1==-3.1416)=3.1415;  %bumps offending phase
                out=rayleigh_test(sph1);
                f.pairoutput(1).rayleigh1=out.p;
            else
                f.pairoutput(1).rayleigh1=NaN;
            end
            A.sph2=sph2;
            if ~isempty(sph2)
                %sph2(sph2==3.1416)=-3.1415;   %bumps 3.1416 spikes forward
                %sph2(sph2==-3.1416)=3.1415;  %bumps 3.1416 spikes forward
                out=rayleigh_test(sph2);
                f.pairoutput(1).rayleigh2=out.p;
            else
                f.pairoutput(1).rayleigh2=NaN;
            end
        else            
            A.index=daytetpair;
            A.time=f.output{1}(r).time;
            A.ac1=sum(ac1,1);
            A.ac2=sum(ac2,1);
            A.c1vsc2=sum(c1vsc2,1);
            A.sph1=sph1;
            if ~isempty(sph1)
                out=rayleigh_test(sph1);
                A.rayleigh1=out.p;
            else
                A.rayleigh1=NaN;
            end
            A.sph2=sph2;
            if ~isempty(sph2)
                out=rayleigh_test(sph2);
                A.rayleigh2=out.p;
            else
                A.rayleigh2=NaN;
            end
            f.pairoutput=horzcat(f.pairoutput,A);
        end
    end
    
end    
end
    
    
if 1
    % collect, score, order, plot cross-correlograms of sig units
    
    p_value = 0.05
    window = 0.15    % window around center point, in sec
    gaussiansd = 5  % in binsize, here since bins are 2 ms, this is 10 ms
    binsize = 2;    % in ms
    
    collectedxcorr = [];
    
    %collect
    windowindices = find(abs(f.pairoutput(1).time) < window);
    timevec = 1000*f.pairoutput(1).time(windowindices);
    for p=1:length(f.pairoutput)
        pair = f.pairoutput(p);
        if (pair.rayleigh1 < p_value) && (pair.rayleigh2 < p_value) && ~isempty(pair.c1vsc2)
            % check if xcorr has at least 100 spikes
            if sum(pair.c1vsc2(windowindices))>=100
                collectedxcorr = [collectedxcorr ; f.pairoutput(p).c1vsc2(windowindices)];
            end
        end
    end
 
    
    % smooth cross-correlogram w/ Gaussian
    for r = 1:size(collectedxcorr)
        smoothedxcorr(r,:) = smoothvect(collectedxcorr(r,:),gaussian(gaussiansd,length(collectedxcorr(r,:))));
    end
    
    %z-score as in Mizuseki--Buzsaki-2009
    zxcorr = bsxfun(@minus,smoothedxcorr,mean(smoothedxcorr,2));
    zxcorr = bsxfun(@rdivide,zxcorr,std(smoothedxcorr,0,2));
   
    %order by max index
    maxindex = [];
    for r=1:size(zxcorr)
        maxindex = [maxindex ; find(max(zxcorr(r,:))==zxcorr(r,:))];
    end
    
    % plot
    
    finalcorr = flipud(sortrows([maxindex zxcorr],1));
    imagesc(timevec,1:size(finalcorr,1),finalcorr(:,2:end),[-3,3])
    hold on
    % plot center line
    plot(zeros(1,9999),1:9999,'Color',[1 1 1],'LineWidth',2)
    % plot maxindices in a black line
    timeindices = [];
    for r=1:length(finalcorr(:,1))
        timeindices = [timeindices ; timevec(finalcorr(r,1))];
    end
    plot(timeindices,1:size(finalcorr),'k','LineWidth',3)
    title([animals{1} ' CA1-CA1 z-scored cross-correlograms (Mizuseki specs)'],'FontSize',14,'FontWeight','bold')
    
    
end



















% plot correlations and fields
if 0
g = gaussian(3, 18);
phbins = [0:pi/12:2*pi];
for a = 1:length(f)
    for e = 1:length(f(a).output)
	for c = 1:length(f(a).output{e})
	    xc = f(a).output{e}(c);
	    if (~isempty(xc.c1vsc2) & (max(xc.ac1) > 100) & (max(xc.ac2) > 100))
		figure(1);
		for i = 1:4
		    subplot(4,1,i);
		    tmpind = find(abs(xc.time) < .5);
		    if (~isempty(xc.c1vsc2{i}))
			plot(xc.time(tmpind), smoothvect(xc.c1vsc2{i}(tmpind), g));
		    end
		    if (i == 1)
			title(sprintf('%d %d - %d %d vs. %d %d', xc.index))
		    end
		end
		figure(2);
		lf1 = xc.lf1.trajdata;
		lf2 = xc.lf2.trajdata;
		pp1 = xc.lf1.phasedist;
		pp2 = xc.lf2.phasedist;
		ph1 = [];
		ph2 = [];
		for i = 1:4
		    subplot(4,2,2*i-1);
		    plot(lf1{i}(:,1), lf1{i}(:,5), 'b');
		    hold on;
		    plot(lf2{i}(:,1), lf2{i}(:,5), 'r');
		    set(gca, 'XLim', [0 180]);
		    subplot(4,2,2*i);
		    plot(pp1{i}.dist, pp1{i}.phase, 'b.');
		    hold on;
		    plot(pp1{i}.dist, pp1{i}.phase+2*pi, 'b.');
		    plot(pp2{i}.dist, pp2{i}.phase, 'r.');
		    plot(pp2{i}.dist, pp2{i}.phase+2*pi, 'r.');
		    ph1 = [ph1 ;pp1{i}.phase];
		    ph2 = [ph2 ;pp2{i}.phase];
		    set(gca, 'YLim', [0 4*pi]);
		    set(gca, 'XLim', [0 180]);
		end
		figure(3);
		subplot(2,1,1)
		c = histc(ph1, phbins);
		bar(phbins, c, 'b');
		subplot(2,1,2)
		c = histc(ph2, phbins);
		bar(phbins, c, 'r');
		'pausing'
		pause
		clf(1);
		clf(2);
	    end
	end
    end
end
end
	




