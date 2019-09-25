%Animal selection
%-----------------------------------------------------
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Dudley','Miles','Conley','Bond','Frank','Ten'};
%animals = {'Dudley'};
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Frank'};
%animals = {'Bond'};
animals = {'Bond'};
%-----------------------------------------------------

runscript = 1;
statictet = 30;   % theta eeg reference for all specified days
     % Chapati: best staticref is tetrode 1. matches Mizuseki for CA1
     % Egypt: best staticref is tetrode 11. matches Mizuseki for CA1
     % Frank: tetrode 30
     % Bond: tetrode 30
     % Corriander: reference is tetrode 24

  

if runscript == 1
    
    clear eegfilter
    clear super;
%Filter creation
%--------------------------------------------------------
epochfilter = [];
%epochfilter{1} = '(isequal($type, ''run'') || isequal($type,''sleep''))';
%epochfilter{1} = '(($exposure > 0 ))';
epochfilter{1} = '(isequal($type, ''run''))';

%ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && ($meanrate < 4))';
%ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && ($meanrate < 4))';

%ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && ($meanrate < 3))';
%ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && ($meanrate < 3))';
ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && isequal($type, ''principal''))';
%ca2cellfilter = '(isequal($area, ''CA2'') && ($numspikes > 100) && isequal($type, ''principal''))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && isequal($type, ''principal''))';

%timefilter = { {'get2dstate', '((abs($velocity) >= 8))'} };
timefilter_ca1 = {{'kk_getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'}, ...
                  {'get2dstate', '((abs($velocity) >= 5))'}, ... 
                  {'gethighthetatimes2','($ntheta >= 1)',[],'tetfilter','(isequal($area, ''CA1'') && ($numcells >= 1))','powerratio1',4,'powerratio2',.3,'mindur',.5}};
            % {'gethighthetatimes2','($ntheta == 1)',30,'powerratio1',4,'powerratio2',.3,'mindur',.5}};
            % {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 5))', 6}, ...
timefilter_ca3 = {{'kk_getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'}, ...
              {'get2dstate', '((abs($velocity) >= 5))'}, ... 
              {'gethighthetatimes2','($ntheta >= 1)',[],'tetfilter','(isequal($area, ''CA3'') && ($numcells >= 1))','powerratio1',4,'powerratio2',.3,'mindur',.5}};            
            

%eegfilter = {'geteegtet', 'theta', 'statictet',statictet};
%ca1tetfilter = '(isequal($area, ''CA1''))';
%eegfilter = {'geteegtet', 'theta', 'maxvar', 1, 'tetfilter',ca1tetfilter};
                   % picks tetrode with largest theta variance 
%eegfilter = {'geteegtet', 'theta', 'statictet', 11};
%eegfilter = {'geteegtet', 'theta', 'statictet',statictet};
eegfilter = {'geteegtet', 'eeggnd', 'sametet',1};

iterator = 'singlecelleeganal';

ca1f = createfilter('animal',animals,'epochs',epochfilter, 'excludetimefilter', timefilter_ca1, 'cells', ca1cellfilter,'eegtetrodes', eegfilter, 'iterator', iterator);
%ca2f = createfilter('animal',animals,'epochs',epochfilter, 'excludetimefilter', timefilter, 'cells', ca2cellfilter,'eegtetrodes', eegfilter, 'iterator', iterator);
ca3f = createfilter('animal',animals,'epochs',epochfilter, 'excludetimefilter', timefilter_ca3, 'cells', ca3cellfilter,'eegtetrodes', eegfilter, 'iterator', iterator);

ca1f = setfilterfunction(ca1f, 'dfakk_sta', {'spikes', 'eeg'},'window',[0.5 0.5]);
ca3f = setfilterfunction(ca3f, 'dfakk_sta', {'spikes', 'eeg'},'window',[0.5 0.5]);

ca1f = runfilter(ca1f);
%ca2f = runfilter(ca2f);
ca3f = runfilter(ca3f);

end

%% Designate hippocampal area to analyze.

for reg=[1 3]
    if reg==1
        region='CA1';
        clr='k';
        caf=ca1f;
    elseif reg==2
        region='CA2';
        clr='g';
        caf=ca2f;
    else
        region='CA3';
        clr='r';
        caf=ca3f;
    end

%% Generate reference string for labelling figures.

if isequal(eegfilter,{'geteegtet', 'eeg', 'sametet', 1})      % local reference
    referencestring = 'local reference';
else                                                            % static reference
    try
    referencestring = sprintf('static reference, tetrode %d',statictet);
    catch
    referencestring = 'maxvar tetrode';
    end
end
    
%% Consolidate single cells' across epochs in a day (.celloutput field)

dummyindex=[];
for i=1:length(caf.output{1})
    dummyindex = [dummyindex;caf.output{1}(i).index];           % collect all day-epoch indices
end
caf.celloutput=struct;

for i=1:size(caf.output{1},2)     % iterate through each cell-epoch
    daytetcell=caf.output{1}(i).index([1 3 4]);
    ind=[];
    while rowfind(daytetcell,dummyindex(:,[1 3 4]))~=0          % collect all rows (epochs)
        ind = [ind rowfind(daytetcell,dummyindex(:,[1 3 4]))];
        dummyindex(rowfind(daytetcell,dummyindex(:,[1 3 4])),:)=[0 0 0 0];
    end
     
    % gather eegwindows into dummy
    dummy=[];
    for r=ind
        dummy=[dummy ; caf.output{1}(r).eegwindows];
    end

    if ~isempty(dummy)
        if ~isfield(caf.celloutput,'index')      % if initial entry
            caf.celloutput(1).index=daytetcell;
            caf.celloutput(1).eegwindows=dummy;
        else            
            % concatenate structs (representing each cell) together
            A.index=daytetcell;
            A.eegwindows=dummy;
            caf.celloutput=horzcat(caf.celloutput,A);
        end
    end
    
end

%% Compute mean phase, modulation depth, phase histograms, Rayleigh test, and plot overall modulation w/ SE.

super{reg}=caf;


%% Plots.


% average of all STAs

means=[];

for i=1:length(caf.celloutput)
    means = [means ; mean(caf.celloutput(i).eegwindows,1)];
end

figure
plot(mean(means,1),'k','LineWidth',2)
title(sprintf('region %d',reg))







end
