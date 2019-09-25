% plots the combined place fields of all place cells for a day (mean normalised rate) to see the
% distribution of placefields on the track

Veqn = '>=0';
minV =  str2num(Veqn(end));
maxstage = 3; % [1 2 3]
minVPF = 3; %cm/sec
minPeakPF = 3;
lessthan=0;
includestates = 6;

%Animal selection
%-----------------------------------------------------
animals = {'I1'};
%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------




%Filter creation
%--------------------------------------------------------
% day filterionno

%days='[1:10]';%,'1:10';
%days = '[1:1]';
days = '[11]';




%% get placefield data
epochtype='Run';
epochfilter = ['isequal($epochtype, ''Run'')'];
%epochfilter{1} = ['isequal($epoch, 2)'];
cellfilter = '(isequal($area, ''CA1'') && ($meanrate <10 ))'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
timefilter = { {'JY_getriptimes','($nripples == 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},{'JY_getlinvelocity', strcat('$velocity > ',num2str(minVPF))}};
f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);
%iterator = 'singlecellanal';
iterator = 'multicellanal';
f = setfilteriterator(f,iterator);
%out=setfilterfunction(f, 'JY_calcopenfieldoccupancy', {'spikes','data'});
%out=runfilter(out);
%outdata=out.output{1,1};

out2=setfilterfunction(f, 'JY_getriprate', {'ripples'});
out2=runfilter(out2);
day=out.epochs{1,1}(1,1);

% collect epoch data
daydata=[];
epochindex=[];
for i=1:size(out.data{1,1},2)
    daydata=[daydata;out.data{1,1}{1,i}];
    epochindex=[epochindex;repmat(out.epochs{1,1}(i,:),size(out.data{1,1}{1,i},1),1)];
end

% plot combined place fields for all cells

for i=1:size(epochindex,1)
    totalimage(:,:,i)=outdata{1,i}.normsmoothedspikerate;
end
    
    totalimage(isnan(totalimage))=0;
    imagedata=mean(totalimage,3);
    imagesc(flipud(imagedata));
    title(sprintf('Combined place field for all cells \n %s for day %s ',...
        animals{1,1}, num2str(day)));
    cd(strcat(out.animal{1,2},'Plot/'));
figurename = sprintf('Combined_plf_%s_%s',animals{1,1},num2str(day));


saveas(gcf, figurename, 'pdf');
close;
clear all;

