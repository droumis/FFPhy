% plots the combined normalised place fields of all places cells active in
% a rippple



Veqn = '>=0'
minV =  str2num(Veqn(end))
maxstage = 3% [1 2 3]
minVPF = 2 %cm/sec
minPeakPF = 3
lessthan=0
includestates = 6

%Animal selection
%-----------------------------------------------------
%animals = {'L2'};

%animals = {'M2','M1','M3','K3','L2','L3','N1'};
animals = {'K3','L3','M2', 'N3', 'N1','L2','M3','M1','P2'};

%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------




%Filter creation
%--------------------------------------------------------
% day filter

days='[1:8]';%,'1:10';



%epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
%epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];
%epochfPF = ['($switchday > 0)'];

%epochfilter{1} = ['isequal($epochtype, ''Run'')'];
epochfilter{1} = ['($epoch <7)'];
%cellfilter = '(isequal($area, ''CA1'') && ($meanrate <7))'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%timefilter = { {'JY_getlinvelocity', '$velocity <0'} };
%timefilter = { {'JY_getriptimes','($nripples ==0)', [], 3,'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}, {'JY_getlinvelocity', '$velocity <0.01'}};
%timefilter2 = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter2 = { {'getriptimes','($nripples > 0)',[],'minthresh',2,'cellfilter',cellfilter,'tetfilter',[1 2 3 4]}};

%timefilter = { {'JY_getlinvelocity', '(($velocity) >= 0))', 6} };
%f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter);
f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter);


%only include cells with placefields
%if minPeakPF>0
%    includecells = calcincludecells(minVPF, minPeakPF, animals, epochfPF, lessthan);
%5    f = excludecellsf(f, includecells);
%end
%-----------------------------------------------------------


%run function- single cells
%--------------------------------------------
iterator = 'JY_singleepochanal';

f = setfilteriterator(f,iterator);

outall=setfilterfunction(f, 'JY_combineoccupancy', {'occupancyoutput'});

outall=runfilter(outall);

%% control group
out={};

for ii=1:4;

out{ii}=cellfun(@(x) x.occupancy, outall(1,ii).output{1,1}, 'Uniformoutput',false);
end

daydata=[];
for jj=0:3:21
    tempday=[];
    for kk=1:4;
        %newdata=vertcat(out{kk}{jj+1:jj+3});
        %newdata=cat(3,newdata(1:end/3,1:end),newdata((end/3)+1:2*end/3,1:end),newdata(2*(end/3)+1:end,1:end));
        
        newdata=vertcat(out{kk}{jj+3});
        
        newdata=newdata./max(newdata(:));
%         newdatamean=mean(newdata(newdata>0));
%         
%         newdatasd=std(newdata(newdata>0));
%         indexzeros=newdata==0;
%         newdata=(newdata-newdatamean)./newdatasd;
%         newdata(indexzeros)=-1;
        
        tempday=cat(3,tempday,newdata);

    end
    daydata=cat(3, daydata,mean(tempday,3));
end

% specify paramters for tight_subplot
figure;
dist=0;
%nrows=size(daydata,2);
nrows=1;
%columns=3; %7 to include sleep sessions
columns=size(daydata,3);
gap_h=0.005;
gap_w=0.005;
marg_h=[0.01 0.1];
marg_w=[0.01 0.2];
ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
ax=1;
cmap=colormap(hot(64));
cmap(1,:)=0;

maxval=0.002

% loop for each day
while ax<(nrows*columns);

    for d=1:size(daydata,3);
                        % get current axis
                        axes(ha(ax));                       
                        % plot image
                        imagedata = daydata(:,:,d);
                         % set gaussian parameters
                        stdv=2;
                        g = gaussian2(stdv,(6*stdv));
                        imagedata = filter2(g,(imagedata));
                        imagedata = flipud(imagedata);
                        %image(imagedata,'CDataMapping', 'scaled');
                        %imagesc(imagedata,[-0.9 -0.7]);
                        imagesc(imagedata,[0 maxval]);
                        
                        %image(imagedata);
        set(gca,'xtick',[],'ytick',[]);
                        set(gca,'PlotBoxAspectRatio',[1 1 1]);
                        ax=ax+1;
                        d=d+1;
                    end
                end
 

% Print title for whole figure

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
    1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 0.7,sprintf('Mean last session occupancy for each day control group'),'HorizontalAlignment'...
   ,'center','VerticalAlignment', 'top','FontSize',24);
step=maxval/6;
colorbarlabel=[0:step:maxval];
hcol=colorbar('YTickLabel', num2str(colorbarlabel(2:end)',2),'Location','EastOutside');
%hcol=colorbar('Location','EastOutside');
cpos=get(hcol,'Position');
cpos(4)=cpos(4)/2;      % Halve the thickness
cpos(1)=cpos(1)-0.05;
cpos(2)=cpos(2)+0.3;   %  Move it down outside the plot
set(hcol,...
    'Position',cpos)

set(gcf,'PaperPositionMode','auto')
set(gcf,'Position',[1145 457 1375 450])
set(gcf,'PaperSize',[20 8.5])

directoryname='/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/';
cd(directoryname);
print(sprintf('lastepochoccupancycontrol'),'-depsc');
close;

%% inactivation group

out={};

for ii=5:9;

out{ii}=cellfun(@(x) x.occupancy, outall(1,ii).output{1,1}, 'Uniformoutput',false);
end

daydata=[];
for jj=0:3:21
    tempday=[];
    for kk=5:9;
       newdata=vertcat(out{kk}{jj+3});
        
        newdata=newdata./max(newdata(:));
%         newdatamean=mean(newdata(newdata>0));
%         
%         newdatasd=std(newdata(newdata>0));
%         indexzeros=newdata==0;
%         newdata=(newdata-newdatamean)./newdatasd;
%         newdata(indexzeros)=-1;
        
        tempday=cat(3,tempday,newdata);
    end
    daydata=cat(3, daydata,mean(tempday,3));
end

% specify paramters for tight_subplot
figure;
dist=0;
%nrows=size(daydata,2);
nrows=1;
%columns=3; %7 to include sleep sessions
columns=size(daydata,3);
gap_h=0.005;
gap_w=0.005;
marg_h=[0.01 0.1];
marg_w=[0.01 0.2];
ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
ax=1;
cmap=colormap(hot(64));
cmap(1,:)=0;

maxval=0.002

% loop for each day
while ax<(nrows*columns);

    for d=1:size(daydata,3);
                        % get current axis
                        axes(ha(ax));                       
                        % plot image
                        imagedata = daydata(:,:,d);
                         % set gaussian parameters
                        stdv=2;
                        g = gaussian2(stdv,(6*stdv));
                        imagedata = filter2(g,(imagedata));
                        imagedata = flipud(imagedata);
                        %image(imagedata,'CDataMapping', 'scaled');
                        %imagesc(imagedata,[-0.9 -0.7]);
                        imagesc(imagedata,[0 maxval]);
                        
                        %image(imagedata);
        set(gca,'xtick',[],'ytick',[]);
                        set(gca,'PlotBoxAspectRatio',[1 1 1]);
                        ax=ax+1;
                        d=d+1;
                    end
                end
 

% Print title for whole figure

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
    1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 0.7,sprintf('Mean last session occupancy for each day muscimol group'),'HorizontalAlignment'...
   ,'center','VerticalAlignment', 'top','FontSize',24);
step=maxval/6;
colorbarlabel=[0:step:maxval];
hcol=colorbar('YTickLabel', num2str(colorbarlabel(2:end)',2),'Location','EastOutside');
%hcol=colorbar('Location','EastOutside');
cpos=get(hcol,'Position');
cpos(4)=cpos(4)/2;      % Halve the thickness
cpos(1)=cpos(1)-0.05;
cpos(2)=cpos(2)+0.3;   %  Move it down outside the plot
set(hcol,...
    'Position',cpos)

set(gcf,'PaperPositionMode','auto')
set(gcf,'Position',[1145 457 1375 450])
set(gcf,'PaperSize',[20 8.5])

directoryname='/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/';
cd(directoryname);
print(sprintf('lastepochoccupancymuscimol'),'-depsc');
close;



