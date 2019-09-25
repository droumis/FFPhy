function [ output_args ] = PlotReroutetraj( animaldir,animalname,day,epoch )
% Plots trajectory of each run

% datadirectory - directory where the output from batchDIOextract is
% located, eg '/data14/jai/F12_/'
% animalname -  name of the animal, eg F1
% days - vector of day or days to be analysed, eg [1], [1:3]
% 


% get the current directory to go back to after the function is done
currentdir = pwd;

% ---- File loading ----
% See if day number needs 0
dsz = '';
   if (day < 10)
      dsz = '0';
   end
   
% Specify data directory and load the file 
%animalname = animaldir(1:end-1);
datadir = '/data14/jai/';

% Converts the day and epoch into text
dayt = num2str(day);
epocht = num2str(epoch);
sfilename =strcat(datadir,animaldir,'/',animalname,'_',dsz,dayt,'.mat');

% Load the mat file


load(sfilename);

figure;

% plot trajectory between times

nrows=ceil(size(Data{1,day}{1,epoch}.Run,1)/5);
 
dist=0;

columns=5;
gap_h=0.01;
gap_w=0.01;
marg_h=[0.01 0.1];
marg_w=[0.01 0.01];


ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);


 
for i=1:size(Data{1,day}{1,epoch}.Run,1);
    tstart=Data{1,day}{1,epoch}.Run(i,3);
    tend=Data{1,day}{1,epoch}.Run(i,4);
    
    [val tstartref]=min(abs((Data{1,day}{1,epoch}.Pos.data(:,1)*10000-tstart)));
    [val tendref]=min(abs(Data{1,day}{1,epoch}.Pos.data(:,1)*10000-tend));

    % get traj between times
    
    traj=[Data{1,day}{1,epoch}.Pos.data(tstartref:tendref,2) Data{1,day}{1,epoch}.Pos.data(tstartref:tendref,3)];
    
    % calculate path traveled
    
    
    for j=1:size(traj,1)-1;
        dt=sqrt((traj(j+1,1)-traj(j,1))^2+(traj(j+1,2)-traj(j,2))^2);
        dist=dist+dt;
        j=j+1;
    end
    
    %subplot(nrows,5,i);
        
    axes(ha(i));
 
    plot(traj(:,1), traj(:,2));
    
    ylimval=get(gca,'ylim');
    xlimval=get(gca,'xlim');
    
    %title(sprintf('Run %d',i),'FontSize',10);
    ylim(ylimval);
    xlim(xlimval);
    %XTickLabel('FontSize',4);
    
    distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.85*(ylimval(2)-ylimval(1))+ylimval(1)];
    titlelabel=[0.6*(xlimval(2)-xlimval(1))+xlimval(1) 0.08*(ylimval(2)-ylimval(1))+ylimval(1)];
    
    
    text(distlabel(1),distlabel(2),sprintf('%s m',num2str(dist/1000,2)),'FontSize',10);
    text(titlelabel(1),titlelabel(2),sprintf('Run %d',i),'FontSize',10);
    
    set(gca,'xtick',[],'ytick',[]);
    
    i=i+1;

end

% Print title for whole figure

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

text(0.5, 0.98,sprintf('%s day %s epoch %s',animalname,dayt,epocht),'HorizontalAlignment'...
,'center','VerticalAlignment', 'top');


end

