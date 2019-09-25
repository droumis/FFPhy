
%%%%%%%%%%%%%%%%%%%%%%%%%
system_dependent memstats
feature('DumpMem')
%%%%%%%%%%%%%%%%%%%%%%%%%

figure;  redimscreen; orient(gcf,'landscape'); hold on;
set(gcf, 'PaperPositionMode', 'auto');
set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold');
    set(0,'defaultaxeslinewidth',2);

%%%%% Plotting etc %%%%%%
    
  plot(mprefr, mepoch,'ko','Linewidth',4,'Markersize',8); hold on;
    sign = find(tt_epochvspre(:,1)==1 );
    if length(sign)>=1,
        plot(mprefr(sign), mepoch(sign),'bo','Linewidth',4,'Markersize',8); hold on;
    end    
    title(['Epoch vs pre']);
    xlabel('Fir. Rate (Hz) in Pre-Whisk period')
    ylabel('Fir. Rate (Hz) in Whisking Epoch')


axis( [0 floor((2*postlength+binplot+sep+1*binsize)/binsize)  0 m+1] );
set(gca,'xtick',[floor([100:100:(postlength+binplot)]/binsize), floor([(postlength+binplot+sep):100:(postlength+binplot+sep+postlength)]/binsize)],...
    'xticklabel',{num2str( [([100:100:(postlength+binplot)]');([0:100:postlength]') ])});
text( floor((1.5*binsize)/binsize), m,['PRE-WHISK'],'FontSize', 16, 'FontWeight','bold');
text( floor((postlength+2*onsettime)/binsize), m,['WHISKING EPOCH'],'FontSize', 20, 'FontWeight','bold');
text( floor((postlength+binplot+1*binsize+sep)/binsize), m,['POST-WHISK'],'FontSize', 16, 'FontWeight','bold');
post1fr=macrofirsumm(nenum).post1(:,1:postlength); Ntrials = size(post1fr,1);
x = sum(reshape(post1fr,size(post1fr,1),binsize,floor(postlength/binsize)),2);
bin_post1fr(nenum,:) = reshape(mean(x,1),floor(postlength/binsize),1);
normbin_post1fr(nenum,:) = bin_post1fr(nenum,:)*(1000/binsize); %% Fir rate in binsize in Hz

m=60; % Hz
arr=0:0.1:m;

%%%% TRIGGERS
plot((floor(postlength/binsize))*ones(1,length(arr)),arr,'r-','Markersize',2,'Linewidth',2);
plot((floor((postlength+binplot)/binsize))*ones(1,length(arr)),arr,'r-','Markersize',2,'Linewidth',2);
plot((floor((postlength+binplot+sep)/binsize))*ones(1,length(arr)),arr,'r-','Markersize',2,'Linewidth',2);
%%% SHADED AREAS
%% Onset
xpoints=floor(postlength/binsize):floor((postlength+onsettime)/binsize);
jbfill(xpoints,(m-1)*ones(size(xpoints)),zeros(size(xpoints)),[0.5,0.5,0.5],[0.5,0.5,0.5],1,0.2);
%% Middle
xpoints=floor((postlength+middle_st)/binsize):floor((postlength+middle_end)/binsize);
jbfill(xpoints,(m-1)*ones(size(xpoints)),zeros(size(xpoints)),[0.5,0.5,0.5],[0.5,0.5,0.5],1,0.2);

    
    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Tip: Use the resshape command to include this in autospikecut
save 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spikes.waveforms = spikes.waveforms';
spikes.waveforms_ch1 = spikes.waveforms_ch1';
spikes.waveforms_ch2 = spikes.waveforms_ch2';
spikes.waveforms_ch3 = spikes.waveforms_ch3';
spikes.waveforms_ch4 = spikes.waveforms_ch4';
spikes


%%%  TO ADD SUMMY FIELDS TO TESTDATA1 FOR SYNCING WITH NEW SORTING GUI AND
%%%  POST ANALYSIS

for i=1:length(unique(spikes.igorstim)), spikes.sweepstim{i}=spikes.sweep(find(spikes.igorstim==i));end

spikes.sweepall=[];
for i=1:size(spikes.sweepstim,2),
    spikes.sweepall=[spikes.sweepall; spikes.sweepstim{i}];
end
spikes.sweepall=sort(spikes.sweepall);

spikes.Alligorstim=[];
for i=1:10
    spikes.Alligorstim=[spikes.Alligorstim randperm(9)];
end
clear i;

cmd=sprintf('swaveforms = spikes.waveforms_ch%d;',i); eval(cmd);


for i=1:100; 
    cmd=sprintf('a%d = 0;',i); eval(cmd);
end
str =['a';'b';'c';'d'];
for i=1:length(str); 
    cmd=sprintf('str(%d) = 0;',i); eval(cmd);
end    


fopen('MIC02.clu.1', 'wt')
fprintf(fid,'%2.0f\n',x);

fget1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load and plot labview acquired data


fn = 'D:\Chronic\Tetrodes\MIC06\Labview\2_16_06_plexon4'

fn=[files{lp_f}];


t1=0; t2=300; Fs=32000;

%data=labview_loadspike4(fn,length(ch),32000,t1,t2);
t1=0; t2=300;
fn = 'spike.dat'; Fs=32000;
data=labview_loadspike4(fn,5,Fs,300,600);  %% NEED TO READ ALL 8 channels
ch=[1:4];
data=data(ch,:); data=data';

figure; hold on;
plot(data(Fs*260:Fs*270,1),'b')
plot(data(Fs*260:Fs*270,2),'r')
plot(data(Fs*260:Fs*270,3),'g')
plot(data(Fs*260:Fs*270,4),'m')

figure; hold on;
subplot(4,1,1);hold on; 
plot(data(Fs*260:Fs*260.5,1),'b'); axis off;
subplot(4,1,2);hold on; 
plot(data(Fs*260:Fs*260.5,2),'r'); axis off;
subplot(4,1,3);hold on; 
plot(data(Fs*260:Fs*260.5,3),'g'); axis off;
subplot(4,1,4);hold on; 
plot(data(Fs*260:Fs*260.5,4),'m'); axis off;



figure; hold on;
plot(data(:,1),'b')
plot(data(:,2),'r')
plot(data(:,3),'g')
plot(data(:,4),'m')

figure; hold on; plot(data(2,2000000:2500000),'b')
plot(data(3,20000:52000),'r')
plot(data(7,1000000:1500000),'g')
plot(data(8,1000000:1500000),'m')

tempthr= [mean(mean(data(2,1:Fs))) mean(std(data(2,1:Fs)))];

tempthr= [mean(mean(data(5*Fs:60*Fs,2))) mean(std(data(5*Fs:60*Fs,2)))];

thrs1= tempthr(1)-3*tempthr(2);
thrs2= tempthr(1)-4*tempthr(2);
thrs3= tempthr(1)-5*tempthr(2);
thrs4= tempthr(1)-6*tempthr(2);

plot(1:10*Fs,-0.08*ones(1,10*Fs),'k')

thr=thrs1;
thr=-0.07;
indices = find( (data(:,1)<=thr) | (data(:,2)<=thr) | (data(:,3)<=thr) | (data(:,4)<=thr));
%indices = find (data(:,3)<=thr); % values over threshold
thr_index = indices([find(diff(indices)>1)+1]); % threshold crossings in +/- direction: point of threshold cross
thr_index =[indices(1); thr_index]; %spikes start:  index of threshold crossing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load and plot matlab acquired data

[data time abstime Events daqinfo]=daqread('MIC05jan31_a010.daq','Channels',[8]);


[data time abstime Events daqinfo]=daqread('MIC07_mar29_0.daq','Channels',[1:4]);
[data time abstime Events daqinfo]=daqread('mar20a.daq','Channels',[2:3]);

figure; hold on; plot(time(1:3000000),data(1:3000000,1),'b')

tempthr= [mean(mean(data(1:Fs,1))) mean(std(data(1:Fs,1)))];
thrs2= tempthr(1)-4.5*tempthr(2);thr=thrs2;ch=1;
thrs1= tempthr(1)-5*tempthr(2);
thrs3= tempthr(1)-6*tempthr(2);
thrs4= tempthr(1)-7*tempthr(2);

plot(1:10,thrs2*ones(1,10),'k')

figure; hold on; plot(time(1:300000),data(1:300000,1),'b')
plot(time(1:300000),data(1:300000,2),'r')
plot(time(1:300000),data(1:300000,3),'g')
plot(time(1:300000),data(1:300000,4),'m')

thr=thrs2;
%indices = find (data(:,ch)<=thr); % values over threshold
indices = find( (data(:,1)<=thr) | (data(:,2)<=thr) | (data(:,3)<=thr) | (data(:,4)<=thr));
thr_index = indices([find(diff(indices)>1)+1]); % threshold crossings in +/- direction: point of threshold cross
thr_index =[indices(1); thr_index]; %spikes start:  index of threshold crossing


samples=32000;

figure; hold on;
for i=1:8, subplot(8,1,i); plot(data1(i,1:32000));axis fill; end

figure; plot(data1(1:32000,5));

figure; plot(data1(1:32000,5:8))


tempthr= [mean(mean(data1(1:32000,5))) mean(std(data1(1:32000,5)))];
thrs2= tempthr(1)-4*tempthr(2);thr=thrs2;ch=1;
thrs1= tempthr(1)-3*tempthr(2);
thrs3= tempthr(1)-2.5*tempthr(2);

figure; hold on; plot(data1(1:320000,5:8));
plot(1:320000,thrs2*ones(1,320000),'r')


figure; hold on; clr=['b' 'r' 'g' 'k' 'c' 'm' 'y' 'b' 'r' 'g' 'k' 'c' 'm' 'y'];
for i=5:8, plot(data1(1:320000,i),clr(i));end
plot(1:320000,thrs2*ones(1,320000),'k');
plot(1:320000,thrs1*ones(1,320000),'r');
plot(1:320000,thrs3*ones(1,320000),'g');

%%%%%%%-------------------------------------------------------

spikes.hierarchy.assigns(find(spikes.hierarchy.assigns==12))=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


c=load('cal.dat');
temp = c(1)*temp.^3+c(2)*temp.^2+c(3)*temp + c(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; pack;
parameters.outliers_a = 1e-11; parameters.kmeans_options.divisions = 5; parameters.reint_out = 1; parameters.tmin = 0.001; parameters.tref=0.002; parameters.cutoff=0.08;

cd 'D:\DATA\CHRONIC\MIC11-G3R1\10_30_06';
load spiketetcut_all_thr-34-part1_tetNF.mat
spikes = ssm_outliers(spikes, parameters.outliers_a);
spikes = ssm_kmeans(spikes, parameters.kmeans_options);
spikes = ssm_energy(spikes);
save spiketetcut_all_thr-34-part1-energy32.mat

clear; pack;
parameters.outliers_a = 1e-11; parameters.kmeans_options.divisions = 5; parameters.reint_out = 1; parameters.tmin = 0.001; parameters.tref=0.002; parameters.cutoff=0.08;

load spiketetcut_all_thr-34-part2_tetNF.mat
spikes = ssm_outliers(spikes, parameters.outliers_a);
spikes = ssm_kmeans(spikes, parameters.kmeans_options);
spikes = ssm_energy(spikes);
save spiketetcut_all_thr-43_energy32.mat


clear; pack;
cd 'D:\DATA\CHRONIC\MIC11-G3R1\10_30_06';
sss_spkcut_labviewcont5('spike.dat',[1:4],32000,4,1);

clear; pack;
cd 'D:\DATA\CHRONIC\MIC11-G3R1\10_30_06\Whisker Data';
d='D:\DATA\CHRONIC\MIC11-G3R1\10_30_06\Whisker Data';
all_whisk_chronic(d,2,0.05);

clear; pack;
cd 'D:\DATA\CHRONIC\MIC11-G3R1\10_28_06\Whisker Data2';
d='D:\DATA\CHRONIC\MIC11-G3R1\10_28_06\Whisker Data2';
all_whisk_chronic(d,2,0.05);


clear; pack;
parameters.outliers_a = 1e-11; parameters.kmeans_options.divisions = 5; parameters.reint_out = 1; parameters.tmin = 0.001; parameters.tref=0.002; parameters.cutoff=0.08;

cd 'D:\DATA\CHRONIC\MIC11-G3R1\10_29_06_2';
load spiketetcut_all_thr-33_tetNF.mat
spikes = ssm_outliers(spikes, parameters.outliers_a);
spikes = ssm_kmeans(spikes, parameters.kmeans_options);
spikes = ssm_energy(spikes);
save spiketetcut_all_thr-30_sort1_clean_energy32.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


jbfill(rng(1):rng(2),lim2*ones(size([rng(1):rng(2)])),0.5*ones(size([rng(1):rng(2)])),[0,1,0],[0,1,0],1,0.8);


%---------------------------------------------------

fnb='b18.dat';fnw='b18_w103.dat'; t1=1; t2=4000;[f1,fb] = frames_chronic2(fnw,fnb,t1,t2); figure; hold on; imagesc(f1);

%---------------------------------------------------------------
files=dir('s*');
for rawtr=1:length(files)    
    w=load(files(rawtr).name);
    if size(w,1)<=2; w=w';end  ;
    tst(rawtr) = w(1,1);   
end


%---------------------------------------------------------------
% SEE CLUSTERSEPARATION.M IN MCLUST FOR IMPORTANT PLOTTING FEATURES

%%%% plot 2 axes with ylabel1 on left and ylabel2 on right
%%%% also, legends for labels

 [AX,H1,H2] = plotyy(bins_In,H_In,bins_Out,H_Out); 
	axes(AX(1)); ylabel('In Cluster');
	axes(AX(2)); ylabel('Out of Cluster');
    
 %%% Writing text
 
 subplot(2,2,1)	
	msgstr = {};
	msgstr{end+1} = ['Cluster ' num2str(iClust)];
	msgstr{end+1} = ' ';
	msgstr{end+1} = sprintf('L-Extra                = %5.3f',CluSep.L_Extra);
	msgstr{end+1} = sprintf('L-Ratio                = %5.4f',CluSep.L_Ratio);
% 	msgstr{end+1} = sprintf('Chi fit (corr)          = %3.2f',C);
	msgstr{end+1} = sprintf('Isolation Distance = %4.1f',CluSep.IsolationDist);
	msgstr{end+1} = ' ';
	msgstr{end+1} = 'Features: ';
	for iF = 1:length(FeaturesToGet)
		msgstr{end+1} = ['   ' FeaturesToGet{iF}];
	end

	axis off;
	h=text(0,0.5, msgstr);
    
    %%%----------------------------------------------------------------
    
    
     savename = strtok(wavefile,'W');
     savename = [savename 'shape'];
     cmd = (['save ''' dn '\' savename ''' shape']); eval(cmd);
     
     
     %%%%%%%%%%
     fid=fopen('TT_tetNF.cut');
assigns=textscan(fid,'%f','headerlines',7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GENERATE VARIABLE NAME 

varname = genvarname(str) 
 
%Read a column header hdr from worksheet trial2 in Excel spreadsheet myproj_apr23:
[data hdr] = xlsread('myproj_apr23.xls', 'trial2');
%Make a variable name from the text of the column header that will not conflict with other names:
v = genvarname(['Column ' hdr{1,3}]);
%Assign data taken from the spreadsheet to the variable in the MATLAB workspace:
eval([v '= data(1:7, 3);']);


%%%%%%%%%%%%%%
%%% SHADED AREAS
    %% Onset
    xpoints=floor(postlength/binsize):floor((postlength+onsettime)/binsize);
    jbfill(xpoints,(m-1)*ones(size(xpoints)),zeros(size(xpoints)),[0.5,0.5,0.5],[0.5,0.5,0.5],1,0.2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


set(0,'defaultaxesfontsize',16);set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',0.5);

set(gca,'Linewidth',0.5); set(gca,'Ticklength',[0.025,0.025]); set(gca,'Fontsize',24,'Fontweight','normal','FontName','Arial');



