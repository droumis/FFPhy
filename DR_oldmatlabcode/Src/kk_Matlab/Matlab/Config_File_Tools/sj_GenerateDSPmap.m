function sj_GenerateDSPmap(tetrodemap, varargin)
% sj_GenerateDSPmap(tetrodemap, options)
% eg. sj_GenerateDSPmap(gen_tetrodemap,'master','hot','slave1','humid','datadir','/data/kenny/kenny19/');
%
% This function creates and writes an nspike config file in the current
% working directory called 'dsplist.txt'.  
% Tetrodemap is a 27-by-n matrix, where n is 1-4, depending on how many
% preamp chips will be used. The 27 numbers in each column assigns a tetrode
% to each channel of a chip. The sequence of the assignemts is shown below:
%
% _________________________________________________________________
% | 1 | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 19 | 21 | 23 | 25 | Gr |
% | 2 | 4 | 6 | 8 | 10| 12 | 14 | 16 | 18 | 20 | 22 | 24 | 26 | 27 |
% 
%So, for example, if the number in index (27,2) of tetrodemap is 14, then
%tetrode 14 has a channel on the pin below ground on the second chip.
%For pins that are unused, assign them to tetrode 0, and they will be
%ignored.
%
%__Options__
%'spectialtetrodes' - any tetrodes that you do not want to view on the
%spike machines, just the EEG machine.  References, for example.
%
%'master' - give the name of the master machine- required
%
%'slave1' - required
%
%'slave2' - if you are using more than 12 full tetrodes, you have the option of using a second slave for the remaining channels. (default '')
%
%'rt' - if you are collecting position data on a rig with an rt machine, the name of the rt machine needs to be declared. 
%If this is not declared, an extra 'camera sync' channel is added to the eeg's. (default '')
%
%'eegexclude' - give a vector of tetrode numbers to exclude from the eeg
%processing on the master machine.  This can be done to free up a DSP for
%spike data.
%
%'datadir' - give directory of where to put data (ending with a slash)- default '/data/' 
%
%'highspikefilter' - the high cutoff for the spike filter (use 6000 or 9000. Default 6000.)
%
%'lowspikefilter' - the low cutoff for the spike filter (use 100, 300, or 600. Default 300.)
%
%Example:
%GenerateDSPmap(mattiasmap,'specialtetrodes',[12 29],'eegexclude',[18:30],'master','drizzle','slave1','rain','slave2','mist','rt','fog','datadir','/data/mkarlsso/')
%

specialtetrodes = [0];
mastermachine = '';
slave1 = '';
slave2 = '';
rtmachine = '';
datadir = '/data';
eegexclude = [0];
lowspikefilter = 300;
highspikefilter = 6000;


%assign any options
for i = 1:length(varargin)
    if (mod(i,2))
        switch varargin{i}
            case 'specialtetrodes'
                specialtetrodes = [varargin{i+1}(:);0];
            case 'master'
                mastermachine = varargin{i+1};
            case 'slave1'
                slave1 = varargin{i+1};
            case 'slave2'
                slave2 = varargin{i+1};
            case 'rt'
                rtmachine = varargin{i+1};
            case 'eegexclude'
                eegexclude = [varargin{i+1}(:);0];
            case 'datadir'
                datadir = varargin{i+1};
            case 'lowspikefilter'
                lowspikefilter = varargin{i+1};
            case 'highspikefilter'
                highspikefilter = varargin{i+1};
            otherwise
                error(['Option ''',varargin{i},''' unknown']);
        end
    end
end
                               
if ~isempty(specialtetrodes)
    %here we assign the default reference tetrode
    defaultref = specialtetrodes(1);
end

%load the dsp info
load dspmap.mat
junkdspchannels =  setdiff([0:126],[dspmap(:); 63]');

out = [];
eegchannels = [];
inputcolumns = size(tetrodemap,2);
inputlength = size(tetrodemap,1);
if (inputlength ~= 27)
    error('Length of input must be 27')
end

if isempty(mastermachine)
    error('You need to give a master machine name');
end
if isempty(slave1)
    error('You need to give a name for for slave1');
end
nslaves = 1;
if ~isempty(slave2)
    nslaves = nslaves+1;
end
if ~isempty(rtmachine)
    nslaves = nslaves+1;
end

junkdspcount = 0;
tetrodemap = tetrodemap(:);
tetrodemap2 = tetrodemap(:);
dspmap2 = dspmap(:);
fulltetrodes = setdiff(tetrodemap(:),specialtetrodes(:));
if ~(sum(rowcount(fulltetrodes,tetrodemap(:)) == 4) == length(fulltetrodes))
    channelcounts = rowcount(fulltetrodes,tetrodemap(:));
    problemtetrodes = fulltetrodes(find(channelcounts ~= 4))';
    problemtetrodecount = channelcounts(find(channelcounts ~= 4));
    for i = 1:length(problemtetrodes)
        if (problemtetrodecount(i) > 4)  %any tetrodes with more than four channels is not allowed  
            error(['Any tetrodes not specified in ''specialtetrodes'' must have no more than four channels defined. This is not the case. Problem tetrode: ',num2str(problemtetrodes(i)),'.']);
        else %if a tetrode has fewer then four channels, we fill in the extra channels with dummy DSP addresses
            tetrodeindex = find(tetrodemap2 == problemtetrodes(i));
            tetrodeindex = tetrodeindex(end);
            numjunktoadd = 4-problemtetrodecount(i);  %number of junk dsp channels that need to be added to complete the tetrode
            tmpdspadd = [];
            tmpchanneladd = [];
            for j = 1:numjunktoadd
                junkdspcount = junkdspcount+1;
                tmpdspadd = [tmpdspadd; junkdspchannels(junkdspcount)];
                tmpchanneladd = [tmpchanneladd; problemtetrodes(i)];
            end
            %insert the junk channels into the tetrode map and the dsp map
            if (tetrodeindex < length(tetrodemap2))
                tetrodemap2 = [tetrodemap2(1:tetrodeindex);tmpchanneladd;tetrodemap2(tetrodeindex+1:end)];
                dspmap2 = [dspmap2(1:tetrodeindex);tmpdspadd;dspmap2(tetrodeindex+1:end)];
            else
                tetrodemap2 = [tetrodemap2(1:tetrodeindex);tmpchanneladd];
                dspmap2 = [dspmap2(1:tetrodeindex);tmpdspadd];
            end
                               
        end
    end
        
end


values = unique(tetrodemap);
channels = tetrodemap2(:);

%for each pin, assign the tetrode number and dsp value
%also, use the first channel of each tetrode as the eeg
for j = 1:length(values)
    tmpval = values(j);
    count = rowcount(tmpval,tetrodemap2(:));
    channels(find(tetrodemap2 == tmpval),2) = [0:count-1]';
    %[rowval, colval] = find(tetrodemap == tmpval);
    tetfind = find(tetrodemap == tmpval);
    tmpdspmap = dspmap(:);
    eegchannels = [eegchannels; [tmpval tmpdspmap(tetfind(1))]];
end

%tmpdspmap = dspmap(:);
channels = [channels dspmap2(1:size(channels,1))];
channels = sortrows(channels,1);

eegchannels = sortrows(eegchannels,1);

%open the output file
FID = fopen('dsplist.txt','w');

%write all the general info
%*********************************************
fprintf(FID,['nslaves	',num2str(nslaves),'\n' , ...
'master	',mastermachine,'	11031\n', ...
'slave	',slave1,'	11026\n']);
if ~isempty(slave2)
    fprintf(FID,['slave	',slave2,'	11028\n']);
end
if ~isempty(rtmachine)
    fprintf(FID,['rtslave  ', rtmachine,'	11030\n']);
end
fprintf(FID,'%% The port numbers must be > 1024 and < 65335\n');
fprintf(FID,['port	11001\n', ...
'port	11002\n', ...
'port	11003\n', ...
'port	11004\n', ...
'port	11005\n', ...
'port	11006\n', ...
'port	11007\n', ...
'port	11008\n', ...
'port	11009\n', ...
'port	11010\n', ...
'port	11011\n', ...
'port	11012\n', ...
'port	11013\n', ...
'port	11014\n', ...
'port	11015\n', ...
'port	11016\n', ...
'port	11017\n', ...
'port	11018\n', ...
'port	11019\n', ...
'port	11020\n', ...
'port	11021\n', ...
'port	11022\n', ...
'port	11023\n', ...
'port	11024\n', ...
'port	11025\n']);

fprintf(FID,['datatype	',mastermachine,'	CONTINUOUS	POSITION	DIGITALIO FSDATA\n', ...	
'datatype	',slave1,'	SPIKE\n']);
if ~isempty(slave2)
    fprintf(FID,['datatype	',slave2,'	SPIKE\n']);
end
if ~isempty(rtmachine)
    fprintf(FID,['datatype	',rtmachine,'	 POSITION\n']);
end
fprintf(FID,['datadirectory	',mastermachine,'	',datadir,'\n', ...
'datadirectory	',slave1,'	',datadir,'\n']);
if ~isempty(slave2)
    fprintf(FID,['datadirectory	',slave2,'	',datadir,'\n']);
end

fprintf(FID,['audio	0	100	1000	0\n', ...
'audio	1	100	0	0\n', ...
'dspclock	1\n', ...
'colorfile	/usr/local/nspike/nspike_rgbcolor\n', ...
'posinput        -1\n', ...
'trackdarkpixels	0\n', ...
'posthresh	230\n', ...
'posxflip	1\n', ...
'posyflip	0\n', ...
'mpegquality	80\n', ...
'mpegslices	2\n', ...
'videocodec      1\n', ...
'videogopsize    0\n' ...
'ndioports	4\n', ...
'dioport	0	output\n', ...
'dioport	1	input\n', ...
'dioport	2	output\n', ...
'dioport	3	output\n', ...
'fsgui stim\n', ...
'fsdata POSITION\n', ...
'fsdata CONTINUOUS 3\n', ...
'eegtracelength	2.000000\n', ...
'calibrationfile	/usr/local/nspike/nspikecalfile\n\n']);

% 'nprograms       4\n' ...
% 'program 0       /home/lorenlab/NSpike/src-user/user_program &\n' ...
% 'program 1       &\n' ...
% 'program 2       &\n' ...
% 'program 3       &\n' ...
% 'usergui stim\n' ...
% 'enable_daq 1\n\n' ...


%***********************************************

%write the dspinfo lines
%****************************************************
numeegs = length(setdiff(tetrodemap(:),eegexclude));
if (numeegs > 16)
    eegdspinfo = [0 1 2; 1 1500 1500]';
else
    eegdspinfo = [0 1; 1 1500]';
end
fprintf(FID,'%%          dsp num    samplingrate     machine assigned to dsp\n');
for i = 1:size(eegdspinfo,1)
    fprintf(FID,'dspinfo     %d              %d               %s\n',eegdspinfo(i,1),eegdspinfo(i,2),mastermachine);
end
numspikedsp = ceil(length(setdiff(channels(:,1),specialtetrodes))/4);
numuseddsp = size(eegdspinfo,1)-1;
for i = 1:numspikedsp
    if (i>3)
        if (nslaves > 1)
            machine = slave2;
        else
            machine = slave1;
        end
    else
        machine = slave1;
    end
    fprintf(FID,'dspinfo     %d              %d              %s\n',i+numuseddsp,30000,machine);
end
fprintf(FID,'\n');
%***********************************************************


%write the dsp channel info for each channel
%************************************************************
fprintf(FID,'%%          electrode chan# dsp_channel\n');
%print the dsp info for each pin into the file
for i = 1:size(channels,1)
    if (channels(i,1) > 0)
        fprintf(FID, 'electmap       %d      %d      %d\n',channels(i,1),channels(i,2),channels(i,3));
    end
end

if isempty(rtmachine) %with no rt machine, we need to add the camera sync channel to the eeg's
    syncchannel = channels(i,1)+1;
    fprintf(FID, 'electmap       %d      %d      %d\n',syncchannel,0,63);
end
%**************************************************************


%print the properties of each eeg channel into the file
%************************************************************
totalchannelcount = 0;
fprintf(FID,'\n%%EEG channels\n');
fprintf(FID,'hostname   %s\n',mastermachine);
for i = 1:size(eegchannels,1)
    if ~ismember(eegchannels(i,1),eegexclude)
        fprintf(FID,'   channel %d\n',totalchannelcount);
        fprintf(FID,'       dspnum      %d\n',floor(totalchannelcount/16)+1);
        fprintf(FID,'       dspchan     %d\n',eegchannels(i,2));
        fprintf(FID,'       number      %d\n',eegchannels(i,1));
        fprintf(FID,'       electchan   %d\n', 0);
        fprintf(FID,'       depth       %d\n', 0);
        fprintf(FID,'       refelect    %d\n', defaultref);
        fprintf(FID,'       refchan     %d\n', 0);
        fprintf(FID,'       thresh      %d\n', 40);
        fprintf(FID,'       maxdispval  %d\n', 2000);
        fprintf(FID,'       filter      %d %d\n', 1,400);
        fprintf(FID,'       color       %d\n', mod(eegchannels(i,1),2)+1);
        totalchannelcount = totalchannelcount + 1;
    end
end
if isempty(rtmachine) %with no rt machine, we need to add the camera sync channel to the eeg's
    fprintf(FID,'   channel %d\n',totalchannelcount);
    fprintf(FID,'       dspnum      %d\n',floor(totalchannelcount/16)+1);
    fprintf(FID,'       dspchan     %d\n',63);
    fprintf(FID,'       number      %d\n',syncchannel);
    fprintf(FID,'       electchan   %d\n', 0);
    fprintf(FID,'       depth       %d\n', 0);
    fprintf(FID,'       refelect    %d\n', 0);
    fprintf(FID,'       refchan     %d\n', 0);
    fprintf(FID,'       thresh      %d\n', 10000);
    fprintf(FID,'       maxdispval  %d\n', 65500);
    fprintf(FID,'       filter      %d %d\n', 1,11000);
    fprintf(FID,'       color       %d\n', mod(eegchannels(i,1),2)+1);
    totalchannelcount = totalchannelcount + 1;
end

%**********************************************************


%skip to the next dsp 
%this line increments the channel counter to beginning of the next set of 16
totalchannelcount = max([16 (floor(((totalchannelcount-1)/16)+1)*16) ]); 


%print the properties of each pin into the file
%*****************************************************
channelcount = 0;
hostswitched = 0;
fprintf(FID,'\n%%Spike channels\n');
fprintf(FID,'hostname   %s\n',slave1);
for i = 1:size(channels,1)
    if ~ismember(channels(i,1), specialtetrodes);
        fprintf(FID,'   channel %d\n',channelcount);
        fprintf(FID,'       dspnum      %d\n',floor(totalchannelcount/16)+1);
        fprintf(FID,'       dspchan     %d\n',channels(i,3));
        fprintf(FID,'       number      %d\n',channels(i,1));
        fprintf(FID,'       electchan   %d\n', channels(i,2));
        fprintf(FID,'       depth       %d\n', 0);
        fprintf(FID,'       refelect    %d\n', defaultref);
        fprintf(FID,'       refchan     %d\n', 0);
        fprintf(FID,'       thresh      %d\n', 40);
        fprintf(FID,'       maxdispval  %d\n', 200);
        fprintf(FID,'       filter      %d %d\n', lowspikefilter,highspikefilter);
        fprintf(FID,'       color       %d\n', mod(channels(i,1),10)+1);
        channelcount = channelcount+1;
        totalchannelcount = totalchannelcount + 1;
        if ((channelcount > 47) & (~hostswitched) & (size(channels,1)>i) & (channels(i+1,1) ~= channels(i,1)) & (nslaves > 1))
            %skip to next dsp            
            %if isempty(slave2)
            %    error('You have too many channels- you must give a name for slave2');
            %end
            totalchannelcount = max([16 (floor(((totalchannelcount-1)/16)+1)*16) ]);
            channelcount = 0;
            hostswitched = 1;
            fprintf(FID,'\n%%Switch to second slave\n');
            fprintf(FID,'hostname   %s\n',slave2);
        end
    end
end
%****************************************************************

fclose(FID);
%open the file for viewing
open('dsplist.txt');
if (totalchannelcount > 128)
    error(['You have assigned more than 128 channels!  This will not work becuase we only have 8 DSPs.  Total channels assigned: ', num2str(totalchannelcount)])
end

disp(['Total channels assigned: ', num2str(totalchannelcount)])


