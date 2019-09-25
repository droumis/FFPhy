function sj_gspikes_to_fetm(gspikefile, nch, elec, Nspike, dopc, domatclust, win)
% sj_gspikes_to_fetm('08_030610_tet7_gspike_tetNF2', 4, '7', 1,1, 0);
% Use the spikes structure in gsort and generate .fet file for KKwik and
% .spk file for Klusters, with PCs included as default option
% and also as a matlcust .mat file and matclust _params.mat file if asked
% for

%%%%%%%%%%%%%%%
if (nargin < 2)
    nch = 4;
end

if (nargin < 3)
    elec = '1';
end

if (nargin < 4)
    Nspike = 1; % when data comes from Nspike
end
if nargin<5,
    dopc=1;
end

if nargin<6,
    domatclust=0;
end
   
if (nargin<7),
    win=0;   % Linux machine (0) or Windows machine for Kkwik
end



%%%%%%%%%%%%%%
load(gspikefile);
Samprange=[4:32];

%%%%%% Time stamps in xx.y ms - resolution of 0.1ms
Fs = spikes.Fs; % in kHz
SamplingRate = Fs*1000; %sampling rate within each waveform
%TimePerStep = (1/SamplingRate)*UnitsPerSec;
%TimePerWindow = TimePerStep*40;
%MAXAMP = 5;     %max voltage (in milivolts)

%%%%%%%%% Get Name for Feature File for Kkwik
fetfile = [gspikefile '.fet.' elec];
spkfile = [gspikefile '.spk.' elec];
matclust_file = [gspikefile 'matc'];
matclust_paramfile = [gspikefile 'matc_params'];

%%%%%%%%% Get PCs - Column 1 to 12
if dopc==1,
    spc = sj_pcasvd_gspikes(spikes.waveforms,4,Samprange);
else
    spc=[];
end

%%%%%%%%% Get Peak Ampl - Column 13-16
if Nspike==1
    tetrodeinfo = getTetrodeConfig2;
    tetnum=str2num(elec);
    thresh = tetrodeinfo(tetnum).thresh;
    triggers = thresh;
    waves=int16(zeros(size(spikes.waveforms_ch1,2),nch,size(spikes.waveforms_ch1,1)));
    for i=1:nch
        cmd=sprintf('waves(:,%d,:)=int16(spikes.waveforms_ch%d\'');',i,i); eval(cmd);
    end
    timestamps=uint32(spikes.fstimes*10);
    filedata.params = parmcalc(timestamps,waves,triggers)';
    ampl = filedata.params(:,1:nch);
    %waves(:,1,:)=int16(spikes.waveforms_ch1');
    %waves(:,2,:)=int16(spikes.waveforms_ch2');
    %waves(:,3,:)=int16(spikes.waveforms_ch3');
    %waves(:,4,:)=int16(spikes.waveforms_ch4');
else
    if nch>1
        for i=1:nch
            cmd=sprintf('ampl(:,i)=abs(max(spikes.waveforms_ch%d,[],2)) + abs(min(spikes.waveforms_ch%d,[],2));',i,i); eval(cmd);
        end
    else
        ampl = abs(max(spikes.waveforms,[],2)) + abs(min(spikes.waveforms,[],2));
    end
    ampl=ampl*1000;  % Multiplier for my old files
end

spc = [spc, ampl];

%%%%%%%%% Kkwik timestamps - Column 17
% Fs*Time in secs
spc = [spc, (Fs.*spikes.swtimes)];

%%%%%%%%% Save ascii file for Klustakwik

nfets = size(spc,2);
if win==0,  % Default- linux machine
    
    if dopc==1
        
        nfets = 17; % its 12+4+1=17
        fet_matrix = zeros(size(spc,1),nfets);
        fet_matrix = spc;
        
        % Meth 1 - use Matlab save
        %save(fetfile,'nfets','spc','-ascii');
        
        %Meth 2 - use fprintf: Note that this assuming that PCs are in there!
        fid = fopen(fetfile, 'wt');
        fprintf(fid, '%2.0f\n', nfets);
        fprintf(fid, '%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.2f %8.2f %8.2f %8.2f %8.2f\n', fet_matrix');
        fclose(fid);
        
    else
        
        nfets=17;
        fet_matrix = zeros(size(spc,1),17);
        fet_matrix_idxs = [13 14 15 16 17];
        fet_matrix(:,fet_matrix_idxs) = spc;
        
        % Meth 1 - use Matlab save
        % save(fetfile,'nfets','fet_matrix','-ascii');
        
        % Meth 2 - use fprintf
        fid = fopen(fetfile, 'wt');
        fprintf(fid, '%2.0f\n', nfets);
        fprintf(fid, '%1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %8.2f %8.2f %8.2f %8.2f %8.2f\n', fet_matrix');
        fclose(fid);
        
    end
    
else % win=1, windows machine
    
    if dopc==1
        
        nfets = 17; % its 12+4+1=17
        fet_matrix = zeros(length(timestamps),nfets);
        fet_matrix = spc;
        
        % Meth 1 - use Matlab save
        %save(fetfile,'nfets','spc','-ascii');
        
        %Meth 2 - use fprintf:
        fid = fopen(fetfile, 'wt');
        fprintf(fid, '%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %6.2f %6.2f %6.2f %6.2f %6.2f\n', fet_matrix');
        fclose(fid);
        
    else
        
        nfets=17;
        fet_matrix = zeros(length(timestamps),17);
        fet_matrix_idxs = [13 14 15 16 17];
        fet_matrix(:,fet_matrix_idxs) = spc;
        
        % Meth 1 - use Matlab save
        % save(fetfile,'nfets','fet_matrix','-ascii');
        
        % Meth 2 - use fprintf
        fid = fopen(fetfile, 'wt');
        fprintf(fid, '%1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %6.2f %6.2f %6.2f %6.2f %6.2f\n', fet_matrix');
        fclose(fid);
        
    end
    
end % if win=0


%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Save binary file for Klusters

if ~exist('waves','var'),
    waves=int16(zeros(size(spikes.waveforms_ch1,2),nch,size(spikes.waveforms_ch1,1)));
    for i=1:nch
        cmd=sprintf('waves(:,%d,:)=int16(spikes.waveforms_ch%d\'');',i,i); eval(cmd);
    end
end

% Save as single column, Samp1Elec1Spk1, Samp1,Elec2Spk1... Samp40Elec4Spk1
wavesli = permute(waves,[2 1 3]);
wavesli = wavesli(:);
fid = fopen(spkfile, 'w'); % Do I need big endian or native or ... Check later
fwrite(fid, wavesli, 'int16'); % Need 16 or 32 bit binary, not default uint8
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matclust format
if domatclust==1,
    
    % 1st file
    if ~exist('timestamps','var'),
        timestamps=uint32(spikes.fstimes*10);
    end
    save(matclust_file,'timestamps','waves');
    
    % 2nd file
    
    if Nspike==1
        filedata.params = [double(timestamps) filedata.params spc(:,1:12)];
        if dopc==1
            filedata.paramnames = {'Time',
                'Channel 1 Max',
                'Channel 2 Max',
                'Channel 3 Max',
                'Channel 4 Max',
                'Max width',
                'Max height change'
                'Wave1PC1',
                'Wave1PC2',
                'Wave1PC3',
                'Wave2PC1',
                'Wave2PC2',
                'Wave2PC3',
                'Wave3PC1',
                'Wave3PC2',
                'Wave3PC3',
                'Wave4PC1',
                'Wave4PC2',
                'Wave4PC3'};
               
        else
            filedata.paramnames = {'Time',
                'Channel 1 Max',
                'Channel 2 Max',
                'Channel 3 Max',
                'Channel 4 Max',
                'Max width',
                'Max height change'};
        end
        
    else
        filedata.params = [double(timestamps) spc];
        if dopc==1
            filedata.paramnames = {'Time',
                'Wave1PC1',
                'Wave1PC2',
                'Wave1PC3',
                'Wave2PC1',
                'Wave2PC2',
                'Wave2PC3',
                'Wave3PC1',
                'Wave3PC2',
                'Wave3PC3',
                'Wave4PC1',
                'Wave4PC2',
                'Wave4PC3',
                'Channel 1 Max',
                'Channel 2 Max',
                'Channel 3 Max',
                'Channel 4 Max',
                'KkwikTime'};
        else
            filedata.paramnames = {'Time',
                'Channel 1 Max',
                'Channel 2 Max',
                'Channel 3 Max',
                'Channel 4 Max',
                'KkwikTime'};
        end
    end
    filedata.filename = gspikefile;
    
    
    
    save(matclust_paramfile,'filedata');
    clear filedata
    
end







