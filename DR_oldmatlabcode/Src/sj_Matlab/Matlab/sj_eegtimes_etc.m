
% For Gideon - getting eeg time base and relating to spikes, etc

% Load dio file
%------------------
DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
load(DIOfile);
stim = DIO{day}{epoch}{16};
% if isempty(stim)
%             stim = DIO{day}{epoch}{15};
%         end
stim_starttime = stim.pulsetimes(:,1)./10; %ms
stim_endtime = stim.pulsetimes(:,2)./10; %ms
stim_length = stim.pulselength;
stim_isi = stim.timesincelast(2:end)./10; %ms

% Load EEG and ripple LFP file
%-------------------------

EEGfile = sprintf('%s/EEG/%seeg%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,tets(1));
load(EEGfile);
e = eeg{day}{epoch}{tets(1)};
t = geteegtimes(e);
pt = stim.pulsetimes ./ 10000; % in s
eind = lookup(pt(:,1), t);
e.samprate=round(e.samprate);

for t=1:length(tets),
    currtet=tets(t);
    ripfile = sprintf('%s/EEG/%sripple%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,currtet);
    load(ripfile);
    ripamp(t,:) = ripple{day}{epoch}{currtet}.data(:,1);
    ripenv(t,:) = ripple{day}{epoch}{currtet}.data(:,3);
end


% Load MU file
%------------------
multifile = sprintf('%s/%smulti%02d.mat', directoryname, prefix, day);
load(multifile);
for t=1:length(tets),
    currtet=tets(t);
    cmd=sprintf('multi%d = multi{day}{epoch}{%d}/10;',currtet,currtet); eval(cmd);
end

%% Set which multiunit firing rate AND which Ripple Power to use
%---------------------------------------------

ripampu=sum(ripamp,1); % Sum across tetrodes if multiple exist
ripenvu=sum(ripenv,1);
multiu=[];
for t=1:length(tets),
    currtet=tets(t);
    cmd=sprintf('curr_multi = multi%d;',currtet); eval(cmd);
    multiu=[multiu;curr_multi];
    % eg. multiu=[multi3; multi4; multi5; multi10; multi11; multi12];
end
tetu=tets;

% Align spikes to stimulation
%---------------------------
for i =1:length(stim_starttime)
    i;
    cnt=cnt+1;
    currstim = stim_starttime(i);
    currspks =  multiu(find( (multiu>=(currstim-pret)) & (multiu<=(currstim+postt)) ));
    currspks = currspks-(currstim-pret);
    histspks = histc(currspks,[0:binsize:pret+postt]);
    stim_spks{cnt}=currspks;
    stim_spkshist(cnt,:) = histspks;
end





