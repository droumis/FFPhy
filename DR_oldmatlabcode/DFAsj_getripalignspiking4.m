function out = DFAsj_getripalignspiking4(index, excludetimes, spikes, ripples, tetinfo, pos, varargin)
% No change from original version. Goes with Ver4 Scripts
% out = DFAsj_getripalignspiking(spike_index, excludeperiods, spikes, ripples, tetinfo, options)

% Called from DFSsj_getripalignspiking
% Use tetinfo and tetfilter passed in, or redefine here to get riptets
% Then use ripples to getriptimes. Use inter-ripple-interval of 1 sec, and use a low-speed criterion.
% Then align spikes to ripples


tetfilter = '';
excludetimes = [];
maxcell = 0;
minstd = 3; 
lowsp_thrs = 4; %cm/sec
highsp_thrs = lowsp_thrs;
dospeed = 0;
minrip=1;

% For ripple trigger
% ------------------
binsize = 10; % ms
pret=550; postt=550; %% Times to plot
push = 500; % either bwin(2) or postt-trim=500. For jittered trigger in background window
trim = 50;
cellcountthresh = 3;  % Can be used to parse ripples
smwin=10; %Smoothing Window - along y-axis for matrix. Carry over from ...getrip4

rwin = [0 200];
bwin = [-500 -300];
push = 500; % either bwin(2) or postt-trim=500. If doing random events. See ... getrip4


for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'tetfilter'
            tetfilter = varargin{option+1};
        case 'excludetimes'
            excludetimes = varargin{option+1};
        case 'minstd'
            minstd = varargin{option+1};
        case 'minrip'
            minrip = varargin{option+1};
        case 'maxcell'
            maxcell = varargin{option+1};
        case 'dospeed'
            dospeed = varargin{option+1};
        case 'lowsp_thrs'
            lowsp_thrs = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

day = index(1);
epoch = index(2);

% Get riptimes
% -------------
if isempty(tetfilter)
    riptimes = sj_getripples_tetinfo(index, ripples, tetinfo, 'tetfilter', '(isequal($descrip, ''riptet''))','minstd',minstd,'minrip',minrip);
else
    riptimes = sj_getripples_tetinfo(index, ripples, tetinfo, 'tetfilter', tetfilter, 'minstd', minstd,'minrip',minrip);
end
% Can opt to have a cellcountthresh for each event as in getpopulationevents2 or  sj_HPexpt_ripalign_singlecell_getrip4
% Not using as of now

% Get triggers as rip starttimes separated by at least 1 sec
% ----------------------------------------------------------
rip_starttime = 1000*riptimes(:,1);  % in ms

% Find ripples separated by atleast a second
% --------------------------------------------
iri = diff(rip_starttime);
keepidx = [1;find(iri>=1000)+1];
rip_starttime = rip_starttime(keepidx);

% Implement speed criterion - Keep in ver4. Try both
% ----------------------------------------
if dospeed
    absvel = abs(pos{day}{epoch}.data(:,5)); % Can also use field 9
    postime = pos{day}{epoch}.data(:,1); % in secs
    pidx = lookup(rip_starttime,postime*1000);
    speed_atrip = absvel(pidx);
    lowsp_idx = find(speed_atrip <= lowsp_thrs);
    highsp_idx = find(speed_atrip > highsp_thrs);
    
    rip_starttime = rip_starttime(lowsp_idx);
end

% Get the spike times
% --------------------
sind = index;
if ~isempty(spikes{sind(1)}{sind(2)}{sind(3)}{sind(4)}.data)
    spikeu = spikes{sind(1)}{sind(2)}{sind(3)}{sind(4)}.data(:,1)*1000;  % in ms
else
    spikeu = [];
end

% Get firing rate of neuron
% -------------------------

totaltime = diff(spikes{sind(1)}{sind(2)}{sind(3)}{sind(4)}.timerange)./10000;
cellfr = length(spikeu)./totaltime;


% Also return cell tag, maybe?
% if exist spikes{sind(1)}{sind(2)}{sind(3)}{sind(4)}.tag
%     celltag = spikes{sind(1)}{sind(2)}{sind(3)}{sind(4)}.tag
% else
%     celltag = 'Und'; % Undefined
%     disp(sprintf('No cell tag for Day %d Ep %d tet %d cell %d',sind(1),sind(2),sind(3),sind(4)));
% end
% 
% if strcmp(celltag(1:3),'CA1')||strcmp(celltag(1:3),'iCA')
%     celltag='CA1Pyr';
% end
% if strcmp(celltag(1:3),'PFC')
%     celltag='PFC';
% end


cntrip=0; nspk = 0; cntbck=0;
% nstd=1: gaussian of length 4. nstd = 2: gaussian of length 7, nstd=3: gaussian of length 10.
nstd = round(binsize*2/binsize); g1 = gaussian(nstd, 3*nstd+1);

% Align to ripples and jittered point in bck window
% ------------------------------------------------
for i=2:length(rip_starttime)-1  % Skip first and last
    i;
    % Align to ripples
    % ------------------------------------
    cntrip=cntrip+1;
    currrip = rip_starttime(i);
    %ripsize_cell(cntrip) = rip_sizes(i); ripsize_multi(cntrip) = rip_sizes(i);
    
    % PFC
    currspks =  spikeu(find( (spikeu>=(currrip-pret)) & (spikeu<=(currrip+postt)) ));
    nspk = nspk + length(currspks);
    currspks = currspks-(currrip);
    histspks = histc(currspks,[-pret:binsize:postt]);
    rip_spks_cell{cntrip}=currspks;
    % Get no of spikes in response window, and back window
    trialResps(cntrip) =  length(find(currspks>=rwin(1) & currspks<=rwin(2)));
    trialResps_bck(cntrip) =  length(find(currspks>=bwin(1) & currspks<=bwin(2)));
    
%     if trialResps(cntrip)~=0,
%         keyboard;
%     end
    
    % Histogram
    % ---------
    %rip_spkshist_cell(cntrip,:) = histspks;
    if isempty(histspks), % Only happens if spikeu is empty
        histspks = zeros(size([-pret:binsize:postt]));
    end
    try
        histspks = smoothvect(histspks, g1);
    catch
        disp('Stopped in DFAsj_getripalignspiking');
        keyboard;
    end
    try
        rip_spkshist_cell(cntrip,:) = histspks;  
    catch
        disp('Stopped in DFAsj_getripalignspiking');
        keyboard;
    end
    
    % THis should not be used - Makse sure you are using bck and not rdm in later codes
    % Align to random jittered triggers in rdm bck window: (trigtime-500ms)+random point in 1-50ms
    % ------------------------------------------------------------------------------------------------
    cntbck=cntbck+1;
    
    %rdmidx = randperm(len_bckwin); % permute len_bckwin 1ms bins
    %rdmidx = rdmidx(1);
    %currtrig = bwin(1)+rdmidx; % This is time in ms. No need for binsize here. Some random time in bck window
    rdmidx = randperm(50); % permute 50 1ms bins
    rdmidx = rdmidx(1);
    currtrig = push-rdmidx; % This is time in ms. No need for binsize here. Some random time in bck window   
    currrip = rip_starttime(i);
    currrip = currrip - currtrig; % Get rdm time behind actual ripple time in given window
 
    currspks =  spikeu(find( (spikeu>=(currrip-pret)) & (spikeu<=(currrip+postt)) ));
    currspks = currspks-(currrip);
    rdm_spks_cell{cntbck}=currspks;
    % Get no of spikes in random response window
    trialResps_rdm(cntrip) =  length(find(currspks>=rwin(1) & currspks<=rwin(2)));
      
%     if trialResps_rdm(cntrip)~=0,
%         keyboard;
%     end
    
    % Histogram
    % ---------
    histspks = histc(currspks,[-pret:binsize:postt]);
  
    if isempty(histspks), % Only happens if spikeu is empty
        histspks = zeros(size([-pret:binsize:postt]));
    end
    rdm_spkshist_cell(cntbck,:) = histspks;
    histspks = smoothvect(histspks, g1);
    rdm_spkshist_cell(cntbck,:) = histspks; 
end


% TRIM ENDS of Histogram and redo bins_resp- Due to Time Smoothing
% -----------------------------------------------
rip_spkshist_cell = rip_spkshist_cell(:,trim/binsize:(pret+postt-trim)/binsize);
% Also do the rdm hist
rdm_spkshist_cell = rdm_spkshist_cell(:,trim/binsize:(pret+postt-trim)/binsize);

% Update pret and postt due to end trims
pret=pret-trim; postt=postt-trim;
timeaxis = -pret:binsize:postt;
bins_center = find(abs(timeaxis)<=100);
bins_resp = find(timeaxis>=rwin(1) & timeaxis<=rwin(2));
bins_bck = find(timeaxis>=bwin(1) & timeaxis<=bwin(2));




% Output
% ------
out.index = sind;
out.rip_spkshist_cell = rip_spkshist_cell;
out.rip_spks_cell = rip_spks_cell;
out.Nspikes = nspk; % Nspks in the entire -pret:postt window across all the events for real ripple events
out.rdm_spkshist_cell = rdm_spkshist_cell;
out.rdm_spks_cell = rdm_spks_cell;
% Nspikes in respective window - summed response
out.trialResps = trialResps;
out.trialResps_rdm = trialResps_rdm;
out.trialResps_bck = trialResps_bck;

% Nspikes summed across trials in response and bckgnd window
out.Nspikes_resp = sum(trialResps);
out.Nspikes_bck = sum(trialResps_bck);

out.pret = pret;
out.postt = postt;
out.binsize = binsize;
out.rwin = rwin;
out.bckwin = bwin;
out.bins_resp  = bins_resp;
out.bins_bck = bins_bck;
out.timeaxis = timeaxis;

% Propoerties
out.cellfr = cellfr;
out.Nrip = cntrip;




