function out = DFAsj_getripmoddata4(index, excludetimes, ripplemod, varargin)
% % Ver4 : Starting 10Feb2014 - Sync codes with everyone
% Do BCK Window instead of Random window. Keep in same structure with "rdm" label
% out = DFAsj_getripalignspiking(spike_index, excludeperiods, spikes, ripples, tetinfo, options)

% Need ripplemod only. 
% If calculating ripplealign from scratch, do as follows
% See DFSsj_getripalignspiking
% Use tetinfo and tetfilter passed in, or redefine here to get riptets
% Then use ripples to getriptimes. Use inter-ripple-interval of 1 sec, and use a low-speed criterion.
% Then align spikes to ripples


tetfilter = '';
excludetimes = [];
maxcell = 0;
minstd = 3;
% lowsp_thrs = 5; %cm/sec
% highsp_thrs = lowsp_thrs;
% dospeed = 0;
% 
% % For ripple trigger
% % ------------------
% binsize = 10; % ms
% pret=550; postt=550; %% Times to plot
% push = 500; % either bwin(2) or postt-trim=500. For jittered trigger in background window
% trim = 50;
% cellcountthresh = 3;  % Can be used to parse ripples
% smwin=10; %Smoothing Window - along y-axis for matrix. Carry over from ...getrip4
% 
% rwin = [0 150];
% bwin = [-500 -300];
% push = 500; % either bwin(2) or postt-trim=500. If doing random events. See ... getrip4


for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'tetfilter'
            tetfilter = varargin{option+1};
        case 'excludetimes'
            excludetimes = varargin{option+1};
        case 'minstd'
            minstd = varargin{option+1};
        case 'maxcell'
            maxcell = varargin{option+1};
        case 'dospeed'
            dospeed = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

day = index(1); epoch = index(2);

% Response window
try
    trialResps1 = ripplemod{index(1)}{index(2)}{index(3)}{index(4)}.trialResps;
catch
    day, epoch, index(3), index(4)
    keyboard;
end
trialResps2 = ripplemod{index(1)}{index(2)}{index(5)}{index(6)}.trialResps;
% IMP - Get no. of co-occurences in window
nsimul = length(find((trialResps1~=0) & (trialResps2~=0)));


% Bck window  - CHANGE TO BCK WINDOW, NOT RDM
trialResps1_rdm = ripplemod{index(1)}{index(2)}{index(3)}{index(4)}.trialResps_bck;
trialResps2_rdm = ripplemod{index(1)}{index(2)}{index(5)}{index(6)}.trialResps_bck;
nsimul_rdm = length(find((trialResps1_rdm~=0) & (trialResps2_rdm~=0)));




% % Coactive Z-score
% % ----------------
% coactivez = coactivezscore(trialResps1,trialResps2);
% coactivez_rdm = coactivezscore(trialResps1_rdm,trialResps2_rdm);
% 
% % Corr Coeff
% % -----------
% [r, p] = corrcoef(trialResps1,trialResps2);
% [r_rdm, p_rdm] = corrcoef(trialResps1_rdm,trialResps2_rdm);
% r = r(1,2); p = p(1,2);
% r_rdm = r_rdm(1,2); p_rdm = p_rdm(1,2);
% 
% 
% % Also do a shuffle for trialResps. 
% % ----------------------------------
% % Does shuffling only once make sense? 
% % Instead, get a p-value from shuffling
% % Cannot shuffle both - only 1
% r_shuf = []; %pshuf = [];
% n = 1000;
% for i = 1:n
%     rorder = randperm(length(trialResps1));
%     trialResps1_shuf = trialResps1(rorder);
%     [rtmp, ptmp] = corrcoef(trialResps1_shuf,trialResps2);
%     r_shuf(i) = rtmp(1,2); %pshuf(i) = ptmp(1,2);
%     coactivez_shuf(i) = coactivezscore(trialResps1_shuf, trialResps2);
% end

% 
% % Get p-value from shuffle. One-tailed test. 
% % FOr one-tailed, use p <0.025, not 0.05.
% % Or else, use absolute values of p, and use p<0.05
% % -------------------------------------------------
% 
% p_shuf_2t = length(find(abs(r_shuf) > abs(r)))./n;
% if r >= 0
%     p_shuf = length(find(r_shuf > r))./n;
% else % r < 0
%     p_shuf = length(find(r_shuf < r))./n;
% end
% 
% coactivez_pshuf_2t = length(find(abs(coactivez_shuf) > abs(coactivez)))./n;
% if coactivez >= 0
%     coactivez_pshuf = length(find(coactivez_shuf > coactivez))./n;
% else 
%     coactivez_pshuf = length(find(coactivez_shuf < coactivez))./n;
% end
% 
% % ----------------------------------
% % Also do a shuffle for trialResps_rdm. 
% % ----------------------------------
% r_rdmshuf = []; %prdmshuf = [];
% n = 1000;
% for i = 1:n
%     rorder = randperm(length(trialResps1_rdm));
%     trialResps1_rdmshuf = trialResps1_rdm(rorder);
%     [rtmp, ptmp] = corrcoef(trialResps1_rdmshuf,trialResps2_rdm);
%     r_rdmshuf(i) = rtmp(1,2); %prdmshuf(i) = ptmp(1,2);
%     coactivez_rdmshuf(i) = coactivezscore(trialResps1_rdmshuf, trialResps2_rdm);
% end
% % Get p-value from shuffle. One-tailed test
% % ------------------------------------------
%  p_rdmshuf_2t = length(find(abs(r_rdmshuf) > abs(r_rdm)))./n;
% if r_rdm >= 0
%     p_rdmshuf = length(find(r_rdmshuf > r_rdm))./n;
% else % r < 0
%     p_rdmshuf = length(find(r_rdmshuf < r_rdm))./n;
% end
% 
% coactivez_prdmshuf_2t = length(find(abs(coactivez_rdmshuf) > abs(coactivez_rdm)))./n;
% if coactivez_rdm >= 0
%     coactivez_prdmshuf = length(find(coactivez_rdmshuf > coactivez_rdm))./n;
% else 
%     coactivez_prdmshuf = length(find(coactivez_rdmshuf < coactivez_rdm))./n;
% end


% % Output
% % ------
out.index = index;
out.nsimul = nsimul;
out.nsimul_rdm = nsimul_rdm;
out.trialResps1 = trialResps1;
out.trialResps2 = trialResps2;
out.trialResps1_rdm = trialResps1_rdm;
out.trialResps2_rdm = trialResps2_rdm;

% out.r = r;
% out.p = p;
% out.p_shuf = p_shuf;
% out.p_shuf_2t = p_shuf_2t;
% out.r_rdm = r_rdm;
% out.p_rdm = p_rdm;
% out.p_rdmshuf = p_rdmshuf;
% out.p_rdmshuf_2t = p_rdmshuf_2t;
% out.coactivez = coactivez;
% out.coactivez_pshuf = coactivez_pshuf;
% out.coactivez_pshuf_2t = coactivez_pshuf_2t;
% out.coactivez_rdm = coactivez_rdm;
% out.coactivez_prdmshuf = coactivez_prdmshuf;
% out.coactivez_prdmshuf_2t = coactivez_prdmshuf_2t;








