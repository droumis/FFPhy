


function out = dfa_lickBoutSpikeCorr(index, excludeperiods, varargin)
% 


fprintf('%d %d %d %d\n',index)
reqData = {'spikes'};
for s = 1:length(reqData)
    if ~any(cellfun(@(x) strcmp(x,reqData{s}), varargin(1:2:end), 'un', 1))
        error('missing required data')
    end
end

pconf = paramconfig;
bin = 0.01; % seconds
sw1 = bin*3; % seconds
sw2 = 0.250;
rmstmax = 0.1;
rmsmincounts = 50;
tmax = .5 ;
if ~isempty(varargin)
    assign(varargin{:});
end

day = index(1);
epoch = index(2);
nt1 = index(3);
cl1 = index(4);
nt2 = index(5);
cl2 = index(6);

try
    t1 = spikes{day}{epoch}{nt1}{cl1}.data(:,1);
    t2 = spikes{day}{epoch}{nt2}{cl2}.data(:,1);
catch
    warning('couldnt load spikes %s d%de%d ntA:%d ntB:%d \n', an, day, ep);
    out = make_blank(index);
    return
end
% Get total time for data
totaleptime = diff(spikes{day}{epoch}{nt1}{cl1}.timerange);
excltime = sum(diff(excludeperiods,[],2));
T = totaleptime - excltime;
fprintf('using %.03f of %.03f sec\n', T, totaleptime);

t1cnd = t1(~logical(isExcluded(t1, excludeperiods)));
t2cnd = t2(~logical(isExcluded(t2, excludeperiods)));

% get xc
xc = spikexcorr(t1cnd,t2cnd, bin, tmax);
% normalize
normxc = xc.c1vsc2 ./ sqrt(xc.nspikes1 * xc.nspikes2);
% smooth
try
    nstd=round(sw1/(xc.time(2) - xc.time(1)));
catch
    if xc.nspikes1 == 0
        emptynt = index(3:4);
    elseif xc.nspikes2 == 0
        emptynt = index(5:6);
    end
    fprintf('no spikes %d %d %d %d\n', day, epoch, emptynt(1),emptynt(2))
    out = make_blank(index);
    return
end
g1 = gaussian(nstd, 2*nstd+1);
smthxc = smoothvect(normxc, g1);
% Get Nevents in raw xc from -200ms to 200ms
% bins = find(abs(xc.time)<=0.2);
%     Neventscorr = sum(xc.c1vsc2(bins));
% Expected probability
p1 = xc.nspikes1/T;
p2 = xc.nspikes2/T; % fr in Hz
expProb = p1*p2; % per sec

% % shuffle
% xcShf = [];
% shiftby = randi([-shuff shuff],1,1000)/1e3; %ms
% for ih = 1:numshuff
%     t1shift = sort(t1cnd+shiftby(ih));
%     xcShf = [xcShf spikexcorr(t1shift, t2cnd, bin, tmax)];
% end
% xcShfxc = cell2mat({xcShf.c1vsc2}');
% xcShfLag0 = xcShfxc(:,knnsearch(xcShf(1).time', [0]));
% xcShfmean = nanmean(xcShfxc);

% compute the excess correlation at 0 lag
exc = excesscorr(xc.time, xc.c1vsc2, xc.nspikes1, xc.nspikes2, sw1, sw2);
% compute RMS
xcrms = xcorrrms(xc.time, xc.c1vsc2, rmstmax, rmsmincounts);

out.index = index;
out.T = T;
out.xc = xc;
out.normxc = normxc;
out.smthxc = smthxc;
out.expProb = expProb;
% out.xcShfmean = xcShfmean;
% out.xcShfLag0 = xcShfLag0;
out.excesscorr = nanmean(exc); % nanmean in case there are two bins equally near lag zero
out.xcrms = xcrms;
end

function out = make_blank(index)
out.index = index;
out.T = [];
out.xc = {};
out.normxc = {};
out.smthxc = {};
out.expProb = {};
% out.xcShfmean = {};
% out.xcShfLag0 = {};
out.excesscorr = {};
out.xcrms = {};
end



