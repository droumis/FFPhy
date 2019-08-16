

function out = make_noiseEvents(lfpstack, varargin)
%{
Inputs
noizevents = mat of indices into lfpstack:.day.epoch.ripStartTimes
%}

noisevents = struct;
    noisevents(1).animal = 'JZ3';
noisevents(1).ripidx = [682:683  777:780 794:795 2924:2959 2963:3000 3013:3038 3054:3077 3083:3110 3157:3172 ...
    3184:3196 3235:3265 3269:3271 3278:3280 3345:3355 3359:3373 3378:3385 3391:3403 ...
    3406:3410 3423:3424 3454:3462 3468:3469 3484:3485 3514:3522 3527:3531 3583:3590 ...
    3599:3614 3627:3634 3717:3749 3760:3770 3811:3819 3909:3921 3927:3929 3933:3975 ...
    5285:5311 5329:5336 5346:5349 5358:5361 5366:5379 5387:5393 5398:5407 3049:3051 ...
    621:623 4001 4022 4025 4072 4258 4121 4281 4152:4581 563 1128 3079 4004 4074 3758 ...
    3759 3706 3707 237:239 3925:3927 3978:3981 3998:4000 4022:4024 4085:4088 4093:4100 ...
    4139:4143 4622:4626 4715:4718 4690:4694 4793:4795 4810:4815 4826:4828 5138:5142 ...
    4106:4108 3989:3993 4981:4983 5004:5019 5162:5165 5172:5174 5478:5503 3074:3081 ...
    3000:3009 3138:3140 3980:3984];
% RUNNING TYPHOON jz3 is cleaned and results are being run right now
noisevents(2).animal = 'JZ2';
noisevents(2).ripidx = [438 716:719];
% RUNNING DERECHO1 jz2 is cleaned and results needs to be rerun
noisevents(3).animal = 'JZ4';
noisevents(3).ripidx = [303 1067 1211 27 165 171 1631 1645 1228 768 1215 1592 1595 1581:1583 1594 1593 1694 1673 1714:1718];
% TO DO jz4 is cleaned and results needs to be reruns
noisevents(4).animal = 'JZ1'; % none
noisevents(4).ripidx = [];
% jz1 is clean.. running plots now
noisevents(5).animal = 'D13';
noisevents(5).ripidx = [];
% RUNNING VIRGA01 trying 10 as ref instead of 11 remake evertthing.. rip detection in progress
noisevents(6).animal = 'D12';
noisevents(6).ripidx = [487];
% RUNNING DERECHO2 trying 10 as mec reference.. remake evertthing.. rip detection in progress
noisevents(7).animal = 'D10';
noisevents(7).ripidx = [];

saveout = 1;
if ~isempty(varargin)
    assign(varargin{:});
end
animals = {lfpstack.animal};
for ian = 1:length(animals)
    animal = animals{ian};
%     lfpanidx = find(strcmp({lfpstack.animal}, animal));
    nzanidx = find(strcmp({noisevents.animal}, animal));
    out(ian).animal = animal;
    out(ian).dims = {'event', {'day', 'epoch', 'time'}};
    r = noisevents(nzanidx).ripidx;
    if isempty(r)
        out(ian).events = [];
        out(ian).perepoch = [];
        continue
    end
    out(ian).events = [lfpstack(ian).day(r) lfpstack(ian).epoch(r) ...
        lfpstack(ian).ripStartTime(r)];
    days = unique(out(ian).events(:,1))';
    for d = days
        dr = out(ian).events(:,1) == d;
        eps = unique(out(ian).events(dr,2))';
        for e = eps
            deridx = ismember(out(ian).events(:,1:2), [d e], 'rows');
            if out(ian).events(deridx,2) ~= e
                error('bad ep');
            elseif out(ian).events(deridx,1) ~= d
                error('bad day')
            end
            out(ian).perepoch{d}{e}{1}.starttime = out(ian).events(deridx,3);
        end
    end
    
end
if saveout
    save_data(out, 'filterframework', 'noiseEvents', 'animpos', 0);
end
end