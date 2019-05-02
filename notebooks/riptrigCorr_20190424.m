



% 
% cellpairfilter = {'allcomb', '(isequal($area, ''mec'') && ($meanrate < 7))', ...
%     '(isequal($area, ''mec'') && ($meanrate < 7))'};
% 
% timefilter = {{'getriptimes', '($nripples >= 0)', [], 'cellfilter', ...
%     '(isequal($area, ''CA1''))'}};

% iterator = 'singlecellanal';




Fp = load_filter_params('riptrigspiking_corr');
Fp.animals = {'D13'};
Fp.days = [1];

% runFilterFramework = 1;

f = createfilter('animal',Fp.animals,'days', Fp.days, 'epochs', Fp.epochfilter, ...
    'cellpairs', Fp.cellpairfilter,'excludetimefilter', Fp.timefilter, ...
    'iterator', Fp.iterator);

f = setfilterfunction(f, Fp.filtfunction, {'spikes', 'linpos'});
f = runfilter(f);
% 
% %% ---------------- Get Data---------------------------------------------------
% % ---------------- Run FIlter -------------------------------------------------
% if runFilterFramework == 1    
%     F = createfilter('animal',Fp.animals,'days',Fp.days,'epochs', ... 
%         Fp.epochfilter, 'cells',Fp.cellfilter, 'excludetime', Fp.timefilter, ...
%         'iterator', Fp.iterator);
%     F = setfilterfunction(F, Fp.filtfunction, {'spikes', ...
%         Fp.eventDataLabel ,'pos','task'}, 'TF',Fp.TF,'window',Fp.window, ...
%         'binsize',Fp.binsize,'frbinsize',Fp.frbinsize,'minthresh', ...
%         Fp.minstdthresh,'maxvelocity',Fp.maxvelocity,'minvelocity', ...
%         Fp.minvelocity, 'consensus_numtets',Fp.consensus_numtets,'welldist', ...
%         Fp.welldist);
%     F = runfilter(F);
%     F.datafilter_params = Fp;
% end