function showEEGana(id)
%function showEEGana(mtm)
%
% Analyze EEG power spectrum.

%id= 'all';
mtm=1;

opt.plot= 'cdf'; opt.legend= 0; opt.size= 1.5*[1 2/3];

plots{1}.outname= 'day1';
plots{1}.plotcol= [[0 0 0]; [1 0 0]];
plots{1}.range= [1 2];

plots{2}.outname= 'day2';
plots{2}.plotcol= [[0 0 0]; [0 1 0]];
plots{2}.range= [3 4];

plots{3}.outname= 'day3';
plots{3}.plotcol= [[0 0 0]; [0 0 1]];
plots{3}.range= [5 6];

sets{1}.selectid= 'pf9-famArm';
sets{1}.novelDay= 1;
sets{1}.label= 'fam1';

sets{2}.selectid= 'pf9-novelArm';
sets{2}.novelDay= 1;
sets{2}.label= 'novel1';

sets{3}.selectid= 'pf9-famArm';
sets{3}.novelDay= 2;
sets{3}.label= 'fam2';

sets{4}.selectid= 'pf9-novelArm';
sets{4}.novelDay= 2;
sets{4}.label= 'novel2';

sets{5}.selectid= 'pf9-famArm';
sets{5}.novelDay= 3;
sets{5}.label= 'fam3';

sets{6}.selectid= 'pf9-novelArm';
sets{6}.novelDay= 3;
sets{6}.label= 'novel3';



newsets= {}; nsets= 0;
for is=1:length(sets)
    if ~isempty(sets{is}); nsets= nsets+1; newsets{nsets}= sets{is}; end
end
sets= newsets;

%vars={'theta', 'width', 'max', 'stdmax'};
%vlabels={'peak freq (Hz)', 'half width (Hz)', 'theta power', 'std peak height'};
vars={'theta', 'width', 'max'};
vlabels={'peak freq (Hz)', 'half width (Hz)', 'theta power'};
%vars={'nwinall', 'nwin', 'fracnwin', 'Twinall', 'Twin', 'fracTwin'};
%vlabels={'# windows', '# windows > 1s', 'fraction of windows > 1s', 'win length(s)', 'win length > 1s (s)', 'time fraction of windows > 1s'};
nvars= length(vars);
%if(mtm)
%    nvars= 4;
%else 
%    nvars= 3;
%end
nvx= ceil(nvars/2);
nvy= 2;

plotcol= {'k', 'r', 'b', 'g'};

setRoot;

for is=1:nsets
    selectid= sets{is}.selectid;

    load(sprintf('%s/data/tetrodes-%s', root, selectid));
    tet= tetrodes;
    ntet= length(tet.rat);
    %    load(sprintf('%s/data/epochs-%s', root, selectid));
    %    ep= epochs;
    %    nep= length(ep.rat);

    if isfield(sets{is}, 'novelDay') & sets{is}.novelDay
        ind= find(sets{is}.novelDay== tet.day);
    else
        ind= [1:ntet]';
    end
%    keyboard

    if strfind(selectid, 'novelArm')
        novelStr= 'novelArm';
    elseif strfind(selectid, 'famArm')
        novelStr= 'famArm';
    else
        novelStr= 'outerArm';
    end
   

    if(mtm)
        anafile= sprintf('%s/work/psana-%s-%s-%s',root,id,selectid, novelStr);
    else
        anafile= sprintf('%s/work/psana-pwelch-%s-%s',root,selectid, novelStr);
    end
    load(anafile);

    for iv=1:nvars
        tmp= ana.(vars{iv})(ind);
        V{iv}{is}= tmp(isfinite(tmp));
%        V{iv}{is}= ana.(vars{iv});
    end
    fprintf(1, '%s: %d/%d EEG tetrodes\n', selectid, sum(isfinite(tmp)),length(ind));
end

for ip=1:length(plots)
    if isempty(plots{ip}); continue; end

    opt.outname= [id '_' plots{ip}.outname];
    opt.plotcol= plots{ip}.plotcol;

    opt.ident= 1;
    opt.ref= 1;
    %opt.plot= 'cdf';
    opt.labels= {};
    for is=1:length(plots{ip}.range)
        labels{is}= sets{plots{ip}.range(is)}.label;
    end

    for iv=1:nvars
    %    subplot(nvx,nvy,iv);
    %    opt.varname= vars{iv};
        opt.varname= vlabels{iv};
        if iv==4; 
            opt.labels= labels;
        else
            opt.labels= {};
        end
        opt.title= ['EEG-' vars{iv}]; opt.label= vlabels{iv};
        opt.plot= 'cdf-p';
        auxCmpDist2({V{iv}{plots{ip}.range}}, {sets{plots{ip}.range}}, opt)
    end
    %myprint(5*[1 2/3], 'EEGPowerAna');
    %    set(gca, 'xlim', [1.5,4])
    %   myprint('mini', 'cdf_EEG-width_day-1')
end
