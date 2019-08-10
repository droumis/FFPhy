function out = compute_tfbstats(tfbvarCont, dm, varargin)
pconf = paramconfig;
saveout = 1;
conds = {'swrvarCont','expvarCont', 'tfbvarCont', 'expvarCatDiff'};
epochEnv = 'wtrack';
if ~isempty(varargin)
    assign(varargin{:});
end

for ian = 1:length(tfbvarCont)
    animal = tfbvarCont(ian).animal;
    tfbstats(ian).animal = animal;
    aninfo = animaldef(animal);
    ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
    ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
    ntrodes = unique(ntrodes(:,3));
    for icorr = 1:length(conds)
        PV = eval(conds{icorr});
        pvanidx = find(strcmp({PV.animal}, animal));
        for iv = 1:length(PV(pvanidx).expvars)
            % get the var vec
            Xvar = PV(pvanidx).dm(:,iv);
            for itfb = 1:size(tfbvarCont(ian).expvars)
                tfb = tfbvarCont(ian).expvars{itfb};
                for nti = 1:length(ntrodes)
                    nt = ntrodes(nti);
                    area = ntinfo{1}{1}{nt}.area;
                    subarea = ntinfo{1}{1}{nt}.subarea;
                    cann = ntinfo{1}{1}{nt}.cannula;
                    ntxy = ntinfo{1}{1}{nt}.ntxy;
                    % get the mean pwr for each tfbox
                    itfbpwr = squeeze(tfbvarCont(ian).dm(:,itfb,nti));

                    tfbstats(ian).dmtype{nti, itfb, end+1} = conds{icorr};
                    tfbstats(ian).expvar{nti, itfb, end+1} = tfb;
                    tmp = fitlm(Xvar, itfbpwr);
                    tfbstats(ian).fitlm{nti, itfb, end+1} = tmp;
                    tfbstats(ian).P{nti, itfb, end+1} = tmp.coefTest;
                    tfbstats(ian).R{nti, itfb, end+1} = tmp.Rsquared.Ordinary;
                    tfbstats(ian).coef{nti, itfb, end+1} = tmp.Coefficients.Estimate(end);
                end
            end
        end
    end
    if saveout
        save_data(tfbstats,[pconf.andef{3},'/tfbstats/'],'tfbstats','filetail', ...
            ['_',epochEnv]); end
end