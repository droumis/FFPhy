%function cmpTheta

setRoot

id= 'pf9';

sets{1}.label= 'novel1';
sets{1}.armId= 'novelArm';
sets{1}.plotcol= [1 0 0];
sets{1}.day= 1;
sets{1}.j= 2;

sets{2}.label= 'fam1';
sets{2}.armId= 'famArm';
sets{2}.plotcol= [0 0 0];
sets{2}.day= 1;
sets{2}.j= 1;


load data_ACTIVATION

load([root '/work/mvel_60s_pf9']);
r={}; s={};
for is=1:length(sets)
    armId= sets{is}.armId;
    if isempty(sets{is}); continue; end

    selectid= [id '-' armId];
    load(sprintf('%s/work/coefAbs_%s', root, sets{is}.label));
    load(sprintf('%s/data/analist-%s', root, selectid));

    cl= collectCellList(selectid);

    nc= length(coefAbs.val);
    coef{is}= coefAbs.val;

    j= sets{is}.j; day= sets{is}.day;

    icr= IC{day}{j};
    ripVar= ACTIVATION{day}{j};
    ix= 0; 
    for ic=1:nc
        i= coefAbs.ind(ic);
        rat= analist.rat{i}; d= analist.cellnum(i,1); e= analist.cellnum(i,2); 
        t= analist.cellnum(i,3); c= analist.cellnum(i,4);
        icell= find(strcmp(cl.rat, rat)' & d==cl.cellnum(:,1) &...
            e==cl.cellnum(:,2) & t==cl.cellnum(:,3) & c==cl.cellnum(:,4));
        ir= find(icr== icell);
        if isempty(ir); continue; end
        ix= ix+1;
        r{is}(ix)= ripVar(ir);
        s{is}(ix)= coef{is}(ic);
    end
end

opt.fitline=1; 

figure
opt.xstr= 'correlation coef.'; opt.ystr= 'reactivation prob.';
opt.outname= 'corr_V_reactivation';
for is=1:length(sets)
    if isempty(sets{is}); continue; end
    opt.plotcol= sets{is}.plotcol;
    linrel(s{is}, r{is}, opt);
    hold on
end
hold off

