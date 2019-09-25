function permTest(selfile)

runnew= 1;
DEBUG= 0;
nReps= 5;
if(DEBUG) nReps= 2; end

traj_pairs= [2 3 0 1];
load(selfile);
%varname= analist.varname;
%switch varname
%case  'MutualInfo-asym'
%    calcVar= 'MutualInfo';
%otherwise
%    error('unknown variable name')
%end
olddir= pwd;

nC= length(analist.rat);

if runnew
    count= 0;
    permstats= {}; realstats= {};
    if(DEBUG) nC= 2; end
    oldrat= {'bla'};
    for iC= 1:nC
        % set up data
        rat= analist.rat{iC};
        num= analist.cellnum(iC,:); d=num(1); e=num(2); t=num(3); c=num(4);
        cd(['/home/chengs/theta/' rat '/' analist.adaptid]);
        if(~strcmp(oldrat,rat)) setLocalOptions; end
        oldrat= rat;
        fmaux.selectid= analist.selectid;
        fmaux.select=[ fmaux.data2dir '/select-' fmaux.selectid];
        currInd= getSelectId(num);
        data= loadData([d e t c]);

        % set up AF 
        global adaptest
        [iniModel, fopts]= getLocalFilterModel(data);

        % set up stats
        global fmaux select adaptest
        load(fmaux.select);
        ana.a= select.a{currInd};
        ana.x= select.x{currInd};
        [sopts, ana]= getLocalStatOpts(data,ana);

        % load real stats
        load(['stats-' fmaux.selectid '.mat']);
        realstats{iC}= stats{d}{e}{t}{c};

        % set up permutation test
        if(DEBUG) fopts.niter= 1; end
        traj= analist.traj(iC,:);
        traj(2)= traj_pairs(traj(1)+1);
        global BlockInd 
        loadVar(fmaux.data2dir, 'BlockInd');
        onTraj= {};
        nB=0;
        for i=1:2
            bind{i}= BlockInd{d}{e}{traj(i)+1};
            lentraj(i)= 0;
            for j=1:length(bind{i});
                nB=nB+1;
                onTraj{nB}= find(data.traj(bind{i}(j,1):bind{i}(j,2))== traj(i)) + bind{i}(j,1)-1 ;
                lentraj(i)= lentraj(i)+ length(onTraj{nB});
            end
        end
        p(iC,:)= lentraj/sum(lentraj);

        % run permuation test
        for iR=1:nReps
            count= count+1;
            monitorProgress(count, nC*nReps);
            % permute passes
            for iB=1:nB
%                if(rand < p(iC,1)) i=1; else i=2; end
                if(rand < 0.5) i=1; else i=2; end %@@ equal probability
                data.traj(onTraj{iB})= traj(i);
            end

            % run adaptive filter algorithm
            model= adaptFilter(data, iniModel, fopts);

            % calc statistics
            permstats{iC}{iR}= calcStatistics(data, model, sopts, ana);
        end % for iR
%        save([olddir '/permstats-' analist.selectid '-' analist.adaptid], 'permstats', 'realstats');
        save([olddir '/permstats-p05-' analist.selectid '-' analist.adaptid], 'permstats', 'realstats');
    end  % for iC
%    p
%        keyboard
else
end % if runnew

cd(olddir)
