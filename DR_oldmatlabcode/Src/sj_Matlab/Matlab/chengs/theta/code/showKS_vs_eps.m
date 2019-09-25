%function showKS_eps(rat, num)
%function showKS_eps(rat, num)
%
% Show KS statistic vs learning rate (eps) for a given cell (num)

makegraphs= 0;

%rat= 'kyl';
%num= [4 4 28 6];
%num= [1     4    25     6];

%rat= 'ter';
%num= [1     4    26     4];

num= [];
rat= [];

adapt_dir= {
'adapt_xp015_t015', 
%'adapt_xp1_t015', 
'adapt_xp2_t015', 
'adapt_xp3_t015', 
'adapt_xp4_t015', 
'adapt_xp5_t015', 
'adapt_xp20_t015', 
%'AS_xp015_t015', 
%'AS_xp1_t015', 
%'AS_xp2_t015', 
%'AS_xp3_t015', 
%'AS_xp4_t015', 
%'AS_xp5_t015', 
%'AS_xp20_t015', 
};
select_id= {
%'placefields4-satellite',
'placefields4-train',
};
rats={'kyl', 'ter', 'sta', 'fel'};

nA= length(adapt_dir);
nS= length(select_id);
olddir= pwd;
global adaptest 
eps_xp= nan*ones(nA,1);
eps_t= nan*ones(nA,1);

if isempty(num)
    % determine number of cells in selection
    ncells= [];
    for iR=1:4
        cd(['/bach4/' rats{iR} '/' adapt_dir{1}]);
        fmaux= {};
        setLocalOptions;
        fmaux.selectid= select_id{1};
        fmaux.select=[ fmaux.data2dir '/select-' fmaux.selectid];
        [d,e,t,c,ncells(iR)]= startCellList;
    end
    nS= 1;
    KS= nan*ones(nA, sum(ncells));
    for iA=1:nA
        for iR=1:4
            cd(['/bach4/' rats{iR} '/' adapt_dir{iA}]);
            fmaux= {};
            setLocalOptions;

            for iS=1:nS
                % select file
                fmaux.selectid= select_id{iS};
                fmaux.select=[ fmaux.data2dir '/select-' fmaux.selectid];

                icell= sum(ncells(1:iR-1))+1;
                [d,e,t,c,nc]= startCellList;
                loadVar('.','adaptest',d);
                eps_xp(iA)= adaptest{d}{e}{t}{c}.filter.eps(1);
                eps_t(iA)= adaptest{d}{e}{t}{c}.filter.eps(end);


                while ~isempty(d)
                    [viol KS(iA,icell)]= auxCalcKS(d,e,t,c, 1, makegraphs);
                    icell= icell+1;
                    [d,e,t,c]= getNextCell;
                end
            end
        end
    end

else
    KS= nan*ones(nA, nS);
    d= num(1); e= num(2); t= num(3); c= num(4);

    for iA=1:nA
        cd(['/bach4/' rat '/' adapt_dir{iA}]);
        fmaux= {};
        setLocalOptions;

        for iS=1:nS
            % select file
            fmaux.selectid= select_id{iS};
            fmaux.select=[ fmaux.data2dir '/select-' fmaux.selectid];

            loadVar('.','adaptest',d);
            eps_xp(iA)= adaptest{d}{e}{t}{c}.filter.eps(1);
            eps_t(iA)= adaptest{d}{e}{t}{c}.filter.eps(end);

            [viol KS(iA,iS)]= auxCalcKS(d,e,t,c, 1, makegraphs);

        end
    end
end



figure
plot(eps_xp, KS, '.-')
mKS=mean(KS,2);
figure
plot(eps_xp, mKS)

minEps= eps_xp(imin(KS));
figure
hist(minEps,20);
nEps= [];
for ie=1:length(eps_xp)
    nEps(ie)= sum(eps_xp(ie)== minEps);
end
bar(eps_xp, nEps)

cd(olddir)
