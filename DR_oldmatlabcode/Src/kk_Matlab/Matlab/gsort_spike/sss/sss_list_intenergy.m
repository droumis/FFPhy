function sss_list_intenergy (spikes,givenclusters, parameters)

%% Shantanu Oct 2007 - do not go back to overcluster: use the 32X32
%% UPDATED interface energy matrix, but since you are choosing the clusters to show, index
%% into them properly, so you do not have to use rawclusts
%% [FOR THIS, YOU HAVE TO NORMALIZE BY RUNNING LOOP OVER NUMPTS, CANNOT USE
%% MATRIX STRUCTURE DIRECTLY UNLESS YOU USE SPARSE LIKE SAMAR DOES
%% IN SSS_AGGREGATE]
%% HAVE NOW IMPLEMENTED THE MATRIX NORMALIZE CALCULATION CORRECTLY HERE, AND
%% IN 2ND METHOD, HAVE RUN LOOP OVER EACH CLUSTER PAIR IN SSS_LIST_INTENERGY2

% waveforms = spikes.waveforms;
% times = spikes.fstimes;
% newplots = spikes.hierarchy.assigns;
% store_times=[];
% clr = ['w' 'y' 'm' 'b' 'c' 'r' 'w' 'y' 'm' 'b' 'c' 'r' 'w' 'y' 'm' 'b' 'c' 'r'];

givenclusters(find(givenclusters==0))=[];

if (isfield(spikes,'parameters'))
    tmin = spikes.parameters.tmin; tref = spikes.parameters.tref; 
elseif (length(parameters.tmin)~=0)
    tmin = parameters.tmin; tref = parameters.tref; 
else
    tmin=0.001; tref=0.002;
end


swtimes=spikes.fstimes/1000;
tree=spikes.hierarchy.tree;
rawassigns= spikes.overcluster.assigns; %rawassigns(find(rawassigns==0))=[];
currassigns=spikes.hierarchy.assigns; %currassigns(find(currassigns==0))=[];

%%%%%%%% Calculations for normalization
if (isfield(spikes.hierarchy,'interface_energyn'))
    energies = spikes.hierarchy.interface_energyn;
    strin='new'; update=1;
else
    energies = spikes.hierarchy.interface_energy;
    strin='RAW' ; update=0;
end
nclu=length(givenclusters);
%ncomb=nclu*(nclu-1)/2;

% if (isfield(spikes.hierarchy,'connect_strength'))
%     connect_strength = spikes.hierarchy.connect_strength;
% else
    currclusts = unique(currassigns); currclusts(find(currclusts==0))=[];
    rawclusts = unique(rawassigns); rawclusts(find(rawclusts==0))=[];
    numclusts =length(rawclusts);
    %numpts = full(sparse(assignments, 1, 1, numclusts, 1));
    for i=1:length(rawclusts), 
        numpts(i) = length(find(currassigns==rawclusts(i))); %% If ant rawclust has been merged, it auto goes to zero size in currassigns
    end
    %for i=1:length(rawclusts), numpts(i) = length(find(rawassigns==rawassigns(i))); end
    numpts(find(numpts==0))=2; %% For clusters which have been merged, avoid getting zero size
    numpts=numpts';

    normalize = ((numpts * numpts') - diag(numpts));    % Off diag: Na*Nb, On diag: Na^2-Na ...
    normalize = normalize - diag(0.5*diag(normalize));  % ... and divide diagonal by 2
    norm_energies = energies ./ normalize;

    self = repmat(diag(norm_energies), [1,numclusts]);
    connect_strength = 2 .* norm_energies ./ (self + self');
    connect_strength = connect_strength .* (1-eye(numclusts));  % diag entries <- 0, so we won't agg clusters with themselves
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=figure;  axis off; 
%redimscreen_halfvert(1)
redimscreen100(101);

count=0;
for n=1:nclu-1
    currclu = givenclusters(n); 
    for cn=n+1:nclu
        count=count+1;
        pairs(count,:)=[currclu,givenclusters(cn)];
        energy(count)=energies(currclu,givenclusters(cn));
        normenergy(count)=norm_energies(currclu,givenclusters(cn));
        conn_str(count)=connect_strength(currclu,givenclusters(cn));
        t1=swtimes(find(currassigns==currclu));
        t2=swtimes(find(currassigns==givenclusters(cn)));
        [allow(count), scores(count,:)] = isiQuality(t1, t2, tmin, 0.010, tref, spikes.Fs);
    end
end

% %% Order by un-normalized energy
% [energyp,idx] = sort(energy); energyp=roundn(energyp,-1);
% pairsp=pairs(idx,:);
% normenergyp=normenergy(idx);  %normenergyp=roundn(normenergyp,-1);
% conn_strp=roundn(conn_str(idx),-2);
% allow=allow(idx);
% scores=roundn(scores(idx,:),-2);

%% Order by Connection Strength
[conn_strp,idx] = sort(conn_str); 
conn_strp=roundn(conn_str(idx),-2); % roundn exists ONLY in mapping toolbox in matlab - Make your own!
conn_strp=roundn(conn_str(idx),-2);
energyp=roundn(energy(idx),-1);
pairsp=pairs(idx,:);
normenergyp=normenergy(idx);  normenergyp=roundn(normenergyp,-1);
allow=allow(idx);
scores=roundn(scores(idx,:),-2);


plotn=1; plotn2=1; plotn3=1;
text(0.35,1.05,['Energies: ' strin ],'FontSize', 16, 'FontWeight','bold');

for comb=count:-1:1
%     text(0,(-0.1+0.04*(comb)),[ num2str(pairsp(comb,1)) ', ' num2str(pairsp(comb,2))],'FontSize', 16, 'FontWeight','bold');
%     text(0.09,(-0.1+0.04*(comb)),['=  ' num2str(energyp(comb))],'FontSize', 16, 'FontWeight','bold');
    if plotn<=28
        text(-0.04,(1.05-0.04*(plotn)),[ num2str(pairsp(comb,1)) ', ' num2str(pairsp(comb,2))],'FontSize', 16, 'FontWeight','bold');
        text(0.05,(1.05-0.04*(plotn)),['=  ' num2str(conn_strp(comb)) ' / '  num2str(normenergyp(comb))  '   //   '   num2str(allow(comb)) ' / '   num2str(scores(comb,3)) ],'FontSize', 16, 'FontWeight','bold');
        plotn=plotn+1;
    elseif plotn2<=28
        text(0.35,(1.05-0.04*(plotn2)),[ num2str(pairsp(comb,1)) ', ' num2str(pairsp(comb,2))],'FontSize', 16, 'FontWeight','bold');
        text(0.35+0.09,(1.05-0.04*(plotn2)),['=   ' num2str(conn_strp(comb)) ' / '  num2str(normenergyp(comb))  '    //   '   num2str(allow(comb)) ' / '   num2str(scores(comb,3)) ],'FontSize', 16, 'FontWeight','bold');
        plotn2=plotn2+1;
    else
        text(0.75,(1.05-0.04*(plotn3)),[ num2str(pairsp(comb,1)) ', ' num2str(pairsp(comb,2))],'FontSize', 16, 'FontWeight','bold');
        text(0.75+0.09,(1.05-0.04*(plotn3)),['=   ' num2str(conn_strp(comb)) ' / '  num2str(normenergyp(comb))  '    //   '   num2str(allow(comb)) ' / '   num2str(scores(comb,3)) ],'FontSize', 16, 'FontWeight','bold');
        plotn3=plotn3+1;
    end
end

plotn

    

