function sss_list_intenergy2 (spikes,givenclusters, parameters)

%% Shantanu Oct 2007 - do not go back to overcluster: use the original
%% interface energy, but since you are choosing the clusters to show, index
%% into them properly, so you do not have to use rawclusts

%% FOR THIS, YOU HAVE TO NORMALIZE BY RUNNING LOOP OVER GivenClusters
%% BEST WAY - AS YOU JUST HAVE TO NORMALIZE FOR YOUR CLUSTER PAIRS THEN -
%% CHECK KLEINFELD AND FEE PAPER

% waveforms = spikes.waveforms;
% times = spikes.fstimes;
% newplots = spikes.hierarchy.assigns;
% store_times=[];
% clr = ['w' 'y' 'm' 'b' 'c' 'r' 'w' 'y' 'm' 'b' 'c' 'r' 'w' 'y' 'm' 'b' 'c' 'r'];

if (isfield(spikes,'parameters'))
    tmin = spikes.parameters.min; tref = spikes.parameters.tref;
elseif (length(parameters.min)~=0)
    tmin = parameters.min; tref = parameters.tref;
else
    tmin=0.001; tref=0.002;
end

givenclusters(find(givenclusters==0))=[];
nclu=length(givenclusters);

tmin = spikes.parameters.min; tref = spikes.parameters.tref; 
swtimes=spikes.fstimes/1000;
tree=spikes.hierarchy.tree;
rawassigns= spikes.overcluster.assigns; %rawassigns(find(rawassigns==0))=[];
currassigns=spikes.hierarchy.assigns; %currassigns(find(currassigns==0))=[];
unqrawassigns=unique(rawassigns); % TO INDEX WITHIN INTERFACE_ENERGY MATRIX WHICH IS 32x32

%%%%%%%% Calculations for normalization
if (isfield(spikes.hierarchy,'interface_energyn'))
    energies = spikes.hierarchy.interface_energyn;
    strin='new'; update=1;
else
    energies = spikes.hierarchy.interface_energy;
    strin='RAW' ; update=0;
end

for i=1:length(givenclusters), 
    numpts(i) = length(find(currassigns==givenclusters(i))); %% Nspikes currently in cluster
    idxwithinmatrix(i)=find(givenclusters(i)==unqrawassigns);
end
numpts(find(numpts==0))=2; %% For clusters which have been merged, avoid getting zero size
numpts=numpts';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=figure;  axis off; 
%redimscreen_halfvert(1)
redimscreen100(101);

count=0;
for n=1:nclu-1
    currclu = givenclusters(n); 
    for cn=n+1:nclu
        count=count+1;
        numpts1=numpts(n); numpts2=numpts(cn);
        first=idxwithinmatrix(n); second=idxwithinmatrix(cn); 
        energythis=energies(first,second);
        
        if n~=cn
%             normalize = ((numpts1 * numpts2) - diag(numpts));    % Off diag: Na*Nb, On diag: Na^2-Na ...
%             normalize = normalize - diag(0.5*diag(normalize));  % ... and divide diagonal by 2
            normalize = ( numpts(n) * numpts(cn)); 
        else
            normalize = ( numpts(n)^2-numpts(n) )/2;
        end       
        norm_energythis = energythis ./ normalize;      
        
        self1en = energies(first,first); self2en = energies(second,second);
        self1norm = ( numpts(n)^2-numpts(n) )/2; self2norm = ( numpts(cn)^2-numpts(cn) )/2;
        self1 = self1en/self1norm;  self2 = self2en/self2norm; 
        
        connect_strength_this = 2*norm_energythis/self1+self2;
        
        pairs(count,:)=[currclu,givenclusters(cn)];
        energy(count)=energythis;
        normenergy(count)=norm_energythis;
        conn_str(count)=connect_strength_this;
        t1=swtimes(find(currassigns==currclu));
        t2=swtimes(find(currassigns==clusters(cn)));
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
[conn_strp,idx] = sort(conn_str); conn_strp=roundn(conn_str(idx),-2);
energyp=roundn(energy(idx),-1);
pairsp=pairs(idx,:);
normenergyp=normenergy(idx);  %normenergyp=roundn(normenergyp,-1);
allow=allow(idx);
scores=roundn(scores(idx,:),-2);


plotn=1; plotn2=1; plotn3=1;
text(0.35,1.05,['Energies: ' strin ],'FontSize', 16, 'FontWeight','bold');

for comb=count:-1:1
%     text(0,(-0.1+0.04*(comb)),[ num2str(pairsp(comb,1)) ', ' num2str(pairsp(comb,2))],'FontSize', 16, 'FontWeight','bold');
%     text(0.09,(-0.1+0.04*(comb)),['=  ' num2str(energyp(comb))],'FontSize', 16, 'FontWeight','bold');
    if plotn<=28
        text(-0.04,(1.05-0.04*(plotn)),[ num2str(pairsp(comb,1)) ', ' num2str(pairsp(comb,2))],'FontSize', 16, 'FontWeight','bold');
        text(0.05,(1.05-0.04*(plotn)),['=  ' num2str(conn_strp(comb)) ' / '  num2str(energyp(comb))  '   //   '   num2str(allow(comb)) ' / '   num2str(scores(comb,3)) ],'FontSize', 16, 'FontWeight','bold');
        plotn=plotn+1;
    elseif plotn2<=28
        text(0.35,(1.05-0.04*(plotn2)),[ num2str(pairsp(comb,1)) ', ' num2str(pairsp(comb,2))],'FontSize', 16, 'FontWeight','bold');
        text(0.35+0.09,(1.05-0.04*(plotn2)),['=   ' num2str(conn_strp(comb)) ' / '  num2str(energyp(comb))  '    //   '   num2str(allow(comb)) ' / '   num2str(scores(comb,3)) ],'FontSize', 16, 'FontWeight','bold');
        plotn2=plotn2+1;
    else
        text(0.75,(1.05-0.04*(plotn3)),[ num2str(pairsp(comb,1)) ', ' num2str(pairsp(comb,2))],'FontSize', 16, 'FontWeight','bold');
        text(0.75+0.09,(1.05-0.04*(plotn3)),['=   ' num2str(conn_strp(comb)) ' / '  num2str(energyp(comb))  '    //   '   num2str(allow(comb)) ' / '   num2str(scores(comb,3)) ],'FontSize', 16, 'FontWeight','bold');
        plotn3=plotn3+1;
    end
end

plotn

    

