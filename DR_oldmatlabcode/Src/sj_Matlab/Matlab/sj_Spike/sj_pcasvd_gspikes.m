function [spc] = sj_pcasvd_gspikes(data, nch, Samprange)

%  Principal Components: Called from sj_gspikestofetm.m
%%% Shantanu: Generate 3PCs from each electrode for Kkwik
% Data here is waveforms
% Samprange: if you want to look at central part of waveform only
% eg: [5:25];
% 04/21/2010


%PCASVD            Principal Components Analysis via (mean-subtracted) SVD.
%   PROJ = PCASVD(DATA), where DATA is an M x N matrix, returns the M x N
%   matrix PROJ where the (M,N)th entry is the projection of the Mth row
%   of DATA onto the Nth eigenvector of the covariance matrix formed from
%   the rows of DATA.
%
%   [PROJ,U,S,V] = PCASVD(DATA) also returns matrices U, S, V such that
%   DATA = U * S * V' and PROJ = DATA * V.
%
%   All of these computations are generally performed taking the mean over
%   all rows of the matrix DATA to be the zero vector.  This is therefore
%   enforced if it is not already the case.

% % Very simple code -- basically just a macro.
% data = detrend(data, 'constant');   % remove mean row
% [u,s,v] = svd(data, 0);             % SVD the data matrix
% proj = data * v;                    % compute (mean-subtracted) pca projections?

if nargin<2,
    nch = 4;
end

if nargin<3,
    Samprange=[];
end

channel_lth = size(data,2)/nch;  %%%128/4=32
if length(Samprange)==0,
    Samprange=1:channel_lth;
end


norm = 1;    % normalize Waveforms (1) or don't normalize (0)
nSpikes =  size(data,1);
I = ones(nSpikes,1);

%% Klustakwik / Klusters / MClust Method

cnt_fet=0;
for ch=1:nch

    w = data(:,channel_lth*(ch-1)+1:channel_lth*ch);
    w = w(:,Samprange);

    if norm
        % normalize waveforms to unit L2 norm (so that only their SHAPE or
        % relative angles but not their length (energy) matters)
        l2norms = sqrt(sum(w.^2,2));
        w = w./l2norms(:,ones(1,length(Samprange)));
    end

    % Other Possible Norms %
    % w = detrend(w, 'constant');

    cv = cov(w);
    sd = sqrt(diag(cv))';        % row std vector
    av = mean(w);                % row mean vector
    pc = sj_wavePCA(cv);            % get PCA eigenvectors (in columns of pc)

    wstd=(w-(I*av))./(I*sd);     % standardize data to zero mean and unit variance
    wpc = wstd*pc;               % project data onto principal component axes

    scores(:,ch) = wpc(:,1); %%% FIRST PC IS SCORES by default

    % RETURN 1st 3 PCs by default
    for npc=1:3,
        cnt_fet=cnt_fet+1;
        spc(:,cnt_fet) = wpc(:,npc);
        %cmd=sprintf('pcadata.pc%d(:,ch) = wpc(:,%d);',npc, npc); eval(cmd);
    end
    
end



