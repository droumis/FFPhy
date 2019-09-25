function [spc] = pcasvd_matspikes(waves, nch, Samprange)

%% Shantanu: Principal Components: Usually called from nmakeparams_pc.m / addpc.m, 
% Samprange: if you want to look at central part of waveform only
% eg: [4:30];

%%
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

%%
if nargin<2,
    nch = size(waves,2);
end

if nargin<3,
    Samprange=[];
end

channel_lth = size(waves,1);  % usually 40, usually 28 after de-jittering
if length(Samprange)==0,
    Samprange=1:channel_lth;
end



%% Very simple code -- basically just a macro.
% data = detrend(data, 'constant');   % remove mean row
% [u,s,v] = svd(data, 0);             % SVD the data matrix
% proj = data * v;                    % compute (mean-subtracted) pca projections?

%% Klustakwik / Klusters Method

norm = 0;    % normalize Waveforms (1) or don't normalize (0)

nSpikes =  size(waves,3);
I = ones(nSpikes,1);

cnt_fet=0;
for ch=1:nch
    
    w = double(squeeze(waves(:,ch,:)))';
    w = w(:,Samprange);
     
    if norm
        % normalize waveforms to unit L2 norm (so that only their SHAPE or
        % relative angles but not their length (energy) matters)
        l2norms = sqrt(sum(w.^2,2));
        w = w./l2norms(:,ones(1,length(Samprange)));
    end
    
    cv = cov(w);
    sd = sqrt(diag(cv))';        % row std vector
    av = mean(w);                % row mean vector
    pc = wavePCA(cv);            % get PCA eigenvectors (in columns of pc)
    
    wstd=(w-(I*av))./(I*sd);     % standardize data to zero mean and unit variance
    wpc = wstd*pc;               % project data onto principal component axes
    
    % RETURN 1st 3 PCs by default
    for npc=1:3,
        cnt_fet=cnt_fet+1;
        spc(:,cnt_fet) = wpc(:,npc);
        %cmd=sprintf('pcadata.pc%d(:,ch) = wpc(:,%d);',npc, npc); eval(cmd);
    end
    
end

