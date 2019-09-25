function [pcadata, scores] = sj_pcasvd(data, nch)

%%% Shantanu: Change to match Mclust definition of PCA


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


norm = 1;    % normalize Waveforms (1) or don't normalize (0)

channel_lth = size(data,2)/nch;  %%%128/4=32
nSpikes =  size(data,1);
I = ones(nSpikes,1);


if nch>=4 %% FOR TETRODE ONLY %%
    for ch=1:nch

        w = data(:,channel_lth*(ch-1)+1:channel_lth*ch);
        if norm
            % normalize waveforms to unit L2 norm (so that only their SHAPE or
            % relative angles but not their length (energy) matters)
            l2norms = sqrt(sum(w.^2,2));
            w = w./l2norms(:,ones(1,channel_lth));
        end

        cv = cov(w);
        sd = sqrt(diag(cv))';        % row std vector
        av = mean(w);                % row mean vector
        pc = wavePCA(cv);            % get PCA eigenvectors (in columns of pc)

        wstd=(w-(I*av))./(I*sd);     % standardize data to zero mean and unit variance
        wpc = wstd*pc;               % project data onto principal component axes

        scores(:,ch) = wpc(:,1); %%% FIRST PC IS SCORES by default

        %%% SAVE FIRST 4 PCS
        for npc=1:4,
            cmd=sprintf('pcadata.pc%d(:,ch) = wpc(:,%d);',npc, npc); eval(cmd);
        end

    end

else %% If 1 ch or stereotrode - use multiple PCs from same channel

    %% SAMARS METHOD
    %data = detrend(data, 'constant');
    %[u,s,v] = svd(data, 0);
    %proj = data * v;

    %% MCLUST METHOD
    %Normalize as above
    if norm
        l2norms = sqrt(sum(data.^2,2));
        w = data./l2norms(:,ones(1,channel_lth));
    end

    cv = cov(w);
    sd = sqrt(diag(cv))';        % row std vector
    av = mean(w);                % row mean vector
    pc = wavePCA(cv);            % get PCA eigenvectors (in columns of pc)

    wstd=(w-(I*av))./(I*sd);     % standardize data to zero mean and unit variance
    wpc = wstd*pc;               % project data onto principal component axes
    %% Save First 4 PCs
    for npc=1:4,
        scores(:,npc) = wpc(:,npc); %%%
        cmd=sprintf('pcadata.pc%d = wpc(:,npc);',npc); eval(cmd);
    end
    
end




