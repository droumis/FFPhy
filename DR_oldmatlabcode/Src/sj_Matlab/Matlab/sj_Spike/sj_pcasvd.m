function [pcadata] = sj_pcasvd(data, nch, pcatype)

%%% PCA- Used in gsort
% Allow possibility of L2 norm as in Mclust

% Samar's PCA Calculation
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

% Samar's PCA calculation
% Very simple code -- basically just a macro.
% data = detrend(data, 'constant');   % remove mean row
% [u,s,v] = svd(data, 0);             % SVD the data matrix
% proj = data * v;                    % compute (mean-subtracted) pca projections?

if nargin<3,
    pcatype=1;
end

norm = 1;    % normalize Waveforms - L2 norm calculation (1) or don't normalize (0)

channel_lth = size(data,2)/nch;  %%%128/4=32
nSpikes =  size(data,1);
I = ones(nSpikes,1);


if pcatype==0
    % Samar's PCA calculation
    % Very simple code -- basically just a macro.
    % For Any Number Of Channels
    cntpc=1;
    for ch=1:nch
        w = data(:,channel_lth*(ch-1)+1:channel_lth*ch);
        w = detrend(w, 'constant');   % remove mean row
        [u,s,v] = svd(w, 0);             % SVD the data matrix
        wpc = w * v;                    % compute (mean-subtracted) pca projections?
        
        %Save PCs for congruency with Klusters - 1st 3 PCs of each channel
        pcadata.scores(:,cntpc) = wpc(:,1);
        pcadata.scores(:,cntpc+1) = wpc(:,2);
        pcadata.scores(:,cntpc+2) = wpc(:,3);
        cntpc=cntpc+3;
        
        %%% SAVE FIRST 4 PCS in structure
        for npc=1:4,
            cmd=sprintf('pcadata.pc%d(:,ch) = wpc(:,%d);',npc, npc); eval(cmd);
        end
        
    end
    
else  % pcatype=1
    
    cntpc=1;
    %if nch>=2 %% FOR Ntrode >=2 %%
    
    for ch=1:nch
        
        w = data(:,channel_lth*(ch-1)+1:channel_lth*ch);
        
        %% MCLUST NORM:
        if norm
            % normalize waveforms to unit L2 norm (so that only their SHAPE or
            % relative angles but not their length (energy) matters)
            l2norms = sqrt(sum(w.^2,2));
            w = w./l2norms(:,ones(1,channel_lth));
        end
        
        %% USUAL NORM %%
        %w = detrend(w, 'constant');
        
        cv = cov(w);
        sd = sqrt(diag(cv))';        % row std vector
        av = mean(w);                % row mean vector
        pc = sj_wavePCA(cv);            % get PCA eigenvectors (in columns of pc)
        
        wstd=(w-(I*av))./(I*sd);     % standardize data to zero mean and unit variance
        wpc = wstd*pc;               % project data onto principal component axes
        
        %         scores(:,ch) = wpc(:,1); %%% FIRST PC IS SCORES by default
        %         scores(:,nch+ch) = wpc(:,2); %%% 2nd pc saved after all 1st PCs
        %         scores(:,2*nch+ch) = wpc(:,3); %%% 3rd pc saved after all 1st and 2nd PCs
        
        %Save PCs for congruency with Klusters - 1st 3 PCs of each channel
        
        pcadata.scores(:,cntpc) = wpc(:,1);
        pcadata.scores(:,cntpc+1) = wpc(:,2);
        pcadata.scores(:,cntpc+2) = wpc(:,3);
        cntpc=cntpc+3;
        
        
        %%% SAVE FIRST 4 PCS in structure
        for npc=1:4,
            cmd=sprintf('pcadata.pc%d(:,ch) = wpc(:,%d);',npc, npc); eval(cmd);
        end
        
    end
    
    %end
    
end




% if nch==1 %% FOR SINGLE ELECTRODE ONLY %%
% %% If 1 ch - use multiple PCs from same channel
%
%     %% SAMARS METHOD
%     %data = detrend(data, 'constant');
%     %[u,s,v] = svd(data, 0);
%     %proj = data * v;
%
%     %% MCLUST METHOD
%     %Normalize as Mclust
%     if norm
%         l2norms = sqrt(sum(data.^2,2));
%         w = data./l2norms(:,ones(1,channel_lth));
%     end
%
%     cv = cov(w);
%     sd = sqrt(diag(cv))';        % row std vector
%     av = mean(w);                % row mean vector
%     pc = sj_wavePCA(cv);            % get PCA eigenvectors (in columns of pc)
%
%     wstd=(w-(I*av))./(I*sd);     % standardize data to zero mean and unit variance
%     wpc = wstd*pc;               % project data onto principal component axes
%     %% Save First 4 PCs
%     for npc=1:4,
%         scores(:,npc) = wpc(:,npc); %%%
%         cmd=sprintf('pcadata.pc%d = wpc(:,npc);',npc); eval(cmd);
%     end
%
% end


