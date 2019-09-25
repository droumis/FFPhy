% MISSPCA: Estimates missing data by multiple
%          regression on principal components
%
%
%     Usage: XEST = misspca(X,C)
%
% XEST = [n x p] final matrix with estimated values for original missing data
%    X = [n x p] matrix with missing data as NaN, - Inf, Inf
%   CC = optional flag; if present and true (=1), principal components are
%        estimated on the basis of the correlation matrix, instead of the
%        covariance matrix

% Written by Joao Alves de Oliveira

function Y = misspca(X,usecorr)
  if (nargin < 2) usecorr = []; end;

  if (isempty(usecorr))
    usecorr = 0;
  end;

  indind = 1:size(X,1);                 % Augment X with final col of obs 
  X = [X indind'];                      %  sequence identifiers
  Y = X;

  [N,P] = size(X);                      
  [R,C] = find(~isfinite(X));           % Find missing values
  [uniqrow,freqrow] = uniquef(R);       % Find rows and freqs of missing values
  [uniqfreq] = uniquef(freqrow);        % Find unique values of freqs

  uniqrow = uniqrow';
  freqrow = freqrow';
  uniqfreq = uniqfreq';

  uniqsort = sort(uniqfreq);            % Sort unique freq values

  keep = 1:N;
  keep(uniqrow) = [];
  X0 = X(keep, 1:P-1);                  % Submatrix with complete data
  XTOT = X(keep, :);                    % Augmented submatric with complete data

  if (~usecorr)                         % PC1 loadings for complete data
     loadings = pcacov(X0,1);             % Covariance matrix
  else
     loadings = pcacorr(X0,1);            % Correlation matrix
  end;

  for i = 1:length(uniqsort)                      % Cycle thru freqs of missing values/obs
    missort = X((uniqrow(freqrow==uniqsort(i))),:); % Isolate obs having curr nmbr missing vals
    [rnew,cnew] = find(~isfinite(missort));    % Find missing values in submatrix

    uniqrnew = uniquef(rnew);
    uniqrnew = uniqrnew';

    rc = [rnew(:) cnew(:)];
    acumind = [];
    loadsprd = [];

    for j = 1:length(uniqrnew)
      spread = rc(rnew == uniqrnew(j),:);
      spread = [spread(:,2)'];
      acumind = [acumind; spread];
      loadsprd = [loadsprd loadings(spread)];
    end;

    loadsum = sums(loadsprd);
    [lsumsort, ord] = sort(loadsum);        % Resorts by the loadings
    tosort = missort(ord, :);               %   of missing variables
    acumind = acumind(ord, :);

    while length(acumind(:)) > 0
      [A,B] = size(acumind);

      index = [1];
      for h=2:A
        if (acumind(h,:) == acumind(1,:));  % Compares indexes of missing variables
          index = [index, h];               %   with the first row and selects 
        end;                                %   when equal
      end;                                

      newtosort = tosort(index,:);
      XTOTpca = [XTOT; newtosort];          % Append the next row(s) with
                                            %   missing observations to be
                                            %   estimated to the "complete" set

      knowobs = XTOT(:,acumind(1,:));       % Column vector(s) of variable(s)
                                            %   that are being estimated, but
                                            %   including only known obs

      vindex = ones(1,P);
      vindex(acumind(1,:)) = vindex(acumind(1,:)) - 1;
      vindex(P) = vindex(P) - 1;
      vindex = find(vindex);
      XTOTpca = XTOTpca(:,vindex);   % delete entire columns for the
                                     % variables with missing cells
      pcr = size(XTOTpca,2);

      if (~usecorr)
        [loadings2,percvar,scores] = pcacov(XTOTpca,pcr);
      else
        [loadings2,percvar,scores] = pcacorr(XTOTpca,pcr);
      end;

      y = size(XTOT,1);
      scorescmpl = scores(1:y,:);        % scores of "complete" dataset
      predind = 1:size(scores,1);
      predind(1:y) = [];
      scorespred = scores(predind,:);    % scores of rows with cells to be
                                       % predicted
      [b,stats,predict] = linregr(scorescmpl,knowobs,[],scorespred);

      newtosort(:,acumind(1,:)) = predict;
      XTOT = [XTOT; newtosort] ;

      for m = 1:length(newtosort(:,1))
        a = find(Y(:,P)==newtosort(m,P));
        Y(a,:) = newtosort(m,:);
      end;

      null = 1:A;
      null(index) = [];

      tosort = tosort(null,:);    % delete rows of missing observations
                                  % already estimated
      acumind = acumind(null,:);  % delete indices for variables already
                                  % estimated
    end; % while
  end;  % for

  Y = Y(:,1:(P-1));

  return;






