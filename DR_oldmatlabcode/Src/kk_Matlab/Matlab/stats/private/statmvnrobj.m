function Objective = statmvnrobj(Data, Design, Param, Covar, Resid, CholCovar)
%STATMVNROBJ Log-likelihood for multivariate normal regression without missing data.

%    Copyright 2006 The MathWorks, Inc.
%    $Revision: 1.1.6.3 $ $Date: 2006/12/15 19:33:09 $

Param = Param(:);

[NumSamples, NumSeries] = size(Data);
NumParams = size(Param,1);

if nargin<6
    [CholCovar, CholState] = chol(Covar);
    if CholState > 0
        error('stats:mvregresslike:NonPosDefCov', ...
            'Covariance SIGMA is not positive-definite.');
    end
elseif rank(CholCovar)<size(Covar,1)
    error('stats:mvregresslike:NonPosDefCov', ...
          'Covariance SIGMA is not positive-definite.');
end

LogTwoPi = log(2.0 * pi);
LogDetCovar = 2.0 * sum(log(diag(CholCovar)));

Count = 0;
Objective = 0.0;

if nargin>=5 && ~isempty(Resid)
    Objective = -0.5 * sum(sum((Resid / CholCovar).^2,2));
    Count = size(Resid,1);
elseif iscell(Design)
    if numel(Design) > 1
        for k = 1:NumSamples
            if ~any(isnan(Data(k,:)))
                Count = Count + 1;
                Resid = CholCovar' \ (Data(k,:)' - Design{k} * Param);
                Objective = Objective - 0.5 * Resid' * Resid;
            end
        end
    else
        for k = 1:NumSamples
            if ~any(isnan(Data(k,:)))
                Count = Count + 1;
                Resid = CholCovar' \ (Data(k,:)' - Design{1} * Param);
                Objective = Objective - 0.5 * Resid' * Resid;
            end
        end
    end
else
    for k = 1:NumSamples
        if ~isnan(Data(k))
            Count = Count + 1;
            Resid = CholCovar' \ (Data(k) - Design(k,:) * Param);
            Objective = Objective - 0.5 * Resid' * Resid;
        end
    end
end

Objective = Objective - 0.5 * Count * (NumSeries * LogTwoPi + LogDetCovar);

if Count < 1
    Objective = NaN;
    warning('stats:mvregresslike:AllNaNData', ...
        'All observations have missing data. Cannot compute log-likelihood.');
end
