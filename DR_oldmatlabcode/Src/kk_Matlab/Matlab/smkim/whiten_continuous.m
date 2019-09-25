function whitened = whiten_continuous(continuous,p);
%WHITEN_CONTINUOUS Fit an autoregressive model to continuous data and filter to normalize frequencies.
%
%   [WHITENED, AR_COEFFS] = WHITEN_CONTINUOUS(CONTINUOUS,P) fits an
%   autoregressive model of order P to the continuous data in CONTINUOUS, and
%   then uses the model coefficients AR_COEFFS as an all-pole IIR prewhitening
%   filter. AR_COEFFS are estimated using the Burg method. If CONTINUOUS is a
%   struct array, then one AR_COEFFS vector is estimated across *all* of the
%   data and applied to each struct element. All elements of CONTINUOUS must
%   have the same subject, recording electrode, region, depth,
%   reference, passband and sampling rate.
%   
%   WHITENED inherits all fields of CONTINUOUS, except that the 'samples'
%   field of WHITENED contains the whitened signal. WHITENED has an additional
%   field, 'autoregressive_model_coefficients'.
%
%   Alternatively, instead of specifying a positive integer P, you can use one
%   of the following string codes to automatically select the model order:
%     'fpe': choose P which minimizes Akaike's final prediction error
%     'aic': choose P which minimizes Akaike's information criterion
%     'bic': choose P which minimizes Schwarz's Bayesian information criterion
%       (similar to AIC, but takes sample size into account so that increasing
%       the amount of data does not increase the model order)
%
%   By default, WHITEN_CONTINUOUS(CONTINUOUS) is the same as
%   WHITEN_CONTINUOUS(CONTINUOUS,'bic')
%
%   References: 
%
%   Akaike H. (1979) A Bayesian extension of the minimum AIC procedure of
%   autoregressive model fitting. _Biometrika_ 66:237-242.
%
%   Percival D.B., Walden A.T. (1993) _Spectral Analysis for Physical
%   Applications: Multitaper and Conventional Univariate Techniques_. (see
%   section 9.10, page 437).
%
%   Mitra P.P., Pesaran B. (1999) Analysis of dynamic brain imaging data.
%   _Biophysical Journal_ 76:691-708.
%
%Depends on:
%   IS_CONTINUOUS (written by smk)
%   ARBURG_CHUNKED (written by smk)
%
%Written by smk 2009 March 1.
%

  if (exist('is_continuous') ~= 2)
    error(['WHITEN_CONTINUOUS depends on m-file IS_CONTINUOUS ' ...
        '(written by smk)']);
  end
  if (exist('arburg_chunked') ~= 2)
    error(['WHITEN_CONTINUOUS depends on m-file ARBURG_CHUNKED ' ...
        '(written by smk)']);
  end

  if ~is_continuous(continuous)
    error(['CONTINUOUS does not appear to be a valid continuous data ' ...
        'struct array']);
  end

  % Verify that all elements of CONTINUOUS share the same subject, electrode,
  % channel, depth, region, reference, passband, sampling rate. (It's not
  % correct to estimate a common AR model across all data if the data are from
  % different recording conditions.)
  if (numel(unique({continuous(:).subject})) > 1) || ...
      (numel(unique([continuous(:).electrode])) > 1) || ...
      (numel(unique([continuous(:).channel])) > 1) || ...
      (numel(unique([continuous(:).depth])) > 1) || ...
      (numel(unique({continuous(:).region})) > 1) || ...
      (numel(unique([continuous(:).Fs])) > 1) || ...
      (numel(unique(cellfun(@(c) c(1),{continuous(:).reference}))) > 1) || ...
      (numel(unique(cellfun(@(c) c(2),{continuous(:).reference}))) > 1) || ...
      (numel(unique(cellfun(@(c) c(1),{continuous(:).passband}))) > 1) || ...
      (numel(unique(cellfun(@(c) c(2),{continuous(:).passband}))) > 1)
    error(['elements of CONTINUOUS struct array are not from the same ' ...
        'recording site']);
  end

  if ischar(p)
    switch p
      case 'fpe'
        % Akaike final prediction error
        criterion = @(n,p,error_variance) error_variance*(n+p+1)/(n-p-1);
      case 'aic'
        % Akaike information criterion
        criterion = @(n,p,error_variance) 2*(p+1) + n*log(error_variance); 
      case 'bic'
        % Bayesian information criterion
        criterion = @(n,p,error_variance) log(n)*(p+1) + n*log(error_variance);
      otherwise
        error('string %s is not recognized as AR model order argument',p);
    end
  elseif ~isnumeric(p) || ~isscalar(p) || ~isreal(p) || ~(round(p) == p) || ...
      ~(p > 0)
    error('AR model order must be a positive integer');
  end

  if (nargin == 1)
    % If no model order is specified, use Schwarz criterion
    p = 'bic';
  end

  % Fit AR model  
  x = cellfun(@double,{continuous(:).samples},'UniformOutput',false);
  if ischar(p)
    n = sum(cellfun(@numel,x));
    last = struct( ...
          'p',{[]}, ...
          'ar_coeffs',{[]}, ...
          'error_variance',{[]}, ...
          'score',{[]});
    next = last;
    last.p = 1;
    [last.ar_coeffs, last.error_variance] = arburg_chunked(x,last.p);
    last.score = criterion(n,last.p,last.error_variance);
    while 1
      next.p = last.p + 1;
      [next.ar_coeffs, next.error_variance] = arburg_chunked(x,next.p);
      next.score = criterion(n,next.p,next.error_variance);
      %{
      disp(sprintf('autoregressive model of order %d has score %f', ...
          next.p,next.score));
      %}
      if (next.score < last.score)
        last = next;
      else
        % We choose model order p when model order (p+1) does not further
        % decrease the score
        p = last.p;
        disp(sprintf('selected AR model of order %d',p));
        ar_coeffs = last.ar_coeffs;
        error_variance = last.error_variance;
        break;
      end
    end
  else
    ar_coeffs = arburg_chunked(x,p);
  end

  whitened = continuous;
  [whitened.autoregressive_model_coefficients] = deal(ar_coeffs);
  % Compute frequency response profile of the filter whose coefficients are the
  % AR model coefficients
  [h,w] = freqz(ar_coeffs);
  h = h/exp(min(log(abs(h))));
  % Compute coefficients for an FIR filter whose magnitude response profile is
  % half (on log-scale, e.g. square root) of the desired magnitude response
  % profile, so that the desired magnitude response profile is achieved when the
  % filter is applied once forward and once backward
  [b,a] = invfreqz(sqrt(h),w,numel(ar_coeffs)-1,0);
  % Construct FIR filter object
  Hd = dfilt.dffir(b);

  whitened = filter_continuous(whitened,Hd);

end % end main function WHITEN_CONTINUOUS

