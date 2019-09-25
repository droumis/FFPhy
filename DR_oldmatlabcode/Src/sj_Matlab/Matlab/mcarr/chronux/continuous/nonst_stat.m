function sigma = nonst_stat(data,A,sumV,params)

% Nonstationarity test - continuous process
%
% Usage:
%
% sigma=nonst_test(data,A,sumV,params)
% Input: 
% Note units have to be consistent. See chronux.m for more information.
%       data (1d array in samples) -- required
%       A   quadratic coefficient matrix - (Compute this separately since
%       the computation is time consuming - [A,sumV]=quadcof(N,NW,order). order
%       has to < 4NW.)
%       sumV   sum of the quadratic inverse basis vectors 
%       params: structure with fields tapers, pad, Fs, fpass, err, trialave
%       -optional
%           tapers (precalculated tapers from dpss, or in the form [NW K] e.g [3 5]) -- optional. If not 
%                                                 specified, use [NW K]=[3 5]
%	        pad		    (padding factor for the FFT) - optional (can take values -1,0,1,2...). 
%                    -1 corresponds to no padding, 0 corresponds to padding
%                    to the next highest power of 2 etc.
%			      	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			      	 to 512 points, if pad=1, we pad to 1024 points etc.
%			      	 Defaults to 0.
%           Fs   (sampling frequency) - optional. Default 1.
% Output:
%       sigma   (nonstationarity index Thomson, 2000) 


if nargin < 1; error('Need data'); end;
if nargin < 2; params=[]; end;

order = length(A);
N = length(data);
%nfft=max(2^(nextpow2(N)+pad),N);
[tapers,pad,Fs]=getparams(params);
tapers=dpsschk(tapers,N,Fs); % check tapers

alpha=zeros(1,order);
for j=1:order
  alpha(j) = trace(squeeze(A(:,:,j))*squeeze(A(:,:,j)));
end;

tmp=mtfftc(data,tapers,N,Fs);
%tmp=mtfftc(data,tapers,nfft,Fs);
sigma = zeros(length(data),1);
% Pbar = sum(abs(tmp).^2,2)./sum(weights.^2,2);
Pbar=mean(abs(tmp).^2,2);
for ii=1:order
  a0=real(sum(tmp'.*(squeeze(A(:,:,ii))*tmp.')))'/alpha(ii);
  sigma=sigma+alpha(ii)*(a0./Pbar-sumV(ii)).^2;
end;

