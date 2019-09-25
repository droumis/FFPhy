% pFOURIER:  Fourier decomposition of a quadratic-spline smoothing of a 2D 
%           boundary that is specified by a series of point coordinates.
%           The decomposition depends on choice of center, which by default is 
%           the centroid.  (Another possible choice would be the center of the 
%           best-fitting circle).
%             Produces plots of radius function, fitted outline, and power 
%           spectrum.
%
%     Usage: [a,b,phase,ampl] = ...
%                 pfourier(crds,{nfc},{start},{center},{nharm},{noplot})
%
%           crds =    [n x 2] matrix of point coordinates.
%           nfc =     optional number of pairs of Fourier coefficients to be 
%                       returned, to a max of n [default = n].
%           start =   optional index (subscript) of starting point 
%                       [default = first point in 'crds'].
%           center =  optional index of center point in 'crds', or 2-element 
%                       vector of center point coordinates [default = centroid].
%           nharm =   optional value indicating the number of Fourier components 
%                       to be plotted for the smoothed outline [default = n].
%           noplot =  optional flag indicating, if true, that plots are not to
%                       be produced.  
%           --------------------------------------------------------------------
%           a,b =   Fourier coefficients (a0,b0), (a1,b1),..., 
%                       (a(nfc-1),b(nfc-1)), where the a's are the coefficients 
%                       of the cosine terms and the b's are the coefficients 
%                       of the sine terms.
%           phase = phase-angle coefficient.
%           ampl =  amplitude coefficient.
%

% RE Strauss, 3/5/00
%   3/13/00 - added phase angle and amplitude coefficients; fixed plot.
%   3/18/00 - replaced cubic spline with quadratic spline; 
%             allowed for specification of start and center points.
%   3/26/00 - fixed problems with display of Fourier function.
%   3/12/01 - plot power function.
%   4/1/02 -  removed scaling code to fixed display of function, which was not working;
%             changed call to quadspline().

function [a,b,phase,ampl] = pfourier(crds,nfc,start,center,nharm,noplot)
  if (nargin < 2) nfc = []; end;
  if (nargin < 3) start = []; end;
  if (nargin < 4) center = []; end;
  if (nargin < 5) nharm = []; end;
  if (nargin < 6) noplot = []; end;

  [n,p] = size(crds);
  if (p~=2)
    error('  FOURIER: 2-dimensional input coordinates only.');
  end;

  if (isempty(nfc))
    nfc = n;
  end;
  if (isempty(noplot))
    noplot = 0;
  end;
  if (isempty(nharm))
    nharm = n;
  end;
  if (isempty(start))
    start = 1;
  end;

  if (start > n)
    error('  FOURIER: start index out of range.');
  end;

  if (nharm < 1)
    hnarm = 1;
  elseif (nharm > n)
    nharm = n;
  end;

  if (isempty(center))                    % Center point
    [a,p,center] = polyarea(crds);
  else
    if (length(center)==1)                  % Index to point in 'crds'
      centerpt = crds(center,:);
      crds(center,:) = [];
      n = n-1;
    else
      if (length(center)~=2)
        error('  FOURIER: invalid center-point specification');
      end;
    end;
  end;

  if (start > 1)                              % Shift around start point
    if (eucl(crds([1,n],:)) < (1e-6)*max(range(crds)))  % Open the polygon
      closed = 1;
      crds(n,:) = [];
      n = n-1;
    end;
    crds = [crds(start:n,:); crds(1:start-1,:)];      % Shift
  end;

  if (eucl(crds([1,n],:)) > (1e-6)*max(range(crds)))  % Close the polygon
    crds = [crds; crds(1,:)];
    n = n+1;
  end;

  scrds = quadspline(crds,1,[],0,256);                % Quadratic spline
  [r,theta,perim] = radiusfn(scrds,[],center,1,1);  % Radius function, by angle
  r(length(r)) = [];                                % Stop short of re-wrap
  fx = fft(r);                                      % Fourier decomposition

  a = zeros(nfc,1);                                 % Allocate output matrices
  b = zeros(nfc,1);

  a(1) = fx(1)/n;                                   % Fourier coefficients
  if (nfc > 1)
    for k = 2:nfc
      a(k) = 2*real(fx(k+1))/n;
      b(k) = -2*imag(fx(k+1))/n;
    end;
  end;
  
  phase = atan(b./a);                     % Phase angles
  ampl = sqrt(a.*a + b.*b);               % Amplitudes

  if (~noplot)                            % Optional plot
    theta = linspace(0,2*pi)';              % One cycle
    rmean = mean(r);

    rp = a(1)*ones(size(theta));             % Accumulate harmonics
    for i = 2:nharm                         
      rp = rp + a(i)*cos(i*theta) + b(i)*sin(i*theta);
    end;

    rp = rp.*rmean./mean(rp);                % Adjust radius to observed value

%     theta = theta + angl(center+[1 0],center,crds(1,:),1);  % Rotate to match crds
%     [x,y] = polarcrd(rp,theta,1);            % Retransform to cartesian crds
%     nscrds = length(x);
%     scrds = [x y] + ones(nscrds,1)*center;
% 
%     rcrds = crds;
%     rcrds(n,:) = [];
%     crd_theta = angl(ones(n-1,1)*(center+[1 0]),center,rcrds,1); % Get angles for crds
%     acrds = polyangl(scrds,center,crd_theta); % Get intersects with Fourier fn
%     dev = mean(rcrds - acrds);
%     scrds = scrds + ones(nscrds,1)*dev;
    
     hold on;
    plot(crds(:,1),crds(:,2),'o');
    plot(center(1),center(2),'kx');
    plot(crds(1,1),crds(1,2),'k*');
    plot(crds(:,1),crds(:,2),'k:');
    plot(scrds(:,1),scrds(:,2),'k');
      hold off;
    
    
    
    sqplot([crds;scrds]);
    
  end;

  return;

