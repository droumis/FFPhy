% CoeffDispersion: Given a set of point coordinates in P dimensions and 
%       information about a P-dimensional grid, determines the coefficient
%       of dispersion (CD=var/mean of numbers of points/cell).  Optionally
%       plots the distribution of points/cell.  For P=2, optionally plots
%       the grid and points.
%         A point sitting on a grid line is considered to be in the higher
%       of the two contiguous grid cells.
%
%     Syntax: [cd,freqs,fmean,fvar] = ...
%                 coeffdispersion(crds,gridcells,{gridmin},{gridmax},{doplots})
%
%         crds =      [N x P] matrix of point coordinates.
%         gridcells = 2-element vector of numbers of grid cells per dimensions.
%         gridmin =   optional 2-element vector of minimum values of grid axes
%                       [default = min(crds)].
%         gridmax =   corresponding vector of maximum values of grid axes
%                       [default = max(crds)].
%         doplots =   optional boolean flag indicating that plots are to be
%                       produced [default = 0].
%         ---------------------------------------------------------------------
%         cd =        coefficient of dispersion.
%         freqs =     vector of frequencies of 0:maxfreq representing the
%                       distribution of points per cell.
%         fmean =     mean frequency of points per cell.
%         fvar =      variance in frequency of points per cell.
%

% RE Strauss, 1/19/03

function [cd,freqs,fmean,fvar] = coeffdispersion(crds,gridcells,gridmin,gridmax,doplots)
  if (nargin < 1) help coeffdispersion; return; end;
  
  if (nargin < 3) gridmin = []; end;
  if (nargin < 4) gridmax = []; end;
  if (nargin < 5) doplots = []; end;
  
  [npts,p] = size(crds);
  if (p<1 | ~isintegr(p))
    error('  CoeffDispersion: invalid number of dimensions of coordinates.');
  end;
  
  if (isempty(gridmin)) gridmin = min(crds)-eps*ones(1,p); end;
  if (isempty(gridmax)) gridmax = max(crds)+eps*ones(1,p); end;
  if (isempty(doplots)) doplots = 0; end;
  
  maxcells = max(gridcells);
  gridlines = NaN*ones(maxcells+1,p);
  cells = zeros(npts,p);
  
  for ip = 1:p                            % Determine gridlines
    ng = gridcells(ip)+1;
    gridlines(1:ng,ip) = linspace(gridmin(ip),gridmax(ip),ng)';
  end;
  
  for in = 1:npts                         % Locate points within grid
    for ip = 1:p                            % Determine subscripts for dimensions
      cells(in,ip) = max(find(gridlines(:,ip)<=crds(in,ip)));
    end;
    cellvals = rowtoval(cells);
  end;
  [x,cellfreqs] = uniquef(cellvals);
  [cellfreqs,freqs] = uniquef(cellfreqs,1);
  cellfreqs = [0; cellfreqs];
  freqs = [prod(gridcells)-sum(freqs);freqs];

  [fmean,fvar] = meanwt(cellfreqs,freqs);
  cd = fvar/fmean;

  if (doplots)
    figure;
    histgramb(cellfreqs,freqs);
    box on;
    puttick(0:max(cellfreqs));
    putxlab('Number of points per grid cell');
    puttext(0.60,0.90,sprintf('Mean = %3.2f',fmean));
    puttext(0.60,0.84,sprintf('Var = %3.2f',fvar));
    puttext(0.60,0.78,sprintf('CD = %3.2f',cd));
    
    if (p==2)
      figure;
      plot(crds(:,1),crds(:,2),'kx');
      v = [gridmin(1) gridmax(1) gridmin(2) gridmax(2)];
      axis(v);
%       puttick('off','off');
      hold on;
      for i = 2:gridcells(1)                % Draw vertical grid lines
        plot([gridlines([i,1]);gridlines([i,1])],[v(3:4)],'k:');
      end;
      for i = 2:gridcells(2)                % Draw horizontal grid lines         
        plot([v(1:2)],[gridlines([i,2]);gridlines([i,2])],'k:');
      end;
      hold off;
    end;
  end;
  
  return;
  