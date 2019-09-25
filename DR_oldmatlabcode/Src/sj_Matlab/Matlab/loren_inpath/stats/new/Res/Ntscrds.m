% NTSCRDS: Converts coordinate points from Rohlf's NTS format (output of 
%           tpsdig.exe) to a matrix of cartesian coordinates.  
%
%     Usage: [crds,sp] = ntscrds(nts_crds,{specs},{fill_crds})
%
%         nts_crds =  vector of NTS coordinates, or matrix of coordinates by row 
%                       (with [0 0] or other values to fill out last row of block).
%         specs =     optional vector of specimen numbers for sets of NTS 
%                       coordinates spanning more than one row.
%         fill_crds = optional 2-element row vector of values for 'fill-in' 
%                       coordinates in last row, to be deleted (default = [0 0]).
%         -------------------------------------------------------------------------
%         crds =      [n x 2] vector of cartesian coordinates.
%         sp =        [n x 1] vector of specimen numbers.
%

% RE Strauss, 9/24/98
%   9/15/99 - modified aspect ratio; added 'imagepro' option.
%   5/5/00 -  allow for multiple specimens.
%   5/24/00 - removed 'imagepro' option; delete fill-in coordinates.

function [crds,sp] = ntscrds(nts_crds,specs,fill_crds)
  if (nargin < 2) specs = []; end;
  if (nargin < 3) fill_crds = []; end;

  if (isempty(specs))
    specs = ones(size(nts_crds,1),1);
  end;
  if (isempty(fill_crds))
    fill_crds = [0 0];
  end;

  uspecs = uniquef(specs);                % Numbers of specimens
  nspecs = length(uspecs);

  i = find(specs == uspecs(1));           % Allocate output matrices
  c = nts_crds(i,:)';                     %   using first specimen
  c = c(:);
  nvals = length(c);
  ncrds = nvals/2;

  if (~isintegr(ncrds))
    error('  NTSCRDS: number of NTS coordinate values must be even.');
  end;

  crds = zeros(ncrds*nspecs,2);
  sp = zeros(ncrds*nspecs,1);

  crds(1:ncrds,:) = [c(1:2:nvals-1) c(2:2:nvals)];
  sp(1:ncrds) = uspecs(1)*ones(ncrds,1);

  b = ncrds+1;
  for is = 2:nspecs                       % Remaining specimens
    i = find(specs == uspecs(is));
    c = nts_crds(i,:)';
    c = c(:);

    e = b+ncrds-1;
    crds(b:e,:) = [c(1:2:nvals-1) c(2:2:nvals)];
    sp(b:e) = uspecs(is)*ones(ncrds,1);
    b = e+1;
  end;

  i = find(crds(:,1)==fill_crds(1) & crds(:,2)==fill_crds(2));
  crds(i,:) = [];
  sp(i,:) = [];

  return;
