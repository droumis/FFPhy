function h = outlier(name, dataset, yl, yh, yle, yhe)

% $Revision: 1.1.8.4 $  $Date: 2004/12/06 16:38:04 $
% Copyright 2003-2004 The MathWorks, Inc.

h = stats.outlier;

if nargin==0
   h.YLow='';
   h.YHigh='';
   h.YLowLessEqual=0;
   h.YHighGreaterEqual=0;
   h.name = '';
   h.dataset = '';
else
   h.YLow=yl;
   h.YHigh=yh;
   h.YLowLessEqual=yle;
   h.YHighGreaterEqual=yhe;
   h.dataset=dataset;
   
   % assumes name is unique
   h.name = name;
end

% add it to the list of outliers
connect(h,dfswitchyard('getoutlierdb'),'up');

list(2) = handle.listener(h,findprop(h,'name'),'PropertyPostSet',...
                          {@updatename,h});
list(1) = handle.listener(h,'ObjectBeingDestroyed', {@cleanup,h});
h.listeners=list;

dfgetset('dirty',true);   % session has changed since last save


%=============================================================================
function updatename(hSrc,event,fit)

dfgetset('dirty',true);   % session has changed since last save


%=============================================================================
function cleanup(hSrc,event,fit)

dfgetset('dirty',true);   % session has changed since last save
