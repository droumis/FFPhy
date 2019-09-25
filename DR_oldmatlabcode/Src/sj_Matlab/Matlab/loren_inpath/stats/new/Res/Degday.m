% DEGDAY: Converts calendar days to degree days.  Interpolates among empirical 
%         temperature measurements using either linear or cubic interpolation.
%
%     Usage: [dd,preddate] = 
%               degday(tempdate,temp,{reftemp},{preddate},{doplots},{interpol})
%
%         tempdate =  [nd x 3] matrix of sequential dates ([y m d]) on which  
%                       temperatures were measured.
%         temp =      [nc x 1] vector of corresponding temperatures.
%         reftemp =   reference temperature at or below which zero growth takes 
%                       place [default = 0].
%         preddate =  [npd x 3] matrix of dates ([y m d]) for which degree-day 
%                       values are to be predicted.  If not passed, then 
%                       'tempdate' is used for prediction.
%         doplots =   boolean flag indicating that plots are (=1) or are not 
%                       (=0) to be produced [default = 0].
%         interpol =  type of interpolation: 'l' for linear, 'c' for cubic 
%                       splines [default = 'c'].
%         ---------------------------------------------------------------------
%         dd =        vector of degree-days corresponding to 'preddate' matrix.
%         preddate =  corresponding Julian dates, accumulated across years.
%

% RE Strauss, 2/25/99
%   9/3/99 - changed plot colors for Matlab v5.

function [dd,preddate] = degday(tempdate,temp,reftemp,preddate,doplots,interpol)
  if (nargin < 2) temp = []; end;
  if (nargin < 3) reftemp = []; end;
  if (nargin < 4) preddate = []; end;
  if (nargin < 5) doplots = []; end;
  if (nargin < 6) interpol = []; end;

  if (isempty(temp))
    error('  DEGDAY: too few input arguments');
  end;
  if (size(temp,1)==1)                    % Convert row vector to col vector
    temp = temp';
  end;

  if (isempty(reftemp))                   % Default input arguments
    reftemp = 0;
  end;
  if (isempty(preddate))
    preddate = tempdate;
  end;
  if (isempty(interpol))
    interpol = 'c';
  end;
  if (isempty(doplots))
    doplots = 0;
  end;

  [tr,tc] = size(tempdate);
  [pr,pc] = size(preddate);
  if (tc~=3 | pc~=3)
    error('  DEGDAY: temperature matrices must have 3 columns');
  end;

  tempdate_save = tempdate;
  preddate_save = preddate;
  len_tempdate = length(tempdate);
  len_preddate = length(preddate);

  tempdate = julian(tempdate,1);          % Convert dates to Julian
  preddate = julian(preddate,1);

  if (preddate(1)<tempdate(1) | preddate(end)>tempdate(end))
    disp('  DEGDAY: dates for which degree-days interpolated must be within range');
    disp('          of dates for which temperatures measured');
    error(' ');
  end;

  temp_save = temp;
  date = [tempdate(1):tempdate(end)];     % Daily interpolation
  if (interpol == 'c')                      % Cubic 
    itemp = interp1(tempdate,temp,date,'cubic');
  elseif (interpol == 'l')                  % Linear 
    itemp = interp1(tempdate,temp,date,'linear');
  end;

  temp = itemp-reftemp;                   % Temperature above reference
  i = find(temp < 0);                     % Below reference become zero
  if (~isempty(i))
    temp(i) = zeros(length(i),1);
  end;
  cumtemp = cumsum(temp);                 % Accumulated degree days

  dd = NaN*ones(size(preddate));
  for i = 1:len_preddate
    j = find(date == preddate(i));
    if (~isempty(j))
      dd(i) = cumtemp(j);
    end;
  end;

  if (doplots)
    figure;
    plot(tempdate,temp_save,'ko',date,itemp,'k');
    hold on;
    plot([tempdate(1) tempdate(end)],[reftemp reftemp],'k--');
    hold off;
    putbnd(tempdate,temp_save);
    putxlab('Julian date');
    putylab('Measured temperature (C)');

    figure;
    plot(preddate,dd,'ko',preddate,dd,'k');
    putbnd(preddate,dd);
    putxlab('Julian date');
    putylab('Accumulated degree-days');
  end;

  return;