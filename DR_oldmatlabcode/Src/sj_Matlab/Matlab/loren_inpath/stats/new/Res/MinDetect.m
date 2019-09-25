% MINDETECT:  Randomization test to determine the minimum sample size needed to 
%             detect an event occurring at a specified relative frequency, at a 
%             specified confidence level.
%
%     Usage: [Nmin,p,c] = mindetect(p0,{conf})
% 
%         p0 =    vector of relative frequencies to be detected.
%         conf =  optional vector of confidence levels [default = 0.95].
%         ----------------------------------------------------------------------
%         Nmin =  column vector of length [length(p0)*length(conf)] containing 
%                   minimum sample size estimates.
%         p =     corresponding column vector of specified relative frequencies.
%         c =     corresponding column vector of confidence levels.
%

% RE Strauss, 5/4/01

function [Nmin,p,c] = mindetect(p0,conf)
  if (nargin < 2) conf = []; end;

  if (isempty(conf))
    conf = 0.95;
  end;

  iter = 10000;                         % Number of iterations for each Nmin

  p0 = p0(:);
  conf = conf(:);

  np0 = length(p0);
  nconf = length(conf);

  Nmin = zeros(np0*nconf,1);            % Allocate output matrices
  p0_out = Nmin;
  conf_out = Nmin;
  
  nout = 0;
  for ip = 1:np0                        % Cycle thru target frequencies
    p0cur = p0(ip);

    pfind1 = mindetectf(p0(ip),1,iter);    % Set initial lower & upper bounds
    pfind100 = mindetectf(p0(ip),100,iter); %   on N


    for ic = 1:nconf                      % Cycle thru confidence levels
      confcur = conf(ic);

      Nlow = 1;                           % Initialize bounds on N
      pfindlow = pfind1;
      Nhigh = 100;
      pfindhigh = pfind100;
      
      while (pfindhigh < confcur)         % Increase upper N if needed
        Nhigh = Nhigh + 100;
        pfindhigh = mindetectf(p0(ip),Nhigh,iter);
      end;

      while ((Nhigh-Nlow) > 1)            % Binary search
        search = 0;
        Nmid = round(0.5*(Nlow+Nhigh));
        pfindmid = mindetectf(p0(ip),Nmid,iter);
        if (pfindmid < confcur)
          Nlow = Nmid;
          pfindlow = pfindmid;
        else
          Nhigh = Nmid;
          pfindhigh = pfindmid;
        end;
      end;

      nout = nout + 1;
      Nmin(nout) = ceil(mean([Nlow,Nhigh]));
      p0_out(nout) = p0cur;
      conf_out(nout) = confcur;
    end;
  end;

  p = p0_out;
  c = conf_out;

  return;


