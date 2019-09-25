% DIVERSF:  Objective function for diversity().  Returns diversity and evenness 
%           indices for a set of assemblages.  Returns NaN's for diversity and 
%           evenness for S<2.
%
%     Usage: [divers,evenness,S] = diversf(freq,kind,kindname,N,use_counts)
%
%         freq =      [S x C] matrix of relative or absolute frequencies for S 
%                       species and C assemblages.
%         kind =      integer 1-5 (as described above) indicating the particular 
%                       diversity measure to be returned.
%         kindname =  'shannon' or 'simpson'.
%         N =         vector (length C) of total numbers of individuals 
%                       (required for the randomization of relative frequencies).
%         use_counts = boolean flag indicating that counts are to be used rather 
%                       than relative frequencies for diversity estimate.
%         -----------------------------------------------------------------------
%         divers =    [1 x C] vector of diversity measures.
%         evenness =  [1 x C] vector of evenness measures.
%         S =         [1 x C] vector of numbers of species.
%

% Krebs, CJ. 1989. Ecological Methodology.  Harper & Row.  Chapter 10.

% RE Strauss, 2/18/00
%   2/24/00 - return NaN for S<2.

function [divers,evenness,S] = diversf(freq,kind,kindname,N,use_counts)
  get_evenness = 0;
  if (nargout > 1)
    get_evenness = 1;
  end;

  C = size(freq,2);                     % Number of assemblages
  divers = zeros(1,C);                  % Allocate output matrices
  evenness = zeros(1,C);
  S = zeros(1,C);

  for ic = 1:C                          % Cycle thru assemblages
    p = freq(:,ic);                       % Isolate current assemblage
    Nic = N(ic);                            %   and sample size

    i = find(p < eps);                    % Remove zero freqs
    if (~isempty(i))
      p(i) = [];
    end;
    S(ic) = length(p);                    % Number of species
    s = S(ic);

    if (s < 2)
      divers(ic) = NaN;
      evenness(ic) = NaN;
    else

      if (~use_counts)                      % Convert to relative frequencies
        p = p./sum(p);             
      end;

      switch (kindname)
        case 'shannon',                       % Shannon estimates
          divers(ic) = -sum(p.*log(p));

          if (get_evenness)
            evenness(ic) = divers(ic)/log(s);
          end;

          if (kind==2)                        % Modified shannon estimates
            divers(ic) = exp(divers(ic));
          end;

        case 'simpson',                       % Simpson estimates
          if (use_counts)
            divers(ic) = sum((p.*(p-1))./(Nic.*(Nic-1)));  % Pielou correction
          else
            divers(ic) = p'*p;
          end;

          if (get_evenness)
            if (use_counts)
              n = prbcount((1/s)*ones(1,s),Nic,Nic,0,1);  % Allocate counts as evenly as possible
              maxD = sum((n.*(n-1))./(Nic.*(Nic-1)));     % Max diversity for this sample size
            else
              maxD = 1./s;                                % Estimate from number of species only
            end;
            evenness(ic) = (1-divers(ic))/(1-maxD);
          end;

          if (kind==4)                        % Modified simpson estimates
            divers(ic) = 1-divers(ic);
          end;
          if (kind==5)
            divers(ic) = 1/divers(ic);
          end;
      end;
    end;
  end;

  return;


