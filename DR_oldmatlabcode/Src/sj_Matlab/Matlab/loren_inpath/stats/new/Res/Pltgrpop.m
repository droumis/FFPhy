% PLTGRPOP: Set option flags for function plotgrps()
%

% RE Strauss, 5/31/99

function [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15] ...
            = pltgrpop(options,default_options)

  max_opt = 15;

  if (size(options,1)>1)        % Transpose col vector to row vector      
    options = options';
  end;

  lenopt = length(options);
  lendef = length(default_options);

  if (lenopt < lendef)          % Append any default options
    options = [options default_options(lenopt+1:lendef)];
    lenopt = lendef;
  end;

  if (lenopt < max_opt)         % Pad options vector to full length
    options = [options, zeros(1,max_opt-lenopt)];
  end;

  p1  = 0; p2  = 0; p3  = 0; p4  = 0; p5  = 0;
  p6  = 0; p7  = 0; p8  = 0; p9  = 0; p10 = 0;
  p11 = 0; p12 = 0; p13 = 0; p14 = 0; p15 = 0;

  if options( 1) p1  = 1; end;
  if options( 2) p2  = 1; end;
  if options( 3) p3  = 1; end;
  if options( 4) p4  = 1; end;
  if options( 5) p5  = 1; end;
  if options( 6) p6  = 1; end;
  if options( 7) p7  = 1; end;
  if options( 8) p8  = 1; end;
  if options( 9) p9  = 1; end;
  if options(10) p10 = 1; end;
  if options(11) p11 = 1; end;
  if options(12) p12 = 1; end;
  if options(13) p13 = 1; end;
  if options(14) p14 = 1; end;
  if options(15) p15 = 1; end;

  return;
