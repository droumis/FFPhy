% BOOTTIME: Computes and displays the elapsed time of a bootstrapped analysis,
%           given the start and stop clock times, determined from calls to
%           clock().
%
%     Usage: boottime(time_in, time_out)
%

% RE Strauss, 4/9/94

function boottime(time_in, time_out)
  et = etime(time_in, time_out);
  min = floor(et/60);
  sec = floor(et-(min*60));
  if (min < 1)
    disp(sprintf('  Bootstrap time: %5.0f seconds\n', sec));
  else
    disp(sprintf('  Bootstrap time: %5.0f minutes %5.0f seconds\n', ...
                 min, sec));
  end;
  return;
