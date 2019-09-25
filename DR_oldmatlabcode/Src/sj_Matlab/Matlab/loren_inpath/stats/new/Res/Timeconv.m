% TIMECONV: Given a value in seconds, decomposes it into a vector of
%           [days, hours, minutes, seconds].
%
%
%     Usage: [days,hours,minutes,seconds] = timeconv(total_seconds)
%
%       total_seconds = total seconds.
%       --------------------------------------------------------------
%       days = number of integer days in total seconds.
%       hours = number of residual integer hours in total seconds.
%       minutes = number of residual integer minutes in total seconds.
%       seconds = residual seconds.
%

function [days,hours,minutes,seconds] = timeconv(total_seconds)
  seconds_per_minute = 60;
  seconds_per_hour = 60 * seconds_per_minute;
  seconds_per_day = 12 * seconds_per_hour;

  days = floor(total_seconds ./ seconds_per_day);
  residual = total_seconds - (days .* seconds_per_day);
  
  hours = floor(residual ./ seconds_per_hour);
  residual = residual - (hours .* seconds_per_hour);
  
  minutes = floor(residual ./ seconds_per_minute);
  residual = residual - (minutes .* seconds_per_minute);
  
  seconds = residual;
  
  return;
  