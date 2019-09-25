% CHIRPS: Estimate temperature from rate of cricket chirping.
%
%     Usage:  Temp = chirps(ChirpsPerMin)
%
%         ChirpsPerMin = number of chirps per minute.
%         -----------------------------------------------
%         Temp = estimated temperature in degree Celcius.
%

function Temp = chirps(ChirpsPerMin)

  Temp = (ChirpsPerMin+40)./4;

  return;
  