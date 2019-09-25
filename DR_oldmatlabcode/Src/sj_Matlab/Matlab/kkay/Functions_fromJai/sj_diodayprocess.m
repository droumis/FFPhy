function sj_diodayprocess(daydirect, animdirect, fileprefix, daynum, varargin)

% function dioprocess(daydirect, animdirect, fileprefix, daynum, varargin)
%
%DAYDIRECT -- folder name where the day's raw data is stored
%ANIMDIRECT -- the path to where the animal's processed data will be stored -- example '/data99/student/Mil'
%FILEPREFIX -- also the first three letters of the animal's name (example 'con'), which will attach to
%              the beginning of the .mat files containing the variables.  I recommend using this.
%              If not, type ''.
%DAYNUM -- the day number for the experiment (starting with 1)

TIMESTAMPRATE = 10000;

%[otherArgs] = procOptions(varargin);

times = gettimes(fullfile(daydirect,'times.mat'));

% Read in DIO file
files = dir(fullfile(daydirect,'*.dio'));
if isempty(files)
   error('No dio file found in %s.',daydirect);
end
if length(files) > 1
   error('Multipled dio files found in %s.',daydirect);
end
fid = fopen(fullfile(daydirect,files(1).name),'r');
dio = textscan(fid,'%d %s');
fclose(fid);
diotimes = dio{1};
diovals = double(str2mat(dio{2})) - str2mat('0');

% Flip DIO bits since every 1-16 is actually 15-0

for step=[0 16 32 48];
    i=1:16;
    diovals(:,i+step)=diovals(:,fliplr(i+step));
end
    
    
active_pins = find(sum(diovals,1) > 0);
%active_pins = [32, 48];
for a = 1:length(active_pins)
   pin = active_pins(a);
   diovector = diovals(:,pin);

   ind = find(diovector,1);
 
   pulse = [];
   if ~isempty(ind)
      % Note - I've noticed both rising and falling edges missed;
      % For simplicity I assume that only falling edges are ever
      % neglected.
      i = 1;
      pulse(i).tstart = diotimes(ind);
      done = 0;
      while ~done
         found = 0;
         k = ind + 1;
         while (found == 0)
            if (k > length(diovector))
               found = 1;
               pulse(i).tend = Inf;
               pulse(i).length = NaN;
               done = 1;
            elseif (diovector(k) == 0)
               found = 1;
               pulse(i).tend = diotimes(k);
               pulse(i).length = double(pulse(i).tend - pulse(i).tstart);
               ind = k + 1;
               % find next pulse start
               if (ind <= length(diovector))
                  while ((diovector(ind) == 0) & (ind < length(diovector)))
                     ind = ind + 1;
                  end
                  if (diovector(ind) == 1)
                     i = i + 1;
                     pulse(i).tstart = diotimes(ind);
                  else
                     done = 1;
                  end
               else
                  done = 1;
               end
            elseif (sum(diovals(k,:) ~= diovals(ind,:)) > 0)
               % Did a different pin change, causing DIO to
               % record?
               k = k + 1;
            else  % This accounts for missed falling edge
               % and only occurs if the same DIO string has
               % occured twice with different time values.
               found = 1;
               pulse(i).tend = diotimes(k);
               pulse(i).length = NaN;
               i = i + 1;
               pulse(i).tstart = diotimes(k);
               ind = k;
            end
         end
      end
   end
   for i = 1:length(pulse)
      pulse(i).ind = i;
   end
   diopulses{pin}.pulsetimes = double([cat(1,pulse.tstart) cat(1,pulse.tend)]);
   diopulses{pin}.timesincelast = [inf;diff(diopulses{pin}.pulsetimes(:,1))];
   if length(diopulses{pin}.timesincelast) > 1
      diopulses{pin}.timeuntilnext = [diopulses{pin}.timesincelast(2:end); inf];
   else
      diopulses{pin}.timeuntilnext = inf;
   end
   diopulses{pin}.pulselength = cat(1,pulse.length);
   diopulses{pin}.pulseind = cat(1,pulse.ind);
end

 
for epoch = 1:length(times)
    if ~isempty(times(epoch).starttime)
        starttime = timetrans({times(epoch).starttime},TIMESTAMPRATE,2);
        endtime = timetrans({times(epoch).endtime},TIMESTAMPRATE,2);

        % read in the raw diode information (front and back positions)
        % this looks through all the .p files and picks out the times that are
        % within the start and end times
        inEpochTimes = find(diotimes >= starttime & diotimes <= endtime);
        rawdio{daynum}{epoch}.times = diotimes(inEpochTimes);
        rawdio{daynum}{epoch}.values = diovals(inEpochTimes);

        for a = active_pins
           j = 0;
           inEpochTimes = find( ...
              (diopulses{a}.pulsetimes(:,1) >= starttime) & ...
              (diopulses{a}.pulsetimes(:,2) <= endtime) );
           DIO{daynum}{epoch}{a}.pulsetimes = diopulses{a}.pulsetimes(inEpochTimes,:);
           DIO{daynum}{epoch}{a}.timesincelast = diopulses{a}.timesincelast(inEpochTimes,:);
           DIO{daynum}{epoch}{a}.pulselength = diopulses{a}.pulselength(inEpochTimes,:);
           DIO{daynum}{epoch}{a}.pulseind = diopulses{a}.pulseind(inEpochTimes,:);
        end
    end
end 

save(fullfile(animdirect,sprintf('%sDIO%02d.mat',fileprefix,daynum)),'rawdio','DIO','diopulses');
% save(fullfile(animdirect,sprintf('%sDIO%02d.mat',fileprefix,daynum)),'rawdio','DIO','diopulses');
 
