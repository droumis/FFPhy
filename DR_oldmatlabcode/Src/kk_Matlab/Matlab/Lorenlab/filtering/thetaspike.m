function thetaspike(directoryname,fileprefix,day, varargin)

spkfn = fullfile(directoryname,sprintf('%sspikes%02d.mat',fileprefix,day));
load(spkfn);

eegdir = fullfile(directoryname,'EEG');

for e = 1:length(spikes{day})
   for t = 1:length(spikes{day}{e})
      if length(spikes{day}{e}{t}) > 0
         eefn = fullfile(eegdir,sprintf('%stheta%02d-%d-%02d.mat',fileprefix,day,e,t));
         clear theta;
         load(eefn);
         Tstart = theta{day}{e}{t}.starttime;
         T_S = theta{day}{e}{t}.samprate;
         for c = 1:length(spikes{day}{e}{t})
            if ~isempty(spikes{day}{e}{t}{c})
               thetadata = nan(size(spikes{day}{e}{t}{c}.data,1),2);
               spk_t = spikes{day}{e}{t}{c}.data(:,1);
               spk_t = round((spk_t - Tstart)*T_S);
               in_range_t = find(spk_t > 0 & spk_t <= size(theta{day}{e}{t}.data,1));
               thetadata(in_range_t,:) = double(theta{day}{e}{t}.data(spk_t(in_range_t),:));
               spikes{day}{e}{t}{c}.data = [spikes{day}{e}{t}{c}.data(:,1:2) thetadata];
            end
         end
      end
   end
end

save(spkfn,'spikes');
