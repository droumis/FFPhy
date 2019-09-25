function [out] = getDIOtheta(index, excludeperiods, theta, DIO, delta, varargin)
% function [out] = getDIOtheta(index, excludeperiods, theta, DIO, delta, varargin)
%
%   index [day epoc tetrode]
%
%   out is [durations]

% assign the options

window = [-0.25 0.5];
pin = 48;
TIMESTAMPRATE = 10000;
filter = '';
[otherArgs] = procOptions(varargin);

for t = 1:size(index,1) % number of tetrodes
   e = theta{index(t,1)}{index(t,2)}{index(t,3)};
   th_tmp = e.data;
   if size(th_tmp,2) == 3
      h_theta = double(th_tmp(:,3)) .* (cos(double(th_tmp(:,2))/10000) + sqrt(-1)*sin(double(th_tmp(:,2))/10000));
   else
      h_theta = double(th_tmp(:,2)) .* (cos(double(th_tmp(:,1))) + sqrt(-1)*sin(double(th_tmp(:,1))));
   end
   pulses = DIO{index(t,1)}{index(t,2)}{pin};

   if ~isempty(pulses.pulsetimes)
      if ~isempty(filter)
         filtresult = evaluatefilter2({pulses}, filter, 'activefield','pulseind');
         pulsetimes = pulses.pulsetimes(find(filtresult{1}(:,2))) / TIMESTAMPRATE;
      else
         pulsetimes = pulses.pulsetimes(:,1) / TIMESTAMPRATE;
      end

      goodPulses = find(~isExcluded(pulsetimes, excludeperiods));

      etimes = e.starttime + [0:length(h_theta)-1]*(1/e.samprate);

      startinds = lookup(pulsetimes(goodPulses),etimes);

      W = window * round(e.samprate);  % convert window from seconds to theta samples
      W = round(W(1)):round(W(2));
      out.data(:,:,t) = nan(length(goodPulses),length(W),1);
      for i = 1:length(goodPulses)
         win = startinds(i) + W;
         out.inds(i,:) = [win(1) win(end)];
         out.times(i,:) = out.inds(i,:) / e.samprate;
         if (win(1) >= 1) && (win(end) <= length(h_theta))
            out.data(i,:,t) = h_theta(win);
         elseif (win(end) > length(h_theta))
            k = find(win > length(h_theta),1);
            win = win(1:k-1);
            out.data(i,1:k-1,t) = h_theta(win);
         else
            k = find(win >= 1,1);
            win = win(k:end);
            out.data(i,k:end,t) = h_theta(win);
         end
         if ~isempty(delta)
            if size(delta{index(t,1)}{index(t,2)}{index(t,3)}.data,2) == 2
               d = mean(delta{index(t,1)}{index(t,2)}{index(t,3)}.data(win,2));
            else
               d = mean(delta{index(t,1)}{index(t,2)}{index(t,3)}.data(win,3));
            end
            out.tdr(i,t) = nanmean(abs(out.data(i,:))) / d;
         end
      end
   % else
      % out(:,:,t) = [];
   end
end

if ~isempty(out)
   out.pulsetimes = pulsetimes(goodPulses);
end
