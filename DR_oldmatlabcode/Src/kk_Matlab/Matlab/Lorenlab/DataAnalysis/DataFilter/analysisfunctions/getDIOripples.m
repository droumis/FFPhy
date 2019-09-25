function [out] = getDIOtheta(index, excludeperiods, ripple, ripples, DIO, varargin)
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
   e = ripple{index(t,1)}{index(t,2)}{index(t,3)};
   th_tmp = e.data;
   if size(th_tmp,2) == 2
      h_ripple = double(th_tmp(:,2)) .* (cos(double(th_tmp(:,1))) + sqrt(-1)*sin(double(th_tmp(:,1))));
   else
      h_ripple = double(th_tmp(:,3)) .* (cos(double(th_tmp(:,2)/10000)) + sqrt(-1)*sin(double(th_tmp(:,2)/10000)));
   end
   pulses = DIO{index(t,1)}{index(t,2)}{pin};

   rips = ripples{index(t,1)}{index(t,2)}{index(t,3)};
   isrip = logical(zeros(length(h_ripple),1));
   for i = 1:length(rips.startind)
      isrip(rips.startind(i):rips.endind(i)) = 1;
   end

   if ~isempty(pulses.pulsetimes)
      if ~isempty(filter)
         filtresult = evaluatefilter2({pulses}, filter, 'activefield','pulseind');
         pulsetimes = pulses.pulsetimes(find(filtresult{1}(:,2))) / TIMESTAMPRATE;
      else
         pulsetimes = pulses.pulsetimes(:,1) / TIMESTAMPRATE;
      end

      goodPulses = find(~isExcluded(pulsetimes, excludeperiods));

      etimes = e.starttime + [0:length(h_ripple)-1]*(1/e.samprate);

      startinds = lookup(pulsetimes(goodPulses),etimes);

      W = window * round(e.samprate);  % convert window from seconds to ripple samples
      W = round(W(1)):round(W(2));
      out.ripeeg(:,:,t) = nan(length(goodPulses),length(W),1);
      for i = 1:length(goodPulses)
         win = startinds(i) + W;
         out.inds(i,:) = [win(1) win(end)];
         out.times(i,:) = out.inds(i,:) / e.samprate;
         out.isrip(i,:,t) = logical(zeros(1,length(W),1));
         if (win(1) >= 1) && (win(end) <= length(h_ripple))
            out.ripeeg(i,:,t) = h_ripple(win);
            out.isrip(i,:,t) = isrip(win);
         elseif (win(end) > length(h_ripple))
            k = find(win > length(h_ripple),1);
            win = win(1:k-1);
            out.ripeeg(i,1:k-1,t) = h_ripple(win);
            out.ripeeg(i,k+length(W),t) = nan;
            out.isrip(i,1:k-1,t) = isrip(win);
         else
            k = find(win >= 1,1);
            win = win(k:end);
            out.ripeeg(i,1:k-1,t) = nan;
            out.ripeeg(i,k:end,t) = h_ripple(win);
            out.isrip(i,k:end,t) = isrip(win);
         end

      end
   % else
      % out(:,:,t) = [];
   end
end

if ~isempty(out)
   out.pulsetimes = pulsetimes(goodPulses);
end
