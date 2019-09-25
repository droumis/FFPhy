function [pp, center_of_mass] = peak_align(waves)
% function pp = peak_align(waves)

OVERSAMPLING = 4;
SNIPPET_LENGTH = 28;
SNIPPET_PRE = 10;

THRESHOLD_IND = 9;
WINDOW = [5,15];
CM_THRESH = 20;


N = size(waves,1);
K = N;
L = length(((WINDOW(1))*OVERSAMPLING) : (WINDOW(2)*OVERSAMPLING));

% generate oversampling matrix
delta = zeros(K,1);
for i = 1:K
  dd = delta;
  dd(i) = 1;
  fw = fft(dd);
  fo = zeros(K*OVERSAMPLING,1);
  fo(1:(K/2)) = fw(1:(K/2));
  fo((end-(K/2-2)):end) = fw(((K/2)+2):end);
  OM(:,i) = OVERSAMPLING*ifft(fo);
end


new_waves = nan(SNIPPET_LENGTH,4,size(waves,3));
% center_of_mass = nan(1,size(waves,3));

wU = nan(N*OVERSAMPLING, size(waves,2));

threshold = min(max(waves(THRESHOLD_IND,:,:),[],2)); % this tells us the approximate threshold

for i = 1:size(waves,3)
  w = waves(:,:,i);
  wU = OM*double(w);

  % look for maximum in post-threshold window across all 4 channels
  sw = wU(((WINDOW(1))*OVERSAMPLING) : (WINDOW(2)*OVERSAMPLING), :);
  sw = sw - CM_THRESH;
  sw = sw .* (sw>0);

  % cm = [1:L]*sw ./ sum(sw); % by component center of mass
  cm = [1:L]*sum(sw,2) / sum(sum(sw)); % sum all components center of mass
  center_of_mass(:,i) = cm;

  cm = round(cm);

  ind_start = max(cm - SNIPPET_PRE*OVERSAMPLING + WINDOW(1),0);
  ind_start = min(ind_start, N*OVERSAMPLING - SNIPPET_LENGTH*OVERSAMPLING - 1);
  inds = ind_start + [1:OVERSAMPLING:SNIPPET_LENGTH*OVERSAMPLING];
  new_waves(:,:,i) = wU(inds,:);
end

pp = new_waves;

