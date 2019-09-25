function [out, names] = peak_align_pca()
% function [out, names] = peak_align()

global clustattrib;

filename = strcat(clustattrib.datafile,'.mat');
load(filename);

OVERSAMPLING = 4;
SNIPPET_LENGTH = 32;
SNIPPET_PRE = 10;
SNIPPET_SPACE = [16,64];

L = diff(SNIPPET_SPACE) + 1;
N = size(waves,1);
K = size(waves,3);

new_waves = nan(SNIPPET_LENGTH,4,K);

wU = nan(N*OVERSAMPLING, size(waves,2));

for i = 1:K
  w = waves(:,:,i);
  for j = 1:4
    fw = fft(double(w(:,j)));
    fo = zeros(N*OVERSAMPLING,1);
    fo(1:(N/2)) = fw(1:(N/2));
    fo((end-(N/2-2)):end) = fw(((N/2)+2):end);
    wU(:,j) = ifft(fo);
  end
  sw = sum(wU(SNIPPET_SPACE(1):SNIPPET_SPACE(2),:),2);
  sw = sw - min(sw);
  cm = round([1:L]*sw / sum(sw));
  ind_start = max(cm - SNIPPET_PRE*OVERSAMPLING + SNIPPET_SPACE(1),0);
  ind_start = min(ind_start, N*OVERSAMPLING - SNIPPET_LENGTH - 1);
  inds = ind_start + [1:OVERSAMPLING:SNIPPET_LENGTH*OVERSAMPLING];
  new_waves(:,:,i) = wU(inds,:);
end

ww = reshape(new_waves,[SNIPPET_LENGTH*4,K]);
[v,d] = eig(ww*ww');
out = zeros(K,9);
for i = 1:9
  out(:,i) = ww'*v(:,end - (i-1));
  names{i} = sprintf('p%d',i);
end
