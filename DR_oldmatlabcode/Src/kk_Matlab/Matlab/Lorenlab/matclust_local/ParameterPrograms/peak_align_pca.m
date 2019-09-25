function [out, names] = peak_align_pca()
% function [out, names] = peak_align()

global clustattrib;

filename = fullfile(clustattrib.currentfilepath,strcat(clustattrib.datafile,'.mat'));
load(filename);
OVERSAMPLING = 4;
SNIPPET_LENGTH = 32;
SNIPPET_PRE = 10;
SNIPPET_SPACE = [16,64];

nTetrodes = 4;

nComponents = 4;

L = diff(SNIPPET_SPACE) + 1;
N = size(waves,1);

new_waves = nan(SNIPPET_LENGTH,4,size(waves,3));

wU = nan(N*OVERSAMPLING, size(waves,2));

for i = 1:size(waves,3)
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

j = 0;
for i = 1:nTetrodes
  ww = squeeze(new_waves(:,i,:));
  [v,d] = eig(ww*ww');
  out(:,j + [1:nComponents]) = ww'*v(:,end - [0:nComponents-1]);
  for k = 1:nComponents
    names{j + k} = sprintf('p%d_%d',i,k);
   % names(j + [1:nComponents]) = {sprintf('p%d_1',i), ...
     % sprintf('p%d_2',i), ...
     % sprintf('p%d_3',i)};
  end
  j = j + nComponents;
end

