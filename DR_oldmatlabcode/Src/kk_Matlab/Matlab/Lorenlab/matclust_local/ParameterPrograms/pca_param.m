function [out, names] = pca_param();

global clustattrib;

filename = strcat(clustattrib.datafile,'.mat');
load(filename);

j = 0;
for i = 1:4
  ww = double(squeeze(waves(:,i,:)));
  [v,d] = eig(ww*ww');
  out(:,j + [1:3]) = ww'*v(:,end - [0:2]);
  names(j + [1:3]) = {sprintf('p%d_1',i), ...
    sprintf('p%d_2',i), ...
    sprintf('p%d_3',i)};
  j = j + 3;
end

