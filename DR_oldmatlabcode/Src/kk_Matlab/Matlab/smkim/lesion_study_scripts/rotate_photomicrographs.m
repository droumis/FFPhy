
comment = '2.43 micrometers/pixel';

%{
cd('/data15/smkim/S48/photomicrographs');
list = dir('*_slide*_section*.tif');
for i = 1:numel(list)
  fname = list(i).name;
  imdata = imread(fname);
  imwrite(permute(imdata,[2 1 3]),[fname(1:end-4) '.png'],'png', ...
      'Comment',comment);
end
  
cd('/data15/smkim/S58/photomicrographs');
list = dir('*_slide*_section*.tif');
for i = 1:numel(list)
  fname = list(i).name;
  imdata = imread(fname);
  imwrite(permute(imdata,[2 1 3]),[fname(1:end-4) '.png'],'png', ...
      'Comment',comment);
end
  
cd('/data15/smkim/S59/photomicrographs');
list = dir('*_slide*_section*.tif');
for i = 1:numel(list)
  fname = list(i).name;
  imdata = imread(fname);
  imwrite(flipdim(flipdim(imdata,1),2),[fname(1:end-4) '.png'],'png', ...
      'Comment',comment);
end
%}

cd('/data15/smkim/S60/photomicrographs');
list = dir('*_slide*_section*.tif');
for i = 1:numel(list)
  fname = list(i).name;
  imdata = imread(fname);
  imwrite(flipdim(imdata,1),[fname(1:end-4) '.png'],'png', ...
      'Comment',comment);
end
  

