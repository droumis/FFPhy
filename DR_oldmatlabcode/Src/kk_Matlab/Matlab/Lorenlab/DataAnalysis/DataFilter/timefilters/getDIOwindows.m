function out = getDIOwindows(animaldir,animalprefix,epochs,window,varargin)
% function out = getDIOwindows(animaldir,animalprefix,epochs,window,varargin)
%   EPOCHS - N by 2 matrix, columns are [day epoch]
%   PIN - which DIO pin to look at
%   WINDOW - [2x1] vector correspond to window around pulse start
%      time to <include>
%   FILTER - (optional) filter on fields in DIO struct

pin = 48;
filter = '';
[otherArgs] = procOptions(varargin);

loaddays = unique(epochs(:,1));
DIO = loaddatastruct(animaldir, animalprefix, 'DIO', loaddays);


if ~isempty(filter)
   includeDIO = evaluatefilter2(DIO,filter,'activeField','pulsetimes');
else
   includeDIO = DIO;
end


E = size(epochs,1);

for i = 1:E
   out{epochs(i,1)}{epochs(i,2)} = [];
   if isempty(DIO)
     continue;
   end
   if ~isempty(includeDIO{epochs(i,1)})
      if ~isempty(includeDIO{epochs(i,1)}{epochs(i,2)})
         pp = includeDIO{epochs(i,1)}{epochs(i,2)}{pin}.pulsetimes;
         out{epochs(i,1)}{epochs(i,2)} = repmat(window,size(pp,1),1) + pp/10000;
      end
   end
end
