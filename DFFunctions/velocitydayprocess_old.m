function velocitydayprocess(directoryname,fileprefix,days,varargin)
% function velocitydayprocess(directoryname,fileprefix,days, varargin)
%
%directoryname - example '/data99/user/animaldatafolder/', a folder 
%            containing processed matlab data for the animal
%
%fileprefix -   animal specific prefix for each datafile (e.g. 'fre')
%
%days -         a vector of experiment day numbers 
%

defaultfilter = 'velocitydayprocess_filter.mat';
quiescence_threshold = 2; % cm/s

[otherArgs] = procOptions(varargin);

eval(['load ', defaultfilter]);

L = length(velocityfilter.kernel);

days = days(:)'; % make days a 1xN

isLinPos = 0;

fprintf('%s:',fileprefix);

for day = days
   fname = fullfile(directoryname,sprintf('%spos%02d.mat',fileprefix,day));
   linpname = fullfile(directoryname,sprintf('%slinpos%02d.mat',fileprefix,day));
   if exist(fname) ~= 2
      warning(sprintf('[velocitydayprocess] No files found for day %d!',day));
      continue;
   end
   load(fname);
   fprintf('(%d)',day);

   if exist(linpname) == 2
     % load(linpname);
     % isLinPos = 1;
   end

   nepochs = length(pos{day});
   for e = 1:nepochs
     if ~isfield(pos{day}{e},'data') || isempty(pos{day}{e}.data)
       continue;
     end
     if length(pos{day}{e}.data) > 3*L
       % smooth epoch velocity
       vv = pos{day}{e}.data(:,5);
       vv = [flipud(vv(1:L,1)); vv; flipud(vv(end-L+1:end,1))];
       v = filtfilt(velocityfilter.kernel,1,vv);
       v = v(L+1:end-L);
       if length(v) ~= length(pos{day}{e}.data)
         error('Bad lengths.');
       end

       if isLinPos && ~isempty(linpos{day}{e})
         lvv = linpos{day}{e}.statematrix.linearVelocity(:,1);
         lvv = [flipud(lvv(1:L,1)); lvv; flipud(lvv(end-L+1:end,1))];
         lv = filtfilt(velocityfilter.kernel,1,lvv);
         lv = lv(L+1:end-L);
       end
     else
       v = pos{day}{e}.data(:,5);
     end
     t = pos{day}{e}.data(:,1);
     % calculate times since threshold crossing
     threshinds = v < quiescence_threshold;
     quiesc_time = zeros(length(v),1);
     k = find(threshinds,1);
     m = 0;
     while ~isempty(k)
       k = k + m;
       m = find(threshinds(k+1:end) == 0,1);
       if isempty(m)
         m = length(threshinds);
       else
         m = m + k;
       end
       quiesc_time(k:m) = t(k:m) - t(k);
       k = find(threshinds(m+1:end),1);
     end

     if isLinPos == 0 || isempty(linpos{day}{e})
       if (size(pos{day}{e}.data,2) == 5) % old format
         pos{day}{e}.data = [pos{day}{e}.data(:,1:5) nan(length(pos{day}{e}.data),5) v quiesc_time];
         pos{day}{e}.descript{11} = sprintf('Smoothed velocity and quiescence time using smoothed velocity, %f cm/s threshold',quiescence_threshold);
         pos{day}{e}.fields = 'time x y dir vel [missing-steve-stuff] smooth-v q-time';
       elseif (size(pos{day}{e}.data,2) == 7) % old format smoothed
         fprintf('[Converting]');
         pos{day}{e}.data = [pos{day}{e}.data(:,1:5) nan(length(pos{day}{e}.data),5) v quiesc_time];
         pos{day}{e}.descript{11} = sprintf('Smoothed velocity and quiescence time using smoothed velocity, %f cm/s threshold',quiescence_threshold);
         pos{day}{e}.fields = 'time x y dir vel [missing-steve-stuff] smooth-v q-time';
       elseif (size(pos{day}{e}.data,2) == 10) % new (steve) format
         pos{day}{e}.data = [pos{day}{e}.data(:,1:10) v quiesc_time];
         pos{day}{e}.descript{11} = sprintf('Smoothed velocity and quiescence time using smoothed velocity, %f cm/s threshold',quiescence_threshold);
         pos{day}{e}.fields = 'time x y dir vel steve-stuff smooth-v q-time';
       else
         error('Unknown pos data format');
       end
     else
       keyboard
       pos{day}{e}.data = [pos{day}{e}.data(:,1:10) v quiesc_time lv];
       pos{day}{e}.descript{11} = sprintf('Smoothed velocity and quiescence time using smoothed velocity, %f cm/s threshold',quiescence_threshold);
       pos{day}{e}.descript{12} = sprintf('Smoothed linearized velocity.');
       pos{day}{e}.fields = 'time x y dir vel steve-stuff smooth-v q-time lin-v';
     end
     fprintf('.');
   end

   % save the resulting file
   save(fname, 'pos');
   clear pos
end
fprintf('\n');
