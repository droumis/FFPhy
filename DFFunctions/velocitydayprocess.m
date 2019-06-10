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
[otherArgs] = procOptions(varargin);

eval(['load ', defaultfilter]);

L = length(velocityfilter.kernel);

days = days(:)'; % make days a 1xN

fprintf('%s:',fileprefix);

for day = days
   fname = fullfile(directoryname,sprintf('%spos%02d.mat',fileprefix,day));
   if exist(fname) ~= 2
      warning(sprintf('[velocitydayprocess] No files found for day %d!',day));
      continue;
   end
   load(fname);
   fprintf('(%d)',day);

   nepochs = length(pos{day});
   for e = 1:nepochs
     if ~isfield(pos{day}{e},'data') || isempty(pos{day}{e}.data)
       continue;
     end

     if length(pos{day}{e}.data) > 3*L
       % smooth epoch velocity
       xx = pos{day}{e}.data(:,2:3);
       clear smooth_xx;
       smooth_xx(:,1) = filtfilt(velocityfilter.kernel,1,xx(:,1));
       smooth_xx(:,2) = filtfilt(velocityfilter.kernel,1,xx(:,2));

       v = dist(smooth_xx(1:end-1,:),smooth_xx(2:end,:)) / ...
         mean(diff(pos{day}{e}.data(:,1)));
       v = [v(1);v];

       % vv = pos{day}{e}.data(:,5);
       % vv = [flipud(vv(1:L,1)); vv; flipud(vv(end-L+1:end,1))];
       % v = filtfilt(velocityfilter.kernel,1,vv);
       % v = v(L+1:end-L);
       % if length(v) ~= length(pos{day}{e}.data)
         % error('Bad lengths.');
       % end
     else
       v = pos{day}{e}.data(:,5);
     end

     if strfind(pos{day}{e}.fields,'time x y dir vel xvel yvel xres yres stopped') 
       % new (steve) format
       pos{day}{e}.data = [pos{day}{e}.data(:,1:10) v];
       pos{day}{e}.descript{2} = sprintf('Smoothed velocity (smoothed pos -> velocity)');
       pos{day}{e}.descript{3} = velocityfilter.descript;
       pos{day}{e}.fields = 'time x y dir vel xvel yvel xres yres stopped smooth-v';
       fprintf('.');
     elseif strfind(pos{day}{e}.fields,'time x y dir vel steve-stuff') 
       % new (steve) format
       pos{day}{e}.data = [pos{day}{e}.data(:,1:10) v];
       pos{day}{e}.descript{2} = sprintf('Smoothed velocity (smoothed pos -> velocity)');
       pos{day}{e}.descript{3} = velocityfilter.descript;
       pos{day}{e}.fields = 'time x y dir vel xvel yvel xres yres stopped smooth-v';
       fprintf('.');
     else
       fprintf('[no steve] (%s)\n',pos{day}{e}.fields);
       continue
       % error('Unknown pos data format');
     end
   end

   % save the resulting file
   save(fname, 'pos');
   clear pos
end
fprintf('\n');
