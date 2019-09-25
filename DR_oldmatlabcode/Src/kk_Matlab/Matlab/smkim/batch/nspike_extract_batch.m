
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
subjects = {'S61'};
masters = cellstr(repmat('drizzle',[length(subjects) 1])); % master hostname
slaves = cellstr(repmat('rain',[length(subjects) 1])); % slave hostname
executable = '/home/smkim/nspike_extract';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:numel(subjects)
  this_subject = subjects{s};
  this_master = masters{s};
  this_slave = slaves{s};
  if (exist([path_prefix '/' this_subject]) ~= 7);  
    error('could not find folder %s/%s',path_prefix,this_subject);
  end 

  fileglob = sprintf('%s/%s/%s_*.dat',path_prefix,this_subject,this_subject);
  l = dir(fileglob);
  if isempty(l)
    error('no .dat files found that match filename pattern %s',fileglob);
  else
    for f = 1:length(l)
      filename = sprintf('%s/%s/%s',path_prefix,this_subject,l(f).name);
      if (system(sprintf('%s -config %s',executable,filename)) < 0)
        error('could not extract -config file from %s',filename);
      end
      if ~isempty(regexp(filename,sprintf('/%s\\S+\\.%s\\.dat$', ...
          this_subject,this_master)))
        disp(sprintf('extracting data from %s',filename));
        if (system(sprintf('%s -cont %s',executable,filename)) < 0)
          error('could not extract -cont records from %s',filename);
        end
        if (system(sprintf('%s -pos64 %s',executable,filename)) < 0)
          error('could not extract -pos64 records from %s',filename);
        end
        if (system(sprintf('%s -dio %s',executable,filename)) < 0)
          error('could not extract -dio records from %s',filename);
        end
      elseif ~isempty(regexp(filename,sprintf('/%s\\S+\\.%s\\.dat$', ...
          this_subject,this_slave)))
        disp(sprintf('extracting data from %s',filename));
        if (system(sprintf('%s -spike %s',executable,filename)) < 0)
          error('could not extract -spike records from %s',filename);
        end
      else
        error('dat file %s does not match an expected hostnames', ...
            filename);
      end
    end
  end

end
clear;

