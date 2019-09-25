
subjects = {'S48'};

for i = 1:numel(subjects)
  subject = subjects{i};

  % load session info for this subject
  load(sprintf('/data14/smkim/%s/%s_session.mat',subject,subject));

  tetfolders = dir(sprintf('/data14/smkim/%s/spike/tetrode*',subject));
  for j = 1:length(tetfolders)
    tetfolder = tetfolders(j).name;
    matclust_files = dir(sprintf('/data14/smkim/%s/spike/%s/*_matclust.mat',subject,tetfolder));
    for k = 1:length(matclust_files)
      matclust_file = sprintf('/data14/smkim/%s/spike/%s/%s',subject,tetfolder, ...
          matclust_files(k).name);

      try
        oldclust = load(matclust_file);
      catch
        error('could not load file %s',matclust_file);
        continue;
      end
    
      try
        day = str2num(regexp(matclust_file,'(?<=_day)\d(?=_)','match','once'));
        assert(isscalar(day) && isreal(day) && (round(day) == day));
      catch
        error('could not parse day from filename %s',matclust_file);
      end
      
      % fix only day 1!
      if (day ~= 1)
        continue;
      end
      
      disp(matclust_file);
      try
        newclust = reset_matclust_struct(oldclust,session([session(:).day] == day));
      catch
        error('RESET_MATCLUST_STRUCT failed on %s',matclust_file);
      end

      newclust = oldclust;
      
      % overwrite the old file. make sure that this is what you really want to do!!!!
      save(matclust_file,'-struct','newclust');

      clear('matclust_file','oldclust','newclust');

    end
  end

  clear('session');
  
end

% reclustered matclust files for S48/spike/:
% tetrode02/S48_day1_tetrode02(070)_matclust.mat
% tetrode05/S48_day1_tetrode05(063)_matclust.mat
% tetrode06/S48_day1_tetrode06(059)_matclust.mat
% tetrode09/S48_day1_tetrode09(091)_matclust.mat
% tetrode11/S48_day1_tetrode11(073)_matclust.mat
% tetrode12/S48_day1_tetrode12(069)_matclust.mat
% tetrode13/S48_day1_tetrode13(095)_matclust.mat


