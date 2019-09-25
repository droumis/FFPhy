function f = singletetrodewindowanal(f)
% function f = singletetrodewindowanal(f)
% Iterator for a filter object.  Calls the function designated in
% f().function.name, after loading the variables designated as strings in
% f().function.loadvariables{:}.  Also the function call appends any
% options in the f().function.options{} cell array.  
% 
% Each function call is for one epoch, and it is assumed that the
% function's first input is a list of indices to the tetrode
% ([day epoch tetrode]). The second input is the center time for
% a time window analsysis. The third input is a list of exclusion
% periods [starttime endtime].  The next inputs are the load
% variables, and the final inputs are the options.  out =
% fname(index, times, excludeperiods, var1, var2, ..., option1,
% option2,...).
%
% The output of the call function can either be numeric, or a structure.
% If the output if numeric, data is appended along the first dimenion for
% the group.

%iterate through all animals
for an = 1:length(f)
   %find all unique epochs to analyze for the current animal
   animaldir = f(an).animal{2};
   animalprefix = f(an).animal{3};

   foptions = f(an).function.options;
   
   %iterate through the epochs within each data group
   for g = 1:length(f(an).epochs)
      E = size(f(an).epochs{g},1);
      if (E == 0) | isempty(f(an).data{g})
        f(an).output{g} = [];
        continue;
      end

      allindices = unique(cat(1,f(an).data{g}{:}),'rows');
      
      for e = 1:E % by epoch
         indices = f(an).data{g}{e};
         numindices = size(indices,1);
         indices = [repmat(f(an).epochs{g}(e,:),[numindices 1]) indices];
         excludeperiods = f(an).excludetime{g}{e};

         %load all the variables that the function requires
         loadstring = [];
         for i = 1:length(f(an).function.loadvariables)
            if sum(strcmp(f(an).function.loadvariables{i}, {'eeg','theta','ripple','wideripple','gamma','delta','hpdata','tdr'}))
               eval([f(an).function.loadvariables{i},' = loadeegstruct(animaldir, animalprefix, f(an).function.loadvariables{i}, f(an).epochs{g}(e,1), f(an).epochs{g}(e,2), f(an).data{g}{e});']);
            elseif strcmp(f(an).function.loadvariables{i}, 'eegfiledir')
              if length(unique(f(an).epochs{g}(:,1))) > 1
                error('One day at a time only for raw eegs.');
              end
              daynum = unique(f(an).epochs{g}(:,1));
              eegfiledir = fullfile(animaldirdata(f(an).animal{1},daynum));
            else
               eval([f(an).function.loadvariables{i},' = loaddatastruct(animaldir, animalprefix, f(an).function.loadvariables{i}, f(an).epochs{g}(e,1));']);
            end
            loadstring = [loadstring, f(an).function.loadvariables{i},','];
         end

         t_ = f(an).times{g};
         if isnumeric(t_)
           t_ = {{t_}};
         end
         if isnumeric(t_{e})
           t_{e} = {t_{e}};
         end
         for fi = 1:length(t_{e})
             t_func = t_{e}{fi};
             eval(['fout = ',f(an).function.name,'(indices,t_func,excludeperiods,', loadstring, 'foptions{:});']);

             %save the function output in the filter variable.  Allows numeric or struct outputs
             if isstruct(fout)
                if length(fout) == 1
                  f(an).output{g}{fi}(e) = fout;
                elseif length(fout) == size(f(an).data{g}{e},1)
                  if e > 1
                    fnames = fieldnames(fout);
                    for c = 1:length(fout)
                      for fn = 1:length(fnames)
                        f(an).output{g}{fi}(c).(fnames{fn}) = [f(an).output{g}{fi}(c).(fnames{fn}), fout(c).(fnames{fn})];
                      end
                    end
                  else
                    f(an).output{g}{fi} = fout;
                  end
                else
                  error('In timewindowanal, analysis function must return a structure of the proper size!');
                end
             elseif isnumeric(fout)
                % fout will either be size length(t_func) x length(indices) or
                %   length(t_func) x W x length(indices)
                % For simplicity, though, all of the data from an
                % animal should have the same length(indices). So
                % we'll create a temp variable of the maximum
                % size and feed fout into it.

                if size(fout,1) == length(t_func) & size(fout,2) ~= size(allindices,1)
                  fout2 = nan(size(fout,1),size(fout,2),size(allindices,1));
                  if ~isempty(fout)
                    fout2(:,:,ismember(allindices,f(an).data{g}{e},'rows')) = fout;
                  end
                else
                  fout2 = nan(size(fout,1),size(allindices,1));
                  if ~isempty(fout)
                    fout2(:,ismember(allindices,f(an).data{g}{e},'rows')) = fout;
                  end
                end

                clear fout;
                if (length(f(an).output) < g)
                   f(an).output{g} = [];
                end
                if (length(t_{e}) > 1) % multiple times filters
                    if (length(f(an).output{g}) < fi)
                       f(an).output{g}{fi} = [];
                    end
                    f(an).output{g}{fi} = stack(f(an).output{g}{fi}, fout2);
                else
                    f(an).output{g} = stack(f(an).output{g}, fout2);
                end

             else
                error(['In calling ', f(an).function.name, ': Function output must be either numeric or a structure']);
             end
         end

      end
         
   end
end
