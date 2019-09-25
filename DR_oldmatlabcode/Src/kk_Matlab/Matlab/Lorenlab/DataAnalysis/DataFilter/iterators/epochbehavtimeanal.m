function f = epochbehavtimeanal(f)
% function f = epochbehavtimeanal(f)
%
% Iterator for a filter object.  Calls the function designated in
% f().function.name, after loading the variables designated as strings in
% f().function.loadvariables{:}.  Also the function call appends any
% options in the f().function.options{} cell array.  
% 
% Each function call is for one epoch, and it is assumed that the
% function's first input is a list of indices to the epoch
% ([day epoch]). The second input is the center time for
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
   totalepochs = [];
   for g = 1:length(f(an).epochs)
      totalepochs = [totalepochs; f(an).epochs{g}];
   end
   totaldays = unique(totalepochs(:,1)); %get all of the days across groups
      
   %load all the variables that the function requires
   loadstring = [];
   for i = 1:length(f(an).function.loadvariables)
     if sum(strcmp(f(an).function.loadvariables{i}, {'eeg','theta','ripple','gamma','delta','hp','tdr'}))
       error('Use eeg window analysis iterator!');
     else
       eval([f(an).function.loadvariables{i},' = loaddatastruct(animaldir, animalprefix, f(an).function.loadvariables{i}, totaldays);']);
     end
     loadstring = [loadstring, f(an).function.loadvariables{i},','];
   end

   foptions = f(an).function.options;
   
   %iterate through the epochs within each data group
   for g = 1:length(f(an).epochs)
      E = size(f(an).epochs{g},1);
      if (E == 0)
        f(an).output{g} = [];
      end

      for e = 1:E
         excludeperiods = f(an).excludetime{g}{e};
         indices = f(an).epochs{g}(e,:);

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
                f(an).output{g}(e,fi) = fout;
             elseif isnumeric(fout)        
                if (length(f(an).output) < g)
                   f(an).output{g} = [];
                end
                if (length(t_{e}) > 1) % multiple times filters
                    if (length(f(an).output{g}) < fi)
                       f(an).output{g}{fi} = [];
                    end
                    f(an).output{g}{fi} = stack(f(an).output{g}{fi}, fout);          
                else
                    f(an).output{g} = stack(f(an).output{g}, fout);          
                end
             else
                error(['In calling ', f(an).function.name, ': Function output must be either numeric or a structure']);
             end
         end

      end
         
   end
end
