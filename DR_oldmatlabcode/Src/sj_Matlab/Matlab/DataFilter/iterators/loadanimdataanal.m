function f = loadanimdataanal(f)
% f = epochbehaveanal(f)
% Iterator for a filter object.  Calls the function designated in
% f().function.name, after loading the variables designated as strings in
% f().function.loadvariables{:}.  Also the function call appends any
% options in the f().function.options{} cell array.  
% 
% Each function call if for an animal. 
% Loads a file in the "direct" folder which has data processed across all days and epochs
% Written for DFSsj_getstimstats_summ
% Shantanu, Oct 2011

% The output of the call function can either be a 1 by N vector, or a structure.
% The outputs are stored in f().output, grouped using the same groupings as
% in the filter.

% Incomplete - use epochbehaveanal instead


%iterate through all animals
for an = 1:length(f)
    %get directory, etc
    animaldir = f(an).animal{2};
    animalprefix = f(an).animal{3};
    
    eval(['fout = ',f(an).function.name,'(animaldir,animalprefix,', loadstring, 'foptions{:});']);
    
    % Function will load animaldir/animalprefix[xxxx], where xxxx is the string in foptions that defines filename 
    
    if isstruct(fout)
                if (isempty(f(an).output) | (length(f(an).output) < g))
                    f(an).output{g}(1) = fout;
                else
                    f(an).output{g}(end+1) = fout;
                end
            elseif isnumeric(fout)
                if ((isempty(f(an).output)) | (length(f(an).output) < g))
                    f(an).output{g} = [];
                end
                if (size(fout,1) == 1)
                    f(an).output{g} = [f(an).output{g}; fout];
                else
                    error(['In calling ', f(an).function.name, ': Numeric function outputs must be 1 by N.  Use a structure output for more complicated outputs']);
                end
            else
                error(['In calling ', f(an).function.name, ': Function output must be either numeric or a structure']);
    end
            
end
    
    
    
    %iterate through the epochs within each data group
    for g = 1:length(f(an).epochs)
        
        for e = 1:size(f(an).epochs{g},1)
            excludeperiods = f(an).excludetime{g}{e};
            tmpindex = [f(an).epochs{g}(e,:)];
            eval(['fout = ',f(an).function.name,'(tmpindex,excludeperiods,', loadstring, 'foptions{:});']);
            if isstruct(fout)
                if (isempty(f(an).output) | (length(f(an).output) < g))
                    f(an).output{g}(1) = fout;
                else
                    f(an).output{g}(end+1) = fout;
                end
            elseif isnumeric(fout)
                if ((isempty(f(an).output)) | (length(f(an).output) < g))
                    f(an).output{g} = [];
                end
                if (size(fout,1) == 1)
                    f(an).output{g} = [f(an).output{g}; fout];
                else
                    error(['In calling ', f(an).function.name, ': Numeric function outputs must be 1 by N.  Use a structure output for more complicated outputs']);
                end
            else
                error(['In calling ', f(an).function.name, ': Function output must be either numeric or a structure']);
            end
            
            
        end
    end
end