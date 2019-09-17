function f = kk_daysinglecellanal(f)
% f = kk_daysinglecellanal(f)

% This function iterates over unique day-tet-cell instead of doing
% epoch-based iteration. So with this iterator, the analysis function receives the

%           i.   the day-tet-cell index of the unit
%           ii.  a vector of epochs the unit was clustered in
%           iii. excludeperiods{epochnum}   i.e. a cell array of excludeperiods by active epoch

% This was written to handle STA analysis, which demanded taking the mean
% and SD of the STA traces over a day and only saving that instead of saving all the STA traces. 
    % dfskk_sta   +   dfakk_getsta

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
        if (~iseegvar(f(an).function.loadvariables{i}))
        eval([f(an).function.loadvariables{i},' = loaddatastruct(animaldir, animalprefix, f(an).function.loadvariables{i}, totaldays);']);
        end
        loadstring = [loadstring, f(an).function.loadvariables{i},','];
    end
    foptions = f(an).function.options;

    % iterate through epoch groupings
    for g = 1:length(f(an).epochs)            
        
        % initialize output
        f(an).output{g} = [];

        % within this epoch grouping, ***** determine all the unique dtc for this animal *****
        % note this is downstream of cellfilter, and thus correctly gets cellfilter target cells
        detc = [];
        dtc = [];
        for ee = 1:size(f(an).epochs{g},1)
            day = f(an).epochs{g}(ee,1) ;
            ep = f(an).epochs{g}(ee,2) ;
            tc_list = f(an).data{g}{ee};            % tet-cell list
                numcells = size(tc_list,1);
            newdetc = [ day*ones(numcells,1)   ep*ones(numcells,1)    tc_list ] ;   % concatenate
            if ~isempty(newdetc)
                detc = [detc ; newdetc];
            end
        end
        dtc = unique(detc(:,[1 3 4]),'rows') ;
        
        % iterate through the dtc
        for cc = 1:size(dtc,1)
        
            epochs = [];                % epochs in which the cell was clustered
            excludeperiods = {};        % store excludeperiods list for EACH active epoch

            % identify row indices of epochs for the day of this cell
            epochrows = find( dtc(cc,1) == f(an).epochs{g}(:,1) )'; 
            
            % iterate through epochs 
            for r = epochrows
                
                epoch_tc = f(an).data{g}{r};   % tet-cell list for this epoch
                epochnum = f(an).epochs{g}(r,2);
                
                % check to see if this cell is in (i.e. was clustered) this epoch -- if not, skip
                if rowfind(dtc(cc,[2 3]) , epoch_tc);
                    epochs = [ epochs   epochnum ];  
                    excludeperiods{epochnum} = f(an).excludetime{g}{r};
                else
                    continue
                end
                
            end
            
            % call analysis function
            eval(['fout = ',f(an).function.name,'( dtc(cc,:), epochs, excludeperiods,', loadstring, 'foptions{:});']);
            
            % concatenate the output
            if isstruct(fout)
                f(an).output{g} = [f(an).output{g}   fout] ;
            else
                error('expecting a struct')
            end
            
        end
    end
    
end
    
end