for a = 1:10
    switch a
        case 1
            animalinfo = {'Bond','/data13/mcarr/Bon/','bon'};
        case 2
            animalinfo = {'Conley','/data13/mcarr/Con/','con'};
        case 3
            animalinfo = {'Corriander','/data13/mcarr/Cor/','Cor'};
        case 4
            animalinfo = {'Dudley','/data13/mcarr/Dud/','dud'};
        case 5
            animalinfo = {'Frank','/data13/mcarr/Fra/','fra'};
        case 6
            animalinfo = {'Miles','/data13/mcarr/Mil/','mil'};
        case 7
            animalinfo = {'Ten','/data13/mcarr/Ten/','ten'};
        case 8
            animalinfo = {'Eight','/data13/mcarr/Eig/','Eig'};
        case 9
            animalinfo = {'Five','/data13/mcarr/Fiv/','Fiv'};
        case 10
            animalinfo = {'Seven','/data13/mcarr/Sev/','Sev'};       
    end
    animals{a}.animdirectory = animalinfo{2};
    animals{a}.fileprefix = animalinfo{3};
end

%% ASSUMES THAT RUNCALCNONREFERENCEEEG HAS BEEN RUN ON ALL DATA!!

%% MULTIPLY ALL EEG FILES BY -1
% for a = 1:length(animals)
%     %get days for each animal
%     load([animals{a}.animdirectory, '/', animals{a}.fileprefix, 'cellinfo'])
%     tempdays = 1:length(cellinfo);
% 
%     %excludedays if empty in cellinfo
%     includedays = zeros(size(tempdays));
%     for d = tempdays
%         if ~isempty(cellinfo{tempdays(d)})
%             includedays(d) = 1;
%         end
%     end
%     
%     days{a} = tempdays(logical(includedays));
%     
%     %flip eeg data
%     flipeegdays(animals{a}.animdirectory, animals{a}.fileprefix, days{a},'nonreference',1)
% end
%% PROCESS RIPPLES (NONEEGREFERENCE)

%number of standard deviations for ripple extraction
for a = 1:length(animals)
    %get days for each animal
    load([animals{a}.animdirectory, '/', animals{a}.fileprefix, 'cellinfo'])
    tempdays = 1:length(cellinfo);

    %excludedays if empty in cellinfo
    includedays = zeros(size(tempdays));
    for d = tempdays
        if ~isempty(cellinfo{tempdays(d)})
            includedays(d) = 1;
        end
    end
    
    days{a} = tempdays(logical(includedays));
    
    %ripple analysis
    rippledayprocess(animals{a}.animdirectory, animals{a}.fileprefix, days{a},'tetfilter','isequal($area,''CA1'') || isequal($area,''CA3'')','nonreference',1,'assignphase',1)
end

%% PROCESS RIPPLES (WITH REFERENCE)

%number of standard deviations for ripple extraction
% nstd = 3; 
% 
% for a = 1:length(animals)
%     %get days for each animal
%     load([animals{a}.animdirectory, '/', animals{a}.fileprefix, 'cellinfo'])
%     tempdays = 1:length(cellinfo);
% 
%     %excludedays if empty in cellinfo
%     includedays = zeros(size(tempdays));
%     for d = tempdays
%         if ~isempty(cellinfo{tempdays(d)})
%             includedays(d) = 1;
%         end
%     end
%     
%     days{a} = tempdays(logical(includedays));
%     
%     %ripple analysis
%     rippledayprocess(animals{a}.animdirectory, animals{a}.fileprefix, days{a},'tetfilter','isequal($area,''CA1'') || isequal($area,''CA3'')')
%     for j=1:length(days{a})
%         d = days{a}(j)
%         extractripples(animals{a}.animdirectory, animals{a}.fileprefix, d, -1, 0.015, nstd)%-1 means all tetrodes processes
%     end
% end

%% ADD CORRECT TRAJECTORY TO LINPOS

for a = 1:length(animals)
    load([animals{a}.animdirectory,'/',animals{a}.fileprefix,'cellinfo'])
    tempdays = 1:length(cellinfo);
    includedays = zeros(size(tempdays));
    for d = tempdays
        if ~isempty(cellinfo{tempdays(d)})
            includedays(d) = 1;
        end
    end
    days = tempdays(logical(includedays));
    trajdayprocess(animals{a}.animdirectory,animals{a}.fileprefix,days)
end


%% PROCESS LOW GAMMA (NONREFERENCE)

for a = 1:length(animals)
    %get days for each animal
    load([animals{a}.animdirectory, '/', animals{a}.fileprefix, 'cellinfo'])
    tempdays = 1:length(cellinfo);
    %excludedays if empty in cellinfo
    includedays = zeros(size(tempdays));
    for d = tempdays
        if ~isempty(cellinfo{tempdays(d)})
            includedays(d) = 1;
        end
    end
    
    days = tempdays(logical(includedays));
    lowgammadayprocess(animals{a}.animdirectory, animals{a}.fileprefix, days,'tetfilter','isequal($area,''CA1'') || isequal($area,''CA3'')','nonreference',1,'assignphase',1)
end

%% PROCESS LOW GAMMA (WITH REFERENCE)
% 
% %number of standard deviations for low gamma extraction
% nstd = 1; 
% 
% for a = 3:length(animals)
%     %get days for each animal
%     load([animals{a}.animdirectory, '/', animals{a}.fileprefix, 'cellinfo'])
%     tempdays = 1:length(cellinfo);
%     %excludedays if empty in cellinfo
%     includedays = zeros(size(tempdays));
%     for d = tempdays
%         if ~isempty(cellinfo{tempdays(d)})
%             includedays(d) = 1;
%         end
%     end
%     
%     days{a} = tempdays(logical(includedays));
%     lowgammadayprocess(animals{a}.animdirectory, animals{a}.fileprefix, days{a},'tetfilter','isequal($area,''CA1'') || isequal($area,''CA3'')')
%        for j=1:length(days{a})
%            d = days{a}(j)
%            extractlowgamma(animals{a}.animdirectory, animals{a}.fileprefix, d, 'isequal($area,''CA1'') | isequal($area, ''CA3'')', nstd)
%        end
% end

%% PROCESS HIGH GAMMA (NONREFERENCE)

for a = 1:length(animals)
    %get days for each animal
    load([animals{a}.animdirectory, '/', animals{a}.fileprefix, 'cellinfo'])
    tempdays = 1:length(cellinfo);
    %excludedays if empty in cellinfo
    includedays = zeros(size(tempdays));
    for d = tempdays
        if ~isempty(cellinfo{tempdays(d)})
            includedays(d) = 1;
        end
    end
    
    days{a} = tempdays(logical(includedays));
    
    highgammadayprocess(animals{a}.animdirectory, animals{a}.fileprefix, days{a},'tetfilter','isequal($area,''CA1'') || isequal($area,''CA3'')','nonreference',1,'assignphase',1)
end

%% PROCESS HIGH GAMMA (WITH REFERENCE)

% %number of standard deviations for high gamma extraction
% nstd = 1; 
% 
% for a = 1:length(animals)
%     %get days for each animal
%     load([animals{a}.animdirectory, '/', animals{a}.fileprefix, 'cellinfo'])
%     tempdays = 1:length(cellinfo);
%     %excludedays if empty in cellinfo
%     includedays = zeros(size(tempdays));
%     for d = tempdays
%         if ~isempty(cellinfo{tempdays(d)})
%             includedays(d) = 1;
%         end
%     end
%     
%     days{a} = tempdays(logical(includedays));
%     
%     highgammadayprocess(animals{a}.animdirectory, animals{a}.fileprefix, days{a},'tetfilter','isequal($area,''CA1'') || isequal($area,''CA3'')')
%     for j=1:length(days{a})
%         d = days{a}(j)
%         extracthighgamma(animals{a}.animdirectory, animals{a}.fileprefix, d, 'isequal($area,''CA1'') | isequal($area,''CA3'')', nstd)
% 	end
% end

% %% PROCESS SLOW RIPPLES
% 
% %number of standard deviations for ripple extraction
% nstd = 3; 
% 
% for a = 1:length(animals)
%     %get days for each animal
%     load([animals{a}.animdirectory, '/', animals{a}.fileprefix, 'cellinfo'])
%     tempdays = 1:length(cellinfo);
% 
%     %excludedays if empty in cellinfo
%     includedays = zeros(size(tempdays));
%     for d = tempdays
%         if ~isempty(cellinfo{tempdays(d)})
%             includedays(d) = 1;
%         end
%     end
%     
%     days{a} = tempdays(logical(includedays));
%     
%     %slow ripple analysis
%     slowrippledayprocess(animals{a}.animdirectory, animals{a}.fileprefix, days{a},'area','CA3')
%     for j=1:length(days{a})
%         d = days{a}(j)
%         extractslowripples(animals{a}.animdirectory, animals{a}.fileprefix, d, 'isequal($area,''CA3'')', 0.015, nstd)
%     end
% end
% 
% %% PROCESS THETA
% 
% for a = 8:length(animals)
%     %get days for each animal
%     load([animals{a}.animdirectory, '/', animals{a}.fileprefix, 'cellinfo'])
%     tempdays = 1:length(cellinfo);
% 
%     %excludedays if empty in cellinfo
%     includedays = zeros(size(tempdays));
%     for d = tempdays
%         if ~isempty(cellinfo{tempdays(d)})
%             includedays(d) = 1;
%         end
%     end
%     
%     days{a} = tempdays(logical(includedays));
%     
%     %theta analysis
%     thetadayprocess(animals{a}.animdirectory, animals{a}.fileprefix, days{a},'downsample',1)
% end
