
% Return an indicator matrix for filter conditions meant to exclude event
% times
% input:
%   datastuct: struct array per animal with day, epoch, ripStartTime fields the length of
%   events. also must have animal and dayeps (length # of unique epochs x 2) fields
% options:
%   excludeSets = cell array of strings matching cases of filter creation

% ex: out = getExclusionFilters(datastruct, 'excludeSets', {'firstwell', 'noise'})

function out = getExclusionFilters(datastruct, varargin)

excludeSets = {'firstwell', 'noise'};
noisewin = 1; %s
noisestd = 15; %std
noisekons = 'ca1noisekons';
if ~isempty(varargin)
    assign(varargin{:})
end

for ian = 1:length(datastruct)
    animal = datastruct(ian).animal;
    dayeps = datastruct(ian).dayeps;
    ripstarts = datastruct(ian).ripStartTime;
    andef = animaldef(animal);
    
    for es = 1:length(excludeSets)
        switch excludeSets{es}
            case 'firstwell'
                firstwell{ian} = getpriortofirstwell(andef{2},andef{3}, dayeps);
                
            case 'noise'
                noise{ian} = loaddatastruct(andef{2}, animal, noisekons);
        end
    end
    
    for irip = 1:length(ripstarts)
        day = datastruct(ian).day(irip);
        epoch = datastruct(ian).epoch(irip);
        for es = 1:length(excludeSets)
            switch excludeSets{es}
                case 'firstwell'
                    % pre first well
                    exclrip = isExcluded(ripstarts(irip), firstwell{ian}{day}{epoch}.prefirst_list);
                    out(ian).firstwell(irip, 1) = exclrip;
                case 'noise'
                    % near noise
                    near = abs(noise{ian}{day}{epoch}{1}.starttime - ripstarts(irip)) <= noisewin;
                    out(ian).noise(irip, 1) = any(noise{ian}{day}{epoch}{1}.maxthresh(near) >= ...
                        noisestd);
            end
        end
    end
    
end



end