for a = 1:7
    switch a
        case 1 %first exposure to first track on recording day 1
            animalinfo = {'Miles', '/data12/mkarlsso/Mil', 'mil'};
        case 2 %first exposure to first track on recording day 1
            animalinfo = {'Ten', '/data12/mkarlsso/Ten', 'ten'};
        case 3 %first exposure to first track on recording day 1
            animalinfo = {'Dudley', '/data12/mkarlsso/Dud', 'dud'};
        case 4 %first exposure to first track on recording day 1
            animalinfo = {'Conley', '/data12/mkarlsso/Con', 'con'};
        case 5
            animalinfo = {'Frank', '/data12/mkarlsso/Fra', 'fra'};
        case 6
            animalinfo = {'Nine', '/data12/mkarlsso/Nin', 'nin'};
        case 7
            animalinfo = {'Bond', '/data12/mkarlsso/Bon', 'bon'}; % slower at first, paused at CP lots
        case 8
            % animalinfo = {'Alex', '/data12/mkarlsso/Ale', 'ale'}; %questionable cuz coded 2 tracks same
            
    end
    animals{a}.animdirectory = animalinfo{2};
    animals{a}.fileprefix = animalinfo{3};
end

nstd = 3; %number of standard deviations for ripple extraction

% allseg = [];
% for a = 1:length(animals)
%     linpos = loaddatastruct(animals{a}.animdirectory, animals{a}.fileprefix, 'linpos');
%     for d=1:length(linpos)
%         for e = 1:length(linpos{d})
%             if ~isempty(linpos{d}{e})
%                 seglength = linpos{d}{e}.segmentInfo.segmentLength(1:2:end);
%                 allseg = [allseg; seglength'];
%             end
%         end
%     end
% end
%

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
    %rippledayprocess(animals{a}.animdirectory, animals{a}.fileprefix, days{a})
    for j=1:length(days{a})
        d = days{a}(j)
        extractripples(animals{a}.animdirectory, animals{a}.fileprefix, d, -1, 0.015, nstd, 'samethreshperday', 0)%-1 means all tetrodes processes
    end
    end












