function f = testexcludetimes(f, mintime)
%f = testexcludetimes(f)
% removes epochs from the filter if entire epoch is excluded by parameters
% in excludetimes
%
% mintime is minimumtime included in each epoch for it to be included in
% analysis

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
    pos  = loaddatastruct(animaldir, animalprefix, 'pos', totaldays);

    %iterate through the epochs within each data group
    for g = 1:length(f(an).epochs)
        keepepochs = [];
        excludeepochs = [];
        for e = 1:size(f(an).epochs{g},1)

            index = [f(an).epochs{g}(e,:) ];% [ d e ]
            excludetimes = f(an).excludetime{g}{e};
            ontime = pos{index(1)}{index(2)}.data(end,1)-pos{index(1)}{index(2)}.data(1,1); %length of total epoch
            if ~isempty(excludetimes)
                totalexclude = sum(excludetimes(:,2) - excludetimes(:,1));
            else
                totalexclude = 0;
            end
            totalontime = ontime-totalexclude;

            if totalontime >= mintime
                keepepochs = [keepepochs; e];
            else %if totalontime <= mintime
                excludeepochs = [excludeepochs; e];
            end

        end
        f(an).excludeepoch{g} = f(an).epochs{g}(excludeepochs,:);
        f(an).epochs{g} = f(an).epochs{g}(keepepochs,:);
        f(an).excludetime{g} = f(an).excludetime{g}(:,keepepochs);
        f(an).data{g} = f(an).data{g}(:,keepepochs);
    end
end

