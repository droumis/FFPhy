function f = setfiltertetrodepairs(f, filterStringArray)
% f = setfiltercells(f, filterStringArray)
% For each epoch in the filter F, this function finds the indices to the
% tetrodes that satisfy each of the filter condtions in the 1st and 2nd
% element of the filterStringArray.  
%
% output list is a double array with one row for each pair and 4 columns:
% dataset epoch tet1 tet2
%
% The syntax for each filterString is defined in EVALUATEFILTER.m. 
% The animal and desired epochs for the filter need to be predefined. 
% Assumes that each animal's data  folder contains a file 'tetinfo.mat' that 
% contains a tetrode structure with information about each tetrode.

if (length(filterStringArray) ~= 2) 
    error('filterStringArray must have 2 elements');
end

for an = 1:length(f)
    if isempty(f(an).animal)
        error(['You must define an animal for the filter before filtering the tetrodes'])
    end
    if isempty(f(an).epochs)
        error(['You must define the desired epochs for the filter before filtering the tetrodes'])
    end
    datadir = f(an).animal{2};
    animalprefix = f(an).animal{3};
    tetinfo = loaddatastruct(datadir,animalprefix,'tetinfo');

    for i = 1:length(f(an).epochs)
        if isempty(f(an).epochs{i})
            f(an).data{i} = [];
        end
        for j = 1:size(f(an).epochs{i},1)
	    for k = 1:length(filterStringArray)
		tmplist{k} = evaluatefilter(tetinfo{f(an).epochs{i}(j,1)}{f(an).epochs{i}(j,2)}, filterStringArray{k});
		tmpind{k} = 1:size(tmplist{k}, 1);
	    end
	    % combine the lists using the option in the first element of the
	    % filterStringArray
	    [i1 i2] = meshgrid(tmpind{1}, tmpind{2});
	    tetlist = [tmplist{1}(i1,:) tmplist{2}(i2,:)];

		% get rid of identical elements
		valid = (sum(tetlist(:,1) == tetlist(:,2), 2) ~= 2);
		tetlist = tetlist(valid,:);

	    if (length(tetlist))
		% check for pairs listed twice because the order was reversed
		tmptetlist = [tetlist(:,2) tetlist(:,1)];
		[c i1 i2] = intersect(tetlist, tmptetlist, 'rows');
		valid = ones(size(tetlist,1), 1);
		% for each member of tmptetlist that is in tetlist, we
		% need to delete the correponding member of tetlist if
		% it's index is larger than the index in tmptetlist
		valid(i2) = (i2 < i1);
		f(an).data{i}{j} = tetlist(logical(valid),:);
	    else
		f(an).data{i}{j} = [];
	    end
        end
    end
end
