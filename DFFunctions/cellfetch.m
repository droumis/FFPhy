
function output = cellfetch(cellinput, field, varargin)
% searchresults = cellfetch(input, field)

%DR updated to add functionality of returning all the tags for that day/epoch/ntrode as a struct 
%if you pass in the optional 'alltags' = 1; default is still 0. You must
%still pass in an empty string for the 'field' field.
%e.g.: tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);

% This function parces a nested cell variable 'input', and assumes nothing about
% the structure of the internal cell tree. It returns a structure with the value of the
% field variable in the field 'field' at the end of each cell branch.  The
% indices of each entry is put in the index field of the output.

alltags = 0;

% process varargin and overwrite default values
if (~isempty(varargin))
    assign(varargin{:});
end

if ~isstr(field)
    error('FIELD must be a string');
end
output.index = [];
output.values = {};
indexcount = 0;
% if ~exist('branchindex');
%     branchindex = [];
% end

if iscell(cellinput)
    for i = 1:length(cellinput)
        if isstruct(cellinput{i})
            indexcount = indexcount + 1;
            output.index(indexcount,1) = i;
            if alltags == 1
                try
                    eval(['output.values{indexcount,1} = cellinput{i};']);
                catch
                    output.values{indexcount,1} = [];
                end
            else
                try
                    eval(['output.values{indexcount,1} = cellinput{i}.',field,';']);
                catch
                    output.values{indexcount,1} = [];
                end
            end
%             if isfield(cellinput{i},field)
%                 output.values{indexcount,1} = getfield(cellinput{i},field);
%             else
%                 output.values{indexcount,1} = [];
%             end
        elseif iscell(cellinput{i})
%             varargin = [varargin, {'branchindex'}, {i}];
            [branch] = cellfetch(cellinput{i}, field, varargin{:});
            numentries = size(branch.index,1);
            branch.index = [repmat(i,[numentries,1]) branch.index];
            output.index(indexcount+1:indexcount+numentries,1:size(branch.index,2)) = branch.index;
            indexcount = indexcount+numentries;
            output.values = [output.values; branch.values];
        end
    end
end
            
