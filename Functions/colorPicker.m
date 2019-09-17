%{
tags is a cell array of key names
returns same size cell array of 3 element rgb arrays
%}

function clrOut = colorPicker(tags, varargin)
colorSet = 'DR1'; subtags = {};
if ~isempty(varargin); assign(varargin{:}); end
if ~isa(tags, 'cell')
    tags = num2cell(tags); end

clrdict = setColors(colorSet);
clrOut = cell(size(tags,1), size(tags,2));
for itagC = 1:size(tags,2)
    for itagR = 1:size(tags,1)
        itag = tags{itagR,itagC};
        if ~isempty(subtags)
            isubtag = subtags{itagR,itagC}; end
        try
            if ~isempty(subtags)
            clrOut{itagR,itagC} =  clrdict{find(strcmp(num2str(itag),clrdict(:,1)) & ...
                strcmp(num2str(isubtag),clrdict(:,2))),3}; 
            else; clrOut{itagR,itagC} =  clrdict{find(strcmp(itag,clrdict(:,1)),1),3}; 
            end
        catch
            fprintf('invalid tag entry at row%d col%d \n', itagR, itagC);
            clrOut{itagR,itagC} = [.9 .9 .9];
        end
    end
end
if length(clrOut) == 1
    clrOut = clrOut{1};
end
end

%add custom color sets as cases below

function clrdict = setColors(colorSet)

switch colorSet
    case 'DR1'
        mecsupf = [0 .6 .2];
        mecdeep = [.47 0 .7];
%         dio = lines(6);
        dio = cbrewer('qual', 'Paired', 6);
        clrdict = ...
            {'ref', 'mec', [.9 .9 .9];
            'ref', 'ca1', [.9 .9 .9];
            'ca1', 'd', [.3 .3 .3];
            'mec', 'supf', mecsupf;
            'mec', 'deep', mecdeep;
            'ca1', 'nca1', [1 1 1];
            'mec', 'nsupf', mecsupf;
            'mec', 'ndeep', mecdeep;
            'well', 'input', [1 0 0];
            'well', 'output', [0 0 1];
            '1', [], dio(1,:);
            '2', [], dio(2,:);
            '3', [], dio(3,:);
            '22', [], dio(4,:);
            '23', [], dio(5,:);
            '24', [], dio(6,:)};
    otherwise
        error('pick an existing colorset or make one')
end
end
% 
%             'ca1', 'all', [.06 .07 .07];
%             'sub', 'd', [.50 .50 .50];
%             'sub', 'all', [.50 .50 .50];
%             'mec', '2', [.23 .63 .70];
%             'mec', '3', [.25 .46 .5];
%             'mec', '5', [.07 .10 .80];
%             'mec', '6', [.13 .07 .53];
%             'mec', 'all', [.13 .07 .53];
%             'por', '2', [.53 .06 .26];
%             'por', '3', [.53 .06 .26];
%             'por', '4', [.53 .06 .26];
%             'por', '5', [.53 .06 .26];
%             'por', '6', [.53 .06 .26];
%             'por', 'all', [.53 .06 .26];
%             'v2l', '2', [.65 .22 .23];
%             'v2l', '3', [.65 .22 .23];
%             'v2l', '4', [.65 .22 .23];
%             'v2l', '5', [.65 .22 .23];
%             'v2l', '6', [.65 .22 .23];
%             'v2l', 'all', [.65 .22 .23];