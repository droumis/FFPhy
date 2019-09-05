

%DR March 2017. set custom colors

function rgbColor = colorPicker(areatags, subareatags, varargin)
colorSet = 'DR1';
if ~isempty(varargin)
    assign(varargin{:});
end

colorlist = setColors(colorSet);
rgbColor = cell(size(areatags,1), size(areatags,2));
for intrpair = 1:size(areatags,2)
    for int = 1:size(areatags,1);
        iareatag = areatags{int,intrpair};
        isubareatag = subareatags{int,intrpair};
        try
            rgbColor{int,intrpair} =  colorlist{find(strcmp(iareatag,colorlist(:,1)) & strcmp(num2str(isubareatag),colorlist(:,2))),3};
        catch
            rgbColor{int,intrpair} = [.9 .9 .9];
        end
    end
end
if length(rgbColor) == 1
    rgbColor = rgbColor{1};
end
end

%add custom color sets as cases below

function colorlist = setColors(colorSet)

switch colorSet
    case 'DR1'
        mecsupf = [0 .6 .2];
        mecdeep = [.47 0 .7];
        colorlist = ...
            {'ref', 'mec', [.9 .9 .9];
            'ref', 'ca1', [.9 .9 .9];
            'ca1', 'd', [.3 .3 .3];
            'mec', 'supf', mecsupf;
            'mec', 'deep', mecdeep;
            'ca1', 'nca1', [1 1 1];
            'mec', 'nsupf', mecsupf;
            'mec', 'ndeep', mecdeep;
            'well', 'input', [1 0 0];
            'well', 'output', [0 0 1]};
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