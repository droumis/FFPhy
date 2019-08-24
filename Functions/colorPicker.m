

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
end

%add custom color sets as cases below

function colorlist = setColors(colorSet)

switch colorSet
    case 'DR1'
        colorlist = ...
            {'ref', 'mec', [.5 .5 .5];
            'ref', 'ca1', [.5 .5 .5];
            'ca1', 'd', [.06 .07 .07];
            'mec', 'supf', [.65 .22 .23]
            'mec', 'deep', [.07 .10 .80];
            'ca1', 'nca1', [1 1 1];
            'mec', 'nsupf', [.65 .22 .23];
            'mec', 'ndeep', [.07 .10 .80]};
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