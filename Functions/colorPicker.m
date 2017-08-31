

%DR March 2017. set custom colors

function rgbColor = colorPicker(colorSet, areatags, subareatags)
colorlist = setColors(colorSet);
rgbColor = [];
for int = 1:length(areatags);
    iareatag = areatags{int};
    isubareatag = subareatags{int};
    rgbColor =  [rgbColor; colorlist{find(strcmp(num2str(iareatag),colorlist(:,1)) & strcmp(num2str(isubareatag),colorlist(:,2))),3}];
end
end

%add custom color sets as cases below

function colorlist = setColors(colorSet)

switch colorSet
    case 'DR1'
        colorlist = ...
            {'ref', ' ', [.31 .46 .19];
            'ca1', 'd', [.06 .07 .07];
            'sub', 'd', [.50 .50 .50];
            'mec', '2', [.23 .63 .70];
            'mec', '3', [.25 .46 .5];
            'mec', '5', [.07 .10 .80];
            'mec', '6', [.13 .07 .53];
            'por', '2', [.53 .06 .26];
            'por', '3', [.53 .06 .26];
            'por', '4', [.53 .06 .26];
            'por', '5', [.53 .06 .26];
            'por', '6', [.53 .06 .26];
            'v2l', '2', [.65 .22 .23];
            'v2l', '3', [.65 .22 .23];
            'v2l', '4', [.65 .22 .23];
            'v2l', '5', [.65 .22 .23];
            'v2l', '6', [.65 .22 .23]};
    otherwise
        error('pick an existing colorset or make one')
end
end