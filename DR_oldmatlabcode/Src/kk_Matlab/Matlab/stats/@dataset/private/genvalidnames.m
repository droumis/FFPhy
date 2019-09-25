function [names,wereModified] = genvalidnames(names)
%GENVALIDNAMES Construct valid identifiers from a list of names.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:31:57 $

wereModified = false;

% Loop over all names and make them valid identifiers.  This is the same code
% that genvarnames uses, without the uniqueness checking.
for k = 1:numel(names)
    name = names{k};

    if ~isvarname(name)
        wereModified = true;

        % Insert x if the first column is non-letter.
        name = regexprep(name,'^\s*+([^A-Za-z])','x$1', 'once');

        % Replace whitespace with camel casing.
        [StartSpaces, afterSpace] = regexp(name,'\S\s+\S');
        name(afterSpace) = upper(name(afterSpace));
        name = regexprep(name,'\s*','');
        if (isempty(name))
            name = 'x';
        end
        % Replace non-word character with its HEXIDECIMAL equivalent
        illegalChars = unique(name(regexp(name,'[^A-Za-z_0-9]')));
        for illegalChar=illegalChars
            if illegalChar <= intmax('uint8')
                width = 2;
            else
                width = 4;
            end
            replace = ['0x' dec2hex(illegalChar,width)];
            name = strrep(name, illegalChar, replace);
        end

        % Prepend keyword with 'x' and camel case.
        if iskeyword(name)
            name = ['x' upper(name(1)) name(2:end)];
        end

        % Truncate name to NAMLENGTHMAX
        name = name(1:min(length(name),namelengthmax));

        names{k} = name;
    end
end
