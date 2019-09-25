function disp(t)
%DISP Display a CLASSREGTREE object.
%   DISP(T) prints the CLASSREGTREE object T.
%
%   See also CLASSREGTREE, CLASSREGTREE/VIEW.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:56:16 $

isLoose = strcmp(get(0,'FormatSpacing'),'loose');
maxWidth = get(0,'CommandWindowSize'); maxWidth = maxWidth(1);

% disp(struct(t));

if (isLoose), fprintf('\n'); end

% Get some information about the whole tree
maxnode = numel(t.node);
nd = 1 + floor(log10(maxnode)); % number of digits for node number
names = t.names;
if isempty(names)
    names = strcat('x',strread(sprintf('%d\n',1:t.npred),'%s\n'));
end
isregression = isequal(t.method,'regression');
fprintf('Decision tree for %s\n',t.method);

% Display information about each node
for j=1:maxnode
    if any(t.children(j,:))
        % branch node
        vnum = t.var(j);
        vname = names{abs(vnum)};
        if vnum>0        % continuous predictor "<" condition
            cond = sprintf('%s<%g',vname,t.cut(j));
        else             % categorical predictor, membership condition
            cats = t.catsplit{t.cut(j),1};
            if isscalar(cats)
                cond = sprintf('%s=%d',vname,cats);
            else
                set = deblank(num2str(cats,'%d '));
                cond = sprintf('%s in {%s}',vname,set);
            end
        end
        fprintf('%*d  if %s then node %d else node %d\n',nd,j,cond,t.children(j,:));
    else
        % terminal node, display fit (regression) or class assignment
        if isregression
            fprintf('%*d  fit = %g\n',nd,j,t.class(j));
        else
            fprintf('%*d  class = %s\n',nd,j,t.classname{t.class(j)});
        end
    end
end
if (isLoose), fprintf('\n'); end