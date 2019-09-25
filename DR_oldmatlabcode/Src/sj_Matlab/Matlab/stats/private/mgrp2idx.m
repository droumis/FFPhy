function [ogroup,glabel,gname,multigroup,maxgroup] = mgrp2idx(group,rows,sep)
%MGRP2IDX Convert multiple grouping variables to index vector
%   [OGROUP,GLABEL,GNAME,MULTIGROUP,MAXGROUP] = MGRP2IDX(GROUP,ROWS)
%   takes the inputs GROUP, ROWS, and SEP.  GROUP is a grouping variable
%   (categorical variable, numeric vector, numeric matrix, string matrix, or
%   cell array of strings) or a cell array of grouping variables.  ROWS is the
%   number of observations.  SEP is a separator for the grouping variable
%   values.
%
%   The output OGROUP is a vector of group indices.  GLABEL is a cell
%   array of group labels, each label consisting of the values of the
%   various grouping variables separated by the characters in SEP.
%   GNAME is a cell array containing one column per grouping variable
%   and one row for each distinct combination of grouping variable
%   values.  MULTIGROUP is 1 if there are multiple grouping variables
%   or 0 if there are not.  MAXGROUP is the number of groups before any
%   unused categories are omitted.

%   Copyright 1993-2006 The MathWorks, Inc.
%   $Revision: 1.4.2.5 $  $Date: 2006/11/11 22:57:34 $

multigroup = (iscell(group) && size(group,1)==1) ||...
    (isnumeric(group) && ~isvector(group) && ~isempty(group));
if (~multigroup)
    [ogroup,gname,maxgroup] = grp2idx(group);
    glabel = gname;
else
    % Group according to each distinct combination of grouping variables
    ngrps = size(group,2);
    grpmat = zeros(rows,ngrps);
    namemat = cell(1,ngrps);

    % Get integer codes and names for each grouping variable
    if iscell(group)
        for j=1:ngrps
            [g,gn] = grp2idx(group{1,j});
            if (size(g,1)~=rows)
                error('stats:mgrp2idx:InputSizeMismatch',...
                      'All grouping variables must have %d rows.',rows);
            end
            grpmat(:,j) = g;
            namemat{1,j} = gn;
        end
    else
        for j=1:ngrps
            [g,gn] = grp2idx(group(:,j));
            grpmat(:,j) = g;
            namemat{1,j} = gn;
        end
    end;

    % Find all unique combinations
    wasnan = any(isnan(grpmat),2);
    grpmat(wasnan,:) = [];
    [urows,ui,uj] = unique(grpmat,'rows');
    % Create a cell array, one col for each grouping variable value
    % and one row for each observation
    ogroup = NaN(size(wasnan));
    ogroup(~wasnan) = uj;
    gname = cell(size(urows));

    for j=1:ngrps
        gn = namemat{1,j};
        gname(:,j) = gn(urows(:,j));
    end

    % Create another cell array of multi-line texts to use as labels
    glabel = cell(size(gname,1),1);
    if (nargin > 2)
        nl = sprintf(sep);
    else
        nl = sprintf('\n');
    end
    fmt = sprintf('%%s%s',nl);
    lnl = length(fmt)-3;        % one less than the length of nl
    for j=1:length(glabel)
        gn = sprintf(fmt, gname{j,:});
        gn(end-lnl:end) = [];
        glabel{j,1} = gn;
    end
    maxgroup = length(glabel);
end
