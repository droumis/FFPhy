function genericApplyToRows = DRrowfun(func, matrix)

%dr 11/2014 version of row-wise function that works on matrices.. doesn't work yet

applyToGivenRow = @(func, matrix) @(row) func(matrix(row, :));
newApplyToRows = @(func, matrix) arrayfun(applyToGivenRow(func, matrix), 1:size(matrix,1), 'UniformOutput', false)';
takeAll = @(x) reshape([x{:}], size(x{1},2), size(x,1))';
genericApplyToRows = @(func, matrix) takeAll(newApplyToRows(func, matrix));