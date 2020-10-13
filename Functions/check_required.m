
function check_required(reqData, vararg)

for s = 1:length(reqData)
    if ~any(cellfun(@(x) strcmp(x,reqData{s}), vararg(1:2:end), 'un', 1))
        error(sprintf('missing data: %s ', reqData{~ismember(reqData,vararg(1:2:end))}));
    end
end
