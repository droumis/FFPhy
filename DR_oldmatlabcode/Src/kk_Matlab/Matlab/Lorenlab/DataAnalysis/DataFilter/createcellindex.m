function out = createcellindex(index)

out = [];
for i = 1:length(index)
    out = [out, '{', num2str(index(i)), '}'];
end
        