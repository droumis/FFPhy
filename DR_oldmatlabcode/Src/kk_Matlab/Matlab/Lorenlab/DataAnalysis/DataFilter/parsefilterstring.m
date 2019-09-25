function [variables, expression] = parsefilterstring(inputstr)

pattern = '\$[\w\d_.]+([^\w\d_]|$)';
variables = regexp(inputstr,pattern,'match');
if isempty(variables)
    error(['Not a valid filter:   ',inputstr]);
end

for i = 1:length(variables)
    variables{i} = regexprep(variables{i},'[^\w\d_.]','');
end
expression = regexprep(inputstr,'\$','structVar.'); %replace all '$' with 'structVar.'
    


