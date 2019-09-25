function savedatastruct(datastruct, animaldir, animalprefix, datatype) 
% savedatastruct(datastruct, animaldir, animalprefix, datatype)
%
% Save the components of a data cell array.  If all of the days have been
% combined into one variable, this saves them individually.  
% variable.  Datatype is a string with the base name of the files (e.g.
% 'spikes'.  

for d = 1:length(datastruct)
    if (~isempty(datastruct{d}))
	eval(['clear ', datatype]);
	eval(sprintf('%s{%d} = datastruct{%d};', datatype, d, d));
	%[datatype, '{', = ', datastruct, '{', num2str(d), '};']);
	eval(sprintf('save %s%s%s%02d.mat %s', animaldir, animalprefix, datatype, d, datatype));
    end
end
