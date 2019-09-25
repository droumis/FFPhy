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
  filename = fullfile(animaldir,sprintf('%s%s%02d.mat',animalprefix,datatype,d));
	eval(sprintf('save %s %s', filename, datatype));
    end
end
