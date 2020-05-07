function assign(varargin)
% assign - make a list of assignments (matlab 5 or higher)
%
%	ASSIGN('VAR1', VAL1, 'VAR2', VAL2, ...) makes the assignments 
%	VAR1 = VAL1; VAR2 = VAL2; ... in the caller's workspace.
%
%	This is most useful when passing in an option list to a
%	function.  Thus in the function which starts:
%		function foo(x,y,varargin)
%		z = 0;
%		assign(varargin{:});
%	the variable z can be given a non-default value by calling the
%	function like so: foo(x,y,'z',4);
startcell = find(cell2mat(cellfun(@(x) ischar(x), varargin, 'UniformOutput', false)),1); %DR. look for first occurance of a string to start the assign in 
vars = {varargin{startcell:2:end}};
vals = {varargin{startcell+1:2:end}};

% use length(vals) not length(vars) in case of odd number of arguments
for i = 1:length(vals)
  assignin('caller', vars{i}, vals{i});
end
