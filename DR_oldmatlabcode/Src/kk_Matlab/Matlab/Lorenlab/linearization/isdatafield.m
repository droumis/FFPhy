% toknum = ISDATAFIELD(fieldstring, fieldname)
%     Checks for fieldname as a token in fieldstring and returns the token
%     number if it exists and zero otherwise
function [f] = isdatafield(fs, fn)


toknum = 0;
f = 0;
rem = fs;
while (~isempty(rem))
	[tok rem] = strtok(rem);
	toknum = toknum + 1;
	if (strcmp(tok, fn))
		f = toknum;
		return;
	end
end
