function [n] = fig2pdf(directory)
% N = fig2pdf(directory)
%	converts all .fig files to pdfs and returns the number of files
%	converted

cwd = pwd;
n = 0;
cd(directory);
d = dir;

for i = 3:length(dir)
    if (~isempty(strfind(d(i).name, 'fig')))
	hgload(d(i).name);
	newname = d(i).name;
	newname(end-2:end) = 'pdf';
	eval(['print -dpdf ', newname]);
	close
	n = n+1;
    end
end

cd(cwd);
