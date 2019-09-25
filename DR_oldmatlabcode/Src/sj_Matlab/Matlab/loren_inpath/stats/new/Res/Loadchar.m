% LOADCHAR: Loads a text (ascii) file of text strings (rows).  All rows (strings) must 
%           be of the same length.  Blanks are removed by Matlab before loading.  
%
%     Usage: s = loadchar('fname',cols)
%
%           'fname' = name of text file, in single quotes.
%           cols =    length of lines (rows) in ascii file.
%           -----------------------------------------------
%           s =       string matrix.
%

% RE Strauss, 9/30/01
%   10/9/02 - renamed from 'loadstr';
%             corrected documentation;
%             changed final 'end' to 'return'.

function s = loadchar(fname,cols)
  fid = fopen(fname);
  s = fscanf(fid,'%s',[cols,inf])';
  fclose(fid);
  return;

