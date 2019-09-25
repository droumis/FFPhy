% LOADSTR: Loads a text (ascii) file of text strings (rows).  All rows (strings) must 
%          be of the same length.  Blanks are removed by Matlab before loading.  
%
%     Usage: s = loadstr('fname',cols)
%
%           'fname' = name of text file, in single quotes.
%           cols =    length of lines (rows) in ascii file.
%           -----------------------------------------------
%           s =       string matrix.
%

function s = loadstr(fname,cols)
  fid = fopen(fname);
  s = fscanf(fid,'%s',[cols,inf])';
  fclose(fid);
  return;

