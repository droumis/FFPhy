% function convertTimestampsTo64(infile,outfile)
function aa = convertTimestampsTo64(infile,outfile)

inf = fopen(infile,'r');

if inf == -1
  error('Cannot open input file %s\n',infile);
end

if exist(outfile,'file')
  error('Output file %s exists!\n',outfile);
end

outf = fopen(outfile,'w');
if outf == -1
  error('Cannot open output file %s\n',outfile);
end

%Proper infile header:
%%BEGINHEADER
% File type:    Binary
% Extraction type:      mpeg file offset
% Fields:       offset (unsigned int)
%%ENDHEADER

ss1 = fgets(inf,1000);
ss2 = fgets(inf,1000);
ss3 = fgets(inf,1000);
ss4 = fgets(inf,1000);
ss5 = fgets(inf,1000);
if ~strcmp(ss1(1:end-1),'%%BEGINHEADER') |  ...
  ~strcmp(ss2(1:end-1),'% File type:    Binary') |  ...
  ~strcmp(ss3(1:end-1),'% Extraction type:      mpeg file offset') |  ...
  ~strcmp(ss4(1:end-1),'% Fields:       offset (unsigned int)') |  ...
  ~strcmp(ss5(1:end-1),'%%ENDHEADER')
  % error('Input file header is not as expected for 32 bit offsets.');
end

fprintf(outf,'%s',ss1);
fprintf(outf,'%s',ss2);
fprintf(outf,'%s',ss3);
fprintf(outf,'% Fields:       offset (unsigned long long)\n',ss4);
fprintf(outf,'%s',ss5);

aa = fread(inf,'uint32');

fwrite(outf,aa,'uint64');

fclose(inf);
fclose(outf);
