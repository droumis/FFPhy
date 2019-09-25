% CORRPRM1: Write current 'best' matrix to file
%

% RE Strauss, 3/10/95

function corrprm1(best);
  [r,c] = size(best);
  fid = fopen('best.out','k');

  for i=1:r
    fprintf(fid,'%8.6f %7.0f', best(i,1:2));
    for j=3:c
      fprintf(fid,'%3.0f', best(i,j));
    end;
    fprintf(fid,'\n');
  end;

  fclose(fid);
  return;

