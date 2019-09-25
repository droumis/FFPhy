% CORRPRM2: Write current status and distribution to files
%

% RE Strauss, 3/15/95

function corrprm2(iter,bestval,p,distrib);
  fid = fopen('iterlog.out','a');    % Append current status
  fprintf(fid,'%7.0f %8.6f', iter,bestval);
  for i=1:length(p)
    fprintf(fid,' %3.0f', p(i));
  end;
  fprintf(fid,'\n');
  fclose(fid);

  save distrib;                       % Save current distribution

  return;

