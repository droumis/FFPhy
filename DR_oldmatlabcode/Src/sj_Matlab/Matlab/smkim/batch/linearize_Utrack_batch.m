
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
subjects = {'S48','S58','S59','S60','S61'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read_rawpos(posfilename,mpegfilename,mpegoffsetfilename,timestampfilename,session)

for s = 1:numel(subjects)
  this_subject = subjects{s};
  try
    load(sprintf('%s/%s/behavior/%s_position.mat',path_prefix, ...
        this_subject,this_subject));
    assert(exist('position') == 1);
  catch
    error('could not load position data for %s',this_subject);
  end

  linpos = linearize_Utrack(position(~strcmp('box',{position(:).environment})));
  save(sprintf('%s/%s/behavior/%s_linpos.mat',path_prefix, ...
      this_subject,this_subject),'linpos');

  clear('position','linpos');

end
