
%clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
subjects = {'S48','S58','S59','S60','S61'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:numel(subjects)
  this_subject = subjects{s};
  load(sprintf('%s/%s/behavior/%s_coordinate_transform.mat', ...
      path_prefix,this_subject,this_subject));

  %{
  load(sprintf('%s/%s/behavior/%s_position_smoothing_params.mat', ...
      path_prefix,this_subject,this_subject));
  %}

  load(sprintf('%s/%s/behavior/%s_linearized_smoothing_params.mat', ...
      path_prefix,this_subject,this_subject));

  filelist = dir(sprintf('%s/%s/behavior/%s_*rawpos.mat', ...
      path_prefix,this_subject,this_subject));
  for i = 1:numel(filelist)
    load(sprintf('%s/%s/behavior/%s', ...
        path_prefix,this_subject,filelist(i).name));
    assert(is_rawpos(rawpos));
    transform_idx = struct_cmp(rawpos,transform, ...
        {'subject','day','epoch','environment'});

    %{
    if (exist('position') == 1)
      position(end+1,1) = estimate_position( ...
          rawpos,transform(transform_idx),params);
    else
      position = estimate_position( ...
          rawpos,transform(transform_idx),params);
    end
    %}

    if ~strcmp(rawpos.environment,'box')
      if (exist('linearized_position') == 1)
        linearized_position(end+1,1) = linearize_Utrack( ...
            rawpos,transform(transform_idx),params);
      else
        linearized_position = linearize_Utrack( ...
            rawpos,transform(transform_idx),params);
      end
    end
    %
  end
  

  %{
  disp(position);
  save(sprintf('%s/%s/behavior/%s_position.mat', ...
      path_prefix,this_subject,this_subject),'position');
  %}

  disp(linearized_position);
  save(sprintf('%s/%s/behavior/%s_linearized_position.mat', ...
      path_prefix,this_subject,this_subject),'linearized_position');
  %
  clear('transform','params','filelist','position','linearized_position');

end

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
% This needs to be run for one subject at a time with subject-specific
% parameter values
this_subject = 'S61';
% All of these parameters remain the same for all subject/sessions except for
% min_front_back_separation, which depends on the length of the head-mounted LED
% boom. Correct values for min_front_back_separation:
% S48: 10 cm
% S58: 10 cm
% S59: 12 cm
% S60: 12 cm
% S61: 12 cm
params = struct( ...
    'front_back_marker_weights' , {[1 0]}         , ...
    'loess_halfwidth'           , {0.5}           , ...
    'num_loess_iterations'      , {3}             , ...
    'roughness_penalty'         , {0.1}           , ...
    'rrm_halfwidths'            , {[0.1 0.1 0.1]} , ...
    'epsilon'                   , {3}             , ...
    'min_stop_duration'         , {1}           , ...
    'min_front_back_separation' , {10}            );

% this function is for S48
%rest_alignment_function = @align_rectangle;

% this function is for S58, S59, S60, S61
rest_alignment_function = @align_cylinder;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (exist([path_prefix '/' this_subject]) ~= 7);  
  error('could not find folder %s/%s',path_prefix,this_subject);
end 
clear('all_rawpos','all_transform');
rawpos_files = dir(sprintf('%s/%s/behavior/%s_*rawpos.mat', ...
    path_prefix,this_subject,this_subject));

for f = 1:numel(rawpos_files)
  clear('rawpos');
  load(sprintf('%s/%s/behavior/%s', ...
      path_prefix,this_subject,rawpos_files(f).name));
  assert(is_rawpos(rawpos) && isscalar(rawpos));
  if (exist('all_rawpos') == 1)
    all_rawpos = [all_rawpos; rawpos];
  else
    all_rawpos = rawpos;
  end
end
assert(is_rawpos(all_rawpos));

environments = {'box','U1','U2'};
for e = 1:numel(environments)
  idx = find(strcmp(environments{e},{all_rawpos(:).environment}));
  if isempty(idx)
    continue;
  end
  switch environments{e}
  case 'box'
    [transform, cp] = rest_alignment_function(all_rawpos(idx(1)));
  case 'U1'
    [transform, cp] = align_Utrack(all_rawpos(idx(1)));
  case 'U2'
    [transform, cp] = align_Utrack(all_rawpos(idx(1)));
  end
  for i = idx
    image(all_rawpos(i).video_frame);
    line(cp(:,1),cp(:,2),'Color','c','Marker','o');
    set(gca,'DataAspectRatio',[1 1 1],'YDir','normal');
    disp(all_rawpos(i));
    pause;
    delete(gcf);
  end
  all_transform(idx) = transform;
end

position = estimate_position(all_rawpos,all_transform,params);
save(sprintf('%s/%s/behavior/%s_position.mat', ...
    path_prefix,this_subject,this_subject),'position');
clear;
%}

