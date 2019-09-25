
% Process LFP on the selected theta reference electrodes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
subjects = {'S48', 'S58', 'S59', 'S60', 'S61'};
% Each element in these two vectors maps one-to-one to a subject; the same theta
% reference is used across all days. These electrodes were selected by hand by
% visual inspection of raw traces.
lfp_electrodes = struct( ...
    'left'  , { 6,  4,  3,  5,  1}, ...
    'right' , {10,  9, 12,  8, 10} );

load('theta_filter.mat');
%load('widebandtheta_filter.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:numel(subjects)

  for day = 0:7

    files_list = dir(sprintf('%s/%s/continuous/%s_day%d_*_lfp.mat', ...
        path_prefix,subjects{s},subjects{s},day));
    if isempty(files_list)
      continue;
    end

    theta = [];
    %slowtheta = [];
    %widebandtheta = [];

    for i = 1:numel(files_list)
      load(sprintf('%s/%s/continuous/%s', ...
          path_prefix,subjects{s},files_list(i).name));

      % Select the desired LFP electrodes
      lfp = lfp(([lfp(:).electrode] == lfp_electrodes(s).left) | ...
          ([lfp(:).electrode] == lfp_electrodes(s).right));

      %
      if isempty(theta)
        theta = filter_continuous(lfp,Hd,'phase');
      else
        theta = [theta; filter_continuous(lfp,Hd,'phase')];
      end
      %
      
      %{
      if isempty(slowtheta)
        slowtheta = filter_continuous(lfp,Hd,'phase');
      else
        slowtheta = [slowtheta; filter_continuous(lfp,Hd,'phase')];
      end
      %}

      %{
      if isempty(widebandtheta)
        widebandtheta = filter_continuous(lfp,Hd);
      else
        widebandtheta = [theta; filter_continuous(lfp,Hd)];
      end
      %}

      clear('lfp');
    end

    %
    disp(theta);
    save(sprintf('%s/%s/continuous/%s_day%d_theta.mat', ...
        path_prefix,subjects{s},subjects{s},day),'theta');
    clear('theta');
    %
    
    %{
    disp(slowtheta);
    save(sprintf('%s/%s/continuous/%s_day%d_slowtheta.mat', ...
        path_prefix,subjects{s},subjects{s},day),'slowtheta');
    clear('slowtheta');
    %}
   
    %{ 
    disp(widebandtheta);
    save(sprintf('%s/%s/continuous/%s_day%d_widebandtheta.mat', ...
        path_prefix,subjects{s},subjects{s},day),'widebandtheta');
    clear('widebandtheta');
    %}

  end
end

