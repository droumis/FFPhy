
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
subjects = {'S48','S58','S59','S60','S61'};

params = struct( ...
    'track_interval'    , [-170, +170], ...
    'minimum_speed'     , 10, ...
    'discretization_timestep'  , uint32(1) );
% Nominal half-length of the track is (22 + 147) = 169 cm. Note that this
% includes end sections of the track where the rat is turning; we will need to
% chop these off later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
linearized_passes = [];
for i = 1:numel(subjects)
  load(sprintf('%s/%s/behavior/%s_linearized_position.mat', ...
      path_prefix,subjects{i},subjects{i}));
  linearized_passes = [ linearized_passes; ...
      find_linearized_passes(linearized_position,params) ];
end
%}

