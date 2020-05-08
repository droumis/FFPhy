function [R, B] = target3d(target, armconfig, tilt, rot, dvintersect, quiet)
% TARGET3D: target a location using a tilted, rotated stereotax
%
% Specify the location of your target, and the angle of approach you would
% like to take to that target. This code generates the translation along
% each of the 3 (possibly non-orthogonal) stereotactic axes needed to reach
% the target.
%
% There are many different conventions in use to describe an angle in 3
% dimensions. Here, we specify the approach angle in terms of its tilt away
% from vertical, and its rotation about the same vertical axis. (See
% INPUTS and EXAMPLES below for specifics).
%
% Note that many possible approach angles are awkward or infeasible using a
% given stereotaxic apparatus and probe, so it is recommended that you
% first 'rough out' the tilt and rotation that you would like to use to
% reach your target, then use this code to generate precise coordinates.
%
% (Protip: you can use a second, untilted and unrotated arm to specify your
% target location in space. See:
%  A. J. Greenshaw, J. Neurosci. Methods 78, 169-172 (1997))
%
% You can also use this method to do a real-world 'sanity check' of the 
% outputs of this function.
%
% ******
% N.B. All coordinates are given as [ML AP DV], relative to Bregma, where:
%          ML (mediolateral): animal right is +ve
%          AP (anteroposterior): anterior is +ve
%          DV (dorsoventral): dorsal (up) is +ve
% (Same order and sign as used by Kopf digital display)
% ******
%
% INPUTS:
%  target = [ML AP DV] atlas coordinates of target location.
%
%  armconfig = 'LR', or 'FB': is the stereotax arm 'knuckle' configured so
%        that, when a probe is located over bregma, the arm tilts 
%        left-right, or front-back?
%
%  tilt = degrees of tilt away from vertical. Must be positive, and less 
%        than 90 degrees! (direction of tilt is specified by rotation)
%
%  rot = degrees of CW rotation of tilted DV axis away from animal nose.
%        (I.e. 0 means arm tilts towards animals nose, 90 means arm tilts
%        to the animal's right, -90 (or 270) means arm tilts to the 
%        animal's left.)
%
% OPTIONAL INPUTS:
%  dvintersect = DV heights (in old coordinates) at which to calculate
%        intersection of approach path. Default = []. Useful for 
%        establishing craniotomy locations in unrotated coordinates before 
%        rotating the stereotax arms.
%
%  quiet = true/false: silence text description of outputs, just return
%        coordinates. Default = false.
%
% OUTPUTS:
%  R = coordinates to use to reach the targeted location, using the 
%        rotated/tilted stereotax arm.
%
%  B = location, in *unrotated* coordinates, of the intersection of the 
%        approach path with a DV plane specified by dvintersect. 
%
% EXAMPLES:
%  To target a location at ML 0.25, AP 0.7, DV -6.6, with a stereotax arm 
%  tilted 20 degrees to the right:
%     R = target3d([0.25 0.7 -6.6], 'LR', 20, 90);
%
%     -->R: +2.49 / +0.70 / -6.12 
%
%  To target a location at ML -2, AP -2, DV -3.1, with a stereotax arm 
%  tilted 30 degrees forwards, and rotated 20 degrees counter-clockwise 
%  (i.e. the arm is leaning forward and towards the animal's left).
%     R = target3d([-2 -2 -3.1], 'FB', 30, -20);
%
%     -->R: -2.78 / +0.63 / -3.58
%
%  Same target, also calculate the location of the intersection of the
%  approach path with a plane slightly below bregma (e.g. for marking a 
%  craniotomy site before tilting the arm:
%     [R B] = target3d([-2 -2 -3.1], 'FB', 30, -20, -0.3);
%
%     -->R: -2.78 / +0.63 / -3.58
%     -->B: -2.55 / -0.48 / -0.30
%
%  Another target, with arm approaching the target from behind and to the
%  right (i.e. arm tilted backwards).
%     R = target3d([1.2 -1.3 -4.4],'FB',25,150);
%
%     -->R: +2.57 / -4.36 / -4.85
%
% Tom Davidson, Stanford University (tjd@alum.mit.edu) 2010-2016

% TODO: Make sterotax arm config an output

%% Check inputs
switch lower(armconfig)
    case 'lr'
        tiltLR = true;
    case 'fb'
        tiltLR = false;
    otherwise
        error('Must provide ''armconfig'': either ''LR'', or ''FB''');
end

if ~exist('dvintersect', 'var')
    dvintersect = [];
end

if ~exist('quiet', 'var')
    quiet = false;
end

% Check target
if target(3)>0
    warning('''target'' specifies a location *above* bregma--use negative values to represent dorsal direction');
end

% Check tilt inputs
if tilt < 0
    error('''tilt'' parameter must be positive. Specify direction of tilting with ''rot'' parameter. See help.');
end

if tilt >90
    error('''tilt'' must be <90 degrees. See help.');
end

if tilt > 40
    warning('''tilt'' greater than 40 degrees is unusual--are you sure?');
end

if mod (tilt * 360,pi) == 0 && tilt ~=0
    warning('You specified ''tilt'' as a fraction of pi. ''tilt'' expects degrees, not radians.');
end

% Check rot inputs
if rot <= -360 || rot >=360
    error('Please specify a rotation ''rot'' between (-360,360) degrees. See help.');
end

if mod (rot * 360,pi) == 0 && rot ~=0
    warning('You specified ''rot'' as a fraction of pi. ''rot'' expects degrees, not radians.');
end

% convert requested rotation to interval [0,360)
rot = mod(rot,360);

% Confirm that the specified stereotax arm configuration is reasonable for
% the requested approach angle.

if tilt ~=0
    
    arm_tilted_frontback = rot >315 || rot <45 || (rot > 135 && rot <225);
    
    if tiltLR
        if arm_tilted_frontback
            error('Stereotax arm configuration is wrong for requested approach: use ''FB'' (front-back) tilt instead');
        end
    else
        if ~arm_tilted_frontback
            error('Stereotax arm configuration is wrong for requested approach: use ''LR'' (left-right) tilt instead');
        end
    end
end

%% Handle 'backward'/'leftward' tilts

% While we use spherical coordinates to specify the angle of approach for
% targeting, we do not actually rotate stereotax arms through 360 degrees
% to achieve this targeting.
%
% For instance, when approaching the target with a 10-degree left tilt
% (tilt = 10, rot = 270), we do not tilt the arm by 10 degrees forward then
% rotate it through 270 degrees. Instead, we just tilt the arm leftwards.
% However this results in ML-positive on the stereotax display being
% towards the animal's right, even though in our rotated coordinates,
% traveling in the positive direction ML should move left.
%
% (Similarly, when the arm is configured for front-back tilt, we implement 
% backwards tilts by tilting the arm backwards, and not by tilting forwards
% then rotating the arm through 180 degrees)
% 
% We correct for both these cases here. (NB, we could equivalently correct 
% for this just by flipping sign of ML output, or rotating x1)

if (~tiltLR && rot>90 && rot<270) || ... % 'backward' tilt
   (tiltLR && rot>180) % 'leftward' tilt
    tilt = -tilt;
    rot = rot-180;
    
    % convert to interval [0,360)
    rot = mod(rot,360);
end

%% Handle tilts towards animal

% Typically, on a stereotax with an arm on either side of the animal the
% user chooses to use the arm on the side that allows the arm to tilt away
% from the animal. When this is not done, 

%% Convert inputs to polar coordinates

% Convert DV vector specified by 'tilt' and 'rot' to polar coordinates as
% used by Matlab's built-in sph2cart (theta, phi, rho): 
%  theta = CCW rotation of DV axis in x0y0 plane away from positive 'X' 
%        (i.e. animal right=0, nose=90)
%  phi = elevation of DV axis from x0y0 plane (i.e. from horizontal) 
theta = sub_deg2rad(-rot+90); % theta on (-3/2 * pi, pi/2]
phi = sub_deg2rad(-tilt+90); % phi on [0 pi/2];

theta = mod(theta,2*pi);
phi = mod(phi,2*pi);

%% Calculate rotation matrix

% Specify x1 (new ML) vector in x0y0z0 space

if tiltLR, % Tilt Left-Right  case:
    % Relative to DV, rotation (theta) is same, tilt (phi) is 90deg CW
    [x01(1) x01(2) x01(3)] = sph2cart(theta, phi-pi/2,1);

else % Tilt Front-Back case
    
    % When knuckle is configured for 'FB', tilt is about the ML axis, so 
    % 'tilt'/phi does not affect the ML axis vector. (i.e. movements in ML 
    % are in the x0y0 plane).
    
    % New ML axis rotation (theta) is 90deg CW relative to the rotation of
    % the DV axis. ML axis not tilted (phi=0);
    [x01(1) x01(2) x01(3)] = sph2cart(theta-pi/2, 0, 1);
        
end

% Specify y1 (new AP) vector in x0y0z0 space 

% (AP stereotax axis can't tilt or rotate, so y01 is still the unit vector 
% along the y-axis)
y01 = [0 1 0];


% Specify z1 (new DV) vector in x0y0z0 space 

% (tilted, rotated as per inputs)
[z01(1) z01(2) z01(3)] = sph2cart(theta, phi,1);


% the inverse of this matrix is a rotation matrix that can be used to get 
% coords in x1y1z1 given target in x0y0z0
rotmat = inv([x01' y01' z01']);


%% Generate coordinates in new reference frame
R = rotmat * target';

% output as row vector
R = R';


%% Calculate intersection of approach vector with horizontal planes

% Optionally calculate intersection of approach vector with requested
% planes parallel to x0y0 plane (Bregma plane), in x0y0z0 coordinates.
% (Note: no dependence on rotation matrix, or on rotated coords x1y1z1)

if ~isempty(dvintersect)
  % radial distance from target to Bregma plane
  Br = [-target(3)+dvintersect] / sin(phi);

  % x0y0 distance from target to Bregma plane
  [B(:,1), B(:,2), B(:,3)] = sph2cart(theta,phi,Br);

  % plus x0y0 distance to target
  B = bsxfun(@plus, B, target);
  
else
  B = [];
end
  

%% Pretty-print outputs

if ~quiet
    disp(sprintf('\nAll coords ML(Right+)/AP(Anterior+)/DV(Dorsal+), as on Kopf stereotax digital display.\n'));
    disp(sprintf('R (target location in new coords): %+0.2f / %+0.2f / %+0.2f \n', R));
    if ~isempty(B),
      disp(sprintf('B (intersection with DV plane in unrotated coords): %+0.2f / %+0.2f / %+0.2f \n', B'));
    end
end



function rad = sub_deg2rad(deg)
% DEG2RAD-converts a matrix of angles in degrees into (-pi pi] radians

rad = mod(deg * (pi/180), 2*pi);

radneg = rad > pi;
rad(radneg) = rad(radneg) - 2*pi;
