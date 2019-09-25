% ASYMAREA: Provides the objective function for minimizing the half-area 
%           asymmetry for a polygon.
%
%       Usage:
%
%           p = parameter vector [p_low,p_high]
%                 p_low =  

% RE Strauss, 3/13/95


  % Find boundary intersections with PC2
  low_intsect = pt(lowpt,:) + p_low*(pt(lowpt,:)-pt(lowpt,:));
  high_intsect = pt(highpt,:) + p_high*(pt(highpt,:)-pt(highpt,:));

  % Given initial boundary points, find asym_area

  config = [low_intsect; high_intsect; pt];   % Concatenate
  newconfig = rotate(1,2,config);   % Rotate config so X-axis is the partition
  newpt = newconfig(3:npt+2,:);

  asym_area = asymarea([p_low,p_high],[lowpt,highpt],pt);


