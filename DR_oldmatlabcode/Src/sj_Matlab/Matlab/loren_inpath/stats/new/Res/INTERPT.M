% INTERPT: Plots a set of interpoint distances among landmarks, and prepares
%           files for plotting in SigmaPlot.  Puts landmark coordinates in
%           MATRIX1.OUT and interpoint line segments in MATRIX2.OUT.
%           Zeros in MATRIX2.OUT must be converted to missing-data values
%           for SigmaPlot.
%
%     Usage: interpt(landmarks, dists)
%
%        landmarks = [N x 2] matrix of X,Y coordinates of landmarks.
%        dists =     [P x 2] matrix of paired indices to landmarks,
%                       indicating interpoint segments (characters).
%

% RE Strauss, 9/20/97
%   8/20/99 - changed plot colors for Matlab v5.
%   1/4/00 -  changed usage of sqplot().

function interpt(landmarks,dists)
   [N,c] = size(landmarks);
   [P,c] = size(dists);

   L = landmarks;
   D = dists;

   min_crds = min(L);            % Make all coordinates positive
   for c=1:2
      if (min_crds(c) < 0)
         L(:,c) = L(:,c) - min_crds(c)+1;
      end;
   end;

   plot(L(:,1),L(:,2),'*r');     % Plot landmarks in red
   sqplot(L);                    % Square plot, equal-unit axes
   hold on;
   for i=1:P                     % Plot line segments in black
      plot(L(D(i,:),1),L(D(i,:),2),'-k');
   end;
   hold off;
   pause;

   % Output landmark coordinates for SigmaPlot.

   tofile1(L);

   % Output forewing line segments for SigmaPlot.  Pairs of coordinates for
   % each line segment are separated by rows of zeros, which must later be
   % changed to missing-data values.

   v=[];
   for i=1:P
      v = [v; L(D(i,1),:); L(D(i,2),:); 0 0];
   end;
   [V,c] = size(v);
   v(V,:) = [];                  % Remove final row of zeros

   tofile2(v);
   return;
