function value = auto_assign(self,index,mask)
% Try to assign the position for video frame index, focusing on pixels selected
% by the logical array mask, using assignments in neighboring frames

  if self.assigned(index)
    value = true;
    return;
  end

  % skip if the frame is too close to the start or end
  if (index > self.data.number_of_frames - self.halfwidth) || ...
      (index < self.halfwidth+1)
    value = false;
    return;
  end

  % look up neighboring frames
  predictors = (index-self.halfwidth):(index+self.halfwidth);
  % select only those that have valid assignments
  predictors = predictors(self.assigned(predictors));
  if isempty(predictors)
    value = false;
    return;
  end

  t_lag = double(self.data.timestamps(predictors)) - ...  
      double(self.data.timestamps(index));
  if (numel(predictors) >= 4)
    % use polyfit to fit local quadratic centered around the frame of interest
    % For velocity estimation, remember that t_lag is expressed in timestamp
    % units (1e4 units = 1 second)
    coeffs_xcenter = polyfit(t_lag, ...
        0.5*(self.xfront(predictors) + self.xback(predictors)),2);
    coeffs_ycenter = polyfit(t_lag, ...
        0.5*(self.yfront(predictors) + self.yback(predictors)),2);
    coeffs_z = polyfit(t_lag,self.z(predictors),2);
    xcenter_proposal = coeffs_xcenter(end);
    ycenter_proposal = coeffs_ycenter(end);
    z_proposal = coeffs_z(end);
    xfront_proposal = xcenter_proposal + 0.5*real(z_proposal);
    xback_proposal = xcenter_proposal - 0.5*real(z_proposal);
    yfront_proposal = ycenter_proposal + 0.5*imag(z_proposal);
    yback_proposal = ycenter_proposal - 0.5*imag(z_proposal);
    % if the rat is moving fast, make the radius of detection bigger
    gain = 1 + coeffs_xcenter(2).^2 + coeffs_ycenter(2).^2;
  else
    % otherwise just take the average of the predictors
    sum_t_lag = sum(t_lag);
    xfront_proposal = mean(self.xfront(predictors));
    yfront_proposal = mean(self.yfront(predictors));
    xback_proposal = mean(self.xback(predictors));
    yback_proposal = mean(self.yback(predictors));
    z_proposal = mean(self.z(predictors));
    gain = 1;
  end

  % avoid auto-assignining if the seeds are too close, because then the
  % auto-assignment will sometimes flip the front and back
  if (abs(z_proposal) < max(self.epsilon))
    value = false;
    return;
  end

  % skip clustering if we can enough neighboring frames to interpolate with
  % confidence
  if (numel(predictors) >= 2*self.halfwidth)
    self.assign_front(index,[xfront_proposal yfront_proposal]);
    self.assign_back(index,[xback_proposal yback_proposal]);
    value = true;
    return;
  end

  % find pixels that are within radius of either the front or back
  criterion = gain * (self.sum_epsilon).^2;
  proximity = ((self.xgrid - xfront_proposal).^2 + ...
      (self.ygrid - yfront_proposal).^2 < criterion) | ...
      ((self.xgrid - xback_proposal).^2 + ...
      (self.ygrid - yback_proposal).^2 < criterion);
  % flag the pixels within this zone which exceed luminance threshold
  cdata = self.data.get_picture(index);
  flagged = proximity & mask;
  if (nnz(flagged) < 2)
    value = false;
    return;
  end

  % perform heuristic, non-convergent k-means clustering on the flagged pixels
  [labels, centroids, sum_dist] = kmeans( ...
      [self.xgrid(flagged(:)) self.ygrid(flagged(:))],2, ...
      'Start',[xfront_proposal, yfront_proposal; ...
      xback_proposal, yback_proposal],'EmptyAction','drop','MaxIter',5);
  if (sum_dist(1)/nnz(labels==1) > self.epsilon(1)) || ...
      (sum_dist(2)/nnz(labels==2) > self.epsilon(2))
    % reject the clustering if some of the flagged pixels are too far away from
    % their putative centroids. if one of the centroids is NaN then this
    % condition will not triger
    value = false;
    return;
  elseif ~any(isnan(centroids(:))) && ...
      ((centroids(1,1) - centroids(2,1)).^2 + ...
      (centroids(1,2) - centroids(2,2)).^2 > max(self.epsilon).^2)
    % update model
    self.assign_front(index,centroids(1,:));
    self.assign_back(index,centroids(2,:));
    value = true;
    return;
  elseif ~isnan(prod(centroids(1,:))) && isnan(prod(centroids(2,:))) && ...
      (nnz(t_lag > 0) > 1) && (nnz(t_lag < 0) > 1) && ...
      (abs(z_proposal) > self.sum_epsilon)
    self.assign_front(index,centroids(1,:));
    self.assign_back(index,centroids(1,:) - [real(z_proposal) imag(z_proposal)]);
    value = true;
    return;
  elseif ~isnan(prod(centroids(2,:))) && isnan(prod(centroids(1,:))) && ...
      (nnz(t_lag > 0) > 1) && (nnz(t_lag < 0) > 1) && ...
      (abs(z_proposal) > self.sum_epsilon)
    self.assign_back(index,centroids(2,:));
    self.assign_front(index,centroids(2,:) + [real(z_proposal) imag(z_proposal)]);
    value = true;
    return;
  else
    value = false;
    return;
  end

  %pause(1);
end

