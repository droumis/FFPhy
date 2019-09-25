function compute_background_luminance(self,n)
% take median of num_frames random images
  if (n > self.number_of_frames)
    error(['num_frames argument must be less than or equal to the number ' ...
        'of frames within the target time range']);
  end

  i = randperm(self.number_of_frames);
  % careful: this can trigger an out-of-memory error
  lum = zeros([self.data.height self.data.width n],'uint8');
  for j = 1:n
    cdata = self.data.get_picture(self.start_index + i(j) - 1);
    lum(:,:,j) = cdata(:,:,1);
  end

  % take the 1% quantile as the background
  self.background_luminance = quantile(lum,0.01,3);

end

