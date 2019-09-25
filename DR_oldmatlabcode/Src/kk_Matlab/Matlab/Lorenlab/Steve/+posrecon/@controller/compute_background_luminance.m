function compute_background_luminance(self,n)
% take median of num_frames random images
  if (n > self.number_of_frames)
    error(['num_frames argument must be less than or equal to the number ' ...
        'of frames within the target time range']);
  end

  if 1
    i = randperm(self.number_of_frames);
    % careful: this can trigger an out-of-memory error
    lum = zeros([self.data.height self.data.width n],'uint8');
    for j = 1:n
      cdata = self.data.get_picture(self.start_index + i(j) - 1);
      lum(:,:,j) = cdata(:,:,1);
    end
  else
    keyboard
    [alli,allj] = ind2sub([self.data.height,self.data.width],[1:self.data.height*self.data.width]);
    % load run1assignments
    inds = find(~isnan(rawpos.xfront) & ~isnan(rawpos.xback));
    i = randperm(length(inds));
    lum = zeros([self.data.height self.data.width n],'single');
    for j = 1:n
      ii = inds(i(j));
      cdata = self.data.get_picture(self.start_index + ii - 1);
      lum(:,:,j) = cdata(:,:,1);

      fr = [round(rawpos.xfront(ii)) round(rawpos.yfront(ii))];
      bk = [round(rawpos.xback(ii)) round(rawpos.yback(ii))];
      for k = 1:length(alli)
        if (abs(alli(k)-fr(1)) < 10) & (abs(allj(k)-fr(2)) < 10)
          lum(allj(k),alli(k),j) = nan;
        end
        if (abs(alli(k)-bk(1)) < 10) & (abs(allj(k)-bk(2)) < 10)
          lum(allj(k),alli(k),j) = nan;
        end
      end
    end
  end

  % take the 1% quantile as the background
  self.background_luminance = quantile(lum,0.01,3);
  % self.background_luminance = quantile(lum,0.1,3);

  self.background_luminance = uint8(self.background_luminance);

end

