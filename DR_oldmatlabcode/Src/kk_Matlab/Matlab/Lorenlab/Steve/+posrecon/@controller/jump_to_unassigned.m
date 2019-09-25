function jump_to_unassigned(self)
% jump to the highest-priority gap in coverage

  % identify frames within the target range that already have assignments
  assigned = self.model.assigned(self.start_index:self.end_index);
  if (nnz(assigned) == self.number_of_frames)
    disp('all video frames in the target range have assignments!');
    return;
  end

  % check the start and end
  if ~assigned(1)
    self.current_index = self.start_index;
    return;
  elseif ~assigned(end)
    self.current_index = self.end_index;
    return;
  end

  % if the start and end have assignments, then find the biggest gap in
  % coverage and go to the (approximate) middle of it
  assigned_idx = find(assigned);
  index_diffs = diff(assigned_idx);
  timestamp_diffs = double(diff( ...
      self.data.timestamps(self.start_index - 1 + find(assigned))));
  gap_locations = find(index_diffs > 1);
  gap_sizes = timestamp_diffs(index_diffs > 1);
  [max_gap, i] = max(gap_sizes);
  self.current_index = self.start_index - 1 + ...
      round(mean(assigned_idx(gap_locations(i) + [0 1])));

end

