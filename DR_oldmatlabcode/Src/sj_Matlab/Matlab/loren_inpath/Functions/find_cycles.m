function [maxima, down, minima, up] = find_cycles(t,x,period)

% dependencies: quantile, find_nearest_mex

  assert( isscalar(period) && isreal(period) && ...
      isfinite(period) && (period > 0) && ...
      isvector(t) && isvector(x) && (numel(x) == numel(t)) && ...
      isreal(t) && all(isfinite(t)) && all(diff(t) > 0) && ...
      (min(t) == t(1)) && (max(t) == t(end)) && ...
      isreal(x) && ~any(isnan(x)) && all(isfinite(x)) );

  % How much do we tolerate variations in the oscillatory period?
  tol = 0.2;

  % Identify all zero-crossings and local extrema. 
  features = struct( ...
      'n'     , [], ...
      't'     , [], ...
      'x'     , [], ...
      'type'  , [], ...
      'zc'    , [], ...
      'ex'    , [] );
  % This isn't exactly mathematically exact, but it's approximately adequate 
  % for a discretely-sampled signal and nicely handles cases when the signal 
  % (or its first-order difference) is exactly zero
  i_up = find(diff(x > 0) > 0);
  i_down = find(diff(x > 0) < 0);
  i_maxima = 1 + find(diff(diff(x) > 0) < 0);
  i_minima = 1 + find(diff(diff(x) > 0) > 0);
  % Zero-crossing times are linearly interpolated whereas local extrema simply 
  % take the sample time. This ensures unique monotonic ordering.
  features.t = [ ...
      (t(i_up).*x(i_up+1) - t(i_up+1) .* x(i_up)) ./ ...
      (x(i_up+1) - x(i_up)); ...
      (t(i_down) .* x(i_down+1) - t(i_down+1) .* x(i_down)) ./ ...
      (x(i_down+1) - x(i_down)); ...   
      t([ i_maxima; i_minima ]) ];
  % Values of the signal at zero-crossings is assumed to be zero, and at
  % local extrema we take the closest sample
  features.x = [ ...
      zeros(size(i_up)); zeros(size(i_down)); ...
      x([ i_maxima; i_minima ]) ];
  % 'type' field contains enum integers
  %   -2: local minimum
  %   -1: downward zero-crossing
  %   +1: upward zero-crossing
  %   +2: local maximum
  features.type = [ +ones(size(i_up)); -ones(size(i_down)); 
      +2*ones(size(i_maxima)); -2*ones(size(i_minima)) ];
  % Sort features in time
  [features.t, idx] = sort(features.t,1,'ascend');
  features.type = features.type(idx);
  features.x = features.x(idx);
  % Number of features
  features.n = numel(features.t);
  % Bit for whether the features is a zero-crossing
  features.zc = (abs(features.type) == 1);
  % Bit for whether the features is a local extremum
  features.ex = (abs(features.type) == 2);
  assert(all(xor(features.zc,features.ex)));
  % 'flag' field contains enumerated integers to indicate whether the feature is
  % to be trusted or rejected.
  %   -1: reject
  %    0: undetermined
  %   +1: trust
  features.flag = zeros(size(features.t));
  % Sanity checks:
  % All times must fall within data range
  assert(all(features.t >= t(1)) & all(features.t <= t(end)));
  % Features must not have coincident times
  assert(all(diff(features.t) > 0));
  % Consecutive features must not be of the same type
  assert(~any(diff(features.type) == 0));
  % Consecutive zero-crossings must be of opposite sign
  assert(~any(diff(features.type(features.zc)) == 0));
  % Consecutive extrema must be of opposite sign
  assert(~any(diff(features.type(features.ex)) == 0));
  for i = 1:features.n
    % For each zero-crossing, the closest preceding extremum (if it exists) must
    % be of the opposite sign and the closest following extremum (if it exists)
    % must be of the same sign
    if features.zc(i)
      last_ex = find(features.ex(1:i),1,'last');
      if ~isempty(last_ex)
        assert(features.type(last_ex) == -2*features.type(i));
      end
      next_ex = i - 1 + find(features.ex(i:end),1,'first');
      if ~isempty(next_ex)
        assert(features.type(next_ex) == +2*features.type(i));
      end
    % For each extremum that is immediately adjacent to a zero-crossing, the
    % sign of that zero-crossing must be consistent with the sign of the
    % extremum and their relative order in time
    else % if features.ex(i)
      if (i > 1) && features.zc(i-1)
        assert(features.type(i-1) == sign(features.type(i)));
      end
      if (i < features.n) && (features.zc(i+1))
        assert(features.type(i+1) == -sign(features.type(i)));
      end
    end
  end
  % Between each pair of consecutive zero-crossings, there must exist at least
  % one extremum of the appropriate sign, and the number of extrema with that
  % sign must always be one greater than the number of extrema with the opposite
  % sign
  idx = find(features.zc);
  for i = 1:(length(idx)-1)
    j = (idx(i)+1):(idx(i+1)-1);
    assert( ...
        ~isempty(j) && ...
        all(features.ex(j)) && ...
        (nnz(features.type(j) == ...
        features.type(idx(i)) - features.type(idx(i+1))) == ...
        1 + nnz(features.type(j) == ...
        features.type(idx(i+1)) - features.type(idx(i)))) ...
        );
  end
  % Successive zero-crossings must be separated by no longer than 2*period
  if any(diff(features.t(abs(features.type) == 1)) > 2*period)
    error(['Successive zero-crossings are separated by longer than twice ' ...
        'the expected period (%f); this can indicate a DC offset, an ' ...
        'incorrectly-filtered signal, or a wrong period argument']);
  end
  % Successive extrema must be separated by no longer than 2*period
  if any(diff(features.t(abs(features.type) == 2)) > 2*period)
    error(['Successive extrema are separated by longer than twice ' ...
        'the expected period (%f); this can indicate a DC offset, an ' ...
        'incorrectly-filtered signal, or a wrong period argument']);
  end

  % We start with the zero-crossings, because these tend to be less numerous and
  % more reliable than the local extrema. (Most of the local extrema correspond
  % to spurious wiggles in the waveform.)
  for i = 2:(features.n - 1)
    if features.zc(i)
    % We flag a zero-crossing and its preceding and following neighbors as
    % trusted if it satisfies the following criteria:
    % (1) It is preceded by at least two zero-crossings and followed by at
    %     least two zero-crossings.
    % (2) The next zero-crossing of the same sign is period +/- tolerance in the
    %     future, and the last zero-crossing of the same sign is period
    %     +/-tolerance in the past.
    % (3) The immediately following zero-crossing of opposite sign is period
    %     +/-tolerance after the immediately preceding zero-crossing of opposite
    %     sign.
    % In other words, we flag a zero-crossing as trusted only if the last/next
    % two zero-crossings correspond exactly to the last/next two consecutive
    % half-cycles of oscillation.
      last_zc = find(features.zc(1:(i-1)),2,'last');
      next_zc = i + find(features.zc((i+1):end),2,'first');
      if (numel(last_zc) == 2) && (numel(next_zc) == 2) && ...
          (features.type(last_zc(1)) == features.type(i)) && ... 
          (features.type(last_zc(2)) == -features.type(i)) && ... 
          (features.type(next_zc(1)) == -features.type(i)) && ...
          (features.type(next_zc(2)) == features.type(i)) && ...
          all(abs([features.t(next_zc(2)) - features.t(i), ...
          features.t(next_zc(1)) - features.t(last_zc(2)), ...
          features.t(i) - features.t(last_zc(1))] - period)/period < tol)
        features.flag([last_zc(2), i, next_zc(1)]) = +1;
      end
    end
  end
  if ~any(features.flag == +1);
    error('Could not find any trustworthy oscillatory cycles');
  end
  % TODO: Improve estimate of period using these trusted zero-crossings

  % Starting from the initial trusted zero-crossings, flag the remaining
  % zero-crossings as trusted or rejected
  while any((features.flag == 0) & features.zc)
    % Forward
    for i = 1:features.n
      if (features.flag(i) ~= 0) || ~features.zc(i)
        continue;
      end
      % Find the last two trusted zero-crossings. These do not necessarily occur
      % in adjacent oscillatory cycles!
      last_zc = find(features.zc(1:i) & (features.flag(1:i) == +1),2,'last');
      % Check whether these two trusted zero-crossings have any intervening
      % undecided zero-crossings between them. If they do, then they are not
      % guaranteed to exactly bound one half-cycle, and we should not rely on
      % them.
      if (numel(last_zc) < 2) || any(features.zc(last_zc(1):last_zc(2)) & ...
          (features.flag(last_zc(1):last_zc(2)) == 0))
        continue;
      end
      % Otherwise, confirm that the signs of the zero-crossings alternate as
      % expected
      assert(features.type(last_zc(1)) == -features.type(last_zc(2)));
      % Find next trusted zero-crossing, if it exists
      next_zc = i - 1 + find( ...
          features.zc(i:end) & (features.flag(i:end) == +1),1,'first');
      % Find all zero-crossings that come after last_zc(2) and before next_zc
      % (if it exists). We expect to capture at least one half-cycle in this
      % time interval, unless we are at the end of the signal or close to
      % next_zc.
      if isempty(next_zc)
        candidates = last_zc(2) + find(features.zc((last_zc(2)+1):end));
      else
        candidates = last_zc(2) + find(features.zc((last_zc(2)+1):(next_zc-1)));
      end
      assert(nnz(i == candidates) == 1);
      % This algorithm works in such a way that all zero-crossings that lie
      % between two trusted zero-crossings have the same flag. Given that we
      % know that (features.flag(i) == 0), we expect that all other
      % zero-crossings on this non-trusted interval should be undecided.
      assert(all(features.flag(candidates) == 0));
      % We also expect that the first candidate has opposite sign as last_zc(2),
      % because of the way that zero-crossings alternate in sign
      assert(features.type(candidates(1)) == -features.type(last_zc(2)));
      if isscalar(candidates)
        % If this is the end of the signal and the candidate is not within
        % tolerance, then we infer that the rest of the cycle has been truncated
        % and mark it as rejected. Otherwise there are two possibilities, either
        % of which requires that we mark this zero-crossing as trusted:
        % (1) we are at the end of the signal, but this final zero-crossing is
        %     within tolerance; or 
        % (2) the next zero-crossing is next_zc. 
        if isempty(next_zc) && (abs(features.t(candidates) - ...
            - period - features.t(last_zc(1)))/period > tol)
          features.flag(candidates) = -1;
        else
          assert( isempty(next_zc) || ...
              (features.type(last_zc(end)) == features.type(next_zc)) );
          features.flag(candidates) = +1;
        end
      else
        % If there exist multiple candidates, try to find one that matches the
        % sign of last_zc(1) and is within tolerance of one period ahead of
        % last_zc(1), such that% there is another later candidate which has the
        % opposite sign and is within tolerance of one period ahead of
        % last_zc(2).
        prediction = features.t(last_zc) + period;
        good = ...
            (features.type(candidates) == features.type(last_zc(1))) & ...
            (abs(features.t(candidates) - prediction(1))/period < tol);
        for j = 1:numel(good)
          good(j) = good(j) && (j < numel(candidates)) && ...
              any(abs(features.t(candidates((j+1):end)) - ...
              prediction(2))/period < tol);
        end
        if any(good)
          % Dividing by good ensures that not-good candidates get a score of
          % infinity
          [dev, best] = min( ...
              abs(features.t(candidates) - prediction(1)) ./ good);
        else



          % TODO: use maxima/minima criteria. we want to find a zero-crossing
          % that separates the waveform between large lobes of opposite sign.



        end
        % Confirm that the match we found has the same sign as last_zc(2) and
        % opposite sign as last_zc(2)
        assert(features.type(candidates(best)) == features.type(last_zc(1)));
        features.flag(candidates(best)) = +1;
        if (best > 1)
          % Reject all candidates that lie between last_zc(2) and
          % candidates(best)
          features.flag(candidates(1:(best-1))) = -1;
        end
      end
    end

    % Backward
    for i = features.n:-1:1
      if (features.flag(i) ~= 0) || ~features.zc(i)
        continue;
      end
      % Find the next two trusted zero-crossings. These do not necessarily occur
      % in adjacent oscillatory cycles!
      next_zc = i - 1 + find(features.zc(i:end) & ...
          (features.flag(i:end) == +1),2,'first');
      % Check whether these two trusted zero-crossings have any intervening
      % undecided zero-crossings between them. If they do, then they are not
      % guaranteed to exactly bound one half-cycle, and we should not rely on
      % them.
      if (numel(next_zc) < 2) || any(features.zc(next_zc(1):next_zc(2)) & ...
          (features.flag(next_zc(1):next_zc(2)) == 0))
        continue;
      end
      % Otherwise, confirm that the signs of the zero-crossings alternate as
      % expected
      assert(features.type(next_zc(1)) == -features.type(next_zc(2)));
      % Find last trusted zero-crossing, if it exists
      last_zc = find( ...
          features.zc(1:i) & (features.flag(1:i) == +1),1,'last');
      % Find all zero-crossings that come after last_zc (if it exists) and
      % before next_zc(1). We expect to capture at least one half-cycle in this
      % time interval, unless we are at the end of the signal or close to
      % next_zc.
      if isempty(last_zc)
        candidates = find(features.zc(1:(next_zc(1)-1)));
      else
        candidates = last_zc + find(features.zc((last_zc+1):(next_zc(1)-1)));
      end
      assert(nnz(i == candidates) == 1);
      % This algorithm works in such a way that all zero-crossings that lie
      % between two trusted zero-crossings have the same flag. Given that we
      % know that (features.flag(i) == 0), we expect that all other
      % zero-crossings on this non-trusted interval should be undecided.
      assert(all(features.flag(candidates) == 0));
      % We also expect that the last candidate has opposite sign as next_zc(1),
      % because of the way that zero-crossings alternate in sign
      assert(features.type(candidates(end)) == -features.type(next_zc(1)));
      if isscalar(candidates)
        % If this is the start of the signal and the candidate is not within
        % tolerance, then we infer that the rest of the cycle has been truncated
        % and mark it as rejected. Otherwise there are two possibilities, either
        % of which requires that we mark this zero-crossing as trusted:
        % (1) we are at the start of the signal, but this initial zero-crossing
        %   is within tolerance; or 
        % (2) the last zero-crossing is last_zc.
        if isempty(last_zc) && (abs(features.t(next_zc(2)) - ...
            period - features.t(candidates))/period > tol)
          features.flag(candidates) = -1;
        else
          assert( isempty(last_zc) || ...
              (features.type(next_zc(1)) == features.type(last_zc)) );
          features.flag(candidates) = +1;
        end
      else
        % If there exist multiple candidates, try to find one that matches the
        % sign of next_zc(2) and is within tolerance of one period behind
        % next_zc(2), such that there is another earlier candidate which has the
        % opposite sign and is within tolerance of one period behind next_zc(1).
        prediction = features.t(next_zc) - period;
        good = ...
            (features.type(candidates) == features.type(next_zc(2))) & ...
            (abs(features.t(candidates) - prediction(2))/period < tol);
        for j = numel(good):-1:1
          good(j) = good(j) && (j > 1) && ...
              any(abs(features.t(candidates(1:(j-1))) - ...
              prediction(1))/period < tol);
        end
        if any(good)
          % Dividing by good ensures that not-good candidates get a score of
          % infinity
          [dev, best] = min( ...
              abs(features.t(candidates) - prediction(2)) ./ good);
        else



          % TODO: use maxima/minima criteria. we want to find a zero-crossing
          % that separates the waveform between large lobes of opposite sign.



        end
        % Confirm that the match we found has the same sign as next_zc(2) and
        % opposite sign as next_zc(1)
        assert(features.type(candidates(best)) == features.type(next_zc(2)));
        features.flag(candidates(best)) = +1;
        if (best < numel(candidates))
          % Reject all candidates that lie between candidates(best) and
          % next_zc(1)
          features.flag(candidates((best+1):end)) = -1;
        end
      end
    end
  end

  while any((features.flag == 0) & features.ex)
    for i = 1:features.n
      if (features.flag(i) ~= 0) || ~features.ex(i)
        continue;
      end
      % Given an undecided extremum, find the two trusted zero-crossings that
      % most tightly bound this local extremum. If there exist no non-trusted
      % zero-crossings on this interval, then then these two trusted
      % zero-crossings reliably bound a half-cycle.
      last_zc = find( ...
          features.zc(1:i) & (features.flag(1:i) == +1),1,'last');
      next_zc = i - 1 + find( ...
          features.zc(i:end) & (features.flag(i:end) == +1),1,'first');
      if isscalar(last_zc) && isscalar(next_zc) && ...
          ~any(features.zc(last_zc:next_zc) & ...
          (features.flag(last_zc:next_zc) == 0))
        % Confirm that these two bounding trusted zero-crossings are of
        % opposite signs.
        assert(features.type(last_zc) == -features.type(next_zc));
        % Find all extrema on this half-cycle (regardless of sign)
        candidates = last_zc - 1 + find(features.ex(last_zc:next_zc));
        assert(nnz(i == candidates) == 1);
        % All extrema within a given half-cycle are flagged at the same time
        % (all or none). Given that we know that (features.flag(i) == 0), we
        % expect that all other extrema within this half-cycle should
        % undecided.
        assert(all(features.flag(candidates) == 0));
        % If there is only one local extremum of the appropriate sign on this
        % half-cycle, then mark it as trusted. If there are multiple
        % candidates, then we compare to a trusted extremum of the same sign
        % in an adjacent cycle (if it exists).
        if isscalar(candidates)
          assert(features.type(candidates) == ...
              features.type(last_zc) - features.type(next_zc)); 
          features.flag(candidates) = +1;
        else
          last_ex = find(features.ex(1:last_zc) & ...
              (features.type(1:last_zc) == ...
              features.type(last_zc) - features.type(next_zc)) & ...
              (features.flag(1:last_zc) == +1),1,'last');
          next_ex = next_zc - 1 + find(features.ex(next_zc:end) & ...
              (features.type(next_zc:end) == ...
              features.type(last_zc) - features.type(next_zc)) & ...
              (features.flag(next_zc:end) == +1),1,'first');
          prediction = mean([ features.t(last_ex) + period, ...
              features.t(next_ex) - period ]);
          if isscalar(prediction) && ~isnan(prediction)
            % Pick the most extreme-valued candidate
            good = (features.type(candidates) == ...
                features.type(last_zc) - features.type(next_zc));
            assert(any(good));
            % We multiply by good before finding max to ensure that only the
            % good candidates are chosen
            [dev, best] = max(sign(features.type(candidates)) .* ...
                features.x(candidates) .* good);
            features.flag(candidates) = -1;
            features.flag(candidates(best)) = +1;
          end
        end
      elseif isscalar(last_zc) && isempty(next_zc) && ~any( ...
          features.zc(last_zc:end) & (features.flag(last_zc:end) == 0))
        % If the waveform is truncated at the end of the signal so that there
        % are no later undecided zero-crossings, then we look behind to the
        % trusted extremum of the same sign that precedes last_zc (if it
        % exists)
        candidates = last_zc - 1 + find(features.ex(last_zc:end));
        assert(nnz(i == candidates) == 1);
        assert(all(features.flag(candidates) == 0));
        % Find the trusted extremum of the same sign in the adjacent cycle
        last_ex = find(features.ex(1:last_zc) & ...
            (features.flag(1:last_zc) == +1) & ...
            (features.type(1:last_zc) == +2*features.type(last_zc)), ...
            1,'last');
        if ~isempty(last_ex) && ~any(features.ex(last_ex:last_zc) & ...
              (features.flag(last_ex:last_zc) == 0))
          prediction = features.t(last_ex) + period;
          % Pick the most extreme-valued candidate among those that are
          % within tolerance and have the correct sign. If none match, then
          % reject all; we presume that the real extremum is after t(end).
          good = ...
              (features.type(candidates) == +2*features.type(last_zc)) & ...
              (abs(features.t(candidates) - prediction)/period < tol);
          if any(good)
            features.flag(candidates(~good)) = -1;
            candidates = candidates(good);
            [dev, best] = max(sign(features.type(candidates) .* ...
                features.x(candidates)));
            features.flag(candidates) = -1;
            features.flag(candidates(best)) = +1;
          else
            features.flag(candidates) = -1;
          end
        end
      elseif isscalar(next_zc) && isempty(last_zc) && ~any( ...
          features.zc(1:next_zc) & (features.flag(1:next_zc) == 0))
        % If the waveform is truncated at the start of the signal so that
        % there are no earlier undecided zero-crossings, then we look ahead to
        % the trusted extremum of the same sign that follows next_zc (if it
        % exists)
        candidates = find(features.ex(1:next_zc));
        assert(nnz(i == candidates) == 1);
        assert(all(features.flag(candidates) == 0));
        % Find the trusted extremum of the same sign in the adjacent cycle
        next_ex = next_zc - 1 + find(features.ex(next_zc:end) & ...
            (features.flag(next_zc:end) == +1) & ...
            (features.type(next_zc:end) == -2*features.type(next_zc)), ...
            1,'first');
        if ~isempty(next_ex) && ~any(features.ex(next_zc:next_ex) & ...
            (features.flag(next_zc:next_ex) == 0))
          prediction = features.t(next_ex) - period;
          % Pick the most extreme-valued candidate among those that are
          % within tolerance and have the correct sign. If none match, then
          % reject all; we presume that the real extremum is before t(1).
          good = ...
              (features.type(candidates) == -2*features.type(next_zc)) & ...
              (abs(features.t(candidates) - prediction)/period < tol);
          if any(good)
            features.flag(candidates(~good)) = -1;
            candidates = candidates(good);
            [dev, best] = max(sign(features.type(candidates) .* ...
                features.x(candidates)));
            features.flag(candidates) = -1;
            features.flag(candidates(best)) = +1;
          else
            features.flag(candidates) = -1;
          end
        end
      end
    end
  end
    
    plot(t,x,'LineWidth',1,'Color','k');
    hold on;
    plot(features.t((features.flag == +1) & features.zc), ...
        features.x((features.flag == +1) & features.zc), ...
        'ro','MarkerSize',6);
    plot(features.t((features.flag == +1) & features.ex), ...
        features.x((features.flag == +1) & features.ex), ...
        'gx','MarkerSize',6);


  % Keep only the trusted features
  features.t = features.t((features.flag == +1)); 
  features.type = features.type((features.flag == +1)); 
  features = rmfield(features,{'flag','n','ex','zc'});

  % Extrapolate as necessary, so that each feature type is represented by a
  % trusted point before the start time and after the end time
  while any( ...
      (features.t(find(features.type == +2,1,'first')) > t(1)) || ...
      (features.t(find(features.type == -1,1,'first')) > t(1)) || ...
      (features.t(find(features.type == -2,1,'first')) > t(1)) || ...
      (features.t(find(features.type == +1,1,'first')) > t(1)) )
    switch features.type(1)
    case +2
      % Maximum is preceded by an upward zero-crossing
      features.t = [features.t(find(features.type == +1,1,'first')) - period; 
          features.t];
      features.type = [+1; features.type];
    case -1
      % Downward zero-crossing is preceded by a maximum
      features.t = [features.t(find(features.type == +2,1,'first')) - period; 
          features.t];
      features.type = [+2; features.type];
    case -2
      % Minimum is preceded by a downward zero-crossing
      features.t = [features.t(find(features.type == -1,1,'first')) - period; 
          features.t];
      features.type = [-1; features.type];
    case +1
      % Upward zero-crossing is preceded by a minimum
      features.t = [features.t(find(features.type == -2,1,'first')) - period; 
          features.t];
      features.type = [-2; features.type];
    end
  end
  while any( ...
      (features.t(find(features.type == +2,1,'last')) > t(end)) || ...
      (features.t(find(features.type == -1,1,'last')) > t(end)) || ...
      (features.t(find(features.type == -2,1,'last')) > t(end)) || ...
      (features.t(find(features.type == +1,1,'last')) > t(end)) )
    switch features.type(end)
    case +2
      % Maximum is followed by an downward zero-crossing
      features.t = [features.t; ...
          features.t(find(features.type == -1,1,'last')) + period];
      features.type = [features.type; -1];
    case -1
      % Downward zero-crossing is followed by a minimum
      features.t = [features.t; ...
          features.t(find(features.type == -2,1,'last')) + period];
      features.type = [features.type; -2];
    case -2
      % Minimum is followed by a upward zero-crossing
      features.t = [features.t; ...
          features.t(find(features.type == +1,1,'last')) + period];
      features.type = [features.type; +1];
    case +1
      % Upward zero-crossing is follows by a maximum
      features.t = [features.t; ...
          features.t(find(features.type == +2,1,'last')) + period]; 
      features.type = [features.type; +2];
    end
  end

  % Confirm that the trusted features occur in repeating interleaved order:
  %     ... maximum, down, minimum, up, ...
  switch features.type(1)
  case +2
    assert( all(features.type(1:4:end) == +2) && ...
        all(features.type(2:4:end) == -1) && ...
        all(features.type(3:4:end) == -2) && ...
        all(features.type(4:4:end) == +1) );
  case -1
    assert( all(features.type(1:4:end) == -1) && ...
        all(features.type(2:4:end) == -2) && ...
        all(features.type(3:4:end) == +1) && ...
        all(features.type(4:4:end) == +2) );
  case -2
    assert( all(features.type(1:4:end) == -2) && ...
        all(features.type(2:4:end) == +1) && ...
        all(features.type(3:4:end) == +2) && ...
        all(features.type(4:4:end) == -1) );
  case +1
    assert( all(features.type(1:4:end) == +1) && ...
        all(features.type(2:4:end) == +2) && ...
        all(features.type(3:4:end) == -1) && ...
        all(features.type(4:4:end) == -2) );
  end

  maxima = features.t(features.type == +2);
  down = features.t(features.type == -1);
  minima = features.t(features.type == -2);
  up = features.t(features.type == +1);
 
end % end function find_cycles


