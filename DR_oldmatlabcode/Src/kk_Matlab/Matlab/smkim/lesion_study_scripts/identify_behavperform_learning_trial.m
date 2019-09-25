
load('behavperform_2');
for s = 1:length(behavperform)

  % inbound
  behavperform(s).lt_in = NaN;
  behavperform(s).ld_in = NaN;
  % scan for stretches of above-chance trials
  abovecriterion_in = find(behavperform(s).inprobcorrect(:,2) > ...
      behavperform(s).criterion);
  splitidx = diff(find(diff([-Inf; ...
      abovecriterion_in - (1:numel(abovecriterion_in))'; Inf])));
  t_in = mat2cell(abovecriterion_in,splitidx,1);
  % iterate through these stretches and find the earliest one which entirely
  % encompasses two full days of trials
  daybounds = [ behavperform(s).dayintrials(1:2:end,1), ...
      behavperform(s).dayintrials(2:2:end,2) ];
  for i = 1:length(t_in)
    days_spanned = find(all(ismember(daybounds,t_in{i}),2));
    if (nnz(diff(find(days_spanned)) == 1) >= 1)
      behavperform(s).lt_in = t_in{i}(1);
      behavperform(s).ld_in = find( ...
          daybounds(:,1) <= behavperform(s).lt_in,1,'last');
      break;
    end
  end
  
  % outbound
  behavperform(s).lt_out = NaN;
  behavperform(s).ld_out = NaN;
  % scan for stretches of above-chance trials
  abovecriterion_out = find(behavperform(s).outprobcorrect(:,2) > ...
      behavperform(s).criterion);
  splitidx = diff(find(diff([-Inf; ...
      abovecriterion_out - (1:numel(abovecriterion_out))'; Inf])));
  t_out = mat2cell(abovecriterion_out,splitidx,1);
  % iterate through these stretches and find the earliest one which entirely
  % encompasses two full days of trials
  daybounds = [ behavperform(s).dayouttrials(1:2:end,1), ...
      behavperform(s).dayouttrials(2:2:end,2) ];
  for i = 1:length(t_out)
    days_spanned = find(all(ismember(daybounds,t_out{i}),2));
    if (nnz(diff(find(days_spanned)) == 1) >= 1)
      behavperform(s).lt_out = t_out{i}(1);
      behavperform(s).ld_out = find( ...
          daybounds(:,1) <= behavperform(s).lt_out,1,'last');
      break;
    end
  end
end


