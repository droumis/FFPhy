function f = combinefilterepochs(f,varargin)
% function fout = combinefilterepochs(f,how,how,...)
% 
% Go through the epochs of f and combine those which share some
% commonality, specified by "how".
% 
% possible values of how: 
%    'days' - combine epochs from the same day
%    'empty' - get rid of empty epochs (leaves at least one) 

if isempty(f)
  return;
end

if nargin < 2
  warning('combinefilterepochs: did not specify how');
  return;
end

for an = 1:length(f)
  for h = 1:length(varargin)
    if strcmp(varargin{h},'days') % also deletes empties
      % combine filter epochs by day
      % first assemble all epochs
      allepochs = cat(1,f(an).epochs{:});
      allepochs = sortrows(allepochs);
      if isempty(allepochs)
        continue;
      end
      uniquedays = unique(allepochs(:,1));
      clear newepochs;
      for d = 1:length(uniquedays)
        newepochs{d} = allepochs(find(allepochs(:,1)==uniquedays(d)),:);
      end
      f(an).epochs = newepochs;
      % reset data and excludetime fields
      f(an).data = [];
      f(an).excludetime = [];
      for g = 1:length(f(an).epochs)
        if ~isempty(f(an).epochs{g})
          f(an).data{g} = [];
          f(an).excludetime{g} = [];
          for e = 1:size(f(an).epochs{g},1)
            f(an).data{g}{e} = [];
            f(an).excludetime{g}{e} = [];
          end
        else
          f(an).data{g} = [];
          f(an).excludetime{g} = [];
        end
      end
    elseif strcmp(varargin{h},'single') | strcmp(varargin{h},'split')
      L = length(f(an).epochs);
      for g = 1:L
        if isempty(f(an).epochs{g})
          continue;
        end
        % -------------------------------
        % here's the meat:
        if size(f(an).epochs{g},1) > 1

          % ee = f(an).epochs{g}(2:end,:);
          % f(an).epochs{g}(2:end,:) = [];
          % for i = 1:size(ee,1)
            % f(an).epochs{end+1} = ee(i,:);
          % end

          i = size(f(an).epochs{g},1)-1;
          while i > 0
            if ... % epochs are adjacent
              (f(an).epochs{g}(i+1,1) == f(an).epochs{g}(i,1)) & ...
                (abs(f(an).epochs{g}(i+1,2) - f(an).epochs{g}(i,2)) == 1)
              i = i - 1;
            else
              f(an).epochs{end+1} = f(an).epochs{g}(i+1:end,:);
              f(an).epochs{g}(i+1:end,:) = [];
              i = i - 1;
            end
          end

        end
        % -------------------------------
      end

      epochs = zeros(length(f(an).epochs),2);
      for g=1:length(f(an).epochs)
        epochs(g,:) = f(an).epochs{g}(1,:);
      end
      [epochs,inds] = sortrows(epochs);
      clear epochs;
      for g=1:length(f(an).epochs)
        epochs{g} = f(an).epochs{inds(g)};
      end
      f(an).epochs = epochs;
      
      % reset data and excludetime fields
      for g = 1:length(f(an).epochs)
        if ~isempty(f(an).epochs{g})
          for e = 1:length(f(an).epochs{g})
            f(an).data{g}{e} = [];
            f(an).excludetime{g}{e} = [];
          end
        else
          f(an).data{g} = [];
          f(an).excludetime{g} = [];
        end
      end
    elseif strcmp(varargin{h},'forcesingle')
      L = length(f(an).epochs);
      for g = 1:L
        if isempty(f(an).epochs{g})
          continue;
        end
        if size(f(an).epochs{g},1) > 1

          ee = f(an).epochs{g}(2:end,:);
          f(an).epochs{g}(2:end,:) = [];
          for i = 1:size(ee,1)
            f(an).epochs{end+1} = ee(i,:);
          end
        end
      end
      
      % reset data and excludetime fields
      for g = 1:length(f(an).epochs)
        if ~isempty(f(an).epochs{g})
          for e = 1:length(f(an).epochs{g})
            f(an).data{g}{e} = [];
            f(an).excludetime{g}{e} = [];
          end
        else
          f(an).data{g} = [];
          f(an).excludetime{g} = [];
        end
      end

    elseif strcmp(varargin{h},'empty')
      e = logical(zeros(length(f(an).epochs),1));
      for g = 1:length(f(an).epochs)
        if isempty(f(an).epochs{g})
          e(g) = 1;
        end
      end
      f(an).epochs(e) = [];
      f(an).data(e) = [];
      f(an).excludetime(e) = [];

    elseif strcmp(varargin{h},'emptytimes')
      e = logical(ones(length(f(an).epochs),1));
      for g = 1:length(f(an).epochs)
        if size(f(an).epochs{g},1) > 1
          error('times option only works for filters with single epoch epochs');
        end
        if iscell(f(an).times{g}) % multiple dio filters
          for ff = 1:length(f(an).times{g}{1})
            if ~isempty(f(an).times{g}{1}{ff})
              e(g) = 0;
            end
          end
        else
          if ~isempty(f(an).times{g})
            e(g) = 0;
          end
        end
      end
      f(an).epochs(e) = [];
      f(an).data(e) = [];
      f(an).excludetime(e) = [];
      f(an).times(e) = [];
      f(an).pulses(e) = [];
    elseif strcmp(varargin{h},'singletimes')
      e = logical(ones(length(f(an).epochs),1));
      for g = 1:length(f(an).epochs)
        if size(f(an).epochs{g},1) > 1
          error('times option only works for filters with single epoch epochs');
        end
        if iscell(f(an).times{g}) % multiple dio filters
          for ff = 1:length(f(an).times{g}{1})
            if size(f(an).times{g}{1}{ff},1)>1
              e(g) = 0;
            end
          end
        else
          if size(f(an).times{g},1)>1
            e(g) = 0;
          end
        end
      end
      f(an).epochs(e) = [];
      f(an).data(e) = [];
      f(an).excludetime(e) = [];
      f(an).times(e) = [];
      f(an).pulses(e) = [];
    end
  end
end

