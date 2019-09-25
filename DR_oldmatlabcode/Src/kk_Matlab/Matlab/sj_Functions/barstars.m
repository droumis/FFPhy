function hh = barstars(cmp, ax, varargin)

mapping = 1:max(cmp(:,2));
[otherArgs] = procOptions(varargin);

xticks = get(ax,'xtick');
N = length(xticks);

if nchoosek(N,2) ~= size(cmp,1)
  error('Mismatch between cmp and x-axis');
end

level = diff(cmp(:,1:2),1,2);
[y,levsort] = sort(level);
cmp = cmp(levsort,:);

h = text(0,0,'*');

starh = get(h,'extent');
delete(h);
starh = starh(4);

tickh = get(ax,'ticklength');
tickh = tickh(1);

ylim = get(ax,'ylim');
basey = ylim(2);

lineheight = 2*tickh + starh;

dx = (xticks(2)-xticks(1))/20;

maxy = ylim(2);
lastline = basey;
past_lines = [];
for i = 1:size(cmp,1)
  if isnan(cmp(i,3))
    continue;
  end
  if cmp(i,3) == 1
      x1 = xticks(round(mapping(cmp(i,1)))) + dx;
      x2 = xticks(round(mapping(cmp(i,2)))) - dx;
      ht = 0;
      hh{i} = sigstar(x1,x2,lastline + ht*lineheight,tickh,sprintf('n.s.'));
  else
        x1 = xticks(round(mapping(cmp(i,1)))) + dx;
        x2 = xticks(round(mapping(cmp(i,2)))) - dx;
        if isempty(past_lines)
            ht = 0;
        else
            overlaps = findoverlaps([x1 x2],past_lines(:,1:2));
            ht = setdiff(past_lines(:,3),past_lines(overlaps==1,3));
            if isempty(ht)
                ht = max(past_lines(:,3)) + 1;
            else
                ht = ht(1);
            end
        end

        if (cmp(i,3) >= 0.01 && cmp(i,3) < 1)
            hh{i} = sigstar(x1,x2,lastline + ht*lineheight,tickh,sprintf('p = %3.2f',cmp(i,3)));
        else
            hh{i} = sigstar(x1,x2,lastline + ht*lineheight,tickh,sprintf('p = %1.0e',cmp(i,3)));
        end
        
  end
  
  maxy = max(maxy, basey + (ht+1)*lineheight);
  past_lines = [past_lines; x1 x2 ht];
end

set(ax,'ylim',[ylim(1) maxy]);

function ovlp = findoverlaps(xx, allxx)

if size(allxx,1) == 0
  ovlp = [];
  return;
end

for i = 1:size(allxx,1)
  if max(xx) < min(allxx(i,:)) | min(xx) > max(allxx(i,:))
    ovlp(i) = 0;
  else
    ovlp(i) = 1;
  end
end
