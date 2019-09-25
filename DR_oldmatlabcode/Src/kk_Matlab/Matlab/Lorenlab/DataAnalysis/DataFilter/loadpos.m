
% load position data
pos_files = dir('*.p');
if isempty(pos_files)
  error('no p files found');
end

postimestamps = [];
pos = [];
smoothpos = [];
for i = 1:length(pos_files)
  [ttmp, ptmp] = getposinfo(pos_files(i).name);
  postimestamps = [postimestamps; ttmp];

  missing = find(ptmp(:,1)<10);
  good = find(ptmp(:,1)>=10);
  fixed = interp1(good,double(ptmp(good,:)),missing,'linear','extrap');
  ptmp2 = zeros(size(ptmp));
  ptmp2(good,:) = double(ptmp(good,:));
  ptmp2(missing,:) = fixed;

  pos = [pos; ptmp2];

  smooth_ptmp = filtfilt(fir1(60,0.005),1,ptmp2);
  smoothpos = [smoothpos; smooth_ptmp];
end

% clear ttmp, ptmp, ptmp2, smooth_ptmp, good, fixed, missing;

% exclude pulses outside of position record
outside = (pulse.timestamp(:,1) > double(postimestamps(end))) | ...
  (pulse.timestamp(:,1) < double(postimestamps(1)));


% mark location of novel track
plot(pos(:,1),pos(:,2),'.');
fprintf('Please mark novel portion of track');
key = waitforbuttonpress;
point1 = get(gca,'CurrentPoint');    % button down detected
rbbox; % returns in figure units
point2 = get(gca,'CurrentPoint');    % button up detected

novel_box(1) = min(point1(1,1),point2(1,1));
novel_box(2) = min(point1(1,2),point2(1,2));
novel_box(3) = max(point1(1,1),point2(1,1));
novel_box(4) = max(point1(1,2),point2(1,2));

novelidx = find( ...
  (pos(:,1) >= novel_box(1)) & ...
  (pos(:,1) < novel_box(3)) & ...
  (pos(:,2) >= novel_box(2)) & ...
  (pos(:,2) < novel_box(4)));

hold on;
plot(pos(novelidx,1),pos(novelidx,2),'r.');
fprintf('\n');

indexes = lookup(pulse.timestamp(:,1),double(postimestamps));

pulse.loc = nan(length(pulse.timestamp),size(pos,2));
pulse.loc(~outside,:) = pos(indexes(~outside),:);

novelpulses = ismember(indexes,novelidx);
novelpulses = find(novelpulses & ~outside);


% mark location of familiar track
fprintf('Please mark familiar portion of track');
key = waitforbuttonpress;
point1 = get(gca,'CurrentPoint');    % button down detected
rbbox; % returns in figure units
point2 = get(gca,'CurrentPoint');    % button up detected

fam_box(1) = min(point1(1,1),point2(1,1));
fam_box(2) = min(point1(1,2),point2(1,2));
fam_box(3) = max(point1(1,1),point2(1,1));
fam_box(4) = max(point1(1,2),point2(1,2));

familiaridx = find( ...
  (pos(:,1) >= fam_box(1)) & ...
  (pos(:,1) < fam_box(3)) & ...
  (pos(:,2) >= fam_box(2)) & ...
  (pos(:,2) < fam_box(4)));

hold on;
plot(pos(familiaridx,1),pos(familiaridx,2),'g.');
fprintf('\n');

fampulses = ismember(indexes,familiaridx);
fampulses = find(fampulses & ~outside);

clear point1 point2 key;

save noveltrials1 novelpulses fampulses novel_box fam_box;
