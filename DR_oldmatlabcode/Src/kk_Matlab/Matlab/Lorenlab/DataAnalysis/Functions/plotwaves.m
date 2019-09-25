function plotwaves(waves,ts,highlight)
% function plotwaves(waves,ts,highlight)

waves = double(waves);

subplot(2,2,1,'align')
seplot(squeeze(waves(:,1,:))','facealpha',1.0, 'vartype', 'std');
subplot(2,2,2,'align')
seplot(squeeze(waves(:,2,:))','facealpha',1.0, 'vartype', 'std');
subplot(2,2,3,'align')
seplot(squeeze(waves(:,3,:))','facealpha',1.0, 'vartype', 'std');
subplot(2,2,4,'align')
seplot(squeeze(waves(:,4,:))','facealpha',1.0, 'vartype', 'std');

aextremes = [inf,-inf];
for i = 1:4
   subplot(2,2,i);
   a = axis;
   aextremes(1) = min(aextremes(1),a(3));
   aextremes(2) = max(aextremes(2),a(4));
end

for i = 1:4
   subplot(2,2,i);
   set(gca,'color','none','tickdir','out')
   box off
   axis([-1 40 aextremes]);
   line([0 0],[a(3) a(3)+100], 'linewidth',4,'color','k'); 
   line([10 40],[a(3)+2 a(3)+2], 'linewidth',4,'color','k');
   axis off
end

if nargin == 3
   if size(ts,1) ~= size(waves,3)
      error('Timestamps and waves do not match in length.')
   end
   inds = find(ismember(ts,highlight));

   if isempty(inds)
      error('Highlight timestamps not found.');
   end

   mn = mean(waves(:,:,inds),3);
%    subplot(2,2,1,'align'); hold on
%    plot(mn(:,1),'r');
%    subplot(2,2,2,'align'); hold on
%    plot(mn(:,2),'r');
%    subplot(2,2,3,'align'); hold on
%    plot(mn(:,3),'r');
%    subplot(2,2,4,'align'); hold on
%    plot(mn(:,4),'r');
   subplot(2,2,1,'align'); hold on
   seplot(squeeze(waves(:,1,inds))','fcolor',{[0.75 0.25 0.25]}, 'facealpha',1.0, 'vartype', 'std');
   subplot(2,2,2,'align'); hold on
   seplot(squeeze(waves(:,2,inds))','fcolor',{[0.75 0.25 0.25]}, 'facealpha',1.0, 'vartype', 'std');
   subplot(2,2,3,'align'); hold on
   seplot(squeeze(waves(:,3,inds))','fcolor',{[0.75 0.25 0.25]}, 'facealpha',1.0, 'vartype', 'std');
   subplot(2,2,4,'align'); hold on
   seplot(squeeze(waves(:,4,inds))','fcolor',{[0.75 0.25 0.25]}, 'facealpha',1.0, 'vartype', 'std');
end
