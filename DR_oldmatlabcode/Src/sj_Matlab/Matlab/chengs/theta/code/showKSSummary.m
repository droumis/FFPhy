function showKSSummary(file)

if nargin < 1; file= 'KStest_summary'; end

load(file)

nana= size(viol,2)/2;
n= size(viol,1);
nviol= sum(viol);
%fprintf(1, '%d (%d) out of %d (%.1f%%/ %.1f%%) within 95%% (99%%)-CI.\n',...
%adapt_dir
%n-nviol, 
100*(n-nviol)/n

vp= zeros(size(viol));
vp(:,1:nana)= viol(:,1:2:2*nana);
vp(:,nana+1:end)= viol(:,2:2:2*nana);

plotsym={'sb', 'sk', 'sr', 'sy', 'sc'};
%cmap= colormap;
%colormap(gray(2));
figure
%imagesc(vind, viol');
imagesc(vp');
axis xy
%colormap(cmap);

print -depsc2 KSsummary
