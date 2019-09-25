
% Version where ripplmod was saved individually for all cells. 
% See updated versons for filter framework

% USe the daved data files - HP_thetamod and HP_ripplemod to plot
% correlations between theta modulation and ripple modulation'

savefig1=0;
savedir = '/data25/sjadhav/HPExpt/ProcessedData/';
thetafile = [savedir 'HP_thetamod_PFC_gather']; area = 'PFC'; clr = 'b'; % PFC
ripplefile = [savedir 'HP_ripplemod_PFC'];

load(ripplefile); % load allripplemod and allripplemod_idx. Has 61 cells as of now
load(thetafile,'allthetamod','allthetamod_idx'); % Has 72 cells as of now

% Match idxs as in xcorrmesaures2

cntcells=0; cnt_mismatch=0;

for i=1:length(allripplemod)
    
    rippleidx = allripplemod_idx(i,:);
    match = rowfind(rippleidx, allthetamod_idx);
    
    if match~=0,
       cntcells = cntcells+1;
       allmod(cntcells).idx = rippleidx;
       % Theta
       allmod(cntcells).sph = allthetamod(match).sph;
       allmod(cntcells).Nspk = allthetamod(match).Nspk;
       allmod(cntcells).kappa = allthetamod(match).kappa;
       allmod(cntcells).modln = allthetamod(match).modln;
       allmod(cntcells).meanphase = allthetamod(match).meanphase;
       allmod(cntcells).prayl = allthetamod(match).prayl;
       % Ripple
       allmod(cntcells).sig_shuf = allripplemod(i).sig_shuf; % Use this to determine significance
       allmod(cntcells).pshuf = allripplemod(i).pshuf;
       allmod(cntcells).D = allripplemod(i).D; % Distance metric
       allmod(cntcells).ripmodln_shuf = allripplemod(i).modln_shuf; % %value of significance
       allmod(cntcells).sig_ttest = allripplemod(i).sig_ttest;
       allmod(cntcells).ripmodln = allripplemod(i).modln; % % change above baseline 
    end

end

% ------------------
% Population Figures
% ------------------
forppr = 0; 
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/HPExpt/Figures/RippleMod/';
summdir = figdir;
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);

if forppr==1
    set(0,'defaultaxesfontsize',16);
    tfont = 18; % title font
    xfont = 16;
    yfont = 16;
else
    set(0,'defaultaxesfontsize',24);
    tfont = 28;
    xfont = 20;
    yfont = 20;
end

% Get data
for i=1:length(allmod)
    % Theta
    allkappas(i) = allmod(i).kappa;
    allmodln(i) = allmod(i).modln;
    allmeanphase(i) = allmod(i).meanphase;
    allprayl(i) = allmod(i).prayl;       % 43/72. ~60% But only 36 have ripmodln defined
    % Ripple
    allripmodln_shuf(i) = allmod(i).ripmodln_shuf;
    allripmodln(i) = allmod(i).ripmodln;
    allD(i) = allmod(i).D;
    allpshuf(i) = allmod(i).pshuf;       % 25/61 are ripple modulated. 41%
    allsigshuf(i) = allmod(i).sig_shuf;
    allsigttest(i) = allmod(i).sig_ttest;
end

sigtheta = find(allprayl<0.05);
sigrip = find(allpshuf<0.05);
allsig = union(sigrip, sigtheta); 

 
 
 
% Theta Modln vs. D
% -------------------------------------
 figure; hold on;
 if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end

 plot(allmodln, allD, 'k.','MarkerSize',18);
 plot(allmodln(sigtheta), allD(sigtheta), 'r.','MarkerSize',18);
 plot(allmodln(sigrip), allD(sigrip), 'co','MarkerSize',12,'LineWidth',2);
 xlabel(['Theta Modln'],'FontSize',xfont,'Fontweight','normal');
 ylabel(['D metric'],'FontSize',yfont,'Fontweight','normal');
 
 [r1,p1] = corrcoef(allmodln,allD)   % p=0.03, r=0.28
 [r1rt,p1rt] = corrcoef(allmodln(allsig),allD(allsig)) % p=0.0056, r=0.28
 [r1t,p1t] = corrcoef(allmodln(sigtheta),allD(sigtheta)) % p=0.016, r=0.40
 [r1r,p1r] = corrcoef(allmodln(sigrip),allD(sigrip)) % p=0.038, r=0.42
 
 title(sprintf('r=%g, p =%g',roundn(r1(1,2),-2),roundn(p1(1,2),-3)));
 
 if savefig1==1,
     figfile = [figdir,area,'_ThetaModlnVsDmetric'];
     print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
 end
 
 
 
 
 % Kappa vs. D
% -------------------------------------
 figure; hold on;
 if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end

 plot(allkappas, allD, 'k.','MarkerSize',18);
 plot(allkappas(sigtheta), allD(sigtheta), 'r.','MarkerSize',18);
 plot(allkappas(sigrip), allD(sigrip), 'co','MarkerSize',12,'LineWidth',2);
 xlabel(['Kappas'],'FontSize',xfont,'Fontweight','normal');
 ylabel(['D metric'],'FontSize',yfont,'Fontweight','normal');
 
 [r2,p2] = corrcoef(allkappas,allD)   % p=0.0083, r=0.33
 [r2rt,p2rt] = corrcoef(allkappas(allsig),allD(allsig)) %p=0.06, r=0.27
 [r2t,p2t] = corrcoef(allkappas(sigtheta),allD(sigtheta)) %p=0.04, r=0.34
 [r2r,p2r] = corrcoef(allkappas(sigrip),allD(sigrip)) %p=0.01, r=0.50
 
 title(sprintf('r=%g, p =%g',roundn(r2(1,2),-2),roundn(p2(1,2),-3)));

  
 if savefig1==1,
     figfile = [figdir,area,'_KappasVsDmetric'];
     print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
 end
 
 

 
 
 
 
%   % Kappas vs. Ripple Shuffle Modulation
% % -------------------------------------
%  figure; hold on;
%  if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end
%  %plot(allkappas, allripmodln_shuf, 'k.','MarkerSize',18);
%  plot(allkappas(sigtheta), allripmodln_shuf(sigtheta), 'r.','MarkerSize',18);
%  plot(allkappas(sigrip), allripmodln_shuf(sigrip), 'co','MarkerSize',12);
%  xlabel(['Kappas'],'FontSize',xfont,'Fontweight','normal');
%  ylabel(['Ripple Shuffle %tile'],'FontSize',yfont,'Fontweight','normal');
%  [r3,p3] = corrcoef(allkappas,allripmodln_shuf)
%  [r3rt,p3rt] = corrcoef(allkappas(allsig),allripmodln_shuf(allsig))
  
  
%  % Theta Modln vs. Ripple Shuffle Modulation
% % -------------------------------------
%  figure; hold on;
%  if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end
% 
%  plot(allmodln, allripmodln_shuf, 'k.','MarkerSize',18);
%  plot(allmodln(sigtheta), allripmodln_shuf(sigtheta), 'r.','MarkerSize',18);
%  plot(allmodln(sigrip), allripmodln_shuf(sigrip), 'co','MarkerSize',12);
%  xlabel(['Theta Modln'],'FontSize',xfont,'Fontweight','normal');
%  ylabel(['Ripple Modln %tile'],'FontSize',yfont,'Fontweight','normal');
%  [r4,p4] = corrcoef(allmodln,allripmodln_shuf)
%  [r4rt,p4rt] = corrcoef(allmodln(allsig),allripmodln_shuf(allsig))
 
 
 
% Kappas vs. Ripple Modulation
% -------------------------------------
% figure; hold on;
% if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end
% 
% %plot(allkappas, allripmodln, 'k.','MarkerSize',18);
% 
% plot(allkappas(sigtheta), allripmodln(sigtheta), 'r.','MarkerSize',18);
% plot(allkappas(sigrip), allripmodln(sigrip), 'co','MarkerSize',12,'LineWidth',2);
% xlabel(['Kappas'],'FontSize',xfont,'Fontweight','normal');
% ylabel(['Ripple Modln'],'FontSize',yfont,'Fontweight','normal');
% 
% if savefig1==1,
%     figfile = [figdir,area,'_KappaVsRipModln'];
%     print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
% end
% 
% [r1,p1] = corrcoef(allkappas,allripmodln)
% [r1rt,p1rt] = corrcoef(allkappas(allsig),allripmodln(allsig))
% [r1t,p1t] = corrcoef(allkappas(sigtheta),allripmodln(sigtheta))
% [r1r,p1r] = corrcoef(allkappas(sigrip),allripmodln(sigrip)) % 0.08


% Theta Modln vs. Ripple Modulation
% -------------------------------------
% figure; hold on;
% if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end
% 
% %plot(allmodln, allripmodln, 'k.','MarkerSize',18);
% plot(allmodln(sigtheta), allripmodln(sigtheta), 'r.','MarkerSize',18);
% plot(allmodln(sigrip), allripmodln(sigrip), 'co','MarkerSize',12,'LineWidth',2);
% xlabel(['Theta Modln'],'FontSize',xfont,'Fontweight','normal');
% ylabel(['Ripple Modln'],'FontSize',yfont,'Fontweight','normal');
% 
% if savefig1==1,
%     figfile = [figdir,area,'_ThetaModlnVsRipModln'];
%     print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
% end
% 
% [r2,p2] = corrcoef(allmodln,allripmodln)   % 0.08
% [r2rt,p2rt] = corrcoef(allmodln(allsig),allripmodln(allsig))
% [r2t,p2t] = corrcoef(allmodln(sigtheta),allripmodln(sigtheta))
% [r2r,p2r] = corrcoef(allmodln(sigrip),allripmodln(sigrip))

 

 
 
 
 








