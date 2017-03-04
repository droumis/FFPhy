






tetinfo = loaddatastruct(directory,fileprefix,'tetinfo',day);
dummy = evaluatefilter(tetinfo,tetfilter);
selected_tets = unique(dummy((dummy(:,1) == day),3))';


%load Ripples
animalinfo = animaldef(fileprefix);
ripout = kk_getconstimes(animalinfo{2},animalinfo{3},[day ep], 'ripplesdons',1,'consensus_numtets',3);
consvec_rip2 = ripout{day}{ep}.cons;
consvectimes_rip2 = ripout{day}{ep}.time;
[periodtimes_rip2 periodripsinds] = dr_vec2list(consvec_rip2,consvectimes_rip2);


%%
%load the spiking rasters
starttime_epenc = postimevec(1);
endtime_epenc = postimevec(end);
%     spikethresh = [0];
clear spktimes_e
for tet = selected_tets;
    load(sprintf('%sparam_nt%d.mat', paramsdir, tet)); %called filedata
    inds_epenc =  (  filedata.params(:,1)  >=  starttime_epenc  )  &  ( filedata.params(:,1) <= endtime_epenc );
    %       inds_thresh = any(filedata.params(:,2:5) > spikethresh, 2) ;
    inds_e = inds_epenc; % & inds_thresh;
    spktimes_e{tet} = filedata.params(inds_e,1)';
end

% ORRRRR ??
spikes = loaddatastruct(adir, aname, 'spikes');

%%
%plot spiking rasters
axes(ha(2))
      row = 0;
      for tet = selected_tets
          row = row + 1;
          usespikes = spktimes_e{tet}(spktimes_e{tet} > posteriortimevec(ripstartind(i))-rippad/round(1/dt) & spktimes_e{tet} < posteriortimevec(ripstartind(i))+rippad/round(1/dt));
          plot([usespikes; usespikes], repmat([row; row-1], size(usespikes)), 'Color', CM(tet,:), 'lineWidth', .1); hold on; %[.2 .5 .6][0 .4 .6]
      end
      xlim([posteriortimevec(ripstartind(i)-rippad) posteriortimevec(ripstartind(i)+rippad)])
      %       pause
      yl = ylim;
      patch([posteriortimevec(ripstartind(i)) posteriortimevec(ripendind(i)) posteriortimevec(ripendind(i)) posteriortimevec(ripstartind(i))], [yl(2) yl(2) yl(1) yl(1)], [1 1 1 1], 'edgecolor','none','FaceColor','k', 'FaceAlpha',.1);
%       set(gca, 'YTick', [1:yl(2)])
      set(gca, 'YTick', [0])
      set(gca, 'YTicklabel', ' ')
      ylabel('Tetrode MU', 'fontweight','bold','fontsize',9);



linpos = loaddatastruct(directory,fileprefix,'linpos',day);

% ORRRRRR
% look at Jai's spike raster ploter... /home/jai/Src/Functions/plotrasters.m
