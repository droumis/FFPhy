

% 2D Place Field Maps with swr responsees

close all
runFilterFramework = 1;
saveFilterOutput = runFilterFramework;
loadFilterOutput = 0;
plotfigs = 1;
savefigs= 1;
pausefigs = 0;

%% ---------------- Data Filters --------------------------
Fp = load_filter_params('occnormfiring_openfield');
Fp.animals = {'JZ4'};
Fp.days = [1:14];
%% ---------------- Run FIlter ---------------------------------------------------
if runFilterFramework == 1
    F = createfilter('animal',Fp.animals,'days', Fp.days,'epochs', ...
        Fp.epochfilter, 'cells', Fp.cellfilter,'excludetime', Fp.timefilter, ...
        'iterator', Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes);
    F = runfilter(F);
end
%% ---------------- Paths ---------------------------------------------------
paths = make_paths(animaldef(lower('Demetris')), Fp.filtfunction, Fp.days, ...
    Fp.epochEnvironment);
%% ---------------- Save Filter Output ---------------------------------------------------
if saveFilterOutput == 1;
    save_filter_output(F, paths.filtOutputDirectory, paths.filename)
end
%% ---------------- Load Filter Output ---------------------------------------------------
if loadFilterOutput == 1;
    F = load_filter_output(paths.filtOutputDirectory, paths.filename);
end
%% ---------------- plot? --------------------------
if ~plotfigs
    return
end
%% ---------------- Plot --------------------------
Pp = load_plotting_params('occnormfiring');
set(0,'DefaultFigureWindowStyle','normal')
for ian = 1:length(Fp.animals); %per anim
    animal = Fp.animals{ian};
    animalinfo = animaldef(lower(animal));
    anID = animalinfo{1,3};
    FFanimdir = animalinfo{2};
    for i = 1:length(F(ian).output{1})
    % ---------- plot per cell --------------------------------
        if savefigs && ~pausefigs; 
            % invisible figs for faster saving
            ifig = figure('Visible','off','units','normalized','position', ...
                Pp.position);
        else
            ifig = figure('units','normalized','position',Pp.position);
        end
        d = num2cell(F(ian).output{1}(i).index);
        [day, epoch, ntrode, cluster] = d{:};
        
        % load pos data
        load(sprintf('%s%s%s%02d.mat',FFanimdir, anID, 'pos', day));
        load(sprintf('%s%s%s.mat',FFanimdir, anID, 'cellinfo'));
        Xstring = 'x-loess';
        Xcol = find(cell2mat(cellfun(@(x) strcmp(x,Xstring), ...
            strsplit(pos{day}{epoch}.fields, ' '), 'un', 0)));
        xpos = pos{day}{epoch}.data(:,Xcol);
        Ystring = 'y-loess';
        Ycol = find(cell2mat(cellfun(@(x) strcmp(x,Ystring), ...
            strsplit(pos{day}{epoch}.fields, ' '), 'un', 0)));
        ypos = pos{day}{epoch}.data(:,Ycol);


        %plot position, with raw spikes overlayed
        subplot(1,4,1)
        plot(xpos, ypos,'Color', [.9 .9 .9], 'LineWidth', 2); 
        hold on;
        xtic = F(ian).output{1}(i).xticks{1};
        ytic = F(ian).output{1}(i).yticks{1};
        spikesmat = F(ian).output{1}(i).spikes{1};
        [spy, spx] = find(rot90(spikesmat,2)>0);
        plot(spx+min(xtic), spy+min(ytic), 'r.', 'MarkerSize', .5);
        axis equal
        axis tight
        xlabel('cm','FontSize',10,'FontWeight','bold','FontName', 'Arial')
        ylabel('cm','FontSize',10,'FontWeight','bold','FontName', 'Arial')
        title('raw pos,spikes','FontSize',12,'FontWeight','bold','FontName', 'Arial')
 
        %plot occnormfiring
        subplot(1,4,2)
        occNormFR = F(ian).output{1}(i).smoothedspikerate{1};
        occNormFR(occNormFR== -1) = 0;
        occNormFR(isnan(occNormFR)) = 0;
        occNormFR = rot90(occNormFR',3);
        s = imagesc(xtic, ytic, occNormFR);
        cmap = colormap(Pp.usecolormap);
        colormap(cmap)
        axis equal
        axis tight
        xlabel('cm','FontSize',10,'FontWeight','bold','FontName', 'Arial')
        ylabel('cm','FontSize',10,'FontWeight','bold','FontName', 'Arial')
        title('occ norm FR','FontSize',12,'FontWeight','bold','FontName', 'Arial')
        
        subplot(1,4,3)
        autcorr2 = conv2(occNormFR, occNormFR , 'full');
        imagesc(xtic, ytic, autcorr2) %, 'CDataMapping', 'scaled'
        axis equal
        axis tight
        xlabel('cm','FontSize',10,'FontWeight','bold','FontName', 'Arial')
        ylabel('cm','FontSize',10,'FontWeight','bold','FontName', 'Arial')
        title('autocorr2D','FontSize',12,'FontWeight','bold','FontName', 'Arial')
        
        subplot(1,4,4)
        efeft = fft2(occNormFR);
        imagesc(xtic, ytic, real(efeft)) %, 'CDataMapping', 'scaled'
        axis equal
        axis tight
        xlabel('cm','FontSize',10,'FontWeight','bold','FontName', 'Arial')
        ylabel('cm','FontSize',10,'FontWeight','bold','FontName', 'Arial')
        title('FFT2D','FontSize',12,'FontWeight','bold','FontName', 'Arial')
        
 
%         s.FaceColor = 'interp';
%         s.FaceAlpha = 'interp';
%         caxis([0 1])
        
%         set(gca,'XColor', [.5 .5 .5], 'FontSize',6,'FontWeight','bold','FontName',...
%                             'Arial')
%         set(gca,'ytick',[], 'yticklabel',[])
%          title(sprintf('%s FR openfield'), 'Color', [.5 .5 .5], 'FontSize',12,...
%             'FontWeight','bold','FontName', 'Arial');
        
%     
%     %find the unique day/tet outputs
%     matInds = cell2mat({F(ian).output{1}.index}');
%     [daytetcells, ~, daytetInds2] = unique(matInds(:,[1 3 4]), 'rows', ...
%         'stable');
% 
%     for icell = 1:length(daytetcells); %for each allepoch-unique cell
%         cellID = daytetcells(icell, 3);
%         icellFoutInds = find(icell == daytetInds2);
%         epochInds = matInds(icellFoutInds,:);
%         numeps = size(epochInds,1);
%         d = num2cell(epochInds(1,:));
%         
% 
%     end
  
    sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
    sprtit = sprintf('%s - %d %d %d %d', animal, day, epoch, ntrode, cluster);
    iStitle = text(.5, .95, {sprtit}, 'Parent', sprtitleax, 'Units', ...
        'normalized');
    set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
        'horizontalAlignment', 'center');
    % ---- pause, save figs ----
    if pausefigs
        pause
    end
    if savefigs
        save_figure(paths.figdirectory, paths.filenamesave, sprtit)
    end
    close all    
    end
end

            
            
%%  
% for i = 1:size(F.output{1},2)
%     fprintf('%d %d %d %d \n',F.output{1}(i).index)
%     try
%         figure
% 
%         subplot(2,2,1)
%         xlabel('xpos')
%         ylabel('ypos')
%         title('smoothed occupancy + spikes pos')
%         hold on
%         occ = flipud(F.output{1}(i).smoothedoccupancy{1});
%         im = image(occ, 'CDataMapping', 'scaled');
%         im.AlphaData = occ > nanmin(nanmin(occ))-2*std(occ(~isnan(occ) & ~isinf(occ)));
%         
%         hold on
%         sp = flipud(F.output{1}(i).spikes{1});
%         image(sp);
%         
%         subplot(2,2,2)
%         hold on
%         sprate = flipud(F.output{1}(i).smoothedspikerate{1});
%         image(real(sprate), 'CDataMapping', 'scaled');
%         
%         subplot(2,2,3)
%         hold on
%         
%         
%         subplot(2,2,4)
%         hold on
%         autcorr2 = conv2(sprate, sprate, 'full');
%         image(autcorr2, 'CDataMapping', 'scaled');
%         
%         break
%     catch
%         continue
%     end
% end
% %%
% position = [.1 .1 .8 .5];
% for i = 1:size(F.output{1},2)
%     try
%         figure
%         sp = F.output{1}(i).smoothedspikerate{1};
%         image(sp)
%         pause
%     catch
%         continue
%     end
% end