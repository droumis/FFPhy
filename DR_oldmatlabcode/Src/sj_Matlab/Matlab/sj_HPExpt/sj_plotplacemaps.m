
% Plot place maps from multiple cells on 1plot to get an idea of population
% First used for PFC


% Similar to combine theta and ripple. Compare run and sleep ripplemodln

% USe the saved data files - HP_ripplemodsleep and HP_ripplemod to plot
% correlations between theta modulation and ripple modulation'

clear;
savefig1=0;
savedir = '/data25/sjadhav/HPExpt/ProcessedData/';
plotep = 2; % Only do epoch 2 for now. Later, you can do both epochs in a day.
clr = {'b','r','g','m','c','y','k','r'};

area = 'PFC';
datafile = [savedir 'HPa_',area,'fields'];
load(datafile);


for an = 1:length(psf)
    for i=1:length(psf(an).output{1}),
        index{an}(i,:)=psf(an).output{1}(i).index;
        alltrajdata{an}{i}=psf(an).output{1}(i).trajdata;
        allmapdata{an}{i}=pmf(an).output{1}(i);
        %allmapdata_sep{an}{i}=pof(an).output{1}(i);
    end
end

% Figures for individual cells

figdir = '/data25/sjadhav/HPExpt/Figures/PlaceFields/IndividualCells/';
summdir = figdir;
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
set(0,'defaultaxesfontsize',25);
tfont = 25;
xfont = 25;
yfont = 25;

figcnt = 0; totalplots = 0; %Count total plots across figures
nsubplots=5; % per figure

for an = 1:length(psf)
    prefix = psf(an).animal{1};
    
    % Get days and epochs for current animal
    days = unique(psf(an).epochs{1}(:,1));
    allepochs = unique(psf(an).epochs{1}(:,2));
    
    for d = 1:length(days)
        day = days(d);
        if day==1,
            allepochs = [4 6];
            plotep=4;
        else
            allepochs = [2 4];
            plotep=2;
        end
        
        
        dayidxs = find(index{an}(:,1)==day);
        % Get tetlist
        tet = unique(index{an}(dayidxs,3));
        
        for elec = 1:length(tet)
            currtet = tet(elec);
            % Get celllist
            daytetidxs = find(index{an}(:,1)==day & index{an}(:,3)==currtet);
            cells = unique(index{an}(daytetidxs,4));
            
            for neuron = 1:length(cells)
                currcell = cells(neuron);
                curridx = find( index{an}(:,1)==day & index{an}(:,2)==plotep & index{an}(:,3)==currtet & index{an}(:,4)==currcell);
                
                % To control plottin
                %if curridx(1)==1 && curridx(2)==5 && curridx(3)==18 && curridx(4)==2

                
                if ~isempty(curridx)
                    trajdata = alltrajdata{an}{curridx};
                    mapdata = allmapdata{an}{curridx}.smoothedspikerate;
                    
                    if mod(totalplots,5)==0 % 5 subplots finished
                        figcnt=figcnt+1;
                    end
                    
                    figure(figcnt); redimscreen; hold on;
                    subplot(2,5,mod(totalplots,5)+1); hold on;
                    imagesc(flipud(mapdata)); %colorbar
                    set(gca,'YLim',[0 100]); %100 for HPa, 110 for HPb
                    set(gca,'XLim',[0 100]);
                    text(65,96,sprintf('%2.1f Hz',max(mapdata(:))),'FontSize',24,'Fontweight','normal');
                    title(sprintf('%s D%d t%d c%d',prefix, day, currtet, currcell),'FontSize',tfont);
                    if mod(totalplots,5)==0
                        xlabel ('X-position (cm)','FontSize',24,'Fontweight','normal');
                        ylabel ('Y-position (cm)','FontSize',24,'Fontweight','normal');
                    end
                    
                    subplot(2,5,mod(totalplots,5)+1+5); hold on;
                    for i=1:length(trajdata),
                        plot(trajdata{i}(:,5),[clr{i} '.-'],'Linewidth',2);
                        %maxrate(i) = max(trajdata{i}(5:end-5,5));
                    end
                    
                    if mod(totalplots,5)==0
                        legend('OutLeft','InLeft','OutRight','InRight');
                        xlabel ('Posn along linear traj (cm)','FontSize',24,'Fontweight','normal');
                        ylabel ('Firing Rate (Hz)','FontSize',24,'Fontweight','normal');
                    end
                    
                    % Update plot number
                    totalplots=totalplots+1;
                    
                end %~isempty curridx
                
                % Saving fig
                if mod(totalplots,5)==0
                    figfile = [figdir,'HPa_PFC_placemaps',num2str(figcnt)];
                    keyboard;
                    %print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
                end
                
                
                
            end % end cell
            
        end % end elec
        
    end % end day
end % end anidx







