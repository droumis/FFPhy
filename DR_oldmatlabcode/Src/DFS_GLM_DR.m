
warning('off','all');
clear;

runscript = 1;
savedata = 1; % save data option - only works if runscript is also on
savefigs=1;
plotstuff = 1; %
plotEachcell = 1;
runnadal =0;

savedir = '/mnt/data25/sjadhav/HPExpt/ProcessedDataDR/';
savefilename = 'HPaGLMtest';
savefile = [savedir savefilename]; area = 'PFC'; %clrunmod = 'r'; clrmod = 'b'; % PFC
figdir = '/mnt/data25/sjadhav/HPExpt/Figures_DR/';

Veqn = '>3';
minV=str2num(Veqn(end));
mintime = 3;
traj = [1:4] ;

% If runscript, run Datafilter and save data
if runscript == 1
    %     for i =  modUnmod;
    %Animal selection
    %-----------------------------------------------------
    %         animals = {'HPa' 'HPb' 'HPc'};
    %         animals = {'HPa' 'HPb' 'HPc'};
    %         animals = {'HPc'};
    animals = {'HPa'};
    %         animals = {'nadal'};
    
    %Filter creation
    %-----------------------------------------------------
    % Epoch filter
    % -------------
    if runnadal == 1;
        dayfilter = '8:17';
    else
        dayfilter = '1:8';
    end
    
    if runnadal == 1;
        runepochfilter = 'isequal($type, ''run'')';
    else
        runepochfilter = 'isequal($environment, ''wtr1'') || isequal($environment, ''wtr2'')';
    end
    
    % %Cell filter
    % %-----------
    placecellfilter = '(strcmp($area, ''PFC'') && ($numspikes > 100))';  % not mod/unmod
    %
    %         if i < 2;
    %                         placecellfilter = '(strcmp($area, ''PFC'') && strcmp($ripmodtag, ''y'') && ($numspikes > 100))';   % Ripple mod
    % %             placecellfilter = '(strcmp($area, ''PFC'') && strcmp($thetamodtag, ''y'') && ($numspikes > 100))';   % theta mod
    %
    %         else
    %                         placecellfilter = '(strcmp($area, ''PFC'') && strcmp($ripmodtag, ''n'') && ($numspikes > 100))'; % Ripple unmod
    % %             placecellfilter = '(strcmp($area, ''PFC'') && strcmp($thetamodtag, ''n'') && ($numspikes > 100))'; % theta unmod
    %         end
    
    % Time filter -
    %%-----------
    
    riptetfilter = '(isequal($descrip, ''riptet''))';
    timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 3))', 6}, {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} }; %DR added velocity filter.. trying to get ride of v high prococc in data
    
    % Iterator
    % --------
    iterator = 'singlecellanal';
    
    % Filter creation
    % ----------------
    spatf = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter,'cells',placecellfilter,'excludetime', timefilter,'iterator',iterator);
    
    %do i need this?
    spatf = testexcludetimes(spatf, mintime); %removes epochs from analysis if all epoch excluded by excludetimes, mintime = 30
    
    disp('Done Filter Creation');
    
    % Set analysis function
    % ----------------------
    %use_____________________
    psf = setfilterfunction(spatf, 'DFAsj_filtercalclinfields_tf',{'spikes', 'linpos'}, 'binsize', 2);
    
    % Run analysis-----------------------
    %         pfs{i} = runfilter(pfs);  % Place Field Stability.. trajectories
    pfs = runfilter(psf);  % Place Field Stability.. trajectories
    %     end
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear runscript  savedata plotstuff cyclemaps savefigs plottrajs runnadal plotEachcell
        save(savefile);
    end
else
    load(savefile);
end % end runscript

if ~exist('savedata')
    return
end

% -------------------------  Filter Format Done -------------------------

% ----------------------------------
% PLOT!

% --------------------------------------------------------------------

for i = 1:length(pfs.output{1,1});
    % outbound
    %generate rate vector: 28:67 ~ center arm
    OUTrightrates = pfs.output{1,1}(1,i).trajdata{1,3}(28:67,5); % OR = 3
    OUTleftrates = pfs.output{1,1}(1,i).trajdata{1,1}(28:67,5); % OL = 1
    INrightrates = pfs.output{1,1}(1,i).trajdata{1,4}(28:67,5); % IR = 4
    INleftrates = pfs.output{1,1}(1,i).trajdata{1,2}(28:67,5); % IL = 2
    
    % concatenate right and left
    xOUT= [OUTrightrates; OUTleftrates];
    xIN= [INrightrates; INleftrates];
    
    %generate logic vector ones(length(right rate)) zeros(length(left rate))
    %and concatenate
    yOUT = [ones(length(OUTrightrates), 1); zeros(length(OUTleftrates),1)];
    yIN = [ones(length(INrightrates), 1); zeros(length(INleftrates),1)];
    
    % fit outbound
    [bOUT,devOUT,statsOUT] = glmfit(xOUT,yOUT,'binomial','link','logit');
    xxOUT = linspace(min(xOUT), max(xOUT),100);
    yfitOUT = glmval(bOUT,xxOUT,'logit');
    
    % fit inbound
    [bIN,devIN,statsIN] = glmfit(xIN,yIN,'binomial','link','logit');
    xxIN = linspace(min(xIN), max(xIN),100);
    yfitIN = glmval(bIN,xxIN,'logit');
    
    %save all results
    bINall{1,i} = bIN;
    devINall{1,i} = devIN;
    statsINall{1,i} = statsIN;
    yfitINall{1,i} = yfitIN;
    statsINallPvalues(1,i) = statsIN.p(1,1);
    
    bOUTall{1,i} = bOUT;
    devOUTall{1,i} = devOUT;
    statsOUTall{1,i} = statsOUT;
    statsOUTallPvalues(1,i) = statsOUT.p(1,1);
    yfitOUTall{1,i} = yfitOUT;
    
    
    if plotEachcell == 1;
        % plot ( green background = p <0.05)
        figure;
        subplot(2,1,1);
        plot(xOUT,yOUT,'o', xxOUT,yfitOUT,'-');
        PvalueOUT = statsOUT.p(1,1);
        str=sprintf('p  = %d', PvalueOUT);
        title({'Probability of outbound right'; str});
        if statsOUT.p(1,1) < pthresh;
            set(gca,'Color',[.85 1 .85]);
        end
        
        subplot(2,1,2);
        plot(xIN,yIN,'o', xxIN,yfitIN,'-');
        PvalueIN = statsIN.p(1,1);
        str=sprintf('p  = %d', PvalueIN);
        title({'Probability of inbound right'; str});
        if statsIN.p(1,1) < pthresh;
            set(gca,'Color',[.85 1 .85]);
        end
        pause
        close
    end
end

sigOUT = sum(statsOUTallPvalues(1,:) < pthresh)
sigIN = sum(statsINallPvalues(1,:) < pthresh)
sigOUTprop = sigOUT./(length(statsOUTallPvalues(1,:)));
sigINprop = sigIN./(length(statsOUTallPvalues(1,:)));

hist(sigINprop)




