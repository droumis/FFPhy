


loadfile = 0;
plotEachcell = 0;
pthresh = 0.05;

%load data in
if loadfile == 1;
    savedir = '/mnt/data25/sjadhav/HPExpt/ProcessedDataDRtemp/';
    savefile = [savedir 'HPb_PFCfields.mat']; % HPa and HPb fields  %%DR
    load(savefile);
end

%%
for i = 1:length(psf.output{1,1}(1,:));
    % outbound
    %generate rate vector: 28:67 ~ center arm
    OUTrightrates = psf.output{1,1}(1,i).trajdata{1,3}(28:67,5); % OR = 3
    OUTleftrates = psf.output{1,1}(1,i).trajdata{1,1}(28:67,5); % OL = 1
    INrightrates = psf.output{1,1}(1,i).trajdata{1,4}(28:67,5); % IR = 4
    INleftrates = psf.output{1,1}(1,i).trajdata{1,2}(28:67,5); % IL = 2
    
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
%%
 
%make this is a loop for each cell.. 
%maybe this can be done as vectorize rather than for loop?



























%% start of doing a t-test.. table  for now until glm is done.
% function avgcenter;
% 
% savedir = '/mnt/data25/sjadhav/HPExpt/ProcessedDataDRtemp/';
% savefile = [savedir 'HPb_PFCfields']; % HPa and HPb fields  %%DR
% load(savefile);
% 
% for cell = 1:length(psf.output(1,:));
%     for traj = 1:4;
%         avgcenter{cell}(:,traj) = avg(psf(1).output{cell}.trajdata{traj}(:,5));
%     end
% end
%%