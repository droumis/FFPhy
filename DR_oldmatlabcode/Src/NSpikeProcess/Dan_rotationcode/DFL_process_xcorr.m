function process_xcorr(directoryname,fileprefix,days, varargin)

lowercasethree = '';

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'lowercasethree'
            lowercasethree = varargin{option+1};
            
    end
end

%epoch numbers hardcoded for convenience
run_epochs = [2 4 6];
sleep_epochs = [1 3 5 7];


days = days(:)';

for day = days
    dsz = '';
    if (day < 10)
        dsz = '0';
    end
    
    eval(['load ',directoryname,fileprefix,'xcorr',dsz,num2str(day), '.mat']);
    
    xcorr_diff = diff(reshape([xcorr_map{day}{run_epochs}],[size(xcorr_map{day}{run_epochs(1)}),length(run_epochs)]),1,3);
    
    for ii = 1:size(xcorr_diff,3)
        figure(1);
        set(gcf,'position',[0 0 800 400]);
        set(gcf,'PaperPositionMode','auto');
        subplot(1,2,ii);
        set(gca,'FontSize',15);
        image(xcorr_diff(:,:,ii)*10);
        title(sprintf('Xcorr Diff Run E%d-E%d (D%d)',run_epochs(ii+1),run_epochs(ii),day));
        xlabel('neuron')
        ylabel('neuron')
        
    end
    [s,mess,messid] = mkdir(sprintf('%sPlot/xcorr/',directoryname));
    print(sprintf('%sPlot/xcorr/%s_xcorr_rundiff_d%d',directoryname,fileprefix,day),'-dpng');

    
    figure(2);
    set(gcf,'position',[0 0 800 400]);
    set(gcf,'PaperPositionMode','auto');
    for ii = 1:2
        subplot(1,2,ii)
        set(gca,'FontSize',15);
        
        image(xcorr_diff(:,:,ii).*xcorr_map{day}{ii}.*10)
        title(sprintf('Xcorr Run E%d-E%d Overlap Sleep E%d (D%d)',run_epochs(ii+1),run_epochs(ii),sleep_epochs(ii+1),day));
        xlabel('neuron')
        ylabel('neuron')
    end
    print(sprintf('%sPlot/xcorr/%s_xcorr_rundiff_sleepoverlap_d%d',directoryname,fileprefix,day),'-dpng');
end