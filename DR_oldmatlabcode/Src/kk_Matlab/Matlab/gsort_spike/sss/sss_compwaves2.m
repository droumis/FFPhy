function sss_compwaves2 (waveforms, times, newplots, arfp)

% compwave2 is a sister code of compwaves.  The way it differs from the original is that
% the citations about window closing is omitted. 
% All rights reserved
% Tansu Celikel

% MODIFIED: SHANTANU JADHAV (FEB 2005)

qlp=1; fstimes=[];
while qlp==1;
    Tlines={'Clusters - Please enter only cluster numbers with ; in between', ['Cluster ' 'color code= combine(0) or no combine(1)']}; 
    defent={'', ''}; % default entries
    infonames= {'Clusters', 'Coloropt'};
    info = inputdlg(Tlines, 'Which clusters you want to plot overlaid', 1, defent); hold on
    
    if ~isempty(info)              %see if user hit cancel
        info = cell2struct(info,infonames);
        clustnum = str2num(info.Clusters);   %convert string to number   
    end
    clr = ['w' 'y' 'm' 'b' 'c' 'r' 'w' 'y' 'm' 'b' 'c' 'r' 'w' 'y' 'm' 'b' 'c' 'r'];
        
    figure(1001); whitebg(1001)
    Tsp=zeros(size(times,1),1);
    
    for lp=1:size(clustnum,1)
        sp_ind=(find(newplots==clustnum(lp)));
        
        if str2num(info.Coloropt)==0; % if units combined 
            temptimes=times(sp_ind);
            tempwaves=waveforms(sp_ind,:);
            
            fstimes=[fstimes;temptimes];
            plot(tempwaves',clr(lp)); hold on;
             
            if lp==size(clustnum,1)
                set (gca, 'XLim',[1 32]); 
                yl=get(gca, 'Ylim');
                dtimes=diff(sort(fstimes/1000));        
                title(['Clusters ' num2str(clustnum') ' | Nspikes=  ' num2str(length(fstimes)) ' | <ISI spikes=  ' num2str(length(find ((dtimes<arfp/1000) & (dtimes >0))))])
            end
            
        else  % if units are NOT combined
            temptimes=times(sp_ind);
            tempwaves=waveforms(sp_ind,:);
            subplot (4,6,[1:18])
            plot(tempwaves',clr(lp)); hold on;
            set (gca, 'XLim',[1 32]); 
            ylabel ('Voltage (uV)');
            
            subplot (4,6,[19:24])
            axis off;
            text(0,1+(-0.2*lp),['Cluster ' num2str(clustnum(lp)) ' (' clr(lp) ')']);
            text(0.2,1+(-0.2*lp),['Nspikes=  ' num2str(length(temptimes))]);
            dtimes=diff(sort(temptimes/1000));  
            text(0.5,1+(-0.2*lp),['<ISI spikes=  ' num2str(length(find ((dtimes<arfp/1000) & (dtimes >0))))]);
               
        end
        clear sp_ind temptimes tempwaves;
    end
    
    ButtonName=questdlg('Would you like to continue raw waveform comparison using COMPWAVES?', ...
        'Ask to continue', ...
        'Yes','No','No');
    
    %keyboard
    if size(ButtonName,2)==3;
        qlp=1; close([1001]); dtimes=[];fstimes=[];Tlines=[];
    else
        qlp=0; close([1001]);
    end
    
end
