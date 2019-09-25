function out = JY_spkxcorr(index, excludetimes, spikes, lag_window, includestates, minV, f, varargin)
% out = DFL_spkxcorr(index, excludetimes, spikes, lag_window, includestates, minV, f, varargin)
% Filter function that loops through all pairs of cells and calculates the
% maximum unnormalized cross correlation for each cell pair. Saves
% correlation pair matrix to file and saves plots.
% Recommended Filter: multicellanal
% Filter Inputs:
%   spikes - spike structure from filter framework
%   lag_window - maximum lag time to look for cross correlation peak
% Options:
%   None
% OUTPUT:
%   []

% Assume day and epoch are all the same
day = index(1,1);
epoch = index(1,2);
 winbin=0.01;
dsz = '';
if (day < 10)
    dsz = '0';
end
filename=strcat(f.animal{2},f.animal{1},'xcorr',dsz,num2str(day),'.mat');
if(exist(filename,'file'));
    load(filename);
end

%initialize variables
ripple_spikes = cell(1,day);
ripple_spikes{1,day} = cell(1,epoch);
ripple_spikes{1,day}{1,epoch} = cell(2,size(excludetimes,1));

%maps to store correlations and lag times to be saved
xcorr_map{day}{epoch} = zeros(size(index,1),size(index,1));
lag_map{day}{epoch} = zeros(size(index,1),size(index,1));

xcorr=[];
time=[];
max_xcorr=[];
convxcorr=[];
unnorm=[];
cellpairindex=[];
nspikesc2=[];
area=[];


%iterates through all cells
for ii = 1:size(index,1)
    
    tetrode_num = index(ii,3);
    neuron_num = index(ii,4);
    
    
    disp(sprintf('epoch %d, tetrode %d, neuron %d',epoch,tetrode_num,neuron_num));
    
    
    %Neural and position data to process
    spike_data = spikes{day}{epoch}{tetrode_num}{neuron_num};
    
    
    
    
    % bail out of current spike if spike_data.data is empty after filtering (this
    % means there are no valid spikes to compute place fields for in this
    % neuron)
    if(isempty(spike_data.data))
        continue;
    end
    
    % Get exclude indices
    exclind = isExcluded(spike_data.data(:,1),excludetimes);
    
    % Remove spikes that are in exclude times
    % spike_data.data(ismember(spike_data.data(:,7), find(~exclind)),:) = [];
    
    spike_data.data(find(exclind),:) = [];
    
    for jj = 1:size(index,1)
        if ii ~= jj
            disp(sprintf('(%d %d)',ii,jj));
            spike_data2 = spikes{day}{epoch}{index(jj,3)}{index(jj,4)};
            exclind = isExcluded(spike_data2.data(:,1),excludetimes);
            spike_data2.data(find(exclind),:) = [];
            
            % bail out of current spike pair if spike data is empty
            if isempty(spike_data.data)
                %disp(sprintf('cell is empty (%d %d %d %d)',index(ii,1),index(ii,2),index(ii,3),index(ii,4)));
                continue;
            end
            if isempty(spike_data2.data)
                %disp(sprintf('cell is empty (%d %d %d %d)',index(jj,1),index(jj,2),index(jj,3),index(jj,4)));
                continue;
            end
            % convert spike times into spike train
            spike_xcorr = spikexcorr(spike_data.data(:,1)', spike_data2.data(:,1)', winbin, lag_window);
            unnorm = [unnorm;spike_xcorr.c1vsc2];
            
            % normalise nspikes2 is the reference
            spike_xcorr.c1vsc2=spike_xcorr.c1vsc2./sqrt(spike_xcorr.nspikes2*spike_xcorr.nspikes1);
            nspikesc2=[nspikesc2;sqrt(spike_xcorr.nspikes2*spike_xcorr.nspikes1)]; % normalise by the number of comparisons made
            
            %             if max(spike_xcorr.c1vsc2)>0
            %
            %                 figure;
            %                 plot(spike_xcorr.time(:)',spike_xcorr.c1vsc2(:)');
            %                 set(gca,'XTick',-0.5:0.1:0.5)
            %
            %                 title(sprintf('epoch %d, tetrode %d, neuron %d',epoch,tetrode_num,neuron_num));
            %                 cd(f.animal{1,2});
            %                 plotdir = dir('Plot');
            %                 if (isempty(plotdir))
            %                     %an a plot folder needs to be created
            %                     !mkdir Plot
            %                 end
            %
            %                 % change to that directory and saves the figure with file name
            %                 % animal_day_epoch
            %                 cd(strcat(f.animal{1,2},'Plot/'));
            %                 figurename = sprintf('hist_%d_%d_%d',epoch,tetrode_num,neuron_num);
            %
            %                 saveas(gcf, figurename, 'pdf');
            %                 close;
            %
            %             end
            
            
            [Max_xcorr, max_i] = max(spike_xcorr.c1vsc2);
            
            xcorr = [xcorr;spike_xcorr.c1vsc2];
            time = [time;spike_xcorr.time];
            max_xcorr = [max_xcorr;Max_xcorr];
            cellpairindex=[cellpairindex;tetrode_num,neuron_num,index(jj,3),index(jj,4)];
            
            % gaussian smoothing
            sigma=1;
            width = round((6*sigma - 1)/2);
            support = (-width:width);
            gaussFilter = exp( -(support).^2 ./ (2*sigma^2) );
            gaussFilter = gaussFilter/ sum(gaussFilter);
            
            win=gausswin(64,1);
            win=win./sum(win);
            %convxcorr{ii,jj}=conv(spike_xcorr.c1vsc2,win);
            convxcorr=[convxcorr;conv(spike_xcorr.c1vsc2,gaussFilter,'same')];
            area=[area;trapz(spike_xcorr.time,conv(spike_xcorr.c1vsc2,gaussFilter,'same'))];
            
            %             xcorr_map{day}{epoch}(ii,jj) = max_xcorr;
            %             lag_map{day}{epoch}(ii,jj) = (max_i-10000)/10000;
        end
    end
    
end


% saving file and figures
% dsz = '';
% if (day < 10)
%     dsz = '0';
% end
% filename=strcat(f.animal{2},f.animal{1},'xcorr',dsz,num2str(day));
% save(filename,'xcorr_map','lag_map');
% figure(day)
% set(gcf,'position',[0 0 1200 550]);
% set(gcf,'PaperPositionMode','x');
% subplot(2,4,epoch);
%
% set(gca,'FontSize',15);
%
% image(xcorr_map{day}{epoch}.*10)
% title(sprintf('Pairwise X-Corr (D%d E%d)',day,epoch));
%
% xlabel('neuron');
% ylabel('neuron');
%
% [s,mess,messid] = mkdir(sprintf('%sPlot/xcorr/',f.animal{2}));
% print(sprintf('%sPlot/xcorr/%s_xcorr_d%d',f.animal{2},f.animal{1},index(1)),'-depsc');

if isempty(xcorr)
    xcorr=zeros(size(index,1)*size(index,1), lag_window/winbin*2);
    unnorm=zeros(size(index,1)*size(index,1), lag_window/winbin*2);
    nspikesc2=zeros(size(index,1)*size(index,1), 1);
    
end


disp('done');
out.xcorr = xcorr;
out.time = time;
out.max_xcorr=max_xcorr;
out.convxcorr = convxcorr;
out.area=area;
out.unnorm=unnorm;
out.cellpairindex=cellpairindex; % columns 1,2 specify cell 1, columns 3,4 specify cell 2
out.normalisenumber=nspikesc2;
