function sss_rawrastersisi (waveforms, times, newassigns, clr, arfp, skip, skipp, stimonset, x_isi, fg, sweepD, nsweeps, plot_clu)
% This function plots raw waverforms of the clusters sorted using Samar&David's spike sorting algorithm.
% In addition to the waves, Raster plots together with ISI histograms are generated.
%
% NOTE: A number of additional codes are provided to draw PSTHs and Autocorrelations for each cluster.  However
% those codes are commented to use later, when needed.
%
% ALSO SEE: SDSS 

% Feldman Lab @ UCSD
% April 2002 - v.1
% March 2003 - v.2
%   What's new in this version?
%   1. In the first version of the program, the maximum number of the units which could be 
%      sorted was limited, considering that the number of the units we can record with the
%      current in vivo set-up can't be substantially large.  
%      In the second version, however, we rule out this limit for the number of units which 
%      could be sorted/plotted considering the fact that following the iontophoretic application
%      (in Dr. Foeller's study)of SR the number of unit which could be recorded should be larger 
%      than the no-drug condition. 
%
% You have any questions?? You can reach me at clklau@yahoo.it  (Tansu CELIKEL)
% MODIFIED BY: SHANTANU JADHAV (shantanu@ucsd.edu)

Gnrun=plot_clu;
%Gnrun= unique (newassigns);

%nrun= unique (newassigns);

s_tms = mod (times,10000); % bin size =0.1 ms
sweeps = (times-s_tms)/10000;
s_tms = round (s_tms/10); % bin size =1 ms


%if length(Gnrun)>10;

for lp_gnrun=1:ceil(length(Gnrun)/5)
    figure ([fg*10+lp_gnrun-1]); whitebg ([fg*10+lp_gnrun-1]);
    redimscreen70s;
    
    if lp_gnrun~=ceil(length(Gnrun)/5)
        nrun=Gnrun((lp_gnrun-1)*5+1:(lp_gnrun-1)*5+5);
    elseif lp_gnrun==ceil(length(Gnrun)/5)
        nrun=Gnrun((lp_gnrun-1)*5+1:length(Gnrun));
    end
    
    for lp_nrun = 1: size (nrun,1)
        
        %%%%%%%% PLOT WAVES
        subplot (size (nrun,1),4,skip(lp_nrun));
        unit = waveforms (find (newassigns == nrun(lp_nrun)), :);
        
        plot (unit', [clr(lp_nrun)]); axis tight;
        
        text (2,(max (max (unit))*0.80), ['#Sp: ' num2str(size (unit,1))])
        ylabel ([ 'C' num2str(nrun (lp_nrun))]);
        if lp_nrun == 1; title ('Raw Waveforms'); elseif lp_nrun == size (nrun,1); xlabel ('Sample (32 pts/1 ms)'); end
        
        %%%%%%%% ISIS
        %isis = diff(sort(times (find (newassigns == nrun(lp_nrun)), :))); % 0.1 ms resolution
        isis = diff(times (find (newassigns == nrun(lp_nrun)), :));
        
        nrcsp = length(find ((isis <= arfp) & (isis>0))); % number of coincident spikes 
                
        text (29,(max (max (unit))*0.8), [ 'CoSp ' num2str(nrcsp)]) % fprint number of coincident spikes
        
        %%%%%%%%%%%%%% RASTERS
        %     % draw PSTH
        %     subplot (size (nrun,1),4,skipp(lp_nrun));
        %     hist (ttime, 700);
        %     axis tight;   
        %     if lp_nrun == 1; title ('PSTH'); elseif lp_nrun == size (nrun,1); xlabel ('Time (ms)'); end
        
        ttime = s_tms (find (newassigns == nrun(lp_nrun)));
        
        
        % draw Raster Plot
        subplot (size (nrun,1),4,skipp(lp_nrun));
        rastersw = sweeps (find (newassigns == nrun(lp_nrun)));
        
        % rastersw = collapse_sweeps(rastersw);   DO NOT COLLAPSE SWEEPS
      
        plot (ttime, rastersw, 'b.')
        %x=sweeps (find (newassigns == nrun(lp_nrun)));
        %plot(ttime, 1:length(x), 'b.')
        ylimit = length(rastersw);
        if ylimit <= 1
            ylimit=5;
        end
            axis ([0 sweepD 1 nsweeps])
        if lp_nrun == 1; title ('Overall'); elseif lp_nrun == size (nrun,1); xlabel ('Time (ms)'); end
        
        
        %%%%%%%%% ONSET RASTER
        % draw Raster Plot but only for the onset response
        subplot (size (nrun,1),4,skipp(lp_nrun)+1);
        
        plot (ttime, rastersw, 'w.');
        axis tight; axis ([stimonset stimonset+50 1 nsweeps])
        if lp_nrun == 1; title ('Onset Response'); elseif lp_nrun == size (nrun,1); xlabel ('Time (ms)'); end
        
        %%%%%%%% ISI HISTOGRAM
        subplot (size (nrun,1),4,skipp(lp_nrun)+2);
        isi_vct = zeros (1, x_isi+1); 
        nrcsp = isis(find (isis <= x_isi)); % find those spikes which have ISI bigger then or equal to 0 but smaller then or equal to x_isi 
        %   disp (['# of entries in nrcsp.................' num2str(length (nrcsp))])
        if length(nrcsp)==0, nrcsp=x_isi; end
        distr= histc (nrcsp,0:x_isi); 
        bar(0:x_isi,distr,'histc'); %axis ([0 x_isi 0 50])
        if lp_nrun == 1; title ('ISI Distributions'); elseif lp_nrun == size (nrun,1); xlabel ('ISI (bin size 1 ms)'); end
        
        
        %%%%%%%%% AUTOCORRELATIONS
        %     % draw the autocorrelations for each unit
        %     subplot (size (nrun,1),4,skipp(lp_nrun)+2);
        %     vct = hist (ttime, 700); plot (xcorr (vct, vct)); 
        %     if lp_nrun == 1; title ('ISI Distributions'); elseif lp_nrun == size (nrun,1); xlabel ('ISI (ms)'); end
        
    end
end 
%end