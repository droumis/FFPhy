function [unit] = ssg_rawwaves (waveforms, times, newassigns, clr, arfp, skip,fg)
% This function plots raw waverforms of the clusters individually and overlaid onto each other.
%
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
%      Therefore in the second version each figure plotted using this function includes 10 figurines
%      each of which plots individual unit, and a figure where all units are plotted on top of each other.
% You have any questions?? You can reach me at clklau@yahoo.it  (Tansu CELIKEL)
% MODIFIED BY: SHANTANU JADHAV (shantanu@ucsd.edu)


Gnrun= unique (newassigns);
Gnrun(find(Gnrun)==0)=[];

%if length(Gnrun)>10;
for lp_gnrun=1:ceil(length(Gnrun)/5)
    
    figure ([fg*10+lp_gnrun-1]); whitebg ([fg*10+lp_gnrun-1]);
    redimscreen80s;
    
    if lp_gnrun~=ceil(length(Gnrun)/5)
        nrun=Gnrun((lp_gnrun-1)*5+1:(lp_gnrun-1)*5+5);
    elseif lp_gnrun==ceil(length(Gnrun)/5)
        nrun=Gnrun((lp_gnrun-1)*5+1:length(Gnrun));
    end
    
    for lp_nrun = 1: size (nrun,1)
        
        % PLOT WAVES
        subplot (size (nrun,1),4,skip(lp_nrun));
        unit = waveforms (find (newassigns == nrun(lp_nrun)), :);
        
        if (size (unit,1) > 50)
            unit_sm = unit(randperm (size (unit,1)),:);
            unit_sm = unit_sm (1:50, :);
            plot (unit_sm', [clr(lp_nrun)]); axis tight;    
            text (2,(max (max (unit))*0.8), ['#Sp: ' num2str(size (unit_sm,1))])
            set (gca, 'XTickLabel', []);
        else
            plot (unit', [clr(lp_nrun)]); axis tight;
            text (2,(max (max (unit))*0.8), ['#Sp: ' num2str(size (unit,1))])
            set (gca, 'XTickLabel', []);
        end
        
        ylabel ([ 'C' num2str(nrun (lp_nrun))]); 
        if lp_nrun == size (nrun,1); xlabel ('Sample (32 pts/1 ms)'); end
        
        isis = diff(sort(times (find (newassigns == nrun(lp_nrun)), :))); % 0.1 ms resolution
        nrcsp = length(find(isis <= arfp)); % number of coincident spikes 
        
        text (29,(max (max (unit))*0.8), ['CoSp ' num2str(nrcsp)])
        
        subplot (2,4,[2 3 4  6 7 8]);
        plot (unit', [clr(lp_nrun)]); axis tight;
        hold on
    end
end
%end
