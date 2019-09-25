%% add back the reference and plot the 2 versions

cd /data19/droumis/bob/bobstim_proc/EEG/

ref = load('bobeeg25-2-02');
ref = ref.eeg{1,25}{1,2}{1,2}.data;
curchan = load('bobeeg25-2-06');
curchan = curchan.eeg{1,25}{1,2}{1,6}.data;

refaddedchan = curchan + ref;

subplot(3,1,1); plot(ref); axis([0 450027 -1000 1000 ]); hold on;
title 'reference (2)'
subplot(3,1,2); plot(curchan, 'Color', 'g'); axis ([ 0 450027 -1000 1000]); hold on;
title 'HC tet (6)'
subplot(3,1,3); plot(refaddedchan, 'Color', 'r'); axis ([ 0 450027 -1000 1000]); hold on;
title 'HC tet, ref added back (2+6)'

%%if they aren't the same length try to use this, but i didn't actually
%%get it to work efficiently

% reflimit = 450027;
% if max(size(curchan)) ~= max(size(ref)); %not of equal length
%     newchan = zeros(reflimit,1);
%     for cutoff = 1:reflimit;
%         newchan(cutoff,1) = curchan(cutoff,1);
%     end
%     refaddedchan = newchan + ref;
%     figure
%     plot(newchan,1); hold on;
% else %equal length
%     refaddedchan = curchan + ref;
%     figure
% 
% end
% plot(refaddedchan, 'Color', 'r');


   
        
        
        
        
        
        
        
        
 