function [initial_width, later_width, fwhm_a, amp_a, peak1, peak2, interpol, allpoints, peaks, hindex_a, swaves, ori_index_a] = new_spikewidths(filename,figopt1, figopt2)
% Load mat file, and get spike widths 
% FWHM: TAKING ACTUAL PEAK-PEAK AMPL AND THEN LOOKING AT FWHM
% Plot spikes if asked for - figopt=1

[recording,cluster]=strtok(filename,'C');
load ([filename]);

if figopt1==1
%PLOT RAW WAVES
%figure;hold on; plot(1:32,swaves,'b.-'); plot (mean (swaves), 'y', 'Linewidth', 4); title([cluster]); hold off;
figure; hold on; errorbar(mean (newwaves),std(newwaves)*(1/sqrt(size(newwaves,1))), 'b.-');  hold off;
end

thr=14; res = 0.03125;   % Threshold crossing at 14, Res=1/32kHz in ms
[nspikes,datapts]=size(swaves);
x = mean(newwaves);

xi=1:48*5;
a=1:5:48*5;
yi = interp1(a,x,xi,'spline');
x=yi;
newres=res/5; thr=a(14);


% % NEW IMPLEMENTATION
% 
% peak = max(x);
% index = min(find(x==peak);
% 
% min_amp = min(x(1:thr));
% trough1 = min(find(x==min_amp);


%---------------------------------------------------------------------
% OLD IMPLEMENTATION BY RUNNING THROUGH THE WHOLE WAVEFORM


% COULD DO THINGS A LOT EASIER BY USING FIND MIN AND MAX- KEEP THIS FOR NOW
% -------------------------------------------------------------------------

% FIND 1st PEAK AFTER THRESHOLD
flag =0; i=0;
while flag ==0,
    if x(thr+i+1) < x(thr+i), 
        peak=x(thr+i);                              % PEAK
        index=thr+i;
        flag =1;
    end,  
    i=i+1; 
end

flag =0; i=0; trough1=1;
while flag==0 & (index-i-1)>0
     if ( x(index-i) <= x(index-i-1) )            % TROUGH1
        trough1=(index-i); 
        flag=1;
     end
i=i+1;
end

flag =0; i=0; trough2=48*5;
while flag==0 & (index+i+1)<=48*5,
    if ( x(index+i) <= x(index+i+1) )     % FIND TROUGH2=PEAK2: IF THERE IS NONE: LAST POINT =48
        trough2 = index+i;
        flag=1;
    end
i=i+1;
end

flag =0; i=0; peak3=48*5;
while flag==0 & (trough2+i+1)<=48*5,
    if ( x(trough2+i+1) <= x(trough2+i) )     % FIND PEAK3: IF THERE IS NONE: LAST POINT =48
        peak3 = trough2+i;
        flag=1;
    end
i=i+1;
end

%-------------------------------

flag=0; i=0; point1 = trough1;
while flag==0 & (index-i-1)>trough1
     if ( x(index-i)>=0 & x(index-i-1)<=0 )   % 1st zero crossing: IF NOT, 1ST POINT = TROUGH1
        point1=(index-i); 
        flag=1;
     end
i=i+1;
end

flag =0; i=0; point2 = trough2;
while flag==0 & (index+i+1)<=trough2,
    if ( x(index+i)>=0 & x(index+i+1)<=0 )     % FIND 2nd crossing: IF NOT, 2nd POINT = TROUGH2
        point2 = index+i;
        flag=1;
    end
i=i+1;
end


flag =0; i=0; point3 = 48*5;
while flag==0 & (trough2+i+1)<=48*5,
    if ( x(trough2+i)<=0 & x(trough2+i+1)>=0 )     % FIND 3rd crossing: IF NOT, 3rd POINT = 48
        point3 = trough2+i;
        flag=1;
    end
i=i+1;
end


% IF AHP CONTAMINATION (MOSTLY IN FS UNITS, THEN TAKE HALF-MAX-PEAK3 AS POINT3)
if (x(trough2)>=0 & point3==48*5)
    
    hmax3 = x(trough2) + (x(peak3)-x(trough2))/2;
    flag =0; i=0; hmax3_id=peak3;
    while flag==0 & (trough2+i+1)<=peak3
        if (x(trough2+i) <= hmax3) &  (x(trough2+i+1) >= hmax3)
            hmax3_id = (trough2+i);
            flag=1;
        end
        i=i+1;
    end
    point3=hmax3_id;
end
    
initial_width = newres*(point2-point1);
later_width = newres*(point3 - point2);
allpoints = [point1 point2 point3];

%----------------------------------------------------------------

% HALF-MAX: MIDPOINT BETWEEN PEAK AND MIN_AMP


% min_amp_a=x(trough1);

 min_amp_a=min([x(trough1) x(trough2)]);

half_max_a = (abs(peak)+abs(min_amp_a))/2; 

hmax_a=peak-half_max_a; hindex_a = [1 1 48*5];

flag =0; i=0;
while flag==0 & (index+i+1)<=48*5,
     if (x(index+i) > hmax_a) &  (x(index+i+1) <= hmax_a)
        hindex_a(3)=(index+i+1);
        flag=1;
    elseif ( x(index+i) <= x(index+i+1) )        % IF YOU GET A TROUGH BEFORE YOU REACH HALF-MAX, TAKE TROUGH AS END POINT
        hindex_a(3)= index+i;
        flag=1;
    end
i=i+1;
end

flag =0; i=0;
while flag==0 & (index-i-1)>0
     if (x(index-i) > hmax_a) &  (x(index-i-1) <= hmax_a)
        hindex_a(2)=(index-i-1);
        flag=1;
     elseif ( x(index-i) <= x(index-i-1) )
        hindex_a(2)=(index-i); 
        flag=1;
     end
i=i+1;
end
%-------------------------------------------

interpol=x;
hindex_a(4)=index; % LOCATION OF PEAK1
ori_index_a=hindex_a/5;
fwhm_a = 2*newres*(hindex_a(3)-hindex_a(2));             % FWHM (factor of 2)
amp_a = (abs(peak)+abs(min_amp_a));

peak1=abs(peak);
peak2=abs(x(trough2));
peaks=[index trough2];

% -------------------------------- OLD OLD OLD --------------------

% peak1=0
% if trough1~=0
%     flag =0; i=0;
%     while flag==0 & (trough1-i-1)>0
%      if ( x(trough1-i) >= x(trough1-i-1) )            
%         peak1=(trough1-i); 
%         flag=1;
%      end
%     i=i+1;
%     end
% end
% 
% fw= newres*(trough-peak1);             % For Full Width (factor of 2)


% ---------------------------------------------------------------------------------




