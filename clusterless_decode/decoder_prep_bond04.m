load('/opt/data13/kkay/Bon/bonpos04.mat');
load('/opt/data13/kkay/Bon/bonspikes04.mat');
load('/opt/data13/kkay/Bon/bonlinpos04.mat');
load('/opt/data13/kkay/Bon/bonripples04.mat');

% preprocess for position

ex=4;
ep=2;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%            Linearization           %%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% indices of segment X
seg1=find(linpos{ex}{ep}.statematrix.segmentIndex==1);  
seg2=find(linpos{ex}{ep}.statematrix.segmentIndex==2);
seg3=find(linpos{ex}{ep}.statematrix.segmentIndex==3);
seg4=find(linpos{ex}{ep}.statematrix.segmentIndex==4);
seg5=find(linpos{ex}{ep}.statematrix.segmentIndex==5);

time=linpos{ex}{ep}.statematrix.time;
lindist=linpos{ex}{ep}.statematrix.lindist;

if 1
    figure
    plot(time(seg1),lindist(seg1),'b.');
    hold on
    plot(time(seg4),lindist(seg4),'k.');
    plot(time(seg5),lindist(seg5),'c.');
    plot(time(seg2),-lindist(seg2),'r.');  % negative for 
    plot(time(seg3),-lindist(seg3),'m.');
    hold off
end

%%% Center arm: Segment 1

posT1=time(seg1);                           % posT1: segment 1 time points
posL1=lindist(seg1)/max(lindist(seg1));     % posL1: segment 1 % dist travelled  

posL1b = [];

for k=1:length(posL1)
    if posT1(k) >= 2460 && posT1(k) <= 2545.681
        posL1b(k)=-posL1(k); %0=-1              % makes negative posL1
    elseif posT1(k)>=2616.645 && posT1(k)<=2649.484
        posL1b(k)=-(posL1(k)+2*(3-posL1(k))); %-5=-6
    elseif posT1(k)>=2726.31 && posT1(k)<=2744.2884
        posL1b(k)=posL1(k)+2*(3-posL1(k)); %5=6
    elseif posT1(k)>=2744.2885 && posT1(k)<=2765.212
        posL1b(k)=-posL1(k); %0=-1
    elseif posT1(k)>=2816.852 && posT1(k)<=2836.1514
        posL1b(k)=-(posL1(k)+2*(3-posL1(k))); %-5=-6
    elseif posT1(k)>=2836.1515 && posT1(k)<=2849.932
        posL1b(k)=-posL1(k); %0=-1
    elseif posT1(k)>=2931.63 && posT1(k)<=2951.1644
        posL1b(k)=posL1(k)+2*(3-posL1(k)); %5=6
    elseif posT1(k)>=2951.1645 && posT1(k)<=2965.08
        posL1b(k)=-posL1(k); %0=-1
    elseif posT1(k)>=3059.031 && posT1(k)<=3077.7624
        posL1b(k)=posL1(k)+2*(3-posL1(k)); %5=6
    elseif posT1(k)>=3077.7625 && posT1(k)<=3091.307
        posL1b(k)=-posL1(k); %0=-1
    elseif posT1(k)>=3124.28 && posT1(k)<=3139.592
        posL1b(k)=-(posL1(k)+2*(3-posL1(k))); %-5=-6
    elseif posT1(k)>=3225.24 && posT1(k)<=3240.9454
        posL1b(k)=posL1(k)+2*(3-posL1(k)); %5=6
    elseif posT1(k)>=3240.9455 && posT1(k)<=3249.698
        posL1b(k)=-posL1(k); %0=-1
    elseif posT1(k)>=3297.737 && posT1(k)<=3315.681
        posL1b(k)=-(posL1(k)+2*(3-posL1(k))); %-5=-6
    elseif posT1(k)>=3385.306 && posT1(k)<=3405
        posL1b(k)=posL1(k)+2*(3-posL1(k)); %5=6
    else
        posL1b(k)=posL1(k);
    end
end


%%% 0-ACE    % OUTER, A > C > E

    % Segment 4: A > C
posT4=time(seg4);                                     % posT4: segment 4 time points
posL4=1+( lindist(seg4)-min(lindist(seg4)) ) / ...    % posL4: segment 4 % dist travelled
        ( max(lindist(seg4))-min(lindist(seg4)) )  ;
posL4b = [];
for k=1:length(posL4)
    if posT4(k)>=2707.462 && posT4(k)<=2765.212
        posL4b(k)=posL4(k)+2*(3-posL4(k)); %3=6
    elseif posT4(k)>=2908.879 && posT4(k)<=2965.08
        posL4b(k)=posL4(k)+2*(3-posL4(k));
    elseif posT4(k)>=3040.688 && posT4(k)<=3091.307
        posL4b(k)=posL4(k)+2*(3-posL4(k));
    elseif posT4(k)>=3193.29 && posT4(k)<=3249.698
        posL4b(k)=posL4(k)+2*(3-posL4(k));
    elseif posT4(k)>=3347.402 && posT4(k)<=3405
        posL4b(k)=posL4(k)+2*(3-posL4(k));
    else
        posL4b(k)=posL4(k);
    end
end

    % Segment 5: C > E
posT5=time(seg5);
posL5=2+(lindist(seg5)-min(lindist(seg5)))/(max(lindist(seg5))-min(lindist(seg5)));
posL5b = [];
for k=1:length(posL5)
    if posT5(k)>=2707.462 && posT5(k)<=2765.212
        posL5b(k)=posL5(k)+2*(3-posL5(k));
    elseif posT5(k)>=2908.879 && posT5(k)<=2965.08
        posL5b(k)=posL5(k)+2*(3-posL5(k));
    elseif posT5(k)>=3040.688 && posT5(k)<=3091.307
        posL5b(k)=posL5(k)+2*(3-posL5(k));
    elseif posT5(k)>=3193.29 && posT5(k)<=3249.698
        posL5b(k)=posL5(k)+2*(3-posL5(k));
    elseif posT5(k)>=3347.402 && posT5(k)<=3405
        posL5b(k)=posL5(k)+2*(3-posL5(k));
    else
        posL5b(k)=posL5(k);
    end
end


%%% 0-ABD

    % Segment 2: A > B
posT2=time(seg2);
posL2=1+(lindist(seg2)-min(lindist(seg2)))/(max(lindist(seg2))-min(lindist(seg2)));
posL2b = [];
for k=1:length(posL2)
    if posT2(k)>=2580.971 && posT2(k)<=2680.835
        posL2b(k)=-(posL2(k)+2*(3-posL2(k))); %-3=-6
    elseif posT2(k)>=2790.035 && posT2(k)<=2849.932
        posL2b(k)=-(posL2(k)+2*(3-posL2(k)));
    elseif posT2(k)>=2858.031 && posT2(k)<=2888.331
        posL2b(k)=-(posL2(k)+2*(3-posL2(k)));
    elseif posT2(k)>=2991.458 && posT2(k)<=3025.768
        posL2b(k)=-(posL2(k)+2*(3-posL2(k)));
    elseif posT2(k)>=3112.713 && posT2(k)<=3152.636
        posL2b(k)=-(posL2(k)+2*(3-posL2(k)));
    elseif posT2(k)>=3274.19 && posT2(k)<=3324.589
        posL2b(k)=-(posL2(k)+2*(3-posL2(k)));
    else
        posL2b(k)=-posL2(k);
    end
end
    

    % Segment 3: B > D
posT3=time(seg3);
posL3=2+(lindist(seg3)-min(lindist(seg3)))/(max(lindist(seg3))-min(lindist(seg3)));
posL3b = [];
for k=1:length(posL3)
    if posT3(k)>=2580.971 && posT3(k)<=2680.835
        posL3b(k)=-(posL3(k)+2*(3-posL3(k)));
    elseif posT3(k)>=2790.035 && posT3(k)<=2849.932
        posL3b(k)=-(posL3(k)+2*(3-posL3(k)));
    elseif posT3(k)>=2858.031 && posT3(k)<=2888.331
        posL3b(k)=-(posL3(k)+2*(3-posL3(k)));
    elseif posT3(k)>=2991.458 && posT3(k)<=3025.768
        posL3b(k)=-(posL3(k)+2*(3-posL3(k)));
    elseif posT3(k)>=3112.713 && posT3(k)<=3152.636
        posL3b(k)=-(posL3(k)+2*(3-posL3(k)));
    elseif posT3(k)>=3274.19 && posT3(k)<=3324.589
        posL3b(k)=-(posL3(k)+2*(3-posL3(k)));
    else
        posL3b(k)=-posL3(k);
    end
end


vecL0 = [posL1b  posL4b  posL5b  posL2b  posL3b];       % concatenated directional ([-6 +6]) assignments
vecT  = [posT1 ; posT4 ; posT5 ; posT2 ; posT3];        % concatenated timestamps of each segment



%%%%%

if 0
    % what unit selection is this?
    selected_tets=[1 2 4 5 7 11 12 13 14 17 18 19 20 29];                                  % Bond tetrodes w/ clusters (but not all of them?)
    a=[1 1  2 4 5 5 5 11 11 11 11 11 12 12 12 13 13 14 14 14 17 18 19 19 20 29];   % unit tetnums
    b=[1 10 4 1 1 2 3 1  2  3  4  5  1  2  3  1  4   1  3 5  1  1  1  2  3  4];    % unit cellnums
elseif 1
    %all CA1 CA3
    selected_tets=[1 2 4 5 7 10 11 12 13 14 17 18 19 20 22 23 27 29];
    a=[1 1 1 1  1 1 2 2   2 4 4 5  5 5 5 7  10 10 10 10 11  11 11 11 11 12  12 12 13 13 13  13 14 14 14 14  14 14 17 18  18 18 19 19  19 20 20 20  20 20 22 22  23 23 27 29];
    b=[1 2 4 8  9 10 1 2  5 1 2 1  2 3 4 1  1 2 3 4 1       2 3 4 5 1       2 3 1 2 3       4 1 2 3 4         5 7 1 1     5 6 1 2        3 1 2 3   6 8 1 2     1 2 1 4];
end

%%%%%

% This processes vecL0 further in the framework

if 1
    figure
    plot(vecL0,'r','linewidth',2);
    hold on
    
end

for i=1:size(vecL0,2)
    if vecT(i,1)>=2540&&vecT(i)<=2550&&vecL0(i)>0.5
        vecL0(i)=-vecL0(i);
    elseif vecT(i)>=2610&&vecT(i)<=2620&&vecL0(i)>0.5
        vecL0(i)=-6+vecL0(i);
    elseif vecT(i)>=2670&&vecT(i)<=2690&&vecL0(i)<-0.5&&vecL0(i)>-2
        vecL0(i)=-vecL0(i);
    elseif vecT(i)>=2670&&vecT(i)<=2690&&vecL0(i)<-4&&vecL0(i)>-6
        vecL0(i)=6+vecL0(i);
    elseif vecT(i)>=2720&&vecT(i)<=2730&&vecL0(i)>0&&vecL0(i)<2
        vecL0(i)=6-vecL0(i);
    elseif vecT(i)>=2720&&vecT(i)<=2730&&vecL0(i)>-2&&vecL0(i)<0
        vecL0(i)=6+vecL0(i);
    elseif vecT(i)>=2760&&vecT(i)<=2770&&vecL0(i)>0.5
        vecL0(i)=-vecL0(i);
    elseif vecT(i)>=2810&&vecT(i)<=2820&&vecL0(i)>0.5
        vecL0(i)=-6+vecL0(i);
    elseif vecT(i)>=2845&&vecT(i)<=2855&&vecL0(i)>0.5&&vecL0(i)<2
        vecL0(i)=-vecL0(i);
    elseif vecT(i)>=2845&&vecT(i)<=2855&&vecL0(i)<-4&&vecL0(i)>-6
        vecL0(i)=-6-vecL0(i);    
    elseif vecT(i)>=2925&&vecT(i)<=2935&&vecL0(i)<-0.5
        vecL0(i)=6+vecL0(i);
    elseif vecT(i)>=2960&&vecT(i)<=2970&&vecL0(i)>0.5
        vecL0(i)=-vecL0(i);
    elseif vecT(i)>=3055&&vecT(i)<=3065&&vecL0(i)<-0.5
        vecL0(i)=6+vecL0(i);
    elseif vecT(i)>=3090&&vecT(i)<=3100&&vecL0(i)>0.5
        vecL0(i)=-vecL0(i);
    elseif vecT(i)>=3120&&vecT(i)<=3130&&vecL0(i)>0.5
        vecL0(i)=-6+vecL0(i);
    elseif vecT(i)>=3145&&vecT(i)<=3155&&vecL0(i)<-0.5&&vecL0(i)>-2
        vecL0(i)=-vecL0(i);
    elseif vecT(i)>=3145&&vecT(i)<=3155&&vecL0(i)<-4&&vecL0(i)>-6
        vecL0(i)=6+vecL0(i);
    elseif vecT(i)>=3220&&vecT(i)<=3230&&vecL0(i)<-0.5
        vecL0(i)=6+vecL0(i);
    elseif vecT(i)>=3245&&vecT(i)<=3255&&vecL0(i)<2&&vecL0(i)>0.5
        vecL0(i)=-vecL0(i);
    elseif vecT(i)>=3245&&vecT(i)<=3255&&vecL0(i)<6&&vecL0(i)>4
        vecL0(i)=-6+vecL0(i);
    elseif vecT(i)>=3290&&vecT(i)<=3300&&vecL0(i)>0.5
        vecL0(i)=-6+vecL0(i);
    elseif vecT(i)>=3320&&vecT(i)<=3330&&vecL0(i)<-0.5&&vecL0(i)>-2
        vecL0(i)=-vecL0(i);
    elseif vecT(i)>=3320&&vecT(i)<=3330&&vecL0(i)<-4&&vecL0(i)>-6
        vecL0(i)=6+vecL0(i);
    elseif vecT(i)>=3380&&vecT(i)<=3390&&vecL0(i)<-0.5
        vecL0(i)=6+vecL0(i);
    end
end

if 0
    plot(vecL0,'k','linewidth',2);
end




%%% pre-processing folding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vecL : "folded" version of vecL0, eliminating directionality
    % while maintaining L vs. R
vecL = [];
for i=1:length(vecL0)
    if vecL0(i)>=3
        vecL(i)=6-vecL0(i);
    elseif vecL0(i)<=-3
        vecL(i)=-6-vecL0(i);
    else
        vecL(i)=vecL0(i);
    end
end

% vecLF : col 1: ascending time; col 2: sorted linearized folded position
vecLF0 = [vecT vecL'];
[d1,d2] = sort(vecLF0(:,1));  
vecLF = vecLF0(d2,:); % 

A=pos{ex}{ep}.data(:,1);                    % time stamps for animal's trajectory
ti=round(A(1)*1000):1:round(A(end)*1000);   % t: 1-ms binned time stamps

