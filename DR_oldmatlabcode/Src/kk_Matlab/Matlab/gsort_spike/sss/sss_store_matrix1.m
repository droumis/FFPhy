function datmat = sss_store_matrix1(files, n)

% sss_store_matrix1({'MIC02d2_1_1.swp';'MIC02d2_2_1.swp';'MIC02d2_3_1.swp';'MIC02d2_4_1.swp';'MIC02d2_5_1.swp';'MIC02d2_6_1.swp';'MIC02d2_7_1.swp';'MIC02d2_8_1.swp'});
% sss_spkecut_1halfms_basic({'CR03p3r2DeL4';'CR03p3r3DeL4';'CR03p3r4DeL4';'CR03p3r5
% DeL4';'CR03p3r1DeL4'}, params, [0 1 0 1 100], 1, [1:1090], {[1:200 601:800]; [201:600 801:1000]; [1001:1090]});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off;

% Create variables to use later in the code
longsp=0; bsize=1; spikestarts=14;
sfreq=32; % sampling frequency in kHz.

recname = strtok(files{1},'.');

for lp_f=1:size(files,1)
    name1 = strtok(files{lp_f},'.');
    data=load ([files{lp_f}]);
    h = waitbar(0,'Reorganizing the data into matrix...');

    data = [-9999 data' -9999];

    cuts = find (data == -9999);
    %datmat = zeros (length(cuts)-1,(cuts(2)-cuts(1))+1); % checked
    datmat = zeros (length(cuts)-1,(cuts(2)-cuts(1))-7); % checked

    for lp_c = 1: (length (cuts)-1)
       % datmat (lp_c,:)= data (cuts(lp_c):cuts(lp_c+1)); %  checked. for an example data set:see test.xls
        datmat (lp_c,:)= data (cuts(lp_c)+7:cuts(lp_c+1)-1); %  TO GET RID OF -9999 and sweep info
        a=1;
    end

    % KEY FOR THE DATMAT
    %   (which has as many rows as the neumber of sweeps and as many columns as the (cuts(2)-cuts(1))+1))
    % 1st column    -9999
    % 2nd           Rate
    % 3rd           Number of samples
    % 4th           Sweep no
    % 5th           Trial no
    % 6th           Repetition no
    % 7th           Stimulus number
    % 8:N           Sweep data -raw-
    % N+1           -9999 Seperator
    
    clear data ndata h lp_c
    datmat_ane=datmat(1:51,:);
    datmat_awake=datmat(52:end,:);
 %   cmd=sprintf('d_ane_%d = datmat_ane;',lp_f); eval(cmd);
 %   cmd=sprintf('d_awake_%d = datmat_awake;',lp_f); eval(cmd);
    
    %cmd=sprintf('save([name], d_ane%d,d_awake%d);',lp_f,lp_f);
    save temp
    anes1=datmat_ane(1:end);
    awake1=datmat_awake(1:50,:);awake1=awake1(:);
    awake2=datmat_awake(51:101,:);awake2=awake2(:);
    
    save ([name1 '_mat'], 'datmat_ane','datmat_awake','anes1','awake1','awake2') ;
    
    anes(lp_f,:)=anes1;
    awake(lp_f,:)=awake2;
    clear datmat_ane datmat_awake anes awake1 awake2
end

save temp
%save ([recname '_allvec' ], 'anes','awake');
save (['vec1'], 'anes','awake')

%save ([recname '_allmat' ], 'd_ane1','d_awake1','d_ane2','d_awake2','d_ane3','d_awake3','d_ane4',...
%    'd_awake4','d_ane5','d_awake5','d_ane6','d_awake6','d_ane7','d_awake7','d_ane8','d_awake8');






