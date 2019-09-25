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
    
    clear data lp_c
    save temp
    datmat=datmat(31:end,:);
    datmat=datmat';
    awake1=datmat(1:end)';
    save ([name1 '_mat'], 'datmat','awake1') ;
    
    awake(lp_f,:)=awake1;
    clear datmat awake1
end

save temp
%save ([recname '_allvec' ], 'anes','awake');
save ('vecdec05','awake')







