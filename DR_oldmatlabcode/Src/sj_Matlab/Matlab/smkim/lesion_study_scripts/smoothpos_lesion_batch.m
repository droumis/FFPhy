function smoothpos_lesion_batch()

% define a Wtrack template and a corresponding image transform 
% for each day; apply to all subjects over both run1 and run2 of that day

subjects = {'M02','M03','M06'};
%subjects = {'M12','M13','M14','M16','M17','M19','M20','M22','M24','M25','M26'};

load('/data15/smkim/smoothpos_default_params.mat');

for s = 1:numel(subjects)
    disp(['processing ' subjects{s}]);
    load([subjects{s} '_Wtrack_registeredpos.mat']);
    for i = 1:10
        for j = 1:2
            try
                smoothedpos{i}{j} = smoothpos_lesion(registeredpos{i}{j},params);
            end
        end
    end
    try
        save([subjects{s} '_Wtrack_smoothedpos.mat'],'smoothedpos');
        clear('smoothedpos');
    end
end

