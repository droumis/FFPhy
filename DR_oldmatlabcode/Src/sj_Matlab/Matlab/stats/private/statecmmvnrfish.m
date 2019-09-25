function Fisher = statecmmvnrfish(Data, Design, Covar, Method, MatrixFormat)
%STATECMMVNRFISH Fisher information for multivariate normal regression model.
%    Fisher information matrix based on current maximum likelihood or
%    least-squares parameter estimates that account for missing data.

%    Copyright 2006 The MathWorks, Inc.
%    $Revision: 1.1.6.3 $ $Date: 2006/12/15 19:33:07 $

[NumSamples, NumSeries] = size(Data);
if iscell(Design)
    if (numel(Design) == 1)
        SingleDesign = true;
    else
        SingleDesign = false;
    end
    NumParams = size(Design{1},2);
else
    SingleDesign = false;
    NumParams = size(Design,2);
end

if ~all(size(Covar) == [NumSeries, NumSeries])
    error('stats:mvregresslike:InconsistentDims', ...
        'The covariance matrix SIGMA has wrong dimensions.');
else
    [CholCovar, CholState] = chol(Covar);
    if CholState > 0
        error('stats:mvregresslike:NonPosDefCov', ...
            'Covariance matrix SIGMA is not positive-definite.');
    end
end

% Step 2 - initialization
if strcmpi(MatrixFormat,'PARAMONLY')
    TotalParams = NumParams;
else
    TotalParams = NumParams + (NumSeries * (NumSeries + 1))/2;
end

Fisher = zeros(TotalParams,TotalParams);

% Step 3 - calculate fisher information matrix (not hessian)

if strcmpi(Method,'FISHER')
    
    % Step 4 - do partials wrt Mean

    Count = 0;
    TestMatrix = zeros(NumParams,NumParams);
    if iscell(Design)
        if SingleDesign
            A = CholCovar' \ Design{1};
            for k = 1:NumSamples
                TestMatrix = TestMatrix + A'*A;
            end
            Count = NumSamples;
        else
            for k = 1:NumSamples
                if ~all(isnan(Data(k,:)))
                    Count = Count + 1;
                    A = CholCovar' \ Design{k};
                    TestMatrix = TestMatrix + A'*A;
                end
            end
        end
    else
        for k = 1:NumSamples
            if ~all(isnan(Data(k,:)))
                Count = Count + 1;
                A = CholCovar' \ Design(k,:);
                TestMatrix = TestMatrix + A'*A;
            end
        end
    end
    TestMatrix = (1.0/Count) .* TestMatrix;
    Fisher(1:NumParams,1:NumParams) = TestMatrix;

    % Step 5 - do partials wrt Covar

    if strcmpi(MatrixFormat,'FULL')

        InvCovar = inv(Covar);
        
        GradC1 = zeros(NumSeries,NumSeries);
        GradC2 = zeros(NumSeries,NumSeries);

        i = NumParams;
        for i1 = 1:NumSeries
            for j1 = 1:i1
                i = i + 1;

                GradC1(i1,j1) = 1.0;                    % do dC/dtheta(i)
                GradC1(j1,i1) = 1.0;

                j = NumParams;
                for i2 = 1:NumSeries
                    for j2 = 1:i2
                        j = j + 1;

                        if (j <= i)
                            GradC2(i2,j2) = 1.0;        % do dC/dtheta(j)
                            GradC2(j2,i2) = 1.0;

                            Temp1 = InvCovar*GradC1;
                            Temp2 = InvCovar*GradC2;

                            Fisher(i,j) = 0.5*trace(Temp1*Temp2);
                            Fisher(j,i) = Fisher(i,j);                        

                            GradC2(i2,j2) = 0.0;        % undo dC/dtheta(j)
                            GradC2(j2,i2) = 0.0;
                        end
                    end
                end

                GradC1(i1,j1) = 0.0;                    % undo dC/dtheta(i)
                GradC1(j1,i1) = 0.0;
            end
        end
    end
    
    Fisher = Count * Fisher;
    
% Step 6 - calculate hessian (not fisher information matrix)

else
    % Step 7 - main loop over data records

    InvCovar = inv(Covar);

    Map = zeros(NumSeries,1);
    Count = 0;

    for kk = 1:NumSamples

        % Step 8 - determine and map available data in current record

        Map(:) = 0;
        Available = 0;
        for ii = 1:NumSeries
            if isnan(Data(kk,ii))
                Map(ii) = 0;
            else
                Map(ii) = 1;
                Available = Available + 1;
            end
        end

        if Available > 0                        % skip over empty records
            Count = Count + 1;

            % Step 9 - construct covariance matrix subarrays

            if iscell(Design)
                if SingleDesign
                    SubDesign = Design{1};
                else
                    SubDesign = Design{kk};
                end
            else
                SubDesign = Design(kk,:);
            end
            SubCovar = Covar;

            if Available < NumSeries
                for ii = NumSeries:-1:1
                    if Map(ii) == 0
                        SubDesign(ii,:) = [];
                        SubCovar(:,ii) = [];
                        SubCovar(ii,:) = [];
                    end
                end
                InvSubCovar = inv(SubCovar);
            else
                InvSubCovar = InvCovar;
            end

            % Step 10 - do partials wrt Mean for current data record

            TempMatrix = SubDesign' * InvSubCovar * SubDesign;
            Fisher(1:NumParams,1:NumParams) = ...
                Fisher(1:NumParams,1:NumParams) + TempMatrix;
            
            % Step 11 - do partials wrt Covar for current data record

            if strcmpi(MatrixFormat,'FULL')

                GradC1 = zeros(Available,Available);
                GradC2 = zeros(Available,Available);

                p1 = 0;
                i = NumParams;
                for i1 = 1:NumSeries
                    if Map(i1) > 0
                        p1 = p1 + 1;
                    end

                    q1 = 0;
                    for j1 = 1:i1
                        i = i + 1;

                        if Map(j1) > 0
                            q1 = q1 + 1;
                        end

                        if (Map(i1) > 0) && (Map(j1) > 0)
                            GradC1(p1,q1) = 1.0;
                            GradC1(q1,p1) = 1.0;
                        end

                        p2 = 0;
                        j = NumParams;
                        for i2 = 1:NumSeries
                            if Map(i2) > 0
                                p2 = p2 + 1;
                            end

                            q2 = 0;
                            for j2 = 1:i2
                                j = j + 1;

                                if Map(j2) > 0
                                    q2 = q2 + 1;
                                end

                                % dC/dtheta(i) = dC/dC(i1,j1)
                                % dC/dtheta(j) = dC/dC(i2,j2)

                                if (j <= i) && (Map(i1) > 0) && (Map(j1) > 0) && (Map(i2) > 0) && (Map(j2) > 0)
                                    GradC2(p2,q2) = 1.0;
                                    GradC2(q2,p2) = 1.0;

                                    Temp1 = InvSubCovar*GradC1;
                                    Temp2 = InvSubCovar*GradC2;

                                    Fisher(i,j) = Fisher(i,j) + 0.5*trace(Temp1*Temp2);
                                    Fisher(j,i) = Fisher(i,j);

                                    GradC2(p2,q2) = 0.0;
                                    GradC2(q2,p2) = 0.0;
                                end
                            end
                        end

                        if (Map(i1) > 0) && (Map(j1) > 0)   % undo dC/dtheta(i)
                            GradC1(p1,q1) = 0.0;
                            GradC1(q1,p1) = 0.0;
                        end
                    end
                end
            end
        end
    end
end
