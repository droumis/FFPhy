%% Run Filter
animals = {'Bond','Frank','Ten'};

%Define epochs
epochfilter = [];
epochfilter{1} = 'isequal($description,''TrackA'')';
epochfilter{2} = 'isequal($description,''TrackB'')';
% Iterator selection
iterator = 'epocheegnonreferenceanal';

%Define time filter
timefilter = {{'get2dstate', '($velocity < 4)'}};


% Tetrode selection
tetrodepairfilter = {'(isequal($area, ''CA1'') & ($maxcell == 1))', '(isequal($area, ''CA3'') & ($maxcell == 1))'};


% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'eegtetrodepairs', tetrodepairfilter, 'iterator', iterator);
f = setfilterfunction(f, 'calcriptheta', {'eeg','ripples','cellinfo'});

f = runfilter(f);

%% PLOT

%Look at theta power for the 400ms before each SWR
ca1 =[]; ca3 = []; ca1ca3 = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            ca1 = [ca1 f(an).output{d}(e).ca1_power];
            ca3 = [ca3 f(an).output{d}(e).ca3_power];
            ca1ca3 = [ca1ca3; f(an).output{d}(e).coherence];
        end
    end
end

% Plot distribution of theta power for the 400 ms before each SWR
figure
x = -5:0.1:5;
plot(x,hist(ca1,x)./length(ca1),'b',x,hist(ca3,x)./length(ca3),'r');
xlabel('Theta power (z-score)')
ylabel('Proportion of SWRs')
legend('CA1','CA3')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_theta_power.pdf', m, d, y);
print('-dpdf', savestring)


%% Determine whether theta power before hand has any relationship to replay pvalue
load '/data13/mcarr/RipplePaper/decodefilterA.mat'
load '/data13/mcarr/RipplePaper/decodefilterB.mat'

valid_ripples = cell(length(decodefilterA),1);

for an = 1:length(decodefilterA)
    for d = 1:length(decodefilterA(an).epochs)
        valid_ripples{an}{d} = cell(size(decodefilterA(an).epochs{d},1),1);
        for e = 1:length(decodefilterA(an).epochs{d})
            tmp_theta1 = f(an).output{d}(e).ca1_power';
            tmp_theta3 = f(an).output{d}(e).ca3_power';
            tmp_thetac = f(an).output{d}(e).coherence;
            tmp_starttime = f(an).output{d}(e).starttime;
            if ~isempty(decodefilterA(an).output{d}(e).pvalue) || ~isempty(decodefilterB(an).output{d}(e).pvalue)
               
                valid_ripples{an}{d}{e}.index = decodefilterA(an).epochs{d}(e,:);
                if ~isempty(decodefilterA(an).output{d}(e).pvalue)
                    
                    valid_theta = lookup(decodefilterA(an).output{d}(e).eventtime(:,1),tmp_starttime);
                    tmpA = [decodefilterA(an).output{d}(e).eventtime(:,1) ...
                        decodefilterA(an).output{d}(e).pvalue...
                        tmp_theta1(valid_theta)...
                        tmp_theta3(valid_theta)...
                        tmp_thetac(valid_theta)];
                else
                    valid_theta = lookup(decodefilterA(an).output{d}(e).eventtime(:,1),tmp_starttime);
                    tmpA = [decodefilterA(an).output{d}(e).eventtime(:,1) ...
                        nan(size(decodefilterA(an).output{d}(e).eventtime(:,1)))...
                        tmp_theta1(valid_theta)...
                        tmp_theta3(valid_theta)...
                        tmp_thetac(valid_theta)];
                end
                if ~isempty(decodefilterB(an).output{d}(e).pvalue)
                    
                    valid_theta = lookup(decodefilterB(an).output{d}(e).eventtime(:,1),tmp_starttime);
                    tmpB = [decodefilterB(an).output{d}(e).eventtime(:,1) ...
                        nan(size(decodefilterB(an).output{d}(e).eventtime(:,1)))...
                        tmp_theta1(valid_theta)...
                        tmp_theta3(valid_theta)...
                        tmp_thetac(valid_theta)];
                else
                    valid_theta = lookup(decodefilterB(an).output{d}(e).eventtime(:,1),tmp_starttime);
                    tmpB = [decodefilterB(an).output{d}(e).eventtime(:,1) ...
                        nan(size(decodefilterB(an).output{d}(e).eventtime(:,1)))...
                        tmp_theta1(valid_theta)...
                        tmp_theta3(valid_theta)...
                        tmp_thetac(valid_theta)];
                end
   
               
                %Find ripples that are valid for both, pick smaller pvalue
                if ~isempty(tmpA) && ~isempty(tmpB)
                    [c ia ib] = intersect(tmpA(:,1),tmpB(:,1));
                    valid_ripples{an}{d}{e}.time = c;
                    valid_ripples{an}{d}{e}.pvalue = min(tmpA(ia,2),tmpB(ib,2));
                    valid_ripples{an}{d}{e}.ca1theta = tmpA(ia,3);
                    valid_ripples{an}{d}{e}.ca3theta = tmpA(ia,4);
                    valid_ripples{an}{d}{e}.thetacoherence = tmpA(ia,5);
                    tmpA(ia,:) = []; tmpB(ib,:) = [];
                    
                    if ~isempty(tmpA) && ~isempty(tmpB)
                        %Find ripples that are valid on for only one
                        valid_ripples{an}{d}{e}.time = [valid_ripples{an}{d}{e}.time; tmpA(:,1); tmpB(:,1)];
                        valid_ripples{an}{d}{e}.pvalue = [valid_ripples{an}{d}{e}.pvalue; tmpA(:,2); tmpB(:,2)];
                        valid_ripples{an}{d}{e}.ca1theta = [valid_ripples{an}{d}{e}.ca1theta; tmpA(:,3); tmpB(:,3)];
                        valid_ripples{an}{d}{e}.ca3theta = [valid_ripples{an}{d}{e}.ca3theta; tmpA(:,4) tmpB(:,4)];
                        valid_ripples{an}{d}{e}.thetacoherence = [valid_ripples{an}{d}{e}.thetacoherence; tmpA(:,5) tmpB(:,5)];
                    end
                elseif isempty(tmpA) && ~isempty(tmpB)
                     valid_ripples{an}{d}{e}.time = tmpB(:,1);
                     valid_ripples{an}{d}{e}.pvalue = tmpB(:,2);
                     valid_ripples{an}{d}{e}.ca1theta = tmpB(:,3);
                     valid_ripples{an}{d}{e}.ca3theta = tmpB(:,4);
                     valid_ripples{an}{d}{e}.thetacoherence = tmpB(:,5);
                     
                 elseif isempty(tmpB) && ~isempty(tmpA)
                     valid_ripples{an}{d}{e}.time = tmpA(:,1);
                     valid_ripples{an}{d}{e}.pvalue = tmpA(:,2);
                     valid_ripples{an}{d}{e}.ca1theta = tmpA(:,3);
                     valid_ripples{an}{d}{e}.ca3theta = tmpA(:,4);
                     valid_ripples{an}{d}{e}.thetacoherence = tmpA(:,5);
                else
                    valid_ripples{an}{d}{e}.time = [];
                    valid_ripples{an}{d}{e}.pvalue = [];
                    valid_ripples{an}{d}{e}.ca1theta = [];
                    valid_ripples{an}{d}{e}.ca3theta = [];
                    valid_ripples{an}{d}{e}.thetacoherence = [];
                end
                
                %Get rid of NaN pvalues
                invalid = isnan(valid_ripples{an}{d}{e}.pvalue);
                valid_ripples{an}{d}{e}.time(invalid) = [];
                valid_ripples{an}{d}{e}.pvalue(invalid) = [];
                valid_ripples{an}{d}{e}.ca1theta(invalid) = [];
                valid_ripples{an}{d}{e}.ca3theta(invalid) = [];
                valid_ripples{an}{d}{e}.thetacoherence(invalid) = [];
                
                %Sort by ripple starttime
                [val idx] = sort(valid_ripples{an}{d}{e}.time);
                valid_ripples{an}{d}{e}.time = valid_ripples{an}{d}{e}.time(idx);
                valid_ripples{an}{d}{e}.pvalue = valid_ripples{an}{d}{e}.pvalue(idx);
                valid_ripples{an}{d}{e}.ca1theta = valid_ripples{an}{d}{e}.ca1theta(idx);
                valid_ripples{an}{d}{e}.ca3theta = valid_ripples{an}{d}{e}.ca3theta(idx);
                valid_ripples{an}{d}{e}.thetacoherence = valid_ripples{an}{d}{e}.thetacoherence(idx);
                
                clear tmpA tmpB ia ib val idx
            end
        end
    end
end

pvalue = [];
ca1theta = [];
ca3theta = [];
thetacoherence =[];
for an =1:length(valid_ripples)
    for d =1:length(valid_ripples{an})
        for e = 1:length(valid_ripples{an}{d})
            if ~isempty(valid_ripples{an}{d}{e})
                if isfield(valid_ripples{an}{d}{e},'ca1theta')
                    pvalue = [pvalue; valid_ripples{an}{d}{e}.pvalue];
                    ca1theta = [ca1theta; valid_ripples{an}{d}{e}.ca1theta];
                    ca3theta = [ca3theta; valid_ripples{an}{d}{e}.ca3theta];
                    thetacoherence = [thetacoherence; valid_ripples{an}{d}{e}.thetacoherence];
                end
            end
        end
    end
end

%400 ms before
%r1 = -0.06, p1 = 0.09
%r3 = 0.01, p3 = 0.8
%r = -0.03, p = 0.44

[r1 p1]=corrcoef(pvalue,ca1theta);
[r3 p3]=corrcoef(pvalue,ca3theta); 
[r p] = corrcoef(pvalue,thetacoherence);
