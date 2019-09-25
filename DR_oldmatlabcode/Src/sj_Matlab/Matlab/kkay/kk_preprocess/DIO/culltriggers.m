
cull = 0;

if cull
%% Using well position, classifies valid DIO triggers into valid (1) or invalid (0), in output structure.

%% specify these if want to later (see below) exclude erroneous points
clear output
leftEnd = 230;   %   X or Y-value of well position on left and right - check by plotting
rightEnd = 230;  %For w-track, leftEnd and rightEnd should be the same since wells will line up with 1 x-value
range = 10;
axis = 2;  % choose 2 if X-axis, choose 3 if Y-axis

for e=epochs
for t = 1:length(times)
    ind = lookup(times(t)/10000,pos{d}{e}.data(:,1));
    dist = abs([leftEnd rightEnd] - pos{d}{e}.data(ind,axis));    
    if dist(1) > range && dist(2) > range             %trigger time is spurious if pos is outside both well ranges
        output(t,:) = [t 0 times(t) ind];
    else
        output(t,:) = [t 1 times(t) ind];
    end
end
end

%% Plot the invalid triggers to visually check.

failInd = find(output(:,2) == 0);               %index of input/output times that do not match well x-value range
for e=epochs
for t=1:length(failInd)   %(length(times)-1)
    figure
    plot(pos{d}{e}.data(:,2),pos{d}{e}.data(:,3),'LineWidth',2,'Color',[0.8 0.8 0.8]);    % all paths in grey
    hold on
    pos_fail=output(failInd(t),4);         
    plot(pos{d}{e}.data(pos_fail,2),pos{d}{e}.data(pos_fail,3),'k*','LineWidth',3); %red dot=start time of interval
            string = sprintf('black/fail %d',output(failInd(t),3));
            title(string)
end
end

%% Remove invalid pulses from original <animalday>DIO.mat variables: diopulses, rawdio, DIO.

failPulses=output(failInd,3);        % list of pulsetimes to weed out

% DIO
for e=epochs
    for p=1:24
        if ~isempty(DIO{d}{e}{p})
            if ~isempty(DIO{d}{e}{p}.pulsetimes)
                q=length(DIO{d}{e}{p}.pulsetimes)
                for i=1:q
                    if sum(DIO{d}{e}{p}.pulsetimes(q-i+1)==failPulses)
                        DIO{d}{e}{p}.pulsetimes(q-i+1,:)=[];
                        DIO{d}{e}{p}.timesincelast(q-i+1)=[];
                        DIO{d}{e}{p}.pulselength(q-i+1)=[];
                        DIO{d}{e}{p}.pulseind(q-i+1)=[];
                    end
                end
                length(DIO{d}{e}{p}.pulsetimes)
            end
        end
        
        %diopulses
        if ~isempty(diopulses{p})
            if ~isempty(diopulses{p}.pulsetimes)
                q=length(diopulses{p}.pulsetimes)
                for j=1:q
                    if sum(diopulses{p}.pulsetimes(q-j+1)==failPulses)
                        diopulses{p}.pulsetimes(q-j+1,:)=[];
                        diopulses{p}.timesincelast(q-j+1)=[];
                        diopulses{p}.pulselength(q-j+1)=[];
                        diopulses{p}.pulseind(q-j+1)=[];
                    end
                end
                length(diopulses{p}.pulsetimes)
            end
        end
        
    end
    
    %rawdio
    if ~isempty(rawdio{d}{e})
        if ~isempty(rawdio{d}{e}.times)
            q=length(rawdio{d}{e}.times)
            for k=1:q
                if sum(rawdio{d}{e}.times(q-k+1)==failPulses)
                    rawdio{d}{e}.times(q-k+1)=[];
                    rawdio{d}{e}.values(q-k+1)=[];
                end
            end
            length(rawdio{d}{e}.times)
        end
    end
    
end




%save('davDIO06fixed.mat','diopulses','rawdio','DIO')


%save('davDIO06.mat','diopulses','rawdio','DIO')
end

% %% Remove erroneous trials
% 
% rewardtimes{d}{e};
% errortimes{d}{e};
% [initSize] = size(errortimes{d}{e});  %switch between reward and error times
% failInd = find(output(:,2) == 0);
% failTime = output(:,3);
% failTime = failTime(failInd);
% for i = 1:length(failTime)
%     [r c] = find(errortimes{d}{e} == failTime(i));
%     errortimes{d}{e}(r,:) = [];
% end;
% [finalSize] = size(errortimes{d}{e});
% errortimes{d}{e};
% fprintf('\n%d rows have been removed\n\n',initSize(1)-finalSize(1));
% 
% end
%     