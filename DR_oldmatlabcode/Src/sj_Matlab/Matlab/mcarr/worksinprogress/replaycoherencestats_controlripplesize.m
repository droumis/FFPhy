%% Run for coherence
load '/data21/mcarr/RipplePaper/decodefilterA.mat'
load '/data21/mcarr/RipplePaper/decodefilterB.mat'
load '/data21/mcarr/RipplePaper/decodefilterA_sleep.mat'
load '/data21/mcarr/RipplePaper/decodefilterB_sleep.mat'

valid_ripples = cell(length(decodefilterA),1);
for an = 1:length(decodefilterA)
    load(sprintf('%s%scellinfo.mat',decodefilterA(an).animal{2},decodefilterA(an).animal{3}))
    for d = 1:length(decodefilterA(an).epochs)
        valid_ripples{an}{d} = cell(size(decodefilterA(an).epochs{d},1),1);
        for e = 1:length(decodefilterA(an).epochs{d})
        	if ~isempty(decodefilterA(an).output{d}(e).pvalue) || ~isempty(decodefilterB(an).output{d}(e).pvalue)
                valid_ripples{an}{d}{e}.index = decodefilterA(an).epochs{d}(e,:);
                if ~isempty(decodefilterA(an).output{d}(e).pvalue)
                    numcells = zeros(size(decodefilterA(an).output{d}(e).eventtime,1),1);
                    numspikes = zeros(size(numcells));
                    for i = 1:length(numcells)
                        numcells(i) = length(unique(decodefilterA(an).output{d}(e).eventdata(i).cellindex));
                        numspikes(i) = length(decodefilterA(an).output{d}(e).eventdata(i).cellindex);
                    end
                    tmpA = [decodefilterA(an).output{d}(e).eventtime(:,1) ...
                        decodefilterA(an).output{d}(e).pvalue...
                        decodefilterA(an).output{d}(e).eventimmobiletime...
                        numcells numspikes];
                else
                    tmpA = [decodefilterA(an).output{d}(e).eventtime(:,1) ...
                        nan(size(decodefilterA(an).output{d}(e).eventtime(:,1)))...
                        decodefilterA(an).output{d}(e).eventimmobiletime...
                        nan(size(decodefilterA(an).output{d}(e).eventtime(:,1)))...
                        nan(size(decodefilterA(an).output{d}(e).eventtime(:,1)))];
                end
                if ~isempty(decodefilterB(an).output{d}(e).pvalue)
                    numcells = zeros(size(decodefilterB(an).output{d}(e).eventtime,1),1);
                    numspikes = zeros(size(numcells));
                    for i = 1:length(numcells)
                        numcells(i) = length(unique(decodefilterB(an).output{d}(e).eventdata(i).cellindex));
                        numspikes(i) = length(decodefilterB(an).output{d}(e).eventdata(i).cellindex);
                    end
                    tmpB = [decodefilterB(an).output{d}(e).eventtime(:,1) ...
                        decodefilterB(an).output{d}(e).pvalue...
                        decodefilterB(an).output{d}(e).eventimmobiletime...
                        numcells numspikes];
                else
                    tmpB = [decodefilterB(an).output{d}(e).eventtime(:,1) ...
                        nan(size(decodefilterB(an).output{d}(e).eventtime(:,1)))...
                        decodefilterB(an).output{d}(e).eventimmobiletime...
                        nan(size(decodefilterB(an).output{d}(e).eventtime(:,1)))...
                        nan(size(decodefilterB(an).output{d}(e).eventtime(:,1)))];
                end
   
                %Determine peak rate for every ripple for decodefilterA
                for event = 1:length(decodefilterA(an).output{d}(e).eventdata)
                    cell_index = decodefilterA(an).output{d}(e).index(decodefilterA(an).output{d}(e).eventdata(event).cellindex,:);
                    tmp = zeros(size(cell_index,1),1);
                    for i = 1:size(cell_index,1)
                        try ~isempty(cellinfo{cell_index(i,1)}{2}{cell_index(i,3)}{cell_index(i,4)}) & ~isempty(cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)});
                           
                            if isfield(cellinfo{cell_index(i,1)}{2}{cell_index(i,3)}{cell_index(i,4)},'peakrate') && isfield(cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)},'peakrate')
                                tmp(i) = max(cellinfo{cell_index(i,1)}{2}{cell_index(i,3)}{cell_index(i,4)}.peakrate,...
                                    cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)}.peakrate);
                            elseif isfield(cellinfo{cell_index(i,1)}{2}{cell_index(i,3)}{cell_index(i,4)},'peakrate') && ~isfield(cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)},'peakrate')
                                tmp(i) = cellinfo{cell_index(i,1)}{2}{cell_index(i,3)}{cell_index(i,4)}.peakrate;
                            elseif isfield(cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)},'peakrate') && ~isfield(cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)},'peakrate')
                                tmp(i) = cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)}.peakrate;
                            else
                                tmp(i) = NaN;
                            end
                        catch
                            try ~isempty(cellinfo{cell_index(i,1)}{2}{cell_index(i,3)}{cell_index(i,4)});
                                if  isfield(cellinfo{cell_index(i,1)}{2}{cell_index(i,3)}{cell_index(i,4)},'peakrate')
                                    tmp(i) = cellinfo{cell_index(i,1)}{2}{cell_index(i,3)}{cell_index(i,4)}.peakrate;
                                end
                            end
                            try ~isempty(cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)});
                                if isfield(cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)},'peakrate')
                                    tmp(i) = cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)}.peakrate;
                                end
                            catch
                                tmp(i) = NaN;
                            end
                        end

                    end
                tmpA(event,6) = sum(tmp>3)./length(tmp);
                end
                %Determine peak rate for every ripple for decodefilterB
                for event = 1:length(decodefilterB(an).output{d}(e).eventdata)
                    cell_index = decodefilterB(an).output{d}(e).index(decodefilterB(an).output{d}(e).eventdata(event).cellindex,:);
                    tmp = zeros(size(cell_index,1),1);
                    for i = 1:size(cell_index,1)
                        try ~isempty(cellinfo{cell_index(i,1)}{2}{cell_index(i,3)}{cell_index(i,4)}) & ~isempty(cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)});
                           
                            if isfield(cellinfo{cell_index(i,1)}{2}{cell_index(i,3)}{cell_index(i,4)},'peakrate') && isfield(cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)},'peakrate')
                                tmp(i) = max(cellinfo{cell_index(i,1)}{2}{cell_index(i,3)}{cell_index(i,4)}.peakrate,...
                                    cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)}.peakrate);
                            elseif isfield(cellinfo{cell_index(i,1)}{2}{cell_index(i,3)}{cell_index(i,4)},'peakrate') && ~isfield(cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)},'peakrate')
                                tmp(i) = cellinfo{cell_index(i,1)}{2}{cell_index(i,3)}{cell_index(i,4)}.peakrate;
                            elseif isfield(cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)},'peakrate') && ~isfield(cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)},'peakrate')
                                tmp(i) = cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)}.peakrate;
                            else
                                tmp(i) = NaN;
                            end
                        catch
                            try ~isempty(cellinfo{cell_index(i,1)}{2}{cell_index(i,3)}{cell_index(i,4)});
                                if  isfield(cellinfo{cell_index(i,1)}{2}{cell_index(i,3)}{cell_index(i,4)},'peakrate')
                                    tmp(i) = cellinfo{cell_index(i,1)}{2}{cell_index(i,3)}{cell_index(i,4)}.peakrate;
                                end
                            end
                            try ~isempty(cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)});
                                if isfield(cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)},'peakrate')
                                    tmp(i) = cellinfo{cell_index(i,1)}{4}{cell_index(i,3)}{cell_index(i,4)}.peakrate;
                                end
                            catch
                                tmp(i) = NaN;
                            end
                        end
                    end
                tmpB(event,6) = sum(tmp>3)./length(tmp);
                end
                
                %Find ripples that are valid for both, pick smaller pvalue
                if ~isempty(tmpA) && ~isempty(tmpB)
                    [c ia ib] = intersect(tmpA(:,1),tmpB(:,1));
                    valid_ripples{an}{d}{e}.time = c;
                    valid_ripples{an}{d}{e}.pvalue = min(tmpA(ia,2),tmpB(ib,2));
                    valid_ripples{an}{d}{e}.eventimmobiletime = tmpA(ia,3);
                    valid_ripples{an}{d}{e}.numcells = tmpA(ia,4);
                    valid_ripples{an}{d}{e}.numspikes = tmpA(ia,5);
                    valid_ripples{an}{d}{e}.placecells = tmpA(ia,6);
                    tmpA(ia,:) = []; tmpB(ib,:) = [];
                    
                    if ~isempty(tmpA) && ~isempty(tmpB)
                        %Find ripples that are valid on for only one
                        valid_ripples{an}{d}{e}.time = [valid_ripples{an}{d}{e}.time; tmpA(:,1); tmpB(:,1)];
                        valid_ripples{an}{d}{e}.pvalue = [valid_ripples{an}{d}{e}.pvalue; tmpA(:,2); tmpB(:,2)];
                        valid_ripples{an}{d}{e}.eventimmobiletime = [valid_ripples{an}{d}{e}.eventimmobiletime; tmpA(:,3); tmpB(:,3)];
                        valid_ripples{an}{d}{e}.numcells = [valid_ripples{an}{d}{e}.numcells; tmpA(:,4) tmpB(:,4)];
                        valid_ripples{an}{d}{e}.numspikes = [valid_ripples{an}{d}{e}.numspikes; tmpA(:,5) tmpB(:,5)];
                        valid_ripples{an}{d}{e}.placecells = [valid_ripples{an}{d}{e}.placecells; tmpA(:,6); tmpB(:,6)];
                    end
                elseif isempty(tmpA) && ~isempty(tmpB)
                     valid_ripples{an}{d}{e}.time = tmpB(:,1);
                     valid_ripples{an}{d}{e}.pvalue = tmpB(:,2);
                     valid_ripples{an}{d}{e}.eventimmobiletime = tmpB(:,3);
                     valid_ripples{an}{d}{e}.numcells = tmpB(:,4);
                     valid_ripples{an}{d}{e}.numspikes = tmpB(:,5);
                     valid_ripples{an}{d}{e}.placecells = tmpB(:,6);
                     
                 elseif isempty(tmpB) && ~isempty(tmpA)
                     valid_ripples{an}{d}{e}.time = tmpA(:,1);
                     valid_ripples{an}{d}{e}.pvalue = tmpA(:,2);
                     valid_ripples{an}{d}{e}.eventimmobiletime = tmpA(:,3);
                     valid_ripples{an}{d}{e}.numcells = tmpA(:,4);
                     valid_ripples{an}{d}{e}.numspikes = tmpA(:,5);
                     valid_ripples{an}{d}{e}.placecells = tmpA(:,6);
                     
                     
                else
                    valid_ripples{an}{d}{e}.time = [];
                    valid_ripples{an}{d}{e}.pvalue = [];
                    valid_ripples{an}{d}{e}.eventimmobiletime = [];
                    valid_ripples{an}{d}{e}.numcells = [];
                    valid_ripples{an}{d}{e}.numspikes = [];
                    valid_ripples{an}{d}{e}.placecells = [];
                end
                
                %Get rid of NaN pvalues
                invalid = isnan(valid_ripples{an}{d}{e}.pvalue);
                valid_ripples{an}{d}{e}.time(invalid) = [];
                valid_ripples{an}{d}{e}.pvalue(invalid) = [];
                valid_ripples{an}{d}{e}.eventimmobiletime(invalid) = [];
                valid_ripples{an}{d}{e}.numcells(invalid) = [];
                valid_ripples{an}{d}{e}.numspikes(invalid) = [];
                valid_ripples{an}{d}{e}.placecells(invalid) = [];   
                
                %Sort by ripple starttime
                [val idx] = sort(valid_ripples{an}{d}{e}.time);
                valid_ripples{an}{d}{e}.time = valid_ripples{an}{d}{e}.time(idx);
                valid_ripples{an}{d}{e}.pvalue = valid_ripples{an}{d}{e}.pvalue(idx);
                valid_ripples{an}{d}{e}.eventimmobiletime = valid_ripples{an}{d}{e}.eventimmobiletime(idx);
                valid_ripples{an}{d}{e}.numcells = valid_ripples{an}{d}{e}.numcells(idx);
                valid_ripples{an}{d}{e}.numspikes = valid_ripples{an}{d}{e}.numspikes(idx);
                valid_ripples{an}{d}{e}.placecells = valid_ripples{an}{d}{e}.placecells(idx);

                clear tmpA tmpB ia ib val idx                
            end
        end
    end
   
end

%Go through for each and load the rip structure
for an = 1:length(valid_ripples)
    for d = 1:length(valid_ripples{an})
        for e = 1:length(valid_ripples{an}{d})
            if ~isempty(valid_ripples{an}{d}{e})
            	day = valid_ripples{an}{d}{e}.index(1);
                epoch = valid_ripples{an}{d}{e}.index(2);
                loadstring = sprintf('%s%srip%02d.mat',decodefilterA(an).animal{2},decodefilterA(an).animal{3},day);
                load(loadstring)
                rip = rip{day}{epoch};
                ind = lookup(valid_ripples{an}{d}{e}.time,rip.starttime);
                %Get rid of ripples that are too close (1 second apart) together
                invalid = [100; diff(rip.starttime)];
                invalid = invalid<1;
                
                valid_ripples{an}{d}{e}.time(invalid(ind)) =[];
                valid_ripples{an}{d}{e}.pvalue(invalid(ind)) = [];
                valid_ripples{an}{d}{e}.eventimmobiletime(invalid(ind)) = [];
                valid_ripples{an}{d}{e}.numcells(invalid(ind)) = [];
                valid_ripples{an}{d}{e}.numspikes(invalid(ind)) = [];
                valid_ripples{an}{d}{e}.placecells(invalid(ind)) = [];
                
                ind = ind(~invalid(ind));
                valid_ripples{an}{d}{e}.eegtime = rip.time;
                valid_ripples{an}{d}{e}.ca1_ca3 = rip.ca1_ca3_gamma_coherence(:,ind);
                valid_ripples{an}{d}{e}.ripple_amp = max(rip.ca1_ripple_power(:,ind));
                valid_ripples{an}{d}{e}.duration = rip.endtime(ind)-rip.starttime(ind);
                
                clear rip day epoch loadstring
            end
        end
    end  
end

c13 = []; pvalue = []; ripple_amp = []; duration = [];
for an =1:length(valid_ripples)
    for d =1:length(valid_ripples{an})
        for e = 1:length(valid_ripples{an}{d})
            if ~isempty(valid_ripples{an}{d}{e})
                if isfield(valid_ripples{an}{d}{e},'ca1_ca3')
                    if size(valid_ripples{an}{d}{e}.pvalue,1) == size(valid_ripples{an}{d}{e}.ca1_ca3,2) 
                        if size(valid_ripples{an}{d}{e}.pvalue,1)== size(valid_ripples{an}{d}{e}.eventimmobiletime,1)
                        	c13 = [c13 valid_ripples{an}{d}{e}.ca1_ca3];
                            pvalue = [pvalue; valid_ripples{an}{d}{e}.pvalue];
                            ripple_amp =[ripple_amp valid_ripples{an}{d}{e}.ripple_amp];
                            duration = [duration; valid_ripples{an}{d}{e}.duration];
                        end
                    end
                end
            end
        end
    end
end

time = valid_ripples{1}{1}{1}.eegtime;

%% CONTROL FOR RIPPLE AMPLITUDE
sig_amp = ripple_amp(pvalue<0.05);non_amp = ripple_amp(pvalue>0.05);
sigc13 = c13(:,pvalue<0.05); nonc13 = c13(:,pvalue>0.05);
invalid = isnan(sig_amp); sig_amp(invalid) = []; sigc13(:,invalid) = [];
invalid = isnan(non_amp); non_amp(invalid) = []; nonc13(:,invalid) = [];
[val ind] = sort(sig_amp);
sig_amp = sig_amp(ind); sigc13 = sigc13(:,ind);
[val ind] = sort(non_amp);
non_amp = non_amp(ind); nonc13 = nonc13(:,ind);

count = 1;
while mean(sig_amp)>mean(non_amp)
    if mod(count,2)
        sig_amp(end) = []; sigc13(:,end) = [];
        count = count+1;
    else
        non_amp(1) = []; nonc13(:,1) = [];
        count = count+1;
    end
end

nboot = 100;

baseline = mean([mean(sigc13(1:6,:)) mean(nonc13(1:6,:))]);
x1 = sigc13 - baseline; x2 = nonc13 - baseline;
N1 = size(x1,2); N2 = size(x2,2);
q1 = zeros(length(time),nboot); q2 = zeros(size(q1));
for b = 1:nboot
    q1(:,b) = mean(x1(:,ceil(N1*rand(N1,1))),2);
    q2(:,b) = mean(x2(:,ceil(N2*rand(N2,1))),2);
end
sig_m = mean(q1,2); non_m = mean(q2,2);
sig_e = std(q1,[],2); non_e = std(q2,[],2);

figure
plot(time,sig_m,'r',time,non_m,'k')
legend([{'Significant replay'},{'Nonsignificant replay'}])
hold on
fill([time time(end:-1:1)],[sig_m+sig_e; sig_m(end:-1:1)-sig_e(end:-1:1)],'r','EdgeColor','none')
fill([time time(end:-1:1)],[non_m+non_e; non_m(end:-1:1)-non_e(end:-1:1)],'k','EdgeColor','none')
set(gca,'xtick',time(6:10:end),'xticklabel',-0.4:0.1:0.4,'ylim',[-0.05 0.15],'xlim',[-0.4 0.4],'ytick',-0.5:0.05:0.7)
ylabel('CA1-CA3 coherence')
xlabel('Time since ripple detection (s)')
box off

%Permutation test for significant > nonsignificant
Nx = length(sig_amp); Ny = length(non_amp);
X = sigc13; Y = nonc13;
Z = cat(2,X,Y);
q = mean(X,2)-mean(Y,2);

nperm = 100;
qperm = zeros(length(q),nperm);
for s = 1:nperm
    [tmp,i] = sort(rand(Nx+Ny,1));
    Zperm = Z(:,i);
    qperm(:,s) = mean(Zperm(:,1:Nx),2) - mean(Zperm(:,Nx+[1:Ny]),2);
end
qp = zeros(size(time));
for i = 1:length(time)
    qp(i) = mean(abs(qperm(i,:))>q(i));
end
%Significant is greater than nonsignificant for 100-230ms following SWR
%detection, p < 0.05 even when controlling for the amplitude SWR

%baseline = 0.58

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1ca3_coherence_controlamplitude.pdf', m, d, y);
print('-dpdf', savestring)

%% CONTROL FOR RIPPLE DURATION
sig_duration = duration(pvalue<0.05)';non_duration = duration(pvalue>0.05)';
sigc13 = c13(:,pvalue<0.05); nonc13 = c13(:,pvalue>0.05);
invalid = isnan(sig_duration);
sig_duration(invalid) = []; sigc13(:,invalid) = [];
invalid = isnan(non_duration);
non_duration(invalid) = []; nonc13(:,invalid) = [];
[val ind] = sort(sig_duration);
sig_duration = sig_duration(ind); sigc13 = sigc13(:,ind);
[val ind] = sort(non_duration);
non_duration = non_duration(ind); nonc13 = nonc13(:,ind);

count = 1;
while mean(sig_duration)>mean(non_duration)
    if mod(count,2)
        sig_duration(end) = []; sigc13(:,end) = [];
        count = count+1;
    else
        non_duration(1) = []; nonc13(:,1) = [];
        count = count+1;
    end
end

nboot = 100;

baseline = mean([mean(sigc13(1:6,:)) mean(nonc13(1:6,:))]);
x1 = sigc13 - baseline; x2 = nonc13 - baseline;
N1 = size(x1,2); N2 = size(x2,2);
q1 = zeros(length(time),nboot); q2 = zeros(size(q1));
for b = 1:nboot
    q1(:,b) = mean(x1(:,ceil(N1*rand(N1,1))),2);
    q2(:,b) = mean(x2(:,ceil(N2*rand(N2,1))),2);
end
sig_m = mean(q1,2); non_m = mean(q2,2);
sig_e = std(q1,[],2); non_e = std(q2,[],2);

figure
plot(time,sig_m,'r',time,non_m,'k')
legend([{'Significant replay'},{'Nonsignificant replay'}])
hold on
fill([time time(end:-1:1)],[sig_m+sig_e; sig_m(end:-1:1)-sig_e(end:-1:1)],'r','EdgeColor','none')
fill([time time(end:-1:1)],[non_m+non_e; non_m(end:-1:1)-non_e(end:-1:1)],'k','EdgeColor','none')
set(gca,'xtick',time(6:10:end),'xticklabel',-0.4:0.1:0.4,'ylim',[-0.05 0.15],'xlim',[-0.4 0.4],'ytick',-0.5:0.05:0.7)
ylabel('CA1-CA3 coherence')
xlabel('Time since ripple detection (s)')
box off

%Permutation test for significant > nonsignificant
Nx = length(sig_amp); Ny = length(non_amp);
X = sigc13; Y = nonc13;
Z = cat(2,X,Y);
q = mean(X,2)-mean(Y,2);

nperm = 100;
qperm = zeros(length(q),nperm);
for s = 1:nperm
    [tmp,i] = sort(rand(Nx+Ny,1));
    Zperm = Z(:,i);
    qperm(:,s) = mean(Zperm(:,1:Nx),2) - mean(Zperm(:,Nx+[1:Ny]),2);
end
qp = zeros(size(time));
for i = 1:length(time)
    qp(i) = mean(abs(qperm(i,:))>q(i));
end
%Significant is greater than nonsignificant for 70-240ms following SWR
%detection, p < 0.05 even when controlling for the amplitude
%of the ripple oscillation

%baseline = 0.58
% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1ca3_coherence_controlduration.pdf', m, d, y);
print('-dpdf', savestring)


%% Run for phase

%Go through for each and load the rip structure
for an = 1:length(valid_ripples)
    for d = 1:length(valid_ripples{an})
        for e = 1:length(valid_ripples{an}{d})
            if ~isempty(valid_ripples{an}{d}{e})
            	day = valid_ripples{an}{d}{e}.index(1);
                epoch = valid_ripples{an}{d}{e}.index(2);
                loadstring = sprintf('%s%sripc%02d.mat',decodefilterA(an).animal{2},decodefilterA(an).animal{3},day);
                load(loadstring)
                ripc = ripc{day}{epoch};
                ind = lookup(valid_ripples{an}{d}{e}.time,ripc.starttime);
                %Get rid of ripples that are too close (1 second apart) together
                invalid = [100; diff(ripc.starttime)];
                invalid = invalid<1;
                
                valid_ripples{an}{d}{e}.ca1_ca3 = ripc.ca1_ca3_phase(:,ind);
                clear ripc day epoch loadstring
            end
        end
    end  
end

c13 = []; pvalue = []; duration = []; ripple_amp = [];
for an =1:length(valid_ripples)
    for d =1:length(valid_ripples{an})
        for e = 1:length(valid_ripples{an}{d})
            if ~isempty(valid_ripples{an}{d}{e})
                if isfield(valid_ripples{an}{d}{e},'ca1_ca3')
                    if size(valid_ripples{an}{d}{e}.pvalue,1) == size(valid_ripples{an}{d}{e}.ca1_ca3,2) 
                        if size(valid_ripples{an}{d}{e}.pvalue,1)== size(valid_ripples{an}{d}{e}.eventimmobiletime,1)
                        	c13 = [c13 valid_ripples{an}{d}{e}.ca1_ca3];
                            pvalue = [pvalue; valid_ripples{an}{d}{e}.pvalue];
                            ripple_amp = [ripple_amp valid_ripples{an}{d}{e}.ripple_amp];
                            duration = [duration; valid_ripples{an}{d}{e}.duration];
                        end
                    end
                end
            end
        end
    end
end

time = valid_ripples{1}{1}{1}.eegtime;

%% CONTROL FOR RIPPLE AMPLITUDE
sig_amp = ripple_amp(pvalue<0.05);non_amp = ripple_amp(pvalue>0.05);
sigc13 = c13(:,pvalue<0.05); nonc13 = c13(:,pvalue>0.05);
invalid = isnan(sig_amp); sig_amp(invalid) = []; sigc13(:,invalid) = [];
invalid = isnan(non_amp); non_amp(invalid) = []; nonc13(:,invalid) = [];
[val ind] = sort(sig_amp);
sig_amp = sig_amp(ind); sigc13 = sigc13(:,ind);
[val ind] = sort(non_amp);
non_amp = non_amp(ind); nonc13 = nonc13(:,ind);

count = 1;
while mean(sig_amp)>mean(non_amp)
    if mod(count,2)
        sig_amp(end) = []; sigc13(:,end) = [];
        count = count+1;
    else
        non_amp(1) = []; nonc13(:,1) = [];
        count = count+1;
    end
end

nboot = 1000;
baseline = zeros(6,1);
for i = 1:6
    [m r] = anglemean([sigc13(i,:) nonc13(i,:)]);
    baseline(i) = r;
end
baseline = mean(baseline);

x1 = sigc13; x2 = nonc13;
N1 = size(x1,2); N2 = size(x2,2);
q1 = zeros(length(time),nboot); q2 = zeros(size(q1));
for b = 1:nboot
    for i = 1:length(time)
        [m r] = anglemean(x1(i,ceil(N1*rand(N1,1))));
        q1(i,b) = r;
        [m r] = anglemean(x2(i,ceil(N2*rand(N2,1))));
        q2(i,b) = r;
    end
end
sig_m = mean(q1,2) - baseline; non_m = mean(q2,2) - baseline;
sig_e = std(q1,[],2); non_e = std(q2,[],2);


figure
plot(time,sig_m,'r',time,non_m,'k')
legend([{'Significant replay'},{'Nonsignificant replay'}])
hold on
fill([time time(end:-1:1)],[sig_m+sig_e; sig_m(end:-1:1)-sig_e(end:-1:1)],'r','EdgeColor','none')
fill([time time(end:-1:1)],[non_m+non_e; non_m(end:-1:1)-non_e(end:-1:1)],'k','EdgeColor','none')
set(gca,'xtick',time(6:10:end),'xticklabel',-0.4:0.1:0.4,'ylim',[-0.05 0.15],'xlim',[-0.4 0.4],'ytick',-0.5:0.05:0.7)
ylabel('CA1-CA3 phase locking')
xlabel('Time since ripple detection (s)')
box off

%Permutation test for significant > nonsignificant
Nx = length(sig_amp); Ny = length(non_amp);
X = sigc13; Y = nonc13;
Z = cat(2,X,Y);
q = zeros(size(time,2),1);
for i = 1:length(time)
    [m rx] = anglemean(X(i,:));
    [m ry] = anglemean(Y(i,:));
    q(i) = rx-ry;
end
nperm = 1000;
qperm = zeros(length(q),nperm);
for s = 1:nperm
    [tmp,i] = sort(rand(Nx+Ny,1));
    Zperm = Z(:,i);
    for i = 1:length(time)
        [m rx] = anglemean(Zperm(i,1:Nx));
        [m ry] = anglemean(Zperm(i,Nx+[1:Ny]));
        qperm(i,s) = rx-ry;
    end
end
qp = zeros(size(time));
for i = 1:length(time)
    qp(i) = mean(abs(qperm(i,:))>q(i));
end

%Significant is greater than nonsignificant for 100-170ms following SWR
%detection, p < 0.05 even when controlling for ripple amplitude

% Baseline = 0.85

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1ca3_phase_controlrippleamp.pdf', m, d, y);
print('-dpdf', savestring)

%% CONTROL FOR RIPPLE DURATION
sig_duration = duration(pvalue<0.05)';non_duration = duration(pvalue>0.05)';
sigc13 = c13(:,pvalue<0.05); nonc13 = c13(:,pvalue>0.05);
invalid = isnan(sig_duration);
sig_duration(invalid) = []; sigc13(:,invalid) = [];
invalid = isnan(non_duration);
non_duration(invalid) = []; nonc13(:,invalid) = [];
[val ind] = sort(sig_duration);
sig_duration = sig_duration(ind); sigc13 = sigc13(:,ind);
[val ind] = sort(non_duration);
non_duration = non_duration(ind); nonc13 = nonc13(:,ind);

count = 1;
while mean(sig_duration)>mean(non_duration)
    if mod(count,2)
        sig_duration(end) = []; sigc13(:,end) = [];
        count = count+1;
    else
        non_duration(1) = []; nonc13(:,1) = [];
        count = count+1;
    end
end

nboot = 100;

baseline = zeros(6,1);
for i = 1:6
    [m r] = anglemean([sigc13(i,:) nonc13(i,:)]);
    baseline(i) = r;
end
baseline = mean(baseline);

x1 = sigc13; x2 = nonc13;
N1 = size(x1,2); N2 = size(x2,2);
q1 = zeros(length(time),nboot); q2 = zeros(size(q1));
for b = 1:nboot
    for i = 1:length(time)
        [m r] = anglemean(x1(i,ceil(N1*rand(N1,1))));
        q1(i,b) = r;
        [m r] = anglemean(x2(i,ceil(N2*rand(N2,1))));
        q2(i,b) = r;
    end
end
sig_m = mean(q1,2) - baseline; non_m = mean(q2,2) - baseline;
sig_e = std(q1,[],2); non_e = std(q2,[],2);


figure
plot(time,sig_m,'r',time,non_m,'k')
legend([{'Significant replay'},{'Nonsignificant replay'}])
hold on
fill([time time(end:-1:1)],[sig_m+sig_e; sig_m(end:-1:1)-sig_e(end:-1:1)],'r','EdgeColor','none')
fill([time time(end:-1:1)],[non_m+non_e; non_m(end:-1:1)-non_e(end:-1:1)],'k','EdgeColor','none')
set(gca,'xtick',time(6:10:end),'xticklabel',-0.4:0.1:0.4,'ylim',[-0.05 0.15],'xlim',[-0.4 0.4],'ytick',-0.5:0.05:0.7)
ylabel('CA1-CA3 phase locking')
xlabel('Time since ripple detection (s)')
box off

%Permutation test for significant > nonsignificant
Nx = length(sig_duration); Ny = length(non_duration);
X = sigc13; Y = nonc13;
Z = cat(2,X,Y);
q = zeros(size(time,2),1);
for i = 1:length(time)
    [m rx] = anglemean(X(i,:));
    [m ry] = anglemean(Y(i,:));
    q(i) = rx-ry;
end
nperm = 1000;
qperm = zeros(length(q),nperm);
for s = 1:nperm
    [tmp,i] = sort(rand(Nx+Ny,1));
    Zperm = Z(:,i);
    for i = 1:length(time)
        [m rx] = anglemean(Zperm(i,1:Nx));
        [m ry] = anglemean(Zperm(i,Nx+[1:Ny]));
        qperm(i,s) = rx-ry;
    end
end
qp = zeros(size(time));
for i = 1:length(time)
    qp(i) = mean(abs(qperm(i,:))>q(i));
end

%Significant is greater than nonsignificant for 80-170ms following SWR
%detection, p < 0.05 even when controlling for the number of spikes in each
%SWR and the number of cells participating

% Baseline = 0.85

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1ca3_phase_controlduration.pdf', m, d, y);
print('-dpdf', savestring)


