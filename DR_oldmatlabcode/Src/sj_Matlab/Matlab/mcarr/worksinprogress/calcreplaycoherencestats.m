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
                    tmpA = [decodefilterA(an).output{d}(e).eventtime(:,1) ...
                        decodefilterA(an).output{d}(e).pvalue...
                        decodefilterA(an).output{d}(e).rvalue...
                        decodefilterA(an).output{d}(e).entropy...
                        decodefilterA(an).output{d}(e).eventimmobiletime];
                else
                    tmpA = [decodefilterA(an).output{d}(e).eventtime(:,1) ...
                        nan(size(decodefilterA(an).output{d}(e).eventtime(:,1)))...
                        nan(size(decodefilterA(an).output{d}(e).eventtime(:,1)))...
                        nan(size(decodefilterA(an).output{d}(e).eventtime(:,1)))...
                        decodefilterA(an).output{d}(e).eventimmobiletime];
                end
                if ~isempty(decodefilterB(an).output{d}(e).pvalue)
                    tmpB = [decodefilterB(an).output{d}(e).eventtime(:,1) ...
                        decodefilterB(an).output{d}(e).pvalue...
                        decodefilterB(an).output{d}(e).rvalue...
                        decodefilterB(an).output{d}(e).entropy...
                        decodefilterB(an).output{d}(e).eventimmobiletime];
                else
                    tmpB = [decodefilterB(an).output{d}(e).eventtime(:,1) ...
                        nan(size(decodefilterB(an).output{d}(e).eventtime(:,1)))...
                        nan(size(decodefilterB(an).output{d}(e).eventtime(:,1)))...
                        nan(size(decodefilterB(an).output{d}(e).eventtime(:,1)))...
                        decodefilterB(an).output{d}(e).eventimmobiletime];
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
                    ind = min(tmpA(ia,2),tmpB(ib,2))==tmpB(ia,2);
                    ind = ind+1;
                    tmpR = [tmpA(ia,3) tmpB(ib,3)];
                    tmpE = [tmpA(ia,4) tmpB(ib,4)];
                    valid_ripples{an}{d}{e}.rvalue = zeros(size(valid_ripples{an}{d}{e}.pvalue));
                    valid_ripples{an}{d}{e}.entropy = zeros(size(valid_ripples{an}{d}{e}.pvalue));
                    for i = 1:length(ind)
                        valid_ripples{an}{d}{e}.rvalue(i) = tmpR(i,ind(i));
                        valid_ripples{an}{d}{e}.entropy(i) = tmpE(i,ind(i));
                    end
                    valid_ripples{an}{d}{e}.eventimmobiletime = tmpA(ia,5);
                    valid_ripples{an}{d}{e}.placecells = tmpA(ia,6);
                    tmpA(ia,:) = []; tmpB(ib,:) = [];
                
                    %Find ripples that are valid on for only one
                    valid_ripples{an}{d}{e}.time = [valid_ripples{an}{d}{e}.time; tmpA(:,1); tmpB(:,1)];
                    valid_ripples{an}{d}{e}.pvalue = [valid_ripples{an}{d}{e}.pvalue; tmpA(:,2); tmpB(:,2)];
                    valid_ripples{an}{d}{e}.rvalue = [valid_ripples{an}{d}{e}.rvalue; tmpA(:,3); tmpB(:,3)];
                    valid_ripples{an}{d}{e}.entropy = [valid_ripples{an}{d}{e}.entropy; tmpA(:,4); tmpB(:,4)];
                    valid_ripples{an}{d}{e}.eventimmobiletime = [valid_ripples{an}{d}{e}.eventimmobiletime; tmpA(:,5); tmpB(:,5)];
                    valid_ripples{an}{d}{e}.placecells = [valid_ripples{an}{d}{e}.placecells; tmpA(:,6); tmpB(:,6)];
                   
                elseif isempty(tmpA) && ~isempty(tmpB)
                     valid_ripples{an}{d}{e}.time = tmpB(:,1);
                     valid_ripples{an}{d}{e}.pvalue = tmpB(:,2);
                     valid_ripples{an}{d}{e}.rvalue = tmpB(:,3);
                     valid_ripples{an}{d}{e}.entropy = tmpB(:,4);
                     valid_ripples{an}{d}{e}.eventimmobiletime = tmpB(:,5);
                     valid_ripples{an}{d}{e}.placecells = tmpB(:,6);
                     
                 elseif isempty(tmpB) && ~isempty(tmpA)
                     valid_ripples{an}{d}{e}.time = tmpA(:,1);
                     valid_ripples{an}{d}{e}.pvalue = tmpA(:,2);
                     valid_ripples{an}{d}{e}.rvalue = tmpA(:,3);
                     valid_ripples{an}{d}{e}.entropy = tmpA(:,4);
                     valid_ripples{an}{d}{e}.eventimmobiletime = tmpA(:,5);
                     valid_ripples{an}{d}{e}.placecells = tmpA(:,6);
                     
                else
                    valid_ripples{an}{d}{e}.time = [];
                    valid_ripples{an}{d}{e}.pvalue = [];
                    valid_ripples{an}{d}{e}.rvalue = [];
                    valid_ripples{an}{d}{e}.entropy = [];
                    valid_ripples{an}{d}{e}.eventimmobiletime = [];
                    valid_ripples{an}{d}{e}.placecells = [];
                end
                
                %Get rid of NaN pvalues
                invalid = isnan(valid_ripples{an}{d}{e}.pvalue);
                valid_ripples{an}{d}{e}.time(invalid) = [];
                valid_ripples{an}{d}{e}.pvalue(invalid) = [];
                valid_ripples{an}{d}{e}.rvalue(invalid) = [];
                valid_ripples{an}{d}{e}.entropy(invalid) = [];
                valid_ripples{an}{d}{e}.eventimmobiletime(invalid) = [];
                valid_ripples{an}{d}{e}.placecells(invalid) = [];   
                
                %Sort by ripple starttime
                [val idx] = sort(valid_ripples{an}{d}{e}.time);
                valid_ripples{an}{d}{e}.time = valid_ripples{an}{d}{e}.time(idx);
                valid_ripples{an}{d}{e}.pvalue = valid_ripples{an}{d}{e}.pvalue(idx);
                valid_ripples{an}{d}{e}.rvalue = valid_ripples{an}{d}{e}.rvalue(idx);
                valid_ripples{an}{d}{e}.entropy = valid_ripples{an}{d}{e}.entropy(idx);
                valid_ripples{an}{d}{e}.eventimmobiletime = valid_ripples{an}{d}{e}.eventimmobiletime(idx);
                valid_ripples{an}{d}{e}.placecells = valid_ripples{an}{d}{e}.placecells(idx);

                clear tmpA tmpB ia ib val idx                
            end
        end
    end

    for d = 1:length(decodefilterA_sleep(an).epochs)
        valid_ripples{an}{end+1} = cell(size(decodefilterA_sleep(an).epochs{d},1),1);
        for e = 1:length(decodefilterA_sleep(an).epochs{d})
            if ~isempty(decodefilterA_sleep(an).output{d}(e).pvalue) || ~isempty(decodefilterB_sleep(an).output{d}(e).pvalue)
                valid_ripples{an}{end}{e}.index = decodefilterA_sleep(an).epochs{d}(e,:);
                if ~isempty(decodefilterA_sleep(an).output{d}(e).pvalue)
                    tmpA = [decodefilterA_sleep(an).output{d}(e).eventtime(:,1) ...
                        decodefilterA_sleep(an).output{d}(e).pvalue...
                        decodefilterA_sleep(an).output{d}(e).rvalue...
                        decodefilterA_sleep(an).output{d}(e).entropy...
                        decodefilterA_sleep(an).output{d}(e).eventimmobiletime];
                else
                    tmpA = [decodefilterA_sleep(an).output{d}(e).eventtime(:,1) ...
                        nan(size(decodefilterA_sleep(an).output{d}(e).eventtime(:,1)))...
                        nan(size(decodefilterA_sleep(an).output{d}(e).eventtime(:,1)))...
                        nan(size(decodefilterA_sleep(an).output{d}(e).eventtime(:,1)))...
                        decodefilterA_sleep(an).output{d}(e).eventimmobiletime];
                end
                if ~isempty(decodefilterB_sleep(an).output{d}(e).pvalue)
                    tmpB = [decodefilterB_sleep(an).output{d}(e).eventtime(:,1) ...
                        decodefilterB_sleep(an).output{d}(e).pvalue...
                        decodefilterB_sleep(an).output{d}(e).rvalue...
                        decodefilterB_sleep(an).output{d}(e).entropy...
                        decodefilterB_sleep(an).output{d}(e).eventimmobiletime];
                else
                    tmpB = [decodefilterB_sleep(an).output{d}(e).eventtime(:,1) ...
                        nan(size(decodefilterB_sleep(an).output{d}(e).eventtime(:,1)))...
                        nan(size(decodefilterB_sleep(an).output{d}(e).eventtime(:,1)))...
                        nan(size(decodefilterB_sleep(an).output{d}(e).eventtime(:,1)))...
                        decodefilterB_sleep(an).output{d}(e).eventimmobiletime];
                end
                
                %Determine peak rate for every ripple for decodefilterA
                for event = 1:length(decodefilterA_sleep(an).output{d}(e).eventdata)
                    cell_index = decodefilterA_sleep(an).output{d}(e).index(decodefilterA_sleep(an).output{d}(e).eventdata(event).cellindex,:);
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
                for event = 1:length(decodefilterB_sleep(an).output{d}(e).eventdata)
                    cell_index = decodefilterB_sleep(an).output{d}(e).index(decodefilterB_sleep(an).output{d}(e).eventdata(event).cellindex,:);
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
                    valid_ripples{an}{end}{e}.time = c;
                    valid_ripples{an}{end}{e}.pvalue = min(tmpA(ia,2),tmpB(ib,2));
                    ind = min(tmpA(ia,2),tmpB(ib,2))==tmpB(ia,2);
                    ind = ind+1;
                    tmpR = [tmpA(ia,3) tmpB(ib,3)];
                    tmpE = [tmpA(ia,4) tmpB(ib,4)];
                    valid_ripples{an}{end}{e}.rvalue = zeros(size(valid_ripples{an}{end}{e}.pvalue));
                    valid_ripples{an}{end}{e}.entropy = zeros(size(valid_ripples{an}{end}{e}.pvalue));
                    for i = 1:length(ind)
                        valid_ripples{an}{end}{e}.rvalue(i) = tmpR(i,ind(i));
                        valid_ripples{an}{end}{e}.entropy(i) = tmpE(i,ind(i));
                    end
                    valid_ripples{an}{end}{e}.eventimmobiletime = tmpA(ia,5);
                    valid_ripples{an}{end}{e}.placecells = tmpA(ia,6);
                    tmpA(ia,:) = []; tmpB(ib,:) = [];
                
                    %Find ripples that are valid on for only one
                    valid_ripples{an}{end}{e}.time = [valid_ripples{an}{end}{e}.time; tmpA(:,1); tmpB(:,1)];
                    valid_ripples{an}{end}{e}.pvalue = [valid_ripples{an}{end}{e}.pvalue; tmpA(:,2); tmpB(:,2)];
                    valid_ripples{an}{end}{e}.rvalue = [valid_ripples{an}{end}{e}.rvalue; tmpA(:,3); tmpB(:,3)];
                    valid_ripples{an}{end}{e}.entropy = [valid_ripples{an}{end}{e}.entropy; tmpA(:,4); tmpB(:,4)];
                    valid_ripples{an}{end}{e}.eventimmobiletime = [valid_ripples{an}{end}{e}.eventimmobiletime; tmpA(:,5); tmpB(:,5)];
                    valid_ripples{an}{end}{e}.placecells = [valid_ripples{an}{end}{e}.placecells; tmpA(:,6); tmpB(:,6)];
                   
                elseif isempty(tmpA) && ~isempty(tmpB)
                     valid_ripples{an}{end}{e}.time = tmpB(:,1);
                     valid_ripples{an}{end}{e}.pvalue = tmpB(:,2);
                     valid_ripples{an}{end}{e}.rvalue = tmpB(:,3);
                     valid_ripples{an}{end}{e}.entropy = tmpB(:,4);
                     valid_ripples{an}{end}{e}.eventimmobiletime = tmpB(:,5);
                     valid_ripples{an}{end}{e}.placecells = tmpB(:,6);
                     
                 elseif isempty(tmpB) && ~isempty(tmpA)
                     valid_ripples{an}{end}{e}.time = tmpA(:,1);
                     valid_ripples{an}{end}{e}.pvalue = tmpA(:,2);
                     valid_ripples{an}{end}{e}.rvalue = tmpA(:,3);
                     valid_ripples{an}{end}{e}.entropy = tmpA(:,4);
                     valid_ripples{an}{end}{e}.eventimmobiletime = tmpA(:,5);
                     valid_ripples{an}{end}{e}.placecells = tmpA(:,6);
                     
                else
                    valid_ripples{an}{end}{e}.time = [];
                    valid_ripples{an}{end}{e}.pvalue = [];
                    valid_ripples{an}{end}{e}.rvalue = [];
                    valid_ripples{an}{end}{e}.entropy = [];
                    valid_ripples{an}{end}{e}.eventimmobiletime = [];
                    valid_ripples{an}{end}{e}.placecells = [];
                end
                
                %Get rid of NaN pvalues
                invalid = isnan(valid_ripples{an}{end}{e}.pvalue);
                valid_ripples{an}{end}{e}.time(invalid) = [];
                valid_ripples{an}{end}{e}.pvalue(invalid) = [];
                valid_ripples{an}{end}{e}.rvalue(invalid) = [];
                valid_ripples{an}{end}{e}.entropy(invalid) = [];
                valid_ripples{an}{end}{e}.eventimmobiletime(invalid) = [];
                valid_ripples{an}{end}{e}.placecells(invalid) = [];   
                
                %Sort by ripple starttime
                [val idx] = sort(valid_ripples{an}{end}{e}.time);
                valid_ripples{an}{end}{e}.time = valid_ripples{an}{end}{e}.time(idx);
                valid_ripples{an}{end}{e}.pvalue = valid_ripples{an}{end}{e}.pvalue(idx);
                valid_ripples{an}{end}{e}.rvalue = valid_ripples{an}{end}{e}.rvalue(idx);
                valid_ripples{an}{end}{e}.entropy = valid_ripples{an}{end}{e}.entropy(idx);
                valid_ripples{an}{end}{e}.eventimmobiletime = valid_ripples{an}{end}{e}.eventimmobiletime(idx);
                valid_ripples{an}{end}{e}.placecells = valid_ripples{an}{end}{e}.placecells(idx);

                clear tmpA tmpB ia ib val idx                
            end
        end
    end
end

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
                valid_ripples{an}{d}{e}.rvalue(invalid(ind)) = [];
                valid_ripples{an}{d}{e}.entropy(invalid(ind)) = [];
                valid_ripples{an}{d}{e}.eventimmobiletime(invalid(ind)) = [];
                valid_ripples{an}{d}{e}.placecells(invalid(ind)) = [];
                
                ind = ind(~invalid(ind));
                valid_ripples{an}{d}{e}.eegtime = rip.time;
                valid_ripples{an}{d}{e}.ca1_power = rip.ca1_gamma_power(:,ind);
                valid_ripples{an}{d}{e}.ca3_power = rip.ca3_gamma_power(:,ind);

                %valid_ripples{an}{d}{e}.ca1_ca1 = rip.ca1_ca1_gamma_coherence(:,ind);
                valid_ripples{an}{d}{e}.ca1_ca3 = rip.ca1_ca3_gamma_coherence(:,ind);
                %valid_ripples{an}{d}{e}.ca3_ca3 = rip.ca3_ca3_gamma_coherence(:,ind);
                
               
                clear rip day epoch loadstring
            end
        end
    end  
end

%% Save valid_ripples: coherence
savestring = '/data21/mcarr/RipplePaper/replaycoherence.mat';
save(savestring,'valid_ripples')

%% Compare significant and non significant replay

%Initialize variables
c13 = [];
pvalue = [];
place = [];
immobile = [];
rvalue = [];
entropy = [];
for an =1:length(valid_ripples)
    for d =1:2%length(valid_ripples{an})
        for e = 1:length(valid_ripples{an}{d})
            if ~isempty(valid_ripples{an}{d}{e})
                if isfield(valid_ripples{an}{d}{e},'ca1_ca3')
                    if size(valid_ripples{an}{d}{e}.pvalue,1) == size(valid_ripples{an}{d}{e}.ca1_ca3,2) 
                        if size(valid_ripples{an}{d}{e}.pvalue,1)== size(valid_ripples{an}{d}{e}.eventimmobiletime,1)
                        	c13 = [c13 valid_ripples{an}{d}{e}.ca1_ca3];
                            pvalue = [pvalue; valid_ripples{an}{d}{e}.pvalue];
                            rvalue = [rvalue; valid_ripples{an}{d}{e}.rvalue];
                            entropy = [entropy; valid_ripples{an}{d}{e}.entropy];
                            immobile = [immobile; valid_ripples{an}{d}{e}.eventimmobiletime];
                            place = [place; valid_ripples{an}{d}{e}.placecells];
                        end
                    end
                end
            end
        end
    end
end

%Compare significant vs. nonsignificant coherence
time = valid_ripples{1}{1}{1}.eegtime;

valid = place>0;
%valid = immobile >60 & place>0.95;
nboot = 100;

baseline = mean(mean(c13(1:6,:)));
c13 = c13-baseline;
x1 = c13(:,valid&pvalue<0.05); x2 = c13(:,valid&pvalue>0.05);
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

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1ca3_coherence_sleep_sig_nonsig_replay.pdf', m, d, y);
print('-dpdf', savestring)

%Mean baseline for sleep sessions: 0.5969
%Mean baseline for run sessions: 0.578
%Immobility is correlated with baseline, Spearman rho =0.3147

%Permutation test for significant > nonsignificant
Nx = sum(valid&pvalue<0.05); Ny = sum(valid&pvalue>0.05);
X = c13(:,valid&pvalue<0.05); Y = c13(:,valid&pvalue>0.05);
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
%Significant is greater than nonsignificant for 100-250ms following SWR
%detection, p < 0.05
%For sleeps, significant and nonsignificant are never different


%% Show pictures of coherence varying as a function of pvalue
time = valid_ripples{1}{1}{1}.eegtime;
subs = [0 0.05 0.1 0.5 1];
p = lookup(pvalue,subs,1);
a13 = zeros(length(subs),length(time));
valid = place>0;

for i = 1:size(c13,1)
    a13(:,i) = accumarray(p(valid),c13(i,valid),[length(subs) 1],@(x) mean(x));
end
g = gaussian(1.5,10);
for i = 1:size(a13,1)
   a13(i,:) = conv(a13(i,:),g,'same'); 
end

figure
imagesc(time,[],a13)
set(gca,'clim',[-0.05 0.15],'ytick',1:length(subs),'yticklabel',subs,'xlim',[-0.4 0.4])
colorbar('ytick',-0.5:0.05:.7)
xlabel('Time since SWR detection (sec)')
ylabel('Replay p-value')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1ca3_coherence_sleep_pvalue_replay.png', m, d, y);
print('-dpng', savestring)

%% Analyze Entropy

subs = prctile(entropy,[25 50 75 100]);
p = lookup(entropy,subs,1);
a13 = zeros(length(subs),length(time));
valid = place>0;

for i = 1:size(c13,1)
    a13(:,i) = accumarray(p(valid),c13(i,valid),[length(subs) 1],@(x) mean(x));
end
g = gaussian(1.5,10);
for i = 1:size(a13,1)
   a13(i,:) = conv(a13(i,:),g,'same'); 
end

figure
imagesc(time,[],a13)
set(gca,'clim',[-0.05 0.15],'ytick',1:length(subs),'yticklabel',subs,'xlim',[-0.4 0.4])
colorbar('ytick',-0.5:0.05:.7)
xlabel('Time since SWR detection (sec)')
ylabel('Rvalue')
box off

[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1ca3_coherence_entropy.png', m, d, y);
print('-dpng', savestring)



%% Analyze Rvlues
subs = prctile(rvalue,[25 50 75 100]);
p = lookup(rvalue,subs,1);
a13 = zeros(length(subs),length(time));
valid = place>0;

for i = 1:size(c13,1)
    a13(:,i) = accumarray(p(valid),c13(i,valid),[length(subs) 1],@(x) mean(x));
end
g = gaussian(1.5,10);
for i = 1:size(a13,1)
   a13(i,:) = conv(a13(i,:),g,'same'); 
end

figure
imagesc(time,[],a13)
set(gca,'clim',[-0.05 0.15],'ytick',1:length(subs),'yticklabel',subs,'xlim',[-0.4 0.4])
colorbar('ytick',-0.5:0.05:.7)
xlabel('Time since SWR detection (sec)')
ylabel('Rvalue')
box off

[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1ca3_coherence_rvalue.png', m, d, y);
print('-dpng', savestring)

%% Run for phase dispersion
load '/data13/mcarr/RipplePaper/decodefilterA.mat'
load '/data13/mcarr/RipplePaper/decodefilterB.mat'
load '/data13/mcarr/RipplePaper/decodefilterA_sleep.mat'
load '/data13/mcarr/RipplePaper/decodefilterB_sleep.mat'

valid_ripples = cell(length(decodefilterA),1);
for an = 1:length(decodefilterA)
    load(sprintf('%s%scellinfo.mat',decodefilterA(an).animal{2},decodefilterA(an).animal{3}))
    for d = 1:length(decodefilterA(an).epochs)
        valid_ripples{an}{d} = cell(size(decodefilterA(an).epochs{d},1),1);
        for e = 1:length(decodefilterA(an).epochs{d})
            if ~isempty(decodefilterA(an).output{d}(e).pvalue) || ~isempty(decodefilterB(an).output{d}(e).pvalue)
                valid_ripples{an}{d}{e}.index = decodefilterA(an).epochs{d}(e,:);
                if ~isempty(decodefilterA(an).output{d}(e).pvalue)
                    tmpA = [decodefilterA(an).output{d}(e).eventtime(:,1) ...
                        decodefilterA(an).output{d}(e).pvalue...
                        decodefilterA(an).output{d}(e).rvalue...
                        decodefilterA(an).output{d}(e).entropy...
                        decodefilterA(an).output{d}(e).eventimmobiletime];
                else
                    tmpA = [decodefilterA(an).output{d}(e).eventtime(:,1) ...
                        nan(size(decodefilterA(an).output{d}(e).eventtime(:,1)))...
                        nan(size(decodefilterA(an).output{d}(e).eventtime(:,1)))...
                        nan(size(decodefilterA(an).output{d}(e).eventtime(:,1)))...
                        decodefilterA(an).output{d}(e).eventimmobiletime];
                end
                if ~isempty(decodefilterB(an).output{d}(e).pvalue)
                    tmpB = [decodefilterB(an).output{d}(e).eventtime(:,1) ...
                        decodefilterB(an).output{d}(e).pvalue...
                        decodefilterB(an).output{d}(e).rvalue...
                        decodefilterB(an).output{d}(e).entropy...
                        decodefilterB(an).output{d}(e).eventimmobiletime];
                else
                    tmpB = [decodefilterB(an).output{d}(e).eventtime(:,1) ...
                        nan(size(decodefilterB(an).output{d}(e).eventtime(:,1)))...
                        nan(size(decodefilterB(an).output{d}(e).eventtime(:,1)))...
                        nan(size(decodefilterB(an).output{d}(e).eventtime(:,1)))...
                        decodefilterB(an).output{d}(e).eventimmobiletime];
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
                    ind = min(tmpA(ia,2),tmpB(ib,2))==tmpB(ia,2);
                    ind = ind+1;
                    tmpR = [tmpA(ia,3) tmpB(ib,3)];
                    tmpE = [tmpA(ia,4) tmpB(ib,4)];
                    valid_ripples{an}{d}{e}.rvalue = zeros(size(valid_ripples{an}{d}{e}.pvalue));
                    valid_ripples{an}{d}{e}.entropy = zeros(size(valid_ripples{an}{d}{e}.pvalue));
                    for i = 1:length(ind)
                        valid_ripples{an}{d}{e}.rvalue(i) = tmpR(i,ind(i));
                        valid_ripples{an}{d}{e}.entropy(i) = tmpE(i,ind(i));
                    end
                    valid_ripples{an}{d}{e}.eventimmobiletime = tmpA(ia,5);
                    valid_ripples{an}{d}{e}.placecells = tmpA(ia,6);
                    tmpA(ia,:) = []; tmpB(ib,:) = [];
                
                    %Find ripples that are valid on for only one
                    valid_ripples{an}{d}{e}.time = [valid_ripples{an}{d}{e}.time; tmpA(:,1); tmpB(:,1)];
                    valid_ripples{an}{d}{e}.pvalue = [valid_ripples{an}{d}{e}.pvalue; tmpA(:,2); tmpB(:,2)];
                    valid_ripples{an}{d}{e}.rvalue = [valid_ripples{an}{d}{e}.rvalue; tmpA(:,3); tmpB(:,3)];
                    valid_ripples{an}{d}{e}.entropy = [valid_ripples{an}{d}{e}.entropy; tmpA(:,4); tmpB(:,4)];
                    valid_ripples{an}{d}{e}.eventimmobiletime = [valid_ripples{an}{d}{e}.eventimmobiletime; tmpA(:,5); tmpB(:,5)];
                    valid_ripples{an}{d}{e}.placecells = [valid_ripples{an}{d}{e}.placecells; tmpA(:,6); tmpB(:,6)];
                   
                elseif isempty(tmpA) && ~isempty(tmpB)
                     valid_ripples{an}{d}{e}.time = tmpB(:,1);
                     valid_ripples{an}{d}{e}.pvalue = tmpB(:,2);
                     valid_ripples{an}{d}{e}.rvalue = tmpB(:,3);
                     valid_ripples{an}{d}{e}.entropy = tmpB(:,4);
                     valid_ripples{an}{d}{e}.eventimmobiletime = tmpB(:,5);
                     valid_ripples{an}{d}{e}.placecells = tmpB(:,6);
                     
                 elseif isempty(tmpB) && ~isempty(tmpA)
                     valid_ripples{an}{d}{e}.time = tmpA(:,1);
                     valid_ripples{an}{d}{e}.pvalue = tmpA(:,2);
                     valid_ripples{an}{d}{e}.rvalue = tmpA(:,3);
                     valid_ripples{an}{d}{e}.entropy = tmpA(:,4);
                     valid_ripples{an}{d}{e}.eventimmobiletime = tmpA(:,5);
                     valid_ripples{an}{d}{e}.placecells = tmpA(:,6);
                     
                else
                    valid_ripples{an}{d}{e}.time = [];
                    valid_ripples{an}{d}{e}.pvalue = [];
                    valid_ripples{an}{d}{e}.rvalue = [];
                    valid_ripples{an}{d}{e}.entropy = [];
                    valid_ripples{an}{d}{e}.eventimmobiletime = [];
                    valid_ripples{an}{d}{e}.placecells = [];
                end
                
                %Get rid of NaN pvalues
                invalid = isnan(valid_ripples{an}{d}{e}.pvalue);
                valid_ripples{an}{d}{e}.time(invalid) = [];
                valid_ripples{an}{d}{e}.pvalue(invalid) = [];
                valid_ripples{an}{d}{e}.rvalue(invalid) = [];
                valid_ripples{an}{d}{e}.entropy(invalid) = [];
                valid_ripples{an}{d}{e}.eventimmobiletime(invalid) = [];
                valid_ripples{an}{d}{e}.placecells(invalid) = [];   
                
                %Sort by ripple starttime
                [val idx] = sort(valid_ripples{an}{d}{e}.time);
                valid_ripples{an}{d}{e}.time = valid_ripples{an}{d}{e}.time(idx);
                valid_ripples{an}{d}{e}.pvalue = valid_ripples{an}{d}{e}.pvalue(idx);
                valid_ripples{an}{d}{e}.rvalue = valid_ripples{an}{d}{e}.rvalue(idx);
                valid_ripples{an}{d}{e}.entropy = valid_ripples{an}{d}{e}.entropy(idx);
                valid_ripples{an}{d}{e}.eventimmobiletime = valid_ripples{an}{d}{e}.eventimmobiletime(idx);
                valid_ripples{an}{d}{e}.placecells = valid_ripples{an}{d}{e}.placecells(idx);

                clear tmpA tmpB ia ib val idx                
            end
        end
    end

    for d = 1:length(decodefilterA_sleep(an).epochs)
        valid_ripples{an}{end+1} = cell(size(decodefilterA_sleep(an).epochs{d},1),1);
        for e = 1:length(decodefilterA_sleep(an).epochs{d})
            if ~isempty(decodefilterA_sleep(an).output{d}(e).pvalue) || ~isempty(decodefilterB_sleep(an).output{d}(e).pvalue)
                valid_ripples{an}{end}{e}.index = decodefilterA_sleep(an).epochs{d}(e,:);
                if ~isempty(decodefilterA_sleep(an).output{d}(e).pvalue)
                    tmpA = [decodefilterA_sleep(an).output{d}(e).eventtime(:,1) ...
                        decodefilterA_sleep(an).output{d}(e).pvalue...
                        decodefilterA_sleep(an).output{d}(e).rvalue...
                        decodefilterA_sleep(an).output{d}(e).entropy...
                        decodefilterA_sleep(an).output{d}(e).eventimmobiletime];
                else
                    tmpA = [decodefilterA_sleep(an).output{d}(e).eventtime(:,1) ...
                        nan(size(decodefilterA_sleep(an).output{d}(e).eventtime(:,1)))...
                        nan(size(decodefilterA_sleep(an).output{d}(e).eventtime(:,1)))...
                        nan(size(decodefilterA_sleep(an).output{d}(e).eventtime(:,1)))...
                        decodefilterA_sleep(an).output{d}(e).eventimmobiletime];
                end
                if ~isempty(decodefilterB_sleep(an).output{d}(e).pvalue)
                    tmpB = [decodefilterB_sleep(an).output{d}(e).eventtime(:,1) ...
                        decodefilterB_sleep(an).output{d}(e).pvalue...
                        decodefilterB_sleep(an).output{d}(e).rvalue...
                        decodefilterB_sleep(an).output{d}(e).entropy...
                        decodefilterB_sleep(an).output{d}(e).eventimmobiletime];
                else
                    tmpB = [decodefilterB_sleep(an).output{d}(e).eventtime(:,1) ...
                        nan(size(decodefilterB_sleep(an).output{d}(e).eventtime(:,1)))...
                        nan(size(decodefilterB_sleep(an).output{d}(e).eventtime(:,1)))...
                        nan(size(decodefilterB_sleep(an).output{d}(e).eventtime(:,1)))...
                        decodefilterB_sleep(an).output{d}(e).eventimmobiletime];
                end
                
                %Determine peak rate for every ripple for decodefilterA
                for event = 1:length(decodefilterA_sleep(an).output{d}(e).eventdata)
                    cell_index = decodefilterA_sleep(an).output{d}(e).index(decodefilterA_sleep(an).output{d}(e).eventdata(event).cellindex,:);
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
                for event = 1:length(decodefilterB_sleep(an).output{d}(e).eventdata)
                    cell_index = decodefilterB_sleep(an).output{d}(e).index(decodefilterB_sleep(an).output{d}(e).eventdata(event).cellindex,:);
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
                    valid_ripples{an}{end}{e}.time = c;
                    valid_ripples{an}{end}{e}.pvalue = min(tmpA(ia,2),tmpB(ib,2));
                    ind = min(tmpA(ia,2),tmpB(ib,2))==tmpB(ia,2);
                    ind = ind+1;
                    tmpR = [tmpA(ia,3) tmpB(ib,3)];
                    tmpE = [tmpA(ia,4) tmpB(ib,4)];
                    valid_ripples{an}{end}{e}.rvalue = zeros(size(valid_ripples{an}{end}{e}.pvalue));
                    valid_ripples{an}{end}{e}.entropy = zeros(size(valid_ripples{an}{end}{e}.pvalue));
                    for i = 1:length(ind)
                        valid_ripples{an}{end}{e}.rvalue(i) = tmpR(i,ind(i));
                        valid_ripples{an}{end}{e}.entropy(i) = tmpE(i,ind(i));
                    end
                    valid_ripples{an}{end}{e}.eventimmobiletime = tmpA(ia,5);
                    valid_ripples{an}{end}{e}.placecells = tmpA(ia,6);
                    tmpA(ia,:) = []; tmpB(ib,:) = [];
                
                    %Find ripples that are valid on for only one
                    valid_ripples{an}{end}{e}.time = [valid_ripples{an}{end}{e}.time; tmpA(:,1); tmpB(:,1)];
                    valid_ripples{an}{end}{e}.pvalue = [valid_ripples{an}{end}{e}.pvalue; tmpA(:,2); tmpB(:,2)];
                    valid_ripples{an}{end}{e}.rvalue = [valid_ripples{an}{end}{e}.rvalue; tmpA(:,3); tmpB(:,3)];
                    valid_ripples{an}{end}{e}.entropy = [valid_ripples{an}{end}{e}.entropy; tmpA(:,4); tmpB(:,4)];
                    valid_ripples{an}{end}{e}.eventimmobiletime = [valid_ripples{an}{end}{e}.eventimmobiletime; tmpA(:,5); tmpB(:,5)];
                    valid_ripples{an}{end}{e}.placecells = [valid_ripples{an}{end}{e}.placecells; tmpA(:,6); tmpB(:,6)];
                   
                elseif isempty(tmpA) && ~isempty(tmpB)
                     valid_ripples{an}{end}{e}.time = tmpB(:,1);
                     valid_ripples{an}{end}{e}.pvalue = tmpB(:,2);
                     valid_ripples{an}{end}{e}.rvalue = tmpB(:,3);
                     valid_ripples{an}{end}{e}.entropy = tmpB(:,4);
                     valid_ripples{an}{end}{e}.eventimmobiletime = tmpB(:,5);
                     valid_ripples{an}{end}{e}.placecells = tmpB(:,6);
                     
                 elseif isempty(tmpB) && ~isempty(tmpA)
                     valid_ripples{an}{end}{e}.time = tmpA(:,1);
                     valid_ripples{an}{end}{e}.pvalue = tmpA(:,2);
                     valid_ripples{an}{end}{e}.rvalue = tmpA(:,3);
                     valid_ripples{an}{end}{e}.entropy = tmpA(:,4);
                     valid_ripples{an}{end}{e}.eventimmobiletime = tmpA(:,5);
                     valid_ripples{an}{end}{e}.placecells = tmpA(:,6);
                     
                else
                    valid_ripples{an}{end}{e}.time = [];
                    valid_ripples{an}{end}{e}.pvalue = [];
                    valid_ripples{an}{end}{e}.rvalue = [];
                    valid_ripples{an}{end}{e}.entropy = [];
                    valid_ripples{an}{end}{e}.eventimmobiletime = [];
                    valid_ripples{an}{end}{e}.placecells = [];
                end
                
                %Get rid of NaN pvalues
                invalid = isnan(valid_ripples{an}{end}{e}.pvalue);
                valid_ripples{an}{end}{e}.time(invalid) = [];
                valid_ripples{an}{end}{e}.pvalue(invalid) = [];
                valid_ripples{an}{end}{e}.rvalue(invalid) = [];
                valid_ripples{an}{end}{e}.entropy(invalid) = [];
                valid_ripples{an}{end}{e}.eventimmobiletime(invalid) = [];
                valid_ripples{an}{end}{e}.placecells(invalid) = [];   
                
                %Sort by ripple starttime
                [val idx] = sort(valid_ripples{an}{end}{e}.time);
                valid_ripples{an}{end}{e}.time = valid_ripples{an}{end}{e}.time(idx);
                valid_ripples{an}{end}{e}.pvalue = valid_ripples{an}{end}{e}.pvalue(idx);
                valid_ripples{an}{end}{e}.rvalue = valid_ripples{an}{end}{e}.rvalue(idx);
                valid_ripples{an}{end}{e}.entropy = valid_ripples{an}{end}{e}.entropy(idx);
                valid_ripples{an}{end}{e}.eventimmobiletime = valid_ripples{an}{end}{e}.eventimmobiletime(idx);
                valid_ripples{an}{end}{e}.placecells = valid_ripples{an}{end}{e}.placecells(idx);

                clear tmpA tmpB ia ib val idx                
            end
        end
    end
end

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
                
                valid_ripples{an}{d}{e}.time(invalid(ind)) =[];
                valid_ripples{an}{d}{e}.pvalue(invalid(ind)) = [];
                valid_ripples{an}{d}{e}.rvalue(invalid(ind)) = [];
                valid_ripples{an}{d}{e}.entropy(invalid(ind)) = [];
                valid_ripples{an}{d}{e}.eventimmobiletime(invalid(ind)) = [];
                valid_ripples{an}{d}{e}.placecells(invalid(ind)) = [];
                
                ind = ind(~invalid(ind));
                valid_ripples{an}{d}{e}.eegtime = ripc.time;
                valid_ripples{an}{d}{e}.ca1_ca3 = ripc.ca1_ca3_phase(:,ind);
                clear ripc day epoch loadstring
            end
        end
    end  
end

%% Save valid_ripples: phase
savestring = '/data13/mcarr/RipplePaper/replayphase.mat';
save(savestring,'valid_ripples')


%Initialize variables
c13 = [];
p = [];
place = [];
immobile = [];
rvalue = [];
entropy = [];
for an =1:length(valid_ripples)
    for d =1:2%length(valid_ripples{an})
        for e = 1:length(valid_ripples{an}{d})
            if ~isempty(valid_ripples{an}{d}{e})
                if isfield(valid_ripples{an}{d}{e},'ca1_ca3')
                    if size(valid_ripples{an}{d}{e}.pvalue,1) == size(valid_ripples{an}{d}{e}.ca1_ca3,2) 
                        if size(valid_ripples{an}{d}{e}.pvalue,1)== size(valid_ripples{an}{d}{e}.eventimmobiletime,1)
                        	c13 = [c13 valid_ripples{an}{d}{e}.ca1_ca3];
                            p = [p; valid_ripples{an}{d}{e}.pvalue];
                            rvalue = [rvalue; valid_ripples{an}{d}{e}.rvalue];
                            entropy = [entropy; valid_ripples{an}{d}{e}.entropy];
                            immobile = [immobile; valid_ripples{an}{d}{e}.eventimmobiletime];
                            place = [place; valid_ripples{an}{d}{e}.placecells];
                        end
                    end
                end
            end
        end
    end
end


%Mean baseline for sleep sessions: 0.8266
%Mean baseline for run sessions: 0.8265
%Immobility is not correlated with baseline

%Run bootstrap resampling to compute confidence intervals
baseline = zeros(6,1);
for i = 1:6
    [m r] = anglemean(c13(i,:));
    baseline(i) = r;
end
baseline = mean(baseline);
valid = place > 0;

nboot = 100;
x1 = c13(:,valid&p<0.05); x2 = c13(:,valid&p>0.05);
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
sig = mean(q1,2)-baseline; non = mean(q2,2)-baseline;
se_s = std(q1,[],2); se_n = std(q2,[],2);

figure
hold on
plot(time,sig ,'r',time,non ,'k')
legend([{'Significant replay'},{'Nonsignificant replay'}],'location','NorthWest')
fill([time time(end:-1:1)],[sig+se_s; sig(end:-1:1)-se_s(end:-1:1)],'r','EdgeColor','none')
fill([time time(end:-1:1)],[non+se_n; non(end:-1:1)-se_n(end:-1:1)],'k','EdgeColor','none')
set(gca,'xtick',time(6:10:end),'xticklabel',-0.4:0.1:0.4,'ylim',[-0.05 0.105],'ytick',-0.05:0.05:15,'xlim',[-0.4 0.4])
ylabel('CA1-CA3 phase locking')
xlabel('Time since ripple detection (s)')
box off

%Permutation test for significant > nonsignificant
Nx = sum(valid&p<0.05); Ny = sum(valid&p>0.05);
X = c13(:,valid&p<0.05); Y = c13(:,valid&p>0.05);
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

%Permutation test indicates difference is significant for run sessions time 50 - 240ms
%Permutation test indicates difference is not significant for sleep sessions

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1ca3_phase_variance_sleep_sig_nonsig_replay.pdf', m, d, y);
print('-dpdf', savestring)

time = valid_ripples{1}{1}{1}.eegtime;
subs = [0 0.05 0.1 0.5 1];
pvalue = lookup(p,subs,1);
a13 = zeros(length(subs),length(time));
valid = place>0;
for i = 1:size(c13,1)
    for j = 1:length(subs)
        [m r] = anglemean(c13(i,pvalue==j&valid));
        a13(j,i) = r-baseline;
    end
end
g = gaussian(1.5,10);
for i = 1:size(a13,1)
   a13(i,:) = conv(a13(i,:),g,'same'); 
end

figure
imagesc(time,[],a13)
set(gca,'clim',[-0.05 0.1],'ytick',1:length(subs),'yticklabel',subs,'xlim',[-0.4 0.4])
colorbar('ytick',-0.05:0.05:.1)
xlabel('Time since SWR detection (sec)')
ylabel('Replay p-value')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1ca3_phaselocking_pvalue_replay_sleep.png', m, d, y);
print('-dpng', savestring)

%% LOOK AT RVALUES
subs = prctile(rvalue,[25 50 75 100]);
pvalue = lookup(rvalue,subs,1);
a13 = zeros(length(subs),length(time));
valid = place>0;

for i = 1:size(c13,1)
    for j = 1:length(subs)
        [m r] = anglemean(c13(i,pvalue==j&valid));
        a13(j,i) = r-baseline;
    end
end
g = gaussian(1.5,10);
for i = 1:size(a13,1)
   a13(i,:) = conv(a13(i,:),g,'same'); 
end

figure
imagesc(time,[],a13)
set(gca,'clim',[-0.05 0.1],'ytick',1:length(subs),'yticklabel',subs,'xlim',[-0.4 0.4])
colorbar('ytick',-0.5:0.025:.7)
xlabel('Time since SWR detection (sec)')
ylabel('Rvalue')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1ca3_phaselocking_rvalue.png', m, d, y);
print('-dpng', savestring)

%% LOOK AT ENTROPY
subs = prctile(entropy,[25 50 75 100]);
pvalue = lookup(entropy,subs,1);
a13 = zeros(length(subs),length(time));
valid = place>0;

for i = 1:size(c13,1)
    for j = 1:length(subs)
        [m r] = anglemean(c13(i,pvalue==j&valid));
        a13(j,i) = r-baseline;
    end
end
g = gaussian(1.5,10);
for i = 1:size(a13,1)
   a13(i,:) = conv(a13(i,:),g,'same'); 
end

figure
imagesc(time,[],a13)
set(gca,'clim',[-0.05 0.15],'ytick',1:length(subs),'yticklabel',subs,'xlim',[-0.4 0.4])
colorbar('ytick',-0.5:0.05:.7)
xlabel('Time since SWR detection (sec)')
ylabel('Entropy')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1ca3_phaselocking_entropy.png', m, d, y);
print('-dpng', savestring)

%% PLOT RVALUES VS. PVALUES
figure
semilogy(rvalue,p,'k.')
xlabel('Rvalue')
ylabel('Pvalue (log scale)')
box off
set(gca,'ytick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0])

%corr(rvalue,p,'type','Spearmam') = -0.78

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_rvalue_vs_pvalue.pdf', m, d, y);
print('-dpdf', savestring)

%% PLOT RVALUES VS. PVALUES
figure
semilogy(entropy,p,'k.')
xlabel('Entropy')
ylabel('Pvalue (log scale)')
box off
set(gca,'ytick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0])

%corr(rvalue,p,'type','Spearman') = 0.22
%% Look at % of place cells active during each candidate event

run_p = []; run_place = [];
sleep_p = []; sleep_place = [];
for an =1:length(valid_ripples)
    for d =1:length(valid_ripples{an})
        for e = 1:length(valid_ripples{an}{d})
            if ~isempty(valid_ripples{an}{d}{e})
                if isfield(valid_ripples{an}{d}{e},'ca1_ca3')
                    if size(valid_ripples{an}{d}{e}.pvalue,1) == size(valid_ripples{an}{d}{e}.ca1_ca3,2) 
                        if size(valid_ripples{an}{d}{e}.pvalue,1)== size(valid_ripples{an}{d}{e}.eventimmobiletime,1)
                        	if d == 1 | d == 2
                                run_p = [run_p; valid_ripples{an}{d}{e}.pvalue];
                                run_place = [run_place; valid_ripples{an}{d}{e}.placecells];
                            else
                                valid = valid_ripples{an}{d}{e}.eventimmobiletime>60;
                                sleep_p = [sleep_p; valid_ripples{an}{d}{e}.pvalue(valid)];
                                sleep_place = [sleep_place; valid_ripples{an}{d}{e}.placecells(valid)];
                            end
                        end
                    end
                end
            end
        end
    end
end
group = [ones(size(run_place)); 2*ones(size(sleep_place)); 3*ones(sum(run_p<0.05),1); 4*ones(sum(sleep_p<0.05),1)];

figure
boxplot([run_place; sleep_place; run_place(run_p<0.05); sleep_place(sleep_p<0.05)],group,'labels',[{'Run'},{'Sleep'},{'Sig Run'},{'Sig Sleep'}],'symbol','b','whisker',0)
set(gca,'ylim',[0 1.1],'xlim',[0.5 4.5],'ytick',0:0.2:1,'yticklabel',0:20:100)
ylabel('Percent place cells active during candidate events')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_percentplacecells_runvssleep.pdf', m, d, y);
print('-dpdf', savestring)