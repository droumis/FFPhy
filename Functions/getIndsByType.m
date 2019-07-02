



function [iEpTetIndsType, DataTypeFields, DataTypeMat] = getIndsByType(indsSetName,eState, ian, varargin)
fullsize = [];
DataTypeFields = '';
dataByDay = [];
if ~isempty(varargin)
    assign(varargin{:});
end
[iEpTetIndsType, DataTypeFields, DataTypeMat]  = setParams(indsSetName,eState, ian, 'DataTypeFields', DataTypeFields, 'fullsize', fullsize, 'dataByDay', dataByDay);
end

%add custom iEpTetIndsType param sets as cases below

function [iEpTetIndsType, DataTypeFields, DataTypeMat] = setParams(indsSetName,eState, ian, varargin)
fullsize = [];
DataTypeMat = [];
DataTypeFields = [];
% DataTypeFields = '';
dataByDay = [];
if ~isempty(varargin)
    assign(varargin{:});
end

switch indsSetName
    case 'performance'
        DataTypeFields = {'all', 'corrOut', 'mistOut', 'outB', 'inB','outB-inB', 'corrOut-mistOut'};
        %         corrOutBInd = find(cell2mat(cellfun(@(x) strcmp(x,'corrOut'), iEpTetDataTypeFields, 'un', false)));
        %         mistOutBInd = find(cell2mat(cellfun(@(x) strcmp(x,'mistOut'), iEpTetDataTypeFields, 'un', false)));
        corrOutBInd = 6;
        mistOutBInd = 7;
        corrOutEventInd = find(eventstate.state(1:end,corrOutBInd) == 1);
        mistOutEventInd = find(eventstate.state(1:end,mistOutBInd) == 1);
        outBEventInd = [corrOutEventInd; mistOutEventInd];
        inBEventInd = setdiff([1:length(eventstate.state(1:end,1))], outBEventInd)';
        iEpTetIndsType = {[1:length(eventstate.state(1:end,1))], corrOutEventInd, mistOutEventInd, outBEventInd, inBEventInd};
%         baselinePhDataIndSet = num2cell(repmat(1, 1, length(DataTypeFields)));  %index into the iEpTetIndsType cell array for which set of indices to use as baseline
    case 'sleep1'
          DataTypeFields = {'sleep'};
          iEpTetIndsType =   {[1:fullsize]};
%         baselinePhDataIndSet = num2cell(repmat(1, 1, length(DataTypeFields)));  %index into the iEpTetIndsType cell array for which set of indices to use as baseline
    case 'sleepByDay'
        unqDays = unique(eventstate.state(:,daycol));
        daysInds = arrayfun(@(x) find(eventstate.state(:,daycol) == x), unqDays, 'un', 0);
        iEpTetIndsType = daysInds';
        DataTypeFields = cellstr(num2str(unqDays))';
    case 'byDay'
        unqDays = unique(dataByDay);
        daysInds = arrayfun(@(x) find(dataByDay == x), unqDays, 'un', 0);
        iEpTetIndsType = daysInds;
        DataTypeFields = [{{'all'}}, {'days'}];
        DataTypeMat = repmat(num2cell(unique(dataByDay, 'stable'))', length(DataTypeFields{1}), 1);
%         baselinePhDataIndSet = num2cell(unqDays');
    case 'performanceByDay'
        daycol =  find(cell2mat(cellfun(@(x) strcmp(x,'day'), eState.fields, 'un', false)));
        corrOutBInd = find(cell2mat(cellfun(@(x) strcmp(x,'correct'), eState.fields, 'un', false)));
        mistOutBInd = find(cell2mat(cellfun(@(x) strcmp(x,'mistake'), eState.fields, 'un', false)));
        
        unqDays = unique(eState.state(:,daycol), 'stable');
        daysInds = arrayfun(@(x) find(eState.state(:,daycol) == x), unqDays, 'un', 0);
        corrData = eState.state(1:end,corrOutBInd);
        dayData = eState.state(1:end,daycol);
        mistData = eState.state(1:end,mistOutBInd);
%         inbData = ixpc.eventstate{ian}.state(1:end,InBInd);
        iEpTetIndsType = [];
        for iday = 1:length(daysInds)
            corrOutEventInd = find(corrData == 1 & dayData == unqDays(iday));
            mistOutEventInd = find(mistData == 1 & dayData == unqDays(iday));
            outBEventInd = find((corrData == 1 | mistData == 1) & dayData == unqDays(iday));
            inBEventInd = setdiff(find(dayData == unqDays(iday)),outBEventInd);
            idaydata = {[1:length(dayData(dayData == unqDays(iday)))]', corrOutEventInd, mistOutEventInd, outBEventInd, inBEventInd};
            iEpTetIndsType = cat(1,iEpTetIndsType, idaydata);
        end
%         DataTypeFields = [{{'all', 'corrOut', 'mistOut', 'outB', 'inB','outB-inB', 'corrOut-mistOut'}}, {'days'}];
%         DataTypeMat = repmat(num2cell(unique(dataByDay, 'stable'))', length(DataTypeFields{1}), 1);
    otherwise
        error('pick an existing param set or make a new one')
end
end