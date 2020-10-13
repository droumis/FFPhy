function [dioOut, dioOutfields] = get_dio_output(DIO, task, day, eps);

dioOutfields = {'day', 'epoch', 'dioNum', 'dioNumOrig', 'startTime', 'endTime'};
dioOut = [];
for e = 1:length(eps)
    epoch = eps(e);
    isinput = cellfun(@(x) isequal(x.input,1), DIO{day}{epoch}, 'un', 1);
    dioID = cellfun(@(x) str2double(regexp(x.original_id,'\d*','Match')), ...
        DIO{day}{epoch}, 'un', 1);
    outputdios = task{day}{epoch}.outputdio;
    outdioIdx = find(all([ismember(dioID, outputdios)' ~isinput'],2));
    for c = 1:length(outdioIdx)
        ch = outdioIdx(c);
        dioTimeCh = double(DIO{day}{epoch}{ch}.times);
        dioValsCh = double(DIO{day}{epoch}{ch}.values);
        % dios are Up (1) to Down (0)
        while ~isempty(dioValsCh) && dioValsCh(1) == 0
            dioValsCh(1) = []; dioTimeCh(1) = [];
        end
        while ~isempty(dioValsCh) && dioValsCh(end) == 1
            dioValsCh(end) = []; dioTimeCh(end) = [];
        end
        if ~isempty(dioValsCh)
            dioVd = [1; find(abs(diff(dioValsCh)))+1];
            diotimes = dioTimeCh(dioVd);
            dioStEnd = [diotimes(1:2:end) diotimes(2:2:end)];
            numDIO = length(dioStEnd(:,1));
            dioOut = [dioOut; repmat(day,numDIO,1) repmat(epoch,numDIO,1) ...
                repmat(ch,numDIO,1) repmat(outputdios(c),numDIO,1) dioStEnd];
        else
            continue
        end
    end
end