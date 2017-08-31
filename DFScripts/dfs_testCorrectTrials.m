




























%% ---------------- Run FIlter ---------------------------------------------------
if runFilterFramework == 1;
    epochfilter =    sprintf('(isequal($type, ''%s'')) && (isequal($environment, ''%s''))',epochType, epochEnvironment);
    iterator = 'multitetrodeanal'; %multitetrodeanal
    %tetfilter: the ntrodeID's that pass this filter get stored into f.eegdata{day}{epoch}[ntrodeIDs]
%     tetfilter = sprintf('(isequal($area,''%s'') || isequal($area,''%s'') || isequal($area,''%s'') || isequal($area,''%s'') || isequal($area,''%s'') || isequal($area,''%s''))', ntAreas{1}, ntAreas{2}, ntAreas{3}, ntAreas{4}, ntAreas{5}, ntAreas{6});
    tetfilter = sprintf('(isequal($area,''%s'')) || ', ntAreas{:});
    yuck = strfind(tetfilter, '||');
    tetfilter = tetfilter(1:yuck(end)-1); %cut off the trailing '||'
    % timefilter{1} = {'get2dstate','($velocity<4)'};
    
    timefilter{1} = {'getWtracktrialstate', '(($correct == 1) && ($outbound == 1) && ($outboundScore > .8))'};
%     timefilter{2} = {'kk_getriptimes','($nripples>=1)',[],'tetfilter',tetfilter,'minthresh',5};
    %---------- save all filtering parameters in workspace into struct
    iF.datafilter = whos;
    %----------F = createfilter('animal', animals, 'days', dayfilter,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);--------
    F = createfilter('animal', animals, 'days', days,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);
    %----------f = setfilteriterator(f, funcname, loadvariables, options)--------
    eval(['F = setfilterfunction(F, [''dfa_'' filtfunction], {[eventSourceArea eventtype],' strjoin(arrayfun(@(x) sprintf('''%s'',', cell2mat(x)), LFPtypes,'UniformOutput',false)) '},''eventtype'', [eventSourceArea eventtype], ''LFPtypes'', LFPtypes);']);
%     eval(['F = setfilterfunction(F, [''dfa_'' filtfunction], {[eventSourceArea eventtype],' strjoin(reshape(repmat(arrayfun(@(x) sprintf('''%s'',', cell2mat(x)), LFPtypes,'UniformOutput',false), 2, 1), 1,length(LFPtypes)*2)) '},''eventtype'', [eventSourceArea eventtype], ''LFPtypes'', LFPtypes);']);
    tic
    F = runfilter(F);
    F(1).filterTimer = toc; F(1).filterTimer
    F(1).worldDateTime = clock;
    F(1).dataFilter = struct('animal', animals, 'days', days,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator, 'filename', filename);
end