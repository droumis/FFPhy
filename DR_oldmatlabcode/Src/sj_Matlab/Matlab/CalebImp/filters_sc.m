

%-----------------------------------------------------------
fprintf('Generating filters...');
%-----------------------------------------------------------
tetrodes = [1:32];
tetrodefilter = '$tetrode>=0';

% (FILTER ORDER IS IMPORTANT B/C OF USE OF SHORT CIRCUIT && !!
epochfilter = 'strcmp($type,''run'') & (abs($exposure)==1)';

animals = {'samson','uriah','tola','nebuchadnezzar'};
runEpochs = createfilter('animal',animals,'epochs',epochfilter);
runEpochs = combinefilterepochs(runEpochs,'single','empty'); % combine same days and clear empty epochs

diopinfilter = {'(strcmp($type,''ca3''))'};
runEpochs = setdiopinfilter(runEpochs,diopinfilter); 
% this should probably change the "data" field to be the pin to
% be used

% (samson, uriah, tola, neb)
SC_channels = [ 11, 5, 7, 8];

% times?
SC_times = {[0.007 0.0095], [0.006 0.0095], [0.009 0.014], [0.010 0.013]};
% SC_times = {[0.0065 0.009], [0.0065 0.009], [0.009 0.013], [0.010 0.013]};
% SC_times = {[0.0065 0.009], [0.007 0.009], [0.009 0.013], [0.010 0.013]};

normSpeedThresh = 5; % cm/s
velCats = [0.125 0.25 0.5 1 2 4 8 16 32 64];
veltick = [0.125 0.25 0.5 1 2 4 8 16 32 64];
% velocity bounds for regression
minV = velCats(1);
maxV = velCats(end);

runEPSP = runEpochs;
for an = 1:length(runEPSP)
  if isempty(runEPSP(an).epochs) || isempty(runEPSP(an).epochs{1})
    continue;
  end
  totaldays = cat(1,runEPSP(an).epochs{:});
  totaldays = unique(totaldays(:,1));
  stimdio = loaddatastruct(runEPSP(an).animal{2}, runEPSP(an).animal{3}, 'stimdio', totaldays);
  for ep_g = 1:length(runEPSP(an).epochs)
    for ep = 1:size(runEPSP(an).epochs{ep_g},1)
      runEPSP(an).times{ep_g}{ep} = stimdio{runEPSP(an).epochs{ep_g}(ep,1)}{runEPSP(an).epochs{ep_g}(ep,2)}{runEPSP(an).pin_index{ep_g}{1}{1}}.pulsetimes(:,1) / 10000;
    end
  end

  stimeeg = loaddatastruct(runEPSP(an).animal{2},runEPSP(an).animal{3},'stimeeg',totaldays);
  for ep_g = 1:length(runEPSP(an).epochs)
    for ep = 1:size(runEPSP(an).epochs{ep_g},1)
      runEPSP(an).epsps{ep_g}{ep} = stimeeg{runEPSP(an).epochs{ep_g}(ep,1)}{runEPSP(an).epochs{ep_g}(ep,2)}{runEPSP(an).pin_index{ep_g}{1}{1}}.data;
      runEPSP(an).t{ep_g}{ep} = stimeeg{runEPSP(an).epochs{ep_g}(ep,1)}{runEPSP(an).epochs{ep_g}(ep,2)}{runEPSP(an).pin_index{ep_g}{1}{1}}.t;
    end
  end
end

runBehav = modifyfilter(runEPSP, 'iterator', 'epochbehavtimeanal', ...
  'filterfunction', {'getbehavstatetime',{'pos','linpos','task','track_regions'}, ...
                     'behavVars',{'smoothvelocity','position'},'t_offset',0});

%-----------------------------------------------------------
fprintf('\nRunning filters...');
%-----------------------------------------------------
fprintf(' [B]'); runBehav = runfilter(runBehav);
fprintf('\n');


global_SC = [];
vv_SC = [];
for an = 1:length(runEPSP)
  if isempty(runEPSP(an).epochs) || isempty(runEPSP(an).epochs{1})
    continue;
  end
  stimeeg = loaddatastruct(runEPSP(an).animal{2},runEPSP(an).animal{3},'stimeeg',totaldays);
  for ep_g = 1:length(runEPSP(an).epochs)
    for ep = 1:size(runEPSP(an).epochs{ep_g},1)
      runEPSP(an).epsp_slope{ep_g}{ep} = ...
        measureEPSPs(runEPSP(an).t{ep_g}{ep},...
        runEPSP(an).epsps{ep_g}{ep}(:,:,SC_channels(an)), SC_times{an});
      vv = runBehav(an).output{ep_g}(:,1);
      epsp = runEPSP(an).epsp_slope{ep_g}{ep};
      norm_epsp = epsp / mean(epsp(vv > normSpeedThresh));
      runEPSP(an).norm_epsp{ep_g}{ep} = norm_epsp;
      goodV = (vv > minV) & (vv < maxV);
      [beta, bint, r, rint, stats] = regress(norm_epsp(goodV) ,[ones(sum(goodV),1) log(vv(goodV))]);
      runEPSP(an).regress_results{ep_g}{ep}.stats = stats;
      runEPSP(an).regress_results{ep_g}{ep}.beta = beta;
      runEPSP(an).regress_results{ep_g}{ep}.bint = bint;
      global_SC = [global_SC; norm_epsp];
      vv_SC = [vv_SC; vv];
      if sum(imag(vv)) ~= 0
        keyboard
      end
      % figure
      % plot(log(vv(vv>minV&vv<maxV)), norm_epsp(vv>minV&vv<maxV),'.');
    end
  end
end
[bSC,bintSC,statsSC] = fastreg(log(vv_SC(vv_SC>minV & vv_SC<maxV)),global_SC(vv_SC>minV & vv_SC<maxV));


