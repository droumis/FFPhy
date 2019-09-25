

%-----------------------------------------------------------
fprintf('Generating filters...');
%-----------------------------------------------------------
tetrodes = [1:32];
tetrodefilter = '$tetrode>=0';

% (FILTER ORDER IS IMPORTANT B/C OF USE OF SHORT CIRCUIT && !!
epochfilter = 'strcmp($type,''run'') & (abs($exposure)==1)';

animals = {'samson','uriah'};
runEpochs = createfilter('animal',animals,'epochs',epochfilter);
runEpochs = combinefilterepochs(runEpochs,'single','empty'); % combine same days and clear empty epochs

diopinfilter = {'(strcmp($type,''ec''))','(strcmp($type,''ca3''))'};
runEpochs = setdiopinfilter(runEpochs,diopinfilter); 
% this should probably change the "data" field to be the pin to
% be used

% (samson, uriah)
EC_channels = [15, 12];
CA3_channels = [11, 5];

% times?
EC_times = {[0.0026 0.0036], [0.0026 0.0038]};
CA3_times = {[0.007 0.0095], [0.007 0.0095]};

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
      for ff = 1:length(diopinfilter)
        runEPSP(an).times{ep_g}{ep}{ff} = stimdio{runEPSP(an).epochs{ep_g}(ep,1)}{runEPSP(an).epochs{ep_g}(ep,2)}{runEPSP(an).pin_index{ep_g}{1}{ff}}.pulsetimes(:,1) / 10000;
      end
    end
  end

  stimeeg = loaddatastruct(runEPSP(an).animal{2},runEPSP(an).animal{3},'stimeeg',totaldays);
  for ep_g = 1:length(runEPSP(an).epochs)
    for ep = 1:size(runEPSP(an).epochs{ep_g},1)
      for ff = 1:length(diopinfilter)
        runEPSP(an).epsps{ep_g}{ep}{ff} = stimeeg{runEPSP(an).epochs{ep_g}(ep,1)}{runEPSP(an).epochs{ep_g}(ep,2)}{runEPSP(an).pin_index{ep_g}{1}{ff}}.data;
        runEPSP(an).t{ep_g}{ep}{ff} = stimeeg{runEPSP(an).epochs{ep_g}(ep,1)}{runEPSP(an).epochs{ep_g}(ep,2)}{runEPSP(an).pin_index{ep_g}{1}{ff}}.t;
      end
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


ff = 1;
global_EC = [];
vv_EC = [];
for an = 1:length(runEPSP)
  if isempty(runEPSP(an).epochs) || isempty(runEPSP(an).epochs{1})
    continue;
  end
  stimeeg = loaddatastruct(runEPSP(an).animal{2},runEPSP(an).animal{3},'stimeeg',totaldays);
  for ep_g = 1:length(runEPSP(an).epochs)
    for ep = 1:size(runEPSP(an).epochs{ep_g},1)
      runEPSP(an).epsp_slope{ep_g}{ep}{ff} = ...
        measureEPSPs(runEPSP(an).t{ep_g}{ep}{ff},...
        runEPSP(an).epsps{ep_g}{ep}{ff}(:,:,EC_channels(an)), EC_times{an});
      vv = runBehav(an).output{ep_g}{ff}(:,1);
      epsp = runEPSP(an).epsp_slope{ep_g}{ep}{ff};
      norm_epsp = epsp / mean(epsp(vv > normSpeedThresh));
      runEPSP(an).norm_epsp{ep_g}{ep}{ff} = norm_epsp;
      goodV = (vv > minV) & (vv < maxV);
      [beta, bint, r, rint, stats] = regress(norm_epsp(goodV) ,[ones(sum(goodV),1) log(vv(goodV))]);
      runEPSP(an).regress_results{ep_g}{ep}{ff}.stats = stats;
      runEPSP(an).regress_results{ep_g}{ep}{ff}.beta = beta;
      runEPSP(an).regress_results{ep_g}{ep}{ff}.bint = bint;
      global_EC = [global_EC; norm_epsp];
      vv_EC = [vv_EC; vv];
      if sum(imag(vv)) ~= 0
        keyboard
      end
    end
  end
end
[bEC,bintEC,statsEC] = fastreg(log(vv_EC(vv_EC>minV & vv_EC<maxV)),global_EC(vv_EC>minV & vv_EC<maxV));

ff = 2;
global_CA3 = [];
vv_CA3 = [];
for an = 1:length(runEPSP)
  if isempty(runEPSP(an).epochs) || isempty(runEPSP(an).epochs{1})
    continue;
  end
  stimeeg = loaddatastruct(runEPSP(an).animal{2},runEPSP(an).animal{3},'stimeeg',totaldays);
  for ep_g = 1:length(runEPSP(an).epochs)
    for ep = 1:size(runEPSP(an).epochs{ep_g},1)
      runEPSP(an).epsp_slope{ep_g}{ep}{ff} = ...
        measureEPSPs(runEPSP(an).t{ep_g}{ep}{ff},...
        runEPSP(an).epsps{ep_g}{ep}{ff}(:,:,CA3_channels(an)), CA3_times{an});
      vv = runBehav(an).output{ep_g}{ff}(:,1);
      epsp = runEPSP(an).epsp_slope{ep_g}{ep}{ff};
      norm_epsp = epsp / mean(epsp(vv > normSpeedThresh));
      runEPSP(an).norm_epsp{ep_g}{ep}{ff} = norm_epsp;
      goodV = (vv > minV) & (vv < maxV);
      [beta, bint, r, rint, stats] = regress(norm_epsp(goodV) ,[ones(sum(goodV),1) log(vv(goodV))]);
      runEPSP(an).regress_results{ep_g}{ep}{ff}.stats = stats;
      runEPSP(an).regress_results{ep_g}{ep}{ff}.beta = beta;
      runEPSP(an).regress_results{ep_g}{ep}{ff}.bint = bint;
      global_CA3 = [global_CA3; norm_epsp];
      vv_CA3 = [vv_CA3; vv];
      if sum(imag(vv)) ~= 0
        keyboard
      end
    end
  end
end
[bCA3,bintCA3,statsCA3] = fastreg(log(vv_CA3(vv_CA3>minV & vv_CA3<maxV)),global_CA3(vv_CA3>minV & vv_CA3<maxV));




