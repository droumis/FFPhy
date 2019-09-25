figure('Units','normalized', ...
    'Name','ORDINATIONS', ...
  'Color',[0.8 0.8 0.8], ...
  'Position',[ 0.01375 0.2 0.0925 0.52 ]);
h1021=uimenu('Label','ORDINATIONS');
uimenu('Parent',h1021, ...
  'Callback','help pcacov', ...
  'Label','pcacov');
uimenu('Parent',h1021, ...
  'Callback','help pcacorr', ...
  'Label','pcacorr');
uimenu('Parent',h1021, ...
  'Callback','help discrim', ...
  'Label','discrim');
uimenu('Parent',h1021, ...
  'Callback','help sizefree', ...
  'Label','sizefree');
uimenu('Parent',h1021, ...
  'Callback','help mahal', ...
  'Label','mahal');
uimenu('Parent',h1021, ...
  'Callback','help classify', ...
  'Label','classify');
%pcacov
uicontrol('Units','normalized', ...
  'BackgroundColor',[0.964705882352941 0.247058823529412 0.105882352941176], ...
  'Callback','[loadings,percvar,scores]= pcacov(x);figure;plotgrps(scores(:,1),scores(:,2),g);title(''PCACOV'')',...
  'Position',[0.156521739130435 0.803664921465969 0.695652173913043 0.120418848167539 ], ...
  'String','pcacov');
%pcacorr
uicontrol('Units','normalized', ...
  'BackgroundColor',[0.964705882352941 0.247058823529412 0.105882352941176], ...
  'Callback','[loadings,percvar,scores]= pcacorr(x);figure;plotgrps(scores(:,1),scores(:,2),g);title(''PCACORR'')',...
  'Position',[ 0.165217391304348 0.659685863874346 0.695652173913043 0.120418848167539 ], ...
  'String','pcacorr');
%discrim
uicontrol('Units','normalized', ...
  'BackgroundColor',[0.964705882352941 0.247058823529412 0.105882352941176], ...
  'Position',[ 0.165217391304348 0.520942408376963 0.704347826086957 0.117801047120419 ], ...
  'Callback','[loadings,percvar,scores]= discrim(x,g,2);figure;plotgrps(scores(:,1),scores(:,2),g);title(''DISCRIM'')',...
  'String','discrim');
%sizefree
uicontrol('Units','normalized', ...
  'BackgroundColor',[0.964705882352941 0.247058823529412 0.105882352941176], ...
  'Position',[ 0.156521739130435 0.37434554973822 0.704347826086957 0.117801047120419 ], ...
  'Callback','[loadings,percvar,scores] = sizefree(x,g,2);figure;plotgrps(scores(:,1),scores(:,2),g);title(''SIZEFREE'')',...
  'String','sizefree');
%mahal
uicontrol('Units','normalized', ...
  'BackgroundColor',[0.964705882352941 0.247058823529412 0.105882352941176], ...
  'Position',[ 0.156521739130435 0.230366492146597 0.704347826086957 0.117801047120419 ], ...
  'Callback','[D2] = mahal(x,g)',...
  'String','mahal');
%classify
uicontrol('Units','normalized', ...
  'BackgroundColor',[0.964705882352941 0.247058823529412 0.105882352941176], ...
  'Position',[ 0.165217391304348 0.0863874345549738 0.704347826086957 0.117801047120419 ], ...
  'Callback','[results,perclass] = classify(x,g)',...
  'String','classify');

