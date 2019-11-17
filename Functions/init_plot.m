function ifig = init_plot(showfigs, varargin)
% initialize plot
% ifig = init_plot(showfigs, varargin)
% DR 19
position = [.1 .1 .5 .5];
if ~isempty(varargin)
    assign(varargin{:})
end
% init fig
close all;
if showfigs
    ifig = figure('units','normalized','position',position, 'color','white', ...
        'InvertHardcopy', 'off');
else    
    ifig = figure('Visible','off','units','normalized','position', position, ...
        'color','white', 'InvertHardcopy', 'off'); 
end