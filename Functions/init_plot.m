function ifig = init_plot(showfigs, position, varargin)
% initialize plot to my defaults
% ifig = init_plot(showfigs, position)
% FFPhy V0.1
% @DR
if nargin == 1
    position = [.1 .1 .5 .5];
end
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