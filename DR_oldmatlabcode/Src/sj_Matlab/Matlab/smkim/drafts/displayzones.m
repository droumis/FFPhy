function fighandle = display_zones(zonesstruct)

% plot the roi polygons on top of the sample frame
for i = 1:length(zonesstruct)
    patch('Faces',1:size(zonesstruct(i).vertices,1), ...
    'Vertices',zonesstruct(i).vertices, 'FaceColor', ...
    [0 0.9 0],'FaceAlpha',0.2,'EdgeColor','none');
end
hold off;
