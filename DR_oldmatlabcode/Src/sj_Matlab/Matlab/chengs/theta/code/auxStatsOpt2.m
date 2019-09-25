function [sopts, ana]= auxStatsOpt2(data,ana,rat,d,e);
%function [sopts, ana]= auxStatsOpt2(data,ana,rat,d,e);

sopts{1}.name= 'LinCorrOffset';
sopts{1}.sizex= 0.1;
sopts{1}.sizey= pi/20;

sopts{5}.name= 'LinCorrShift';
sopts{5}.sizex= 0.1;
sopts{5}.sizey= pi/20;

sopts{2}.name= 'CircMeasures';
sopts{2}.sizex= 0.1;
sopts{2}.sizey= pi/20;

%sopts{3}.name= 'MutualInfo';
%sopts{3}.nbinx= 200;
%sopts{3}.nbiny= 200;

sopts{3}.name= 'Moments2d';
sopts{3}.sizex= 0.1;
sopts{3}.sizey= pi/20;
sopts{3}.circy= 1;


sopts{4}.name=  'Integral2d';
sopts{4}.dx= 0.1;
sopts{4}.dy= pi/20;

for i=1:length(ana.a)
    ana.a{i}.phase= [0;2*pi];
    traj= ana.a{i}.traj;
    if isempty(ana.x) | (length(ana.x)<i) | ~isfield(ana.x{i},'time')
        % calc stats at every second of occupancy
        ana.x{i}.time= getTimes('occ', rat, d, e, ana.a{i}.traj);
    end
end

