function h= auxPlotMeanDev(x, y, opt, indev)
%function h= auxPlotMeanDev(x, y, opt, indev)
%
% y{i}
% opt.nmin
% opt.plotcol
% opt.center= {'mean', 'median'}
% opt.dev= {'se', 'std', 'none'}
% opt.xrange
% opt.yrange

ny= length(y);

if nargin<3 opt=[]; end
if ~isfield(opt, 'nmin'); opt.nmin= 1; end
if ~isfield(opt, 'plotcol'); opt.plotcol= jet(3); end
if ~isfield(opt, 'center'); opt.center= 'mean'; end
if ~isfield(opt, 'dev'); opt.dev= 'se'; end
if nargin>=4; opt.dev='input'; end

fillcol= rgb2hsv(opt.plotcol);
fillcol(:,2)= .4;
fillcol= hsv2rgb(fillcol);
ind= find(all(fillcol==0,2));
fillcol(ind,:)= .7*ones(length(ind),3);
clear ind

hold on


fh= [];
for iy=1:ny
    nval= sum(isfinite(y{iy}));
    nval(find(nval==0))= nan;
%    stdev= nanstd(y{iy},0,1);
    stdev= nanstd(y{iy});

    switch opt.center
    case 'mean'
%        m{iy}= nanmean(y{iy},1);
        m{iy}= nanmean(y{iy});
    case 'single'
        m{iy}= y{iy};
    otherwise
        error('unknown opt.center');
    end

    switch opt.dev
    case 'std'
        dev= stdev;
    case 'se'
        dev= stdev./sqrt(nval);
    case 'none'
        dev= [];
    case 'input'
        dev= indev{iy};
    otherwise
        error('unknown opt.dev');
    end

    if ~isempty(dev)
        ind= find(~isfinite(m{iy}) | ~isfinite(dev) | nval < opt.nmin);
    else
        ind= find(~isfinite(m{iy}) | nval < opt.nmin);
    end
    m{iy}(ind)= nan;

    if ~isempty(dev)
        ind= find(isfinite(m{iy}));
        [lo hi]= findcontiguous(ind);
        for i=1:length(lo)
            ind= lo(i):hi(i);
            fh(end+1)= fill([x(ind) fliplr(x(ind))], ...
                [m{iy}(ind)-dev(ind), fliplr(m{iy}(ind)+dev(ind))], fillcol(iy,:));
            set(fh(end), 'EdgeColor', fillcol(iy,:));
        end
%keyboard
    end
end


for iy=1:ny
    hp(iy)= plot(x, m{iy});
    set(hp(iy), 'Color', opt.plotcol(iy,:));
end

if isfield(opt, 'xrange'); set(gca, 'xlim', opt.xrange); end
if isfield(opt, 'yrange'); set(gca, 'ylim', opt.yrange); end

%    keyboard

h= [fh, hp];
hold off
