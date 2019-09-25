function bin= auxBin2d(bin, varargin)
% bin.nx bin.ny indicate the edges of the bins (not their centers!)
%
% The last bins contain count datapoints that are outside the valid range.

if nargin<1; error('need to specify at least the bin structure'); end

switch length(varargin)
case 0
    sprintf('initializing bin structure');
    bin.nx= length(bin.x);
    bin.ny= length(bin.y);
    bin.z= zeros(bin.nx,bin.ny);
case 1
    switch varargin{1}
    case 'show'
%        imagesc(dens.z);

        figure
%        x= (bin.x(1:end-1)+bin.x(2:end))/2;
%        y= (bin.y(1:end-1)+bin.y(2:end))/2;
        x= (bin.x(1:end-1)+bin.x(2:end))/2;
        y= (bin.y(1:end-1)+bin.y(2:end))/2;
        imagesc(x, y, log10(bin.z(1:end-1, 1:end-1)));
%        imagesc(x, y, bin.z(1:end-1, 1:end-1));

        set(gca, 'xlim', [3 50])
        set(gca, 'ylim', [-1 50])
        axis xy
        axis square
        h= colorbar
        set(h, 'ytick', [1 2 3]);
        set(h, 'yticklabel', [10^1 10^2 10^3]);
        xlabel('reference tetrode');
        ylabel('other tetrodes');

        nout= sum(bin.z(1:end,end))+ sum(bin.z(end,1:end-1));
        N= sum(sum(bin.z));
        fprintf(1, '%d/ %d (=%.1f%%) datapoints outside limits.\n', nout, N, nout/N*100);
        saveas(gcf, 'rip_coh_sum')

        figure
        subplot(2,1,1)
        mx= sum(bin.z(1:end-1, 1:end-1),2);
        bar(x,mx);
        xlabel('reference tetrode');
        ylabel('count');
        set(gca, 'xlim', [3 50])
        set(gca, 'ylim', [0 8000])

        subplot(2,1,2)
        my= sum(bin.z(1:end-1, 1:end-1));
        bar(y,my);
        xlabel('other tetrodes');
        ylabel('count');
        set(gca, 'xlim', [-1 50])
        set(gca, 'ylim', [0 8000])

        
        keyboard
    otherwise
        error(['unknown command ' varargin{1}]);
    end
case 2
%    adding datapoints
    x= varargin{1}; y= varargin{2};
    n= length(x);
    if n~=length(y); error('x and y have different length'); end
    for i=1:n
        ix= max(find(bin.x<=x(i)));
        if isempty(ix) | x(i)>=bin.x(end); ix= bin.nx; end
        iy= max(find(bin.y<=y(i)));
        if isempty(iy) | y(i)>=bin.y(end); iy= bin.ny; end
        bin.z(ix,iy)= bin.z(ix,iy)+1;
    end

    otherwise
    error('no defined actions');
end
