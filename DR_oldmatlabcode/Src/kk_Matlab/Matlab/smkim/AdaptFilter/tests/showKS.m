function showKS(file_id)
%function showKS(file_id)

if nargin<1; file_id= ''; end
color= true;

filename= 'dynamic';

% show AF results
load([filename '_data']);
load([filename file_id '_result']);

sorted= sort(result.rescaled_isi);

KSdist95 = 1.36;
KSdist99 = 1.63;

N= length(sorted);
if (color)
    figure
    % Find KS statistic in uniform domain
    plot(sorted, ([1:N]-.5)/N, 'b', 'LineWidth', 2);  
    hold on
    plot(0:.01:1,0:.01:1, 'g', 'LineWidth', 2); 
%    plot(0:.01:1, [0:.01:1]+KSdist/sqrt(N), 'r', 0:.01:1,[0:.01:1]-KSdist/sqrt(N), 'r' ); 
    plot(0:.01:1, [0:.01:1]+KSdist95/sqrt(N), 'r', 0:.01:1,[0:.01:1]-KSdist95/sqrt(N), 'r' ); 
    plot(0:.01:1, [0:.01:1]+KSdist99/sqrt(N), 'r', 0:.01:1,[0:.01:1]-KSdist99/sqrt(N), 'r' ); 
    hold off
else
    plot(sorted, ([1:N]-.5)/N, 'k', 'LineWidth', 2.0);  
    hold on
    h = plot(0:.01:1,0:.01:1, 'LineWidth', 1); 
    set(h, 'Color', [.4 .4 .4]);
    h = plot(0:.01:1, [0:.01:1]+KSdist/sqrt(N), '--', 0:.01:1,[0:.01:1]-KSdist/sqrt(N), '--' ); 
    set(h(1), 'Color', [.4 .4 .4]);
    set(h(2), 'Color', [.4 .4 .4]);
    hold off

    axis([0 1 0 1]);
end


if (sum( abs([1:N]'/N-sorted) > KSdist95/sqrt(N)))
    viol(1)= 1;
    fprintf(1, '  outside 95%%-CI. \t');
else
    viol(1)= 0;
    fprintf(1, '  ok! inside 95%%-CI. \t');
end
if (sum( abs([1:N]'/N-sorted) > KSdist99/sqrt(N)))
    viol(2)= 1;
    fprintf(1, '  outside 99%%-CI.\n');
else
    viol(2)= 0;
    fprintf(1, '  ok! inside 99%%-CI.\n');
end

