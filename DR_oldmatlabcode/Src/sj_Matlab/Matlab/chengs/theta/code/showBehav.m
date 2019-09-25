function showBehav
%function showBehav
%  
%  Show behavior

setRoot;
epfile= sprintf('%s/data/epochs-%s', root, selectid);
load(epfile);
ep= epochs;
nep= length(ep.rat);

for ie=1:nep
    rat= ep.rat{ie}; d= ep.num(ie,1); e= ep.num(ie,2);
    fprintf(1, '%s [%d %d], %d/%d\n',rat,d,e,ie,nep);

    % load behavior data
    if ie==1 | ~strcmp(rat, ep.rat{ie-1}) | d~= ep.num(ie-1,1)
        bfile= sprintf('%s/%s/data2/behavdata%.2d.mat',root,rat,d);
        load(bfile);
    end
    bd= behavdata{d}{e};
    hp= plot(bd.xpos([1:1000:end]), bd.ypos([1:1000:end]), '.');
    set(hp, 'Color', .5*[1 1 1]);
    axis ij
    if isfield(bd, 'ripple')
        ind= find(bd.ripple);
        fprintf(1, '  %d ripple events.\n', sum(diff(ind)>1)+1);
        hold on
        hp= plot(bd.xpos(ind), bd.ypos(ind), 'r.');
        hold off
%        keyboard
    else
        fprintf(1, '  no ripple field\n');
    end
    pause
end
