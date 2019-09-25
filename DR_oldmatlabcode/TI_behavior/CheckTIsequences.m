







randhome = load('/home/droumis/MATLAB/TI_behavior/T03/ph2_setrnd/16_T03_ph2_setrnd.txt');

i =1;

ab(i) = length(randhome(randhome(1:150) == 1));
    bc(i) = length(randhome(randhome(1:150) == 2));
    de(i) = length(randhome(randhome(1:150) == 3));
    ef(i) = length(randhome(randhome(1:150) == 4));
    CD(i) = length(randhome(randhome(1:150) == 5));
    bd(i) = length(randhome(randhome(1:150) == 6));
    be(i) = length(randhome(randhome(1:150) == 7));
    ce(i) = length(randhome(randhome(1:150) == 8));
    af(i) = length(randhome(randhome(1:150) == 9));
    

    figure
    bar([1 2 3 4 5 6 7 8 9], [mean(ab) mean(bc) mean(de) mean(ef) mean(CD) mean(bd) mean(be) mean(ce) mean(af)]);  hold on;
    errorbar2([1 2 3 4 5 6 7 8 9], [mean(ab) mean(bc) mean(de) mean(ef) mean(CD) mean(bd) mean(be) mean(ce) mean(af) ],  [stderr(ab) stderr(bc) stderr(de) stderr(ef) stderr(CD) stderr(bd) stderr(be) stderr(ce) stderr(af)] , 0.3, 'k')
