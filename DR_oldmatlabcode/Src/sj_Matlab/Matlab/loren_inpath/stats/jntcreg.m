% JNTCREG.M
% Updated for ver 5.0
% Example Joint Conf Region  uses modelcr1.m curvefit.m jointcr.m,sumsscr.m
% example from es413/512 assig4

%x=[1 4 7 10 13 16]';
%y=[185 618 847 901 1000 955]';
x=[4 4 4 16 16 16 ]';
y=[634 641 618 934 966 1011]';
disp(' enter guess for 2 params e.g. 900 .2 ')
%pin=readv(2)
pin=[900  .2] % init guess

[pfinal,options,error,jac]=curvefit('modelcr1',pin,x,y);
pfinal
ymodel=y+error;
% calc regression data
[std,varresid,r2,cor,vcv,varinf]=regdata(pfinal,ymodel,y,jac);

% calc joint sum of squares region

jointcr(x,y,pfinal,1,'sumsscr',std,0.05);
