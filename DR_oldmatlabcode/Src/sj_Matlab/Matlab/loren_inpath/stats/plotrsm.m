function z=plotrsm(xo,b0,b,B,lohi,ngc)
%
%function z= plotrsm(xo,b0,b,B,lohi,{ngc})
% Plot combined mesh +  contour plots of RSM models
% Hard Copy: print -dljet3
%     xo = vector of optimum factor values with ZERO 0 for 2 factors
%	    which will form x and y axes. Example xo=[1.3, 0, 3.4, 0]
%	    xo(2) and xo(4) used for x and y axes.xo(1),xo(3) fixed.
%	    Height z calculated from RSM model.
%	    Usually xo= Optimum x vector from optrsm or smaxrsm
%     b0,b,B = Parameter Set from leasrsm(1) B usually lower diag
%     lohi   = Constraint matrix (n*2) for each xo(i). Used to calc
%		  ploting range for x and y axis.
%     ngc = [number==grid points ,number==contour lines],default [20,5]
%     z (output) can be re mesh(ed) or contour(ed) with different angles,etc
%	   A. Jutan. UWO 1993
%

if nargin == 5 , npts=20;nc=5;else, npts=ngc(1);nc=ngc(2);end
 [nlohi,junk]=size(lohi);
 if length(b) ~= nlohi ,error('lohi must have #rows = Length(xo)'),end;
 xo=xo(:)'; % force row vector
 % check if B lower diag make full BF
 nx=length(xo);
if B(nx,nx) == 0
      c=B+B';
      d=diag(diag(c));
      BF=c-d/2; % full B matrix
      else
      BF=B;
end
% find zeros in xo vector(corresponds to x and y axis elements)
i=find(xo == 0);
if length(i)~=2,error('only two varibles in xo vector(0"s) can be plotted'),end
% use lohi matrix to calculate ranges to plot
xr=lohi(i(1),:);yr=lohi(i(2),:);
% set up x,y grid with npts points
xgrid=linspace(xr(1),xr(2),npts);
ygrid=linspace(yr(1),yr(2),npts);
for j=1:npts
for k=1:npts
      xo(i)=[xgrid(j),ygrid(k)]; % replace Zero's in xo with grid points
      z(k,j)=b0 + xo*b'  + xo*BF*xo'; % RSM Model height =z axis ,k= yaxis
end
end
% plot Surface and Contour
%subplot(211);
%mesh(z,[-10,25]);
%subplot(212);
%cs= contour(z,nc,xgrid,ygrid);
%clabel(cs);
%xlabel(['factor ',int2str(i(1))]); ylabel(['factor ',int2str(i(2))]);
%grid
%pause
%subplot(111)

% plot combined surface and contour plot
figure(1)
cs=contour(xgrid,ygrid,z,nc);grid;% calculate labels using user grid
clabel(cs); % put labels on contours
figure(2)
meshclabel(xgrid,ygrid,z);% plot mesh+ contour plot
view([-10,25])
xlabel(['factor ',int2str(i(1))]); ylabel(['factor ',int2str(i(2))]);
%set(gca,'GridLineStyle',':'); % set linestyle of grid to ..