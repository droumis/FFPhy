function showGauss(model, fig, style, Tmax)

if nargin < 2
    fig= 1;
end

if nargin < 4
    Tmax= length(model.mx);
else
    if Tmax > length(model.mx)
	warning('requested more data than available in file');
    end
end

figure(fig);

subplot(2,3,1);
plot(model.mx(1:Tmax), style);
hold on
title('mx');

subplot(2,3,4);
plot(model.my(1:Tmax), style);
hold on
title('my');

subplot(2,3,2);
plot(model.Sx(1:Tmax), style);
hold on
title('Sx');

subplot(2,3,5);
plot(model.Sy(1:Tmax), style);
hold on
title('Sy');

subplot(2,3,3);
plot(2/pi*atan(model.r(1:Tmax)), style);
hold on
title('r');

subplot(2,3,6);
plot(exp(model.alpha(1:Tmax)), style);
hold on
title('exp(alpha)');

