function showMix(model, fig, style, Tmax)

if nargin < 2
    fig= 1;
end

if nargin < 4
    Tmax= length(model.GaussComp{1}.mx);
else
    if Tmax > length(model.GaussComp{1}.mx)
	warning('requested more data than available in file');
    end
end

nMix= length(model.GaussComp);

for f=1:nMix
    figure(fig+f-1);

    subplot(2,3,1);
    plot(model.GaussComp{f}.mx(1:Tmax), style);
    hold on
    title('mx');
    
    subplot(2,3,4);
    plot(model.GaussComp{f}.my(1:Tmax), style);
    hold on
    title('my');
    
    subplot(2,3,2);
    plot(model.GaussComp{f}.Sx(1:Tmax), style);
    hold on
    title('Sx');
    
    subplot(2,3,5);
    plot(model.GaussComp{f}.Sy(1:Tmax), style);
    hold on
    title('Sy');
    
    subplot(2,3,3);
    plot(2/pi*atan(model.GaussComp{f}.r(1:Tmax)), style);
    hold on
    title('r');
    
    subplot(2,3,6);
    plot(model.GaussComp{f}.alpha(1:Tmax), style);
    hold on
    title('alpha');
end
