out{1} = control_gammapower_spikemodulation('/data13/mcarr/Bon/','bon',3:10);
out{2} = control_gammapower_spikemodulation('/data13/mcarr/Fra/','fra',2:12);
out{3} = control_gammapower_spikemodulation('/data13/mcarr/Ten/','ten',1:7);

bin = -pi:pi/4:pi;
rip1 = []; rip3 = []; norip1 = []; norip3 = [];
for an = 1:length(out)
    for d = 1:length(out{an})
        for e = 1:length(out{an}{d})
            if ~mod(e,2)
                if ~isempty(out{an}{d}{e}.ripple_ca1_mod)
                    rip1 = [rip1; out{an}{d}{e}.ripple_ca1_mod];
                    rip3 = [rip3; out{an}{d}{e}.ripple_ca3_mod];
                    norip1 = [norip1; out{an}{d}{e}.no_ripple_ca1_mod];
                    norip3 = [norip3; out{an}{d}{e}.no_ripple_ca3_mod];
                end
            end
        end
    end
end
rip1(isnan(rip1)) =[]; rip3(isnan(rip3))=[]; norip1(isnan(norip1))=[]; norip3(isnan(norip3))=[];
rip1(rip1<-pi) = -pi; rip1(rip1>pi) = pi; norip1(norip1<-pi) = -pi; norip1(norip1>pi) = pi;
rip3(rip3<-pi) = -pi; rip3(rip3>pi) = pi; norip3(norip3<-pi) = -pi; norip3(norip3>pi) = pi;

h1 = histc(rip1,bin); h3= histc(rip3,bin); h1_norip = histc(norip1,bin); h3_norip = histc(norip3,bin);

h1 = [h1(1:end-1)./length(rip1); h1(1:end-1)./length(rip1); h1(1)./length(rip1)];
h3 = [h3(1:end-1)./length(rip3); h3(1:end-1)./length(rip3); h3(1)./length(rip3)];
h1_norip = [h1_norip(1:end-1)./length(norip1); h1_norip(1:end-1)./length(norip1); h1_norip(1)./length(norip1)];
h3_norip = [h3_norip(1:end-1)./length(norip3); h3_norip(1:end-1)./length(norip3); h3_norip(1)./length(norip3)];

x = [bin(1:end-1) 2*pi+bin(1:end)];

figure
plot(x,h1,'r',x,h1_norip,'k')
set(gca,'xlim',x([1 end]),'ylim',[0.1 0.16],'xtick',x(1:2:end),'xticklabel',x(1:2:end)*180/pi)
legend([{'Ripple'},{'No Ripple'}])
xlabel('Gamma Phase')
ylabel('Proportion of spikes')

figure
plot(x,h3,'r',x,h3_norip,'k')
set(gca,'xlim',x([1 end]),'ylim',[0.1 0.16],'xtick',x(1:2:end),'xticklabel',x(1:2:end)*180/pi)
legend([{'Ripple'},{'No Ripple'}])
xlabel('Gamma Phase')
ylabel('Proportion of spikes')

%Rayleigh test on gam3 and gam1 to test for significant modulation:

%RIPPLE
%   gam1: theta = -0.125, p = 0.001
%   gam3: theta = -0.8536, p = 0.001

%NO RIPPLE
%   gam1: theta = -0.464, p >0.15
%   gam3: theta = -1.07,  p <1e-5


%% COMPARE DEPTH OF MODULATION DURING SWRS AND POWER MATCHED GAMMA TIMES

q = [(max(h1)-min(h1))./(max(h1)+min(h1)) (max(h3)-min(h3))./(max(h3)+min(h3))...
    (max(h1_norip)-min(h1_norip))./(max(h1_norip)+min(h1_norip)) (max(h3_norip)-min(h3_norip))./(max(h3_norip)+min(h3_norip)) ];

nboot = 1000;
Q = nan(length(nboot),4);
for s = 1:nboot
    boot1 = rip1(ceil(length(rip1)*rand(length(rip1),1)));
    boot3 = rip3(ceil(length(rip3)*rand(length(rip3),1)));
    q1 = histc(boot1,bin);  q3 = histc(boot3,bin);
    q1 = [q1(1:end-1)./length(rip1); q1(1:end-1)./length(rip1); q1(1)./length(rip1)];
    q3 = [q3(1:end-1)./length(rip3); q3(1:end-1)./length(rip3); q3(1)./length(rip3)];
    
    Q(s,[1 2]) = [(max(q1)-min(q1))./(max(q1)+min(q1)) (max(q3)-min(q3))./(max(q3)+min(q3))];
    
    boot1 = norip1(ceil(length(norip1)*rand(length(norip1),1)));
    boot3 = norip3(ceil(length(norip3)*rand(length(norip3),1)));
    q1 = histc(boot1,bin);  q3 = histc(boot3,bin);
    q1 = [q1(1:end-1)./length(norip1); q1(1:end-1)./length(norip1); q1(1)./length(norip1)];
    q3 = [q3(1:end-1)./length(norip3); q3(1:end-1)./length(norip3); q3(1)./length(norip3)];
    
    Q(s,[3 4]) = [(max(q1)-min(q1))./(max(q1)+min(q1)) (max(q3)-min(q3))./(max(q3)+min(q3))];
    
end

%CA1:   sum(Q(:,1)>mean(Q(:,3)))./nboot
% Ripple modulation is greater than nonripple times p<0.001

%CA3:   sum(Q(:,2)>mean(Q(:,4)))./nboot
% No significant difference

figure
bar([1 2 4 5],mean(Q(:,[2 4 1 3])),1,'b')
hold on
errorbar2([1 2 4 5],mean(Q(:,[2 4 1 3])),std(Q(:,[2 4 1 3])),'k')
set(gca,'xtick',[1 2 4 5],'xticklabel',[{'CA3 SWR'},{'CA3 No SWR'},{'CA1 SWR'},{'CA1 No SWR'}])
ylabel('Depth of modulation')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_spikemodulation_controlforgammapower.pdf', m, d, y);
print('-dpdf', savestring)
