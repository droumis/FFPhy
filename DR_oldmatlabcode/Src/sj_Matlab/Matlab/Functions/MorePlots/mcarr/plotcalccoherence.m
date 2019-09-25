% Day 1 Epoch 2 coherence per pair
figure;
subplot(2,1,1); hold on;
d1e2 = zeros(size(f.output{1}(1).coherence(1:250)));
d1e2int = zeros(1,25);
avint = 0;

for i = 1:25
    
    % Compute Average
    plot(f.output{1}(i).frequency(1:250), f.output{1}(i).coherence(1:250), 'x');
    d1e2 = (d1e2 + f.output{1}(i).coherence(1:250));
    
    % Compute Integral under curve
    timestep = f.output{1}(1).frequency(2)-f.output{1}(1).frequency(1);
    int = 0;
    
    for j = 1:250
        int = int + f.output{1}(i).coherence(j);
    end
    
    d1e2int(i) = int;
    avint = avint +int;
    
end

d1e2 = d1e2/25;
avint = avint/25;
plot(f.output{1}(1).frequency(1:250), d1e2, 'r');
subplot(2,1,2);
bar(d1e2int); hold on
bar(26,avint,'r')

%__________________________________________________________________________
% Day 1 Epoch 2 coherence at frequencies in Hz

figure;
subplot(2,1,1); hold on;
d1e2 = zeros(size(f.output{1}(1).coherence(1:250)));
Hz = [1 25 50 75 100 125 150 175 200 225 250];
d1e2int = zeros(25,10);
avint = zeros(25,1);
timestep = f.output{1}(1).frequency(2) - f.output{1}(1).frequency(1);

for i = 1:25
    
    % Compute Average
    plot(f.output{1}(i).frequency(1:250), f.output{1}(i).coherence(1:250), 'x');
    d1e2 = (d1e2 + f.output{1}(i).coherence(1:250));
    
    % Compute Integral under curve
    
    for j = 1:10
        for k = Hz(j):Hz(j+1)
            d1e2int(i,j) = d1e2int(i,j) + f.output{1}(i).coherence(k)*timestep;
        end
    end
    
end

d1e2 = d1e2/25;
avint = sum(d1e2int)/25;
plot(f.output{1}(1).frequency(1:250), d1e2, 'r');
subplot(2,1,2);
plot(Hz(2:11),d1e2int', 'x'); hold on
plot(Hz(2:11), avint)

%__________________________________________________________________________
% Day 9 Epoch 2

% Compute Average
figure;
subplot(2,1,1);hold on;

d9e2 = zeros(size(f.output{1}(300).coherence(1:250)));
d9e2int = zeros(1,25);
avint = 0;

for i = 276:300
    
    plot(f.output{1}(i).frequency(1:250), f.output{1}(i).coherence(1:250), 'x');
    d9e2 = (d9e2 + f.output{1}(i).coherence(1:250));
    
    % Compute Integral under curve
    timestep = f.output{1}(300).frequency(2)-f.output{1}(300).frequency(1);
    int = 0;
    
    for j = 1:250
        int = int + f.output{1}(i).coherence(j);
    end
    
    d9e2int(i-275) = int;
    avint = avint +int;
    
end

d9e2 = d9e2/25;
avint = avint/25;
plot(f.output{1}(300).frequency(1:250), d9e2, 'r');
subplot(2,1,2);
bar(d9e2int); hold on
bar(26,avint,'r')

%__________________________________________________________________________
% Day 9 Epoch 2 coherence at frequencies in Hz

figure;
subplot(2,1,1); hold on;
d9e2 = zeros(size(f.output{1}(300).coherence(1:250)));
d9e2int = zeros(25,10);
Hz = [1 25 50 75 100 125 150 175 200 225 250];
avint = zeros(25,1);
timestep = f.output{1}(300).frequency(2) - f.output{1}(300).frequency(1);

for i = 276:300
    
    % Compute Average
    plot(f.output{1}(i).frequency(1:250), f.output{1}(i).coherence(1:250), 'x');
    d9e2 = (d9e2 + f.output{1}(i).coherence(1:250));
    
    % Compute Integral under curve
    
    for j = 1:10
        for k = Hz(j):Hz(j+1)
            d9e2int(i-275,j) = d1e2int(i-275,j) + f.output{1}(i).coherence(k)*timestep;
        end
    end
    
end

d9e2 = d9e2/25;
avint = sum(d9e2int)/25;
plot(f.output{1}(1).frequency(1:250), d9e2, 'r');
subplot(2,1,2);
plot(Hz(2:11),d9e2int', 'x'); hold on
plot(Hz(2:11), avint)

%__________________________________________________________________________


