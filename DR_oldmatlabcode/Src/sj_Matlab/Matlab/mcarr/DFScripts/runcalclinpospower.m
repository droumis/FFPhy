% Animal selection
animals = {'Five','Seven','Eight','Ten'};

% Epoch selection
epochfilter = [];
epochfilter{1} = '($exposure <= 3)';
epochfilter{2} = '($exposure >= 8)';

% Tetrode selection
tetfilter = '(isequal($area, ''CA1'') & ($maxcell == 1))';

% Iterator selection
iterator = 'eeganal';

%% Low Gamma, Choice point

% Time selection
timefilter = {{'getlinstate', '($traj == 1 | $traj == 3)', 6}, {'getlindist', '$lindist>45 & $lindist<75'}};

% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'eegtetrodes',tetfilter,'iterator', iterator);
v = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter);

f = setfilterfunction(f, 'calclinpospower', {'eeg'},'fpass',[30 55]);
v = setfilteriterator(v, 'singleepochanal');
v = setfilterfunction(v, 'calclinvelocitysegment', {'linpos'},'minO',0.5,'maxO',10);

f = runfilter(f);
v = runfilter(v);

lgcp = f;
cpv = v;

clear f

%% High Gamma, Choice point

% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'eegtetrodes',tetfilter,'iterator', iterator);

f = setfilterfunction(f, 'calclinpospower', {'eeg'},'fpass',[65 90]);

f = runfilter(f);

hgcp = f;

clear f

%% Low Gamma, Not Choice Point

% Time selection
timefilter = {{'getlinstate', '($traj == 2 | $traj == 4)', 6}, {'getlindist', '$lindist>45 & $lindist<75'}};

% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'eegtetrodes',tetfilter,'iterator', iterator);
v = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter);

f = setfilterfunction(f, 'calclinpospower', {'eeg'},'fpass',[30 55]);
v = setfilteriterator(v, 'singleepochanal');
v = setfilterfunction(v, 'calclinvelocitysegment', {'linpos'},'minO',0.5,'maxO',10);

f = runfilter(f);
v = runfilter(v);

lgin = f;
inv = v;

clear f

%% High Gamma, Not Choice Point

% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'eegtetrodes',tetfilter,'iterator', iterator);

f = setfilterfunction(f, 'calclinpospower', {'eeg'},'fpass',[65 90]);

f = runfilter(f);

hgin = f;

clear f

%% Combine Data

Conditions = [{lgcp} {hgcp} {lgin} {hgin}];
% Combine data across animals and across Novel and Familiar conditions

value = {};
novel = {};

for c = 1:length(Conditions)

    f = Conditions(c);
    f = f{1};
    
    %Combine for novel and familiar
    e1 = [];
    e2 = [];

    
    for an = 1:length(f)
        for g = 1:length(f(an).output{1})
        e1 = [e1 f(an).output{1}(g).zPower];
        end
        for g = 1:length(f(an).output{2})
        e2 = [e2 f(an).output{2}(g).zPower];
        end
    end

    Novel = [ones(length(e1),1); 2*ones(length(e2),1)];
    
    value{c} = [e1 e2];
    novel{c} = Novel;
end

% Combine data across Choice point and inbound conditions
LowGamma = [value{1}'; value{3}'];
LowGammaNovel = [novel{1}; novel{3}];
LowGammaChoice = [ones(length(value{1}),1); 2*ones(length(value{3}),1)];

HighGamma = [value{2}'; value{4}'];
HighGammaNovel = [novel{2}; novel{4}];
HighGammaChoice = [ones(length(value{2}),1); 2*ones(length(value{4}),1)];

% Combine data across Low and High Gamma
LinPower = [LowGamma; HighGamma];
Novelty = [LowGammaNovel; HighGammaNovel];
ChoiceP = [LowGammaChoice; HighGammaChoice];
Gamma = [ones(length(LowGamma),1); 2*ones(length(HighGamma),1)];

%% Combine Velocity
Conditions = [{cpv} {inv}];

%Combine data across animals and across Novel and Familiar conditions

vvalue = {};
vnovel = {};

for c = 1:length(Conditions)
        
    Novel = [ones(length(e1),1); 2*ones(length(e2),1)];
    
    v = Conditions(c);
    v = v{1};
    
    %Combine for novel and familiar
    e1 = [];
    e2 = [];
    
    for an = 1:length(v)
        for g = 1:length(v(an).output{1})
        e1 = [e1; v(an).output{1}(g).vout];
        end
        for g = 1:length(v(an).output{2})
        e2 = [e2; v(an).output{2}(g).vout];
        end
    end

    Novel = [ones(length(e1),1); 2*ones(length(e2),1)];
    
    vvalue{c} = [e1; e2];
    vnovel{c} = Novel;
end

%Combine data across Choice point and inbound conditions
Velocity = [vvalue{1}; vvalue{2}];
Novel = [vnovel{1}; vnovel{2}];
Choice = [ones(length(vvalue{1}),1); 2*ones(length(vvalue{2}),1)];

%% Subsample velocity so distributions same for both novelty and choice

inb = Velocity;
inb(Choice==1) = NaN;
outb = Velocity;
outb(Choice==2) = NaN;

inb(inb < 10 | inb > 40) = NaN;
outb(outb < 10 | outb > 40) = NaN;
tempoutb = outb;
tempinb = inb(~isnan(inb));
subsample = zeros(length(tempinb),1);

for i = 1:length(tempinb)
    temp = abs(tempoutb - tempinb(i));
    subsample(i) = (find(temp == min(temp)));
    tempoutb(find(temp == min(temp))) = NaN;
end
outb(:) =NaN;
outb(subsample)=1;
outb(isnan(outb))=0;
inb(~isnan(inb))=1;
inb(isnan(inb))=0;


nov = Velocity;
nov(Novel==2) = NaN;
fam = Velocity;
fam(Novel==1) = NaN;

nov(nov < 10 | nov > 40) = NaN;
fam(fam < 10 | fam > 40) = NaN;
tempfam = fam;
tempnov = nov(~isnan(nov));
subsample = zeros(length(tempnov),1);

for i = 1:length(tempnov)
    temp = abs(tempfam - tempnov(i));
    subsample(i) = (find(temp == min(temp)));
    tempfam(find(temp == min(temp))) = NaN;
end
fam(:) =NaN;
fam(subsample)=1;
fam(isnan(fam))=0;
nov(~isnan(nov))=1;
nov(isnan(nov))=0;

LinPower = LinPower((logical([nov;nov]) | logical([fam;fam])) & (logical([outb;outb]) | logical([inb;inb])));
Novelty = Novelty((logical([nov;nov]) | logical([fam;fam])) & (logical([outb;outb]) | logical([inb;inb])));
ChoiceP = ChoiceP((logical([nov;nov]) | logical([fam;fam])) & (logical([outb;outb]) | logical([inb;inb])));
Gamma = Gamma((logical([nov;nov]) | logical([fam;fam])) & (logical([outb;outb]) | logical([inb;inb])));

%% Run Anova
varnames = {'Novelty','Choice','Gamma'};
[p,table,stats] = anovan(LinPower,{Novelty,ChoiceP,Gamma},'varnames',varnames,'model','full');

%% Plot Comparisons

figure
multcompare(stats,'dimension',1);
pause
multcompare(stats,'dimension',2);
pause
multcompare(stats,'dimension',3);
pause
multcompare(stats,'dimension',[1 2]);
pause
multcompare(stats,'dimension',[1 3]);
pause
multcompare(stats,'dimension',[2 3]);
pause
multcompare(stats,'dimension',[1 2 3]);