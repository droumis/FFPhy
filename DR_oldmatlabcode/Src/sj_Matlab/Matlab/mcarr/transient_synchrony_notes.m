% Analyses for Transient Synchrony

%% Figure 1
    % Figure 1 stays as is.
    
%% Figure 2
    % A: runcalcripstats.m
    % B: runcalcripstats.m
    % C: runcalcripstats_new.m
        % changed threshold to 5 std
        % Anova:
            % Time 0-400ms,peak > baseline at p<1e-5
            % Time -100ms > baseline at p<0.01
	% D: runcalcripstats_new.m
        % changed threshold to 5 std
        % Anova:
            % Time 0-400ms,peak > baseline at p<1e-5
	% E: plotreplay_example.m
        % Replotted on different time scale
    % F: runcalcripstats_new.m
        % changed threshold to 5 std
        % Anova:
            % Time 300ms > baseline at p<0.01
            % Time 100, 200ms, peak > baseline at p<1e-5
    % G: runcalcripstats_new.m
        % changed threshold to 5 std
        % Anova:
            % Time 100, 200ms,peak > baseline at p<0.001
    % H: run_gamma_ripple_glm.m
        % changed threshold to 5 std
        % 69% of sessions significant for CA3, 89% for CA1

%% Figure 3
    % A: runcalcripstats.m
        % changed threshold to 5 std
    % B: runcalcripstats.m
        % changed threshold to 5 std
        % Anova
            % 100-200ms >baseline at p<1e-5
            % 0,300-400ms > baseline at p<0.001
            % -100ms > baseline at p < 0.05
            %mean baseline for phase locking: 0.82
    % C: runripcoherencestats.m
        % changed threshold to 5 std
        % Anova
            % 100ms >baseline at p<1e-5
            % 0 > baseline at p<0.01
            % baseline: 0.698
    % D: runcalcripstats.m
        % changed threshold to 5 std
    % E: runcalcripstats_new.m
        % Anova
            % Time 0-400ms,peak > baseline at p<1e-5
            % baseline = 0.58
    % F: runripcoherencestats.m
            % changed threshold to 5 std
            % Anova
            	% Time 0-400ms > baseline at p<1e-5
                % baseline: 0.545
    % G: runcalcripstats.m
        % changed threshold to 5 std
        % Anova
            % CA1 - CA1:
                % 100ms >baseline at p<1e-5
                % 0ms > baseline at p<0.001
                % 200-300ms > baseline at p < 0.05
                %mean baseline for phase locking: 0.91
            % CA3-CA3:
                % Time 100, 200ms,peak > baseline at p<0.001
                % 100ms >baseline at p<1e-5
                % 0ms,200-400ms > baseline at p<0.05
                %mean baseline for phase locking: 0.98
    % H: runcalcripstats_new.m
        % changed threshold to 5 std
        % Anova:
            % CA1-CA1
                % Time 0-100ms, peak > baseline at p<1e-5
                % Time 200ms > baseline at p<0.05
                % baseline = 0.69
            % CA3-CA3
                % Time 0-400ms,peak > baseline at p<1e-5
                % baseline = 0.78
    % I: run_gamma_ripple_glm.m
        % changed threshold to 5 std
        % 48% sessions have a significant relationship

%% Figure 4
    % A: runcalcspikingmodinripples.m
        % changed threshold to 3 std
        %Rayleigh test on gam3 and gam1 to test for significant modulation:
            %   gam3: theta = -0, p = 0.01
            %   gam1: theta = 1.74, p = 1e-5

        %Run permutation test to determine whether CA3 comes before CA1
            % p <1e-5


    % B: runcalcspikingmodinripples.m
        %Bootstrap test to determine if modulation is greater during SWRs than 500ms preceding
            %CA1:   p<1e-5
            %CA3:   p>0.5

        %Rayleigh test on gam3_preceding and gam1_preceding to test for significant modulation:
            %   gam3_preceding: theta = -0.8, p = 0.001
            %   gam1_preceding: theta = 0.1, p = 0.001

    % C:

    % D: runcalcspikingmodinripples.m
        %Rayleigh test on gam3_first and gam1_first to test for significant modulation:
            %   gam3_first: theta = -0.24, p = 0.001
            %   gam1_first: theta = 0.6, p = 0.001

        %Run permutation test to determine whether CA3 comes before CA1
            % p <0.01

    %Additional statistics:
        
        %Run permutation test to determine whether preceding comes before all
            %CA1: p < 0.001
            %CA3: p < 0.001

        %Run permutation test to determine whether preceding comes before first
            %CA1: p > 0.3
            %CA3: p < 0.03

        %Run permutation test to determine whether first comes before all
            %CA1: p < 0.03
            %CA3: p > 0.4

 
%% Figure 5

    % A:
    % B:
    % C: paircorr.m
        %FOR RUN: gamma is significantly correlated with distance between
        %placefield peaks for last three groups, p<0.01 p<0.001 p<0.001

        %time is significantly correlated with distance between placefield peaks
        %for last two groups, p<0.001

        %gamma is significantly more correlated with time for place fields with
        %distances far apart, p<0.001

        % 75% of cells that are closer than 24cm fire within the same cycle
%% Figure 6

%% Figure 7

    % A: runcalcripstats_new.m
        % changed threshold to 5 std
        % Anova
            % -100-400ms,peak > baseline at p<1e-5

    % B: runcalcripstats_new.m
        % changed threshold to 5std
        % Anova
            % 0-400ms,peak > baseline at p<1e-5
            
    % C: run_gamma_ripple_glm.m
    	% changed threshold to 5 std
        % 70% of sessions significant for CA3, 89% for CA1
        
    % D: runcalcripstats.m
        % changed threshold to 5 std
        % Anova
            % None significantly different
            %mean baseline: 0.87
	% E: runcalcripstats_new.m
        % changed threshold to 5 std
        % Anova
            % Time 200-400ms > baseline at p<0.05
            % Time 100ms, peak > baseline at p<1e-5
            % baseline = 0.62
% F:
% G:
    % H: run_gamma_ripple_glm.m
        % changed threshold to 5 std
        % 74% of sessions significant
    % I: runcalcspiking modinripples.m
        % changed threshold to 5 std
        % X spikes from X CA1 neurons.
        % X spikes from X CA3 neurons.
    
        %Rayleigh test on gam3 and gam1 to test for significant modulation:
            %   gam3: theta = -0.49, p < 0.01
            %   gam1: theta = 2.36,  p <0.01

        %Run permutation test to determine whether CA3 comes before CA1
            % p > 0.2
    
        %Run permutation test to determine whether preceding comes before all
            %CA1: p > 0.1
            %CA3: p<0.05
        
    	%Run permutation test to determine whether preceding comes before first
            %CA1: p > 0.1
            %CA3: p > 0.5
        
    	%Run permutation test to determine whether first comes before all
            %CA1: p > 0.5
            %CA3: p < 0.05
        
    	%Run bootstrap test to determine whether significant difference in
        %angle mean between awake and quiescence
            %CA1: p> 0.1
            %CA3: p< 0.05
        
% J:

        
%% Supplemental Figure 1
    % Unchanged
    
%% Supplemental Figure 2
    % runcalchighgammapower.m
        % changed threshold to 5std
        % Main effect of time and gamma, interaction between two with slow 
        % gamma significantly larger than fast gamma, p<1e-5

%% Supplemental Figure 3

%% New Figures
    % High gamma coherence during SWRs
        % runcalchighgammapower.m
            % Using 3 std threshold, no significant change in CA1-CA3 high
            % gamma coherence during SWRs
                % baseline: 0.54
    % Ripple coherence during SWRs
        % runcalchighgammapower.m
            % Using 3std threshold, significant decrease in CA1-CA3
            % coherence in the 150-250Hz band during SWRs
            % kruskalwallisanova
                % Time -200,-100,0, 100ms < baseline at p<0.05
                % baseline = 0.51
    % Distribution of SWR frequency
        % runcalcinstantaneousfrequency.m
            %3 std SWRs
            % 2% of SWRs have frequencies that are best described as <20Hz, the mean
            % SWR frequency is 29.6 +- 0.17 and there is no evidence of a bimodal
            % distribution. It seems that during SWRs there is a single frequency
            % closer to 30Hz rather than two sepearable 20Hz and 50Hz oscillations
            % median frequency = 29.3Hz, 25th prctile = 25.6Hz, 75th prctile = 33Hz
    % Distribution of SWR duration
        % runcalcrippleduration.m
            % No SWR was shorter than 0.05 seconds suggesting that 20Hz is
            % a reasonable cutoff.
    % Cross frequency coupling
        % runcalccrossfrequencycoupling.m
            %Example modulation for one session and distribution of peak phase
    % Gamma tail
        % runcalcgammatail.m
            %On average, gamma returns to baseline 1ms after ripple power returns to
            %baseline. This effect was not significantly different for all animals
            %by signrank (p<0.001, p<0.05, p>0.1) or ttest.

%% Additional Notes

% There is no real difference between the novel an familiar track beyond
% the enhancement of SWRs as previously reported.
