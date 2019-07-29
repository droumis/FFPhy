
% now that all the prelim power var corr tf maps are plotted.. 
% choose nts, var, tfbox
% make scatter with fit
% get corrcoef r2 and pval per anim

% aside from these continuous corrs.. i need to  also look at
% outbound-inbound and rewarded-unrewarded.
% these are more inline with the experimental conditions described in the
% book.. since i should already have the mean power tfmaps for each condition.. 
% now i want the diff tfmaps, then choose nts, var, tfbox
% once i choose the tfbox.. i can get a distribution for each animal, condition,
% and then test if the distributions are different
% what tests to use? outside of matlab?

% didn't i already do the diff power/phase sig testing? i did.. i just need
% to redo that and actually do the permutation testing again..
% start with D10, D12

%%%%%
%%%%% Results part 1A:
%%%% need to pick nt,var,tfbox for these
% 1. ca1, mec superficial, deep vs overall performance
% 2. ca1, mec superficial, deep vs speed
% 3. ca1, mec superficial, deep vs since epoch start
%%still need to do these:
% 4. ca1, mec superficial, deep vs lindist from well (indicative of both
% task phase and reward well proximity)
% 5. ca1, mec superficial, deep vs ca1 ripple peak z power

%%%% need to run the diff meanpower + shuffle for these
% 4. ca1, mec superficial, deep :: outbound vs inbound
% 5. ca1, mec superficial, deep :: rewarded vs unrewarded

%%%% need to pick nt,var,tfboxA, tfboxB and run these
% 6. MEC post ripple high frequency power vs preripple low freq power 
% --- mec 'context' effects ripple-induced activity

% 7. CA1 ripple frequency power vs preripple MEC low freq power
% --- mec 'context' effects ca1 ripples
% this could go either way.. and would interact with the interpretation of
% the previous thing..

%%%%% 
%%%%% Results part 2
% need to find ripple chains to asses how many i have from each animal
% to what extent do i need to show that two tetrodes are in different
% modules?

%{
mice
anna used 30-50Hz for slow gamma
Analysis of SWRs was restricted to periods of extended immobility, 
whenmouse velocity was <1 cm/s for 30 s or more.
- SO i likely will see a lot more theta and less slow gamma than her
    since i have a < 4cm thresh without a immobility time duration threshold
- Quantification of slowgamma power reflects the averageZscored power over a 
    30–50 Hz frequencyband at 1–100 ms after ripple detection, averaged over SWRs.
- Baseline coherence was quantified as the average coherence overthe 30–50 Hz 
    frequency band 400–300 ms before SWR detection, peakcoherence was calculated 
    0–100 ms after SWR detection, and the SWR-induced coherence increase was the 
    difference between baseline and peakcoherence.
- All data were normally distributed as shown by Lillieforstest, and variances 
    between groups were similar as shown by Bartlett’s test.For these reasons, 
    we used two-tailed paired and unpaired t tests and one-and two-way ANOVA. 
    Post hoc testing was done with Holm-Sidak correctionfor multiple comparisons

maggie
- We found that in addition to the expected increase in ripplepower, 
    there was a substantial increase in a 20–50 Hz slowgamma band in both CA3 and CA1. 
    There was also an increaselow frequency power (<20 Hz) in CA1, but not in CA3 
... so i should expect like 0-50Hz activity during ca1 swrs..

oliva
- rips are a slow (5–15 Hz) sharp-wave, ‘ripple’’ oscillation (150–200 Hz), 
    and a slow ‘‘gamma’’ oscillation (20–40 Hz). sharp wave (40–-120 ms) (8.3 - 25 Hz)
    - they later in intro cite ripple range as  (120-250 Hz) and then (120–200 Hz)
    - in results, cite sharp wave as  2–15 Hz. then later 5–10 Hz

- We show that increased po-wer in the 20–40 Hz band does not reflect an entrain-ment 
    of CA1 and CA3 neurons at gamma frequencybut the power envelope of overlapping ripples.
- Recent reports described that SPW-Rs are phase-coupledwith a power spectral peak in the 
    slow gamma band(20–40 Hz). From this relationship, the authors assumed that rip-ple 
    power is modulated by the phase of a genuine gamma oscil-lation and associated with a 
    higher fidelity replay of past experi-ences and of place cell trajectories 
    (Carr et al., 2012; Pfeiffer andFoster, 2015). 
%}

%{
stories:

PORTAL STATES
% the relevance of the response in superficial MEC could be the ongoing
representation of current location.. and that if there's a positive modulation in
activity preferentially away from vs close to the reward well, then
it could be argued that the replay while the animal is in a navigational
state combines with the ongoing representation.. thereby creating these
'portal states' that fuse that location of the ripple with state transition
ability to enter a remote state.. 
- would this mean that that cell representating a remote state is then more
likely to firing that the swr location, even outside of swrs?
- i think the effect would be huge because of the speed and ypos corr.. 
- swrs are forging these portal states

- i could show decoding examples of remote ca1 replay while the mec keeps
firing in a way predicted by it's grid fields while the animal is
navigating, vs a lack of response or even a reversal of response to ripples
while the animal is at the reward wells

% would i have to show that mec continues to represent the current location
during replay that BOTH includes local posterior density and remote
posterior density?

- ripple propagation efficacy.. i would think that transient gamma or ripple band would be
the proxy for spike propagation as opposed to the low frequency, which
would be the pre or post 'context'


what is expected in MEC during a ripple?
- superficial MEC ecx firing pre ripple only for longer ripples, not
shorter
- 
%}











