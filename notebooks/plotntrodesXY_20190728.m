


%{
Part 1
Hypothesis: TFcontext, expVar --> MECripple propagation (FastGamma, Ripple, MUspikesSWR)

*** ASSIGN /nt area, IDENTIFY TFcontext, expVar
*** can i get help getting the replay scores?
*** the results plotted are like unit tests for the noisedebugging... so i
won't really know when to stop noisedebugging unless i can test.. so i
should do the TFmapZ cXYplots first then noise debugging

#1NOW
=== cXYplots: plot ntrodes at relative cannula position =====
    - fix whatever is happening to mislabel the ntrode numbers, as apparent with D10
    - add ntrode relative XY plot position to tetinfo struct based on each animals cannula
    - plot the ntrodes in the xy subfigure position specified in the tetinfostruct

=== TFmapZ cXYplots - combine all into one designmatrix and one call per @type
= output: /@type [time freq ::]

wtrack    
    #2NOW
    - ::expvarCat /nt @mean @ITPC $time all, outbound, inbound, rewarded, unrewarded,
                                    proximalWell, distalWell, wtrackEp2, wtrackEp4, 
                                    ripnum, day, exposure *********************************
    - ::expvarCatDiff /nt @dmean @dITPC $condition outbound-inbound, rewarded-unrewarded, 
                                            proximalWell-distalWell,wtrackEp2-wtrackEp4****
    - ::expvarCont /nt @corr $swr hd, speed, xpos, ypos, learn, perform, timeDay, timeEp
    - ::tfbvarCont /nt @corr $swr {pre prepre swr post postpost, rip gammas beta theta}
    - ::swrvarCont /nt @corr $swr std, duration, MUspikesSWR************************************
    
    
    #5NEXT
    - ::swrvarCont /nt @corr $swr CA1replaySpeed, CA1replayCoverage, CA1replayFidelity
    - ::swrvarCat /nt @corr $swr CA1reverse, CA1forward, CA1remote, CA1local
    - ::swrvarCatDiff /nt @corr $swr CA1reverse-forward, CA1remote-local
    
#4NEXT
sleep
    - ::expvarCat /nt @mean $time all, sleepEp1, sleepEp3
    - ::expvarCatDiff /nt @dmean $condition sleepEp1-sleepEp3
    - ::expvarCont /nt @corr $swr timeDay, timeEpoch
    - ::tfbvarCont /nt @corr $swr pre, prepre, swr, post, postpost,rip,gammas,beta,theta
    - ::swrvarCont /nt @corr $swr std, duration, MUspikesSWR
    - ::swrvarCont /nt @corr $swr CA1replaySpeed, CA1replayCoverage, CA1replayFidelity
    - ::swrvarCat /nt @corr $swr CA1reverse, CA1forward
    - ::swrvarCatDiff /nt @corr $swr CA1reverse-forward

#3NOW
=== Noise debugging
... check referencing, noise 
    - timeXripHeatRast /nt /frq (1:400, 250:300)
    - timeXfreqPwrMean /nt

LATER
=== Network Measures
    - the designmat is now per ntrode, i.e. and has several measures per other ntrode /nt
    - each TFmapZ result is a confusion matrix of TFmaps instead of /nt it's /nt*nt
    - instead of cXYplots, organize the plots into circle
    - is mecDVCorr Corr w any exp, tf, or swr Var?
    - is mecCa1Corr Corr w any exp, tf, or swr Var?
%}

%{
LATER LATER
Part 2
Hypothesis: ripChains multiscale unlockedGate ++> ripchainDVsnowball, multiscaleCompletion
Discussion: multiscaling is conducive for or the same thing as consolidation; baking new 
            epsiodes into the caserole graph form of the neocortex)
*** ASSIGN DV
*** DETECT ripchains
*** MEASURE ripchainDVsnowball vs ripchainDVpingpong xs awakeSleep
    - 
%}















