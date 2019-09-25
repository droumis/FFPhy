#!/bin/bash
## some useful commands
#       forallrats('discretizeSpikeTimes(0.006)', 'placefields4', 'data2'); \
#       runall('discrete_4ms_xp4_t005'); \
#        runall('disc_12ms_jit_5ms_xp4_t005'); \
#       runall; \

#ssh -x aurora "LD_PRELOAD=/home/chengs/lib/libstdc++.so.6.0.9:/home/chengs/lib/libgcc_s.so.1 nohup nice -n 18 \
#    /usr/local/matlab2007b/bin/matlab -nosplash -nojvm -r  \"\
#        runall('gen_longISI_x', 'gen', 16, [0.001   3e-4    0.002], 0); \
#    quit\" \
#    </dev/null >> out.aur 2>&1  &"

#ssh -x aurora "LD_PRELOAD=/home/chengs/lib/libstdc++.so.6.0.9:/home/chengs/lib/libgcc_s.so.1 nohup nice -n 18 \
#    /usr/local/matlab2007b/bin/matlab -nosplash -nojvm -r  \"\
#        runall('gen_longISI_x', 'gen', [3 10 18 30 100], 0, 0); \
#    quit\" \
#    </dev/null >> out.aur2 2>&1  &"

#ssh -x blizzard "LD_PRELOAD=/home/chengs/lib/libstdc++.so.6.0.9:/home/chengs/lib/libgcc_s.so.1 nohup nice -n 18 \
#    /usr/local/matlab2007b/bin/matlab -nosplash -nojvm -r  \"\
#        runall('gen_longISI_x', 'gen', [0.03 0.1 0.3] , 600, 0); \
#    quit\" \
#    </dev/null >> out.bli 2>&1  &"

#ssh -x blizzard "LD_PRELOAD=/home/chengs/lib/libstdc++.so.6.0.9:/home/chengs/lib/libgcc_s.so.1 nohup nice -n 18 \
#    /usr/local/matlab2007b/bin/matlab -nosplash -nojvm -r  \"\
#        runall('gen_longISI_x', 'gen', 12, [3e-4    0.001   0.002   0.003   0.004 ] , 0); \
#    quit\" \
#    </dev/null >> out.bli2 2>&1  &"

#ssh -x cyclone "LD_PRELOAD=/home/chengs/lib/libstdc++.so.6.0.9:/home/chengs/lib/libgcc_s.so.1 nohup nice -n 18 \
#    /usr/local/matlab2007b/bin/matlab -nosplash -nojvm -r  \"\
#        runall('gen_longISI_x', 'gen', [1 2 3 4 6 8 10 14 ] , .004 , 0); \
#    quit\" \
#    </dev/null >> out.cyc  2>&1  &"

#ssh -x cyclone "LD_PRELOAD=/home/chengs/lib/libstdc++.so.6.0.9:/home/chengs/lib/libgcc_s.so.1 nohup nice -n 18 \
#    /usr/local/matlab2007b/bin/matlab -nosplash -nojvm -r  \"\
#        runall('gen_longISI_x', 'gen', [0.01 0.03 0.1 0.3 1 3 4 6] , .02 , 0); \
#    quit\" \
#    </dev/null >> out.cyc2  2>&1  &"


# ssh -x sunshine "nohup nice -n 18 \
#     /usr/notlocal/matlab2007b/bin/matlab -nosplash -nojvm -r  \"\
#        runall('gen_longISI_x', 'gen', 5, .03, 0); \
#     quit\" \
#     </dev/null >> out.sun 2>&1  &"

ssh -x tornado "nohup nice -n 18 \
    /usr/local/matlab2007b/bin/matlab -nosplash -nojvm -r  \"\
        runall('gen_longISI_x', 'gen', [14 16 18], 6e-4, 0); \
    quit\" \
    </dev/null >> out.tor  2>&1  &" 


#ssh -x tsunami "nohup nice -n 18 \
#    /usr/local/matlab2007b/bin/matlab -nosplash -nojvm -r  \"\
#    quit\" \
#    </dev/null >> out.tsu  2>&1  &" 

#ssh -x tsunami "nohup nice -n 18 \
#    /usr/local/matlab2007b/bin/matlab -nosplash -nojvm -r  \"\
#    quit\" \
#    </dev/null >> out.tsu2  2>&1  &" 


ssh -x virga "nohup nice -n 18 \
    /usr/local/matlab2007b/bin/matlab -nosplash -nojvm -r  \"\
        runall('gen_longISI_x', 'gen', [12 10 8], 6e-4, 0); \
    quit\" \
    </dev/null >> out.vir  2>&1  &" 


#        runall('nonripple_xp', 'pf9-fam', , );\
#        runall('gen_longISI_x', 'gen', , , 0); \

########## more ##########
