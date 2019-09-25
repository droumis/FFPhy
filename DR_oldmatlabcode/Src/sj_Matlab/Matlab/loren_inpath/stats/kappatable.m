% enter the values for the kappa table from _Circular Statistics in Biology_ 
% by Edward Batshelet, 1981, Academic Press


% set the values for n
kappan = [(5:15)' ; (20:5:50)' ; 100 ; 150 ; 200 ; 1e100];
kappar = [.1:.05:.95];

%kappa is length(n) X length(r)
kappa(1,:) =  [0    0    0    0    0    0    0    .15  .67   .94 1.18 1.41 1.68 2.0  2.44 3.1  4.39 8.33];
kappa(2,:) =  [0    0    0    0    0    0    0    .56  .83  1.04 1.25 1.48 1.74 2.07 2.51 3.20 4.54 8.66];
kappa(3,:) =  [0    0    0    0    0    0    .38  .69  .9   1.1  1.3  1.52 1.78 2.11 2.56 3.27 4.65 8.89];
kappa(4,:) =  [0    0    0    0    0    0    .53  .76  .95  1.13 1.33 1.58 1.81 2.15 2.6  3.32 4.73 9.06];
kappa(5,:) =  [0    0    0    0    0    .31  .61  .8   .98  1.16 1.35 1.57 1.84 2.17 2.63 3.36 4.79 9.19];
kappa(6,:) =  [0    0    0    0    0    .42  .65  .83  1.00 1.18 1.37 1.59 1.86 2.19 2.66 3.39 4.84 9.30];
kappa(7,:) =  [0    0    0    0    0    .48  .69  .85  1.02 1.19 1.38 1.61 1.87 2.21 2.68 3.42 4.89 9.39];
kappa(8,:) =  [0    0    0    0    .23  .53  .71  .87  1.03 1.20 1.40 1.62 1.88 2.22 2.69 3.44 4.92 9.46];
kappa(9,:) =  [0    0    0    0    .32  .56  .73  .88  1.04 1.21 1.41 1.63 1.89 2.23 2.71 3.46 4.95 9.53];  
kappa(10,:) = [0    0    0    0    .37  .58  .74  .89  1.05 1.22 1.41 1.63 1.90 2.24 2.72 3.47 4.98 9.58];  
kappa(11,:) = [0    0    0    0    .41  .60  .75  .90  1.06 1.23 1.42 1.64 1.91 2.25 2.73 3.49 5.00 9.63];  
kappa(12,:) = [0    0    0    .30  .50  .65  .79  .93  1.09 1.26 1.45 1.67 1.94 2.28 2.76 3.53 5.07 9.79];  
kappa(13,:) = [0    0    0    .38  .54  .67  .81  .95  1.10 1.27 1.46 1.68 1.95 2.30 2.79 3.56 5.12 9.88];  
kappa(14,:) = [0    0    .22  .42  .56  .69  .82  .96  1.11 1.28 1.47 1.69 1.96 2.31 2.80 3.58 5.15 9.95];  
kappa(15,:) = [0    0    .27  .44  .57  .70  .83  .97  1.12 1.29 1.48 1.70 1.97 2.32 2.81 3.60 5.17 9.99];  
kappa(16,:) = [0    0    .31  .45  .58  .70  .83  .97  1.12 1.29 1.48 1.70 1.98 2.33 2.82 3.61 5.19 10.03];  
kappa(17,:) = [0    .04  .33  .46  .58  .71  .84  .98  1.13 1.30 1.49 1.71 1.98 2.33 2.82 3.62 5.20 10.06];  
kappa(18,:) = [0    .14  .34  .47  .59  .71  .84  .98  1.13 1.30 1.49 1.71 1.98 2.34 2.83 3.62 5.21 10.08];  
kappa(19,:) = [0    .26  .38  .49  .61  .73  .86  1.00 1.15 1.31 1.50 1.73 2.00 2.35 2.85 3.65 5.26 10.18];  
kappa(20,:) = [.18  .28  .39  .50  .62  .74  .86  1.00 1.15 1.32 1.51 1.73 2.00 2.36 2.86 3.66 5.27 10.21];  
kappa(21,:) = [.19  .29  .40  .51  .62  .74  .87  1.00 1.15 1.32 1.51 1.73 2.01 2.36 2.86 3.67 5.28 10.22];  
kappa(22,:) = [.20  .30  .41  .52  .62  .75  .87  1.01 1.16 1.33 1.52 1.74 2.01 2.37 2.87 3.68 5.31 10.27];  

save /home/loren/matlab/stats/kappa kappa*
