
% This script saves MI data files (day-epoch-chan/tet) in a 'MIDATA' folder in the animal's
% animdirectory. 

% Specify each state in statespec. -- if the file already exists, the analysis function will append the
% latest calculation for that tetrode to the file.


setparameters = 1;
runscript = 0;
    if runscript
        
        selection_phasetet = 0;  %  0: Phase tet is same tetrode
                                 %  1, 2, 3, 4 : other options refer to
                                 %  STAreftet (cc, CA2, CA3, DG)
    end
plot_adet = 1;
if plot_adet
    manual_adet = [14 7 4 19];  [1 8 6 24]; [1 9 3 13];
    animal_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander',...
                        'Iceland','Justice','Kapital','Laplace','Eight','Ten','Conley','Miles'};
    MI_scale = [0 .002];
    xminmax = [1 12];
    yminmax = [20 250];
    maxphasefreq = 12;
    states_toplot = [1 2 3 ]; %:2 %:2;
    lastentry = 1;   % plots just the latest data calculated
end

if runscript
        %target_det = [ 11 * ones(14,1)   7 * ones(14,1)    15 ; ... 
        %               8 * ones(14,1)    7 * ones(14,1)    15   ]; %Kapital  -- sleep
       animals_torun = {'Frank','Bond','Corriander','Ten','Conley','Miles'};
       epochfilter =  epochmaker('runlinear_all') ;
       
       % tetinfo 
       target_det =   [];
       stareftet_only = 0;
       tetfilter = '(isequal($area, ''CA3'') && ($numcells >= 2))';
       
                        %[9 5 13 ; 9 5 18];
       
                        %[ 9 6 9 ; 9 6 8 ; 9 6 2 ; 9 6 14 ; 9 6 11 ;9 6 19 ; 9 6 12; 9 6 13; 9 6 18 ; ...
                        %9 7 9 ; 9 7 8 ; 9 7 2 ; 9 7 14 ; 9 7 11 ;9 7 19 ; 9 7 12; 9 7 13; 9 7 18 ; ...
                        %11 9 9 ; 11 9 8 ; 11 9 2 ; 11 9 14 ; 11 9 11 ; 11 9 19 ; 11 9 12; 11 9 13; 11 9 18 ; ...
                        %];  %Government  -- SLEEP 
                    
                        %[ 9 3 9 ; 9 3 8 ; 9 3 2 ; 9 3 14 ; 9 3 11 ;9 3 19 ; 9 3 12; 9 3 13; 9 3 18 ; ...
                        %9 2 9 ; 9 2 8 ; 9 2 2 ; 9 2 14 ; 9 2 11 ;9 2 19 ; 9 2 12; 9 2 13; 9 2 18 ; ...
                        %];  %Government  -- RUN     
                        %[ 9 5 9 ; 9 5 8 ; 9 5 20 ; 9 5 2 ; 9 5 14 ; 9 5 11 ;9 5 19 ; 9 5 12];      %Government  -- sleep     
       
       
                        %[ 15    7    14 ; ... 
                        %15    7    16 ; ... 
                        %15    8    14 ;
                        %15    9    14 ]; %Laplace  -- sleep     
       
       
                      %[ 8    7    15 ; ... 
                      % 11    7    15   ]; %Kapital  -- sleep     
        
        
                       %11 * ones(18,4) 2 * ones(18,4)  (4:18)' ; ... 
                       %11 * ones(18,4) 4 * ones(18,4)  (4:18)' ; ... 
                       %7 * ones(18,4) 2 * ones(18,4)  (4:18)' ; ... 
                       %7 * ones(18,4) 6 * ones(18,4)  (4:18)'  ]; %Kapital
        
                      % Justice  
                      %[ 11 * ones(18,1)   2 * ones(18,1)  (1:18)' ; ...   
                      % 11 * ones(18,1)   8 * ones(18,1)  (1:18)' ; ...
                      % 10 * ones(18,1)   2 * ones(18,1)  (1:18)' ; ...
                      % 10 * ones(18,1)   8 * ones(18,1)  (1:18)' ; ]; 
            %  [ 15 * ones(18,1)   2 * ones(18,1)  (1:18)' ]; Laplace
            %  [ 15 * ones(18,1)   9 * ones(18,1)  (1:18)' ];   % Laplace
            %  [ 8 * ones(18,1) 2 * ones(18,1) (1:18)' ]; Kapital
        
        %tetfilter = '' ;  %'isequal($area, ''CA1'')';  'isequal($area, ''DG'')';
        %tetfilter = '(isequal($area, ''DG'') && ($numcells >= 1))';
        
    end





% pre-processing parameters
if setparameters
% states
    timefilterscript
    %statespec_descript = 'state1: vel4 // state2: vel4, nogf_ca3dg, nogl_ca3dg // state3: immobile, norip ';
    clear statespec
    %statespec_descript = 'state1: vel4 // state2: vel20 // state3: immobile3, norip // state4: gf2_ca3dg ';
    %statespec_descript = 'state1: vel4 // state2: vel20 // state3: immobile3, norip // state4: wides2_ca3b // state5: gf2_ca3b // state6: gl2_ca3b ';
    
    if 1
        statespec_descript = 'state1: vel20 // state2: immobile2   ';
        statespec{1} = { vel20 } ;
        statespec{2} = { vel4, velunder20, norip }; 
        statespec{3} = { immobile2 } ;                                   
    elseif 1
        statespec_descript = 'state1: sleepc NREM // state2: sleepc, REM  ';
        statespec{1} = { sleepc, noREM } ;
        statespec{2} = { sleepc, REM } ;
    end
 
% version 1: original
    %PhaseFreqVector=0:2:50;
    %AmpFreqVector=10:5:300;
    %PhaseFreq_BandWidth=4;
    %AmpFreq_BandWidth=10;
% version 2: low freq combination
    %PhaseFreqVector=0:.5:12;
    %AmpFreqVector=10:5:250;
    %PhaseFreq_BandWidth=2;
    %AmpFreq_BandWidth=5;
% version 3: expansive full-rangeope
    PhaseFreqVector=0:0.5:50;
    AmpFreqVector=10:5:250;
    PhaseFreq_BandWidth=2;
    AmpFreq_BandWidth=5;
end


if runscript
    
    iterator = 'eeganal';
    
    for aa = 1:length(animals_torun)
        
        animal = animals_torun{aa};
        
        timefilter = {};  % done via statespec, within analysis function
        
        f = createfilter('animal',animal,'epochs',epochfilter,'excludetime',timefilter,'iterator', iterator,'eegtetrodes',tetfilter);
        f = setfilterfunction(f,'dfakk_MI',{'eeg'},target_det,animal,statespec,statespec_descript,...
                                'PhaseFreqVector',PhaseFreqVector,'AmpFreqVector',AmpFreqVector,...
                                'PhaseFreq_BandWidth',PhaseFreq_BandWidth,'AmpFreq_BandWidth',AmpFreq_BandWidth,...
                                'selection_phasetet',selection_phasetet,'stareftet_only',stareftet_only);
        f = runfilter(f);
        
    end
    
end


%% Plot results by tetrode across all days and epochs.

% iterate across tetrodes

% first collect all epoch indices [d e t]

if plot_adet
    
    for tt = 1:size(manual_adet,1)
        
        an = manual_adet(tt,1);
        day = manual_adet(tt,2);
        ep = manual_adet(tt,3);
        tet = manual_adet(tt,4);
        
        % load tetinfo
        animalname = animal_order{an};
            animalinfo = animaldef(animalname);
            tetinfo = loaddatastruct(animalinfo{2},animalinfo{3},'tetinfo');
            if isfield(tetinfo{day}{ep}{tet},'area')
                area = tetinfo{day}{ep}{tet}.area;
                if isfield( tetinfo{day}{ep}{tet},'subarea')
                    subarea = tetinfo{day}{ep}{tet}.subarea;
                end
            else
                area = 'nolabel';
                subarea = 'nolabel';
            end
                
                %areascript
                %regionscript

        % load MI data
        cd([animalinfo{2} 'MIDATA_new/'])
            filename = dir(sprintf('*MI-%d-%d-%d.mat',day,ep,tet));
            load([filename.name],'out')
                    
        numentries = length(out);
        if lastentry
             entry_toplot = numentries;
        else
            entry_toplot = 1:numentries;
        end
        for ee = entry_toplot  % plot 
            
            K = figure('units','normalized','outerposition',[0 0 .15 .9]);
            
            entry = out(ee);
            
            statespec_descript = [entry.statespec_descript ];
            durationstring = ['(' num2str(round(entry.duration)) ')'];            
            numstates = length(entry.statespec) ;
            
            for state = states_toplot
                
                subplot(3,1,state)
                freqs = entry.frequency_phase + entry.phasefreq_bandwidth/2;
                amplitudes = entry.frequency_amplitude + entry.ampfreq_bandwidth/2;
                
                pind = lookup(maxphasefreq,entry.frequency_phase);
                
                contourf(freqs(1:pind),...
                    amplitudes, ...
                    entry.COMOD{state}(1:pind,:)',30,'lines','none')
                %colorbar
                caxis(MI_scale);
                set(gca,'YDir','normal');
                colormap('jet')
                %colormap(flipud(colormap))
                % title
                %titlestring = statestring(entry.statespec{state});
                title(sprintf('state %d (%d)',state, round(entry.duration(state))),'fontweight','bold','fontsize',18);
                 set(gca,'fontsize',16,'fontweight','bold')
                 
                 xlim(xminmax);
                 ylim(yminmax)
            end
            
                % supertitle
        [~,title_handle] = suplabel([num2str(entry.index) ' (' area ') --- ' statespec_descript],'t');
        set(title_handle,'Fontsize',12,'FontWeight','bold','Color','k')
            
        end

    end
end
    
    
    
if 0
            flag = 0;
            figure
            subplot_no = 1;
            
            % collect the tetrode's epochs
            tetepochs = indices(indices(:,3)==tet,:);
            numepochs = size(tetepochs,1);
            
            for ep = 1:numepochs
                if subplot_no > subplot_max
                    subplot_no = 1;
                    figure
                end
                subplot(1,1,subplot_no)
                ind = rowfind(tetepochs(ep,:),indices);
                contourf(f{state}.output{1}(ind).frequency_phase +f{state}.output{1}(end).phasefreq_bandwidth/2,...
                    f{state}.output{1}(ind).frequency_amplitude+f{state}.output{1}(end).ampfreq_bandwidth/2, ...
                    f{state}.output{1}(ind).comodulogram',30,'lines','none')
                colorbar
                caxis(MI_scale);
                set(gca,'YDir','normal');
                colormap('gray')
                title([f{state}.animal{1} ' ' num2str(tetepochs(ep,1)) '  ' num2str(tetepochs(ep,2)) ' ' num2str(tetepochs(ep,3))])
                subplot_no = subplot_no + 1;
            end
            
end





%cd('/opt/data13/kkay/__WG/MI')







