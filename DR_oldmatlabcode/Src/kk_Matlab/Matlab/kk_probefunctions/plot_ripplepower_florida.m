%% plots ripple mean power and variance across sites

animalprefix='flo';
sessions=1;
epochs=14;
channels=1:28;
sitemap=[[1 3 5 NaN 9 11 NaN]' [2 4 6 8 10 12 NaN]' [15 17 19 21 23 NaN NaN]' [16 NaN 20 22 24 26 NaN]'];
sitemap=[[1 3 5 NaN 9 11 NaN]' [2 4 6 8 10 12 NaN]' [15 17 19 21 23 NaN NaN]' [16 NaN 20 22 24 26 NaN]'];

% first load filtered eeg (for ripples)

ripple=loadeegstruct('/data12/kkay/Flo/',animalprefix,'ripple',sessions,epochs,channels);

% next 

for s=sessions

    for e=epochs
        
        rms=zeros(size(sitemap)); 
        variance=zeros(size(sitemap));
    
        for c=channels
            site=find(c==sitemap);
            rms(site)=sqrt(mean(ripple{s}{e}{c}.data(:,1).^2));
            variance(site)=var(double(ripple{s}{e}{c}.data(:,1)));
        end
        
        figure

        subplot(1,2,1)
            imagesc(rms)   ;  colormap(bone)    ; colorbar ; title(e) ;
        subplot(1,2,2)     
            imagesc(variance) ;  colormap(bone)    ; colorbar ; title(e) ;

    end
    
end



% plot ripple # to check

for s=sessions
        figure
        ripcount=zeros(size(sitemap));
        for h=1:2%size(sitemap,2)                         % shank
            for z=1:length(sitemap(:,h))
                channel=sitemap(z,h);
                if channel~=0 && ~isnan(channel)
                    sum=0;
                    for e=13
                        sum=sum+length(ripples{s}{e}{channel}.startind);
                    end
                    ripcount(z,h)=sum;
                else
                    ripcount(z,h)=0;
                end
            end
        
        subplot(1,size(sitemap,2),h)
        bar(ripcount(:,h));                % z position
        ylim([0 max(max(ripcount))]);
        string=sprintf('%s shank %d',animalprefix,h);
        title(string)
        
        end
        
        figure
        imagesc(ripcount,[0 max(max(ripcount))]); colormap('bone') ; colorbar
        title(sprintf('%s ripple count, epoch %d',animalprefix,e));
        
end














