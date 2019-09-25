%% plots ripple mean power and variance across sites

animalprefix='mit';
sessions=13;
epochs=1:2;
channels=1:31;
sitemap=[[1 2 0 3:15]' (16:31)'];

% first load filtered eeg (for ripples)

ripple=loadeegstruct('/data13/anna/Mitt/',animalprefix,'ripple',sessions,epochs,channels);


% plot power
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
            imagesc(rms)   ;  colormap(bone)    ; colorbar ; title(animalprefix) ;
        subplot(1,2,2)     
            imagesc(variance) ;  colormap(bone)    ; colorbar ; title(animalprefix) ;
    end
end

% plot ripple # to check

for s=sessions
        figure
        ripcount=zeros(size(sitemap));
        for h=1:size(sitemap,2)                         % shank
            for z=1:length(sitemap(:,h))
                channel=sitemap(z,h);
                if channel~=0
                    sum=0;
                    for e=1
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
        title(sprintf('%s ripple count',animalprefix));
        
end















