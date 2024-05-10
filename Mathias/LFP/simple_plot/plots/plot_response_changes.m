function plot_response_changes(psth_M, neurons_M, psth_Y, neurons_Y, bins, center, norm)
frf_M = [];
pk_M = [];
conti = 0;
    
if norm  
 %normalize on max
    for i = 1:length(psth_M)
        for j = neurons_M
            conti = conti+1;
            psth_M{i}(j,:) = psth_M{i}(j,:)./max(psth_M{i}(j,:));
        end
    end
    for i = 1:length(psth_Y)
        for j = neurons_Y
            conti = conti+1;
            psth_Y{i}(j,:) = psth_Y{i}(j,:)./max(psth_Y{i}(j,:));
            frb = mean(psth_Y{i}(j,bins<=-0.002));
        end
    end
end
if center 
%center data
    for j = neurons_M
        meancenter = [];
        for i = 1:length(psth_M)
            meancenter = [meancenter psth_M{i}(j,:)];
        end
        meancenter = mean(meancenter);
        for i = 1:length(psth_M)
            psth_M{i}(j,:) = psth_M{i}(j,:)-meancenter;
        end
    end
    for j = neurons_Y
        meancenter = [];
        for i = 1:length(psth_Y)
            meancenter = [meancenter psth_Y{i}(j,:)];
        end
        meancenter = mean(meancenter);
        for i = 1:length(psth_Y)
            psth_Y{i}(j,:) = psth_Y{i}(j,:)-meancenter;
        end
    end
end

        
for i = 1:length(psth_M)
    for j = neurons_M
        
        frb = mean(psth_M{i}(j,bins<=-0.002));
        
        conti = conti+1;
        %find offset of response
        idx = find(bins>=0.002);
        idxoff = find(psth_M{i}(j,idx)>=mean(psth_M{i}(j,bins<=-0.002))); 

        fra = mean(psth_M{i}(j,idx(idxoff)));
        if ~isempty(idxoff)
            [pks, locs] = findpeaks(psth_M{i}(j,idx(idxoff)));
            [~,imax] = max(pks);

            if ~isempty(imax)
                pk_M(i,conti) = bins(idx(idxoff(locs(imax))));
                off_M(i,conti) = bins(idx(idxoff(end)));
                if fra > frb
                    frf_M(i,conti) = fra/frb;
                else
                    frf_M(i,conti) = frb/fra;
                end
            else
                pk_M(i,conti) = NaN;
                off_M(i,conti) = NaN;
                frf_M(i,conti) = NaN;
            end  
        else
           pk_M(i,conti) = NaN;
           off_M(i,conti) = NaN;
           frf_M(i,conti) = NaN; 
        end
    end
end

frf_Y = [];
pk_Y = [];
conti = 0;
for i = 1:length(psth_Y)
    for j = neurons_Y
        conti = conti+1;

        frb = mean(psth_Y{i}(j,bins<=-0.002));
        %find offset of response
        idx = find(bins>=0.002);
        idxoff = find(psth_Y{i}(j,idx)>=mean(psth_Y{i}(j,bins<=-0.002))); 

        fra = mean(psth_Y{i}(j,idx(idxoff)));
        if ~isempty(idxoff)
            [pks, locs] = findpeaks(psth_Y{i}(j,idx(idxoff)));
            [~,imax] = max(pks);
            if ~isempty(imax)
                pk_Y(i,conti) = bins(idx(idxoff(locs(imax))));
                off_Y(i,conti) = bins(idx(idxoff(end)));

                if fra > frb
                    frf_Y(i,conti) = fra/frb;
                else
                    frf_Y(i,conti) = frb/fra;
                end
            else
                pk_Y(i,conti) = NaN;
                off_Y(i,conti) = NaN;
                frf_Y(i,conti) = NaN;
            end   
        else
            pk_Y(i,conti) = NaN;
            off_Y(i,conti) = NaN;
            frf_Y(i,conti) = NaN;
        end
    end
end

frf_M(isinf(frf_M)) = 1;
frf_Y(isinf(frf_Y)) = 1;

figure
subplot(3,2,1)
errorbar(1:length(psth_M), nanmean(frf_M,2), std(frf_M,0,2,'omitnan')/sqrt(length(neurons_M)),'ko-')
xlim([0.5 length(psth_M)+0.5])
set(gca, 'Xtick', 1:length(psth_M), 'XtickLabel', {'adaptation','learning'})
xlabel('Condition')
ylabel('Relative FR')
title('Masters')

subplot(3,2,2)
errorbar(1:length(psth_Y), nanmean(frf_Y,2), std(frf_Y,0,2,'omitnan')/sqrt(length(neurons_Y)),'ko-')
xlim([0.5 length(psth_Y)+0.5])
set(gca, 'Xtick', 1:length(psth_Y),'XtickLabel', {'adaptation','learning'})
xlabel('Condition')
ylabel('Relative FR')
title('Yokes')


subplot(3,2,3)
errorbar(1:length(psth_M), nanmean(pk_M,2), std(pk_M,0,2,'omitnan')/sqrt(length(neurons_M)),'ko-')
xlim([0.5 length(psth_M)+0.5])
set(gca, 'Xtick', 1:length(psth_M),'XtickLabel', {'adaptation','learning'})
xlabel('Condition')
ylabel('Peak amp')
title('Masters')

subplot(3,2,4)
errorbar(1:length(psth_Y), nanmean(pk_Y,2), std(pk_Y,0,2,'omitnan')/sqrt(length(neurons_Y)),'ko-')
xlim([0.5 length(psth_Y)+0.5])
set(gca, 'Xtick', 1:length(psth_Y),'XtickLabel', {'adaptation','learning'})
xlabel('Condition')
ylabel('peak amp')
title('Yokes')


subplot(3,2,5)
errorbar(1:length(psth_M), nanmean(off_M,2), std(off_M,0,2,'omitnan')/sqrt(length(neurons_M)),'ko-')
xlim([0.5 length(psth_M)+0.5])
set(gca, 'Xtick', 1:length(psth_M),'XtickLabel', {'adaptation','learning'})
xlabel('Condition')
ylabel('resp length')
title('Masters')

subplot(3,2,6)
errorbar(1:length(psth_Y), nanmean(off_Y,2), std(off_Y,0,2,'omitnan')/sqrt(length(neurons_Y)),'ko-')
xlim([0.5 length(psth_Y)+0.5])
set(gca, 'Xtick', 1:length(psth_Y),'XtickLabel', {'adaptation','learning'})
xlabel('Condition')
ylabel('resp length')
title('Yokes')
end


