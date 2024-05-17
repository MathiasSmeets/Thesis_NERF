function [FplotOut, psdmean_out, Pout, fout] = lfp_comparison_stim(lfp, events)

fs = 30000;
% par.fc = 50;
% par.fs = fs;
% par.type = 'notch';
% par.wd = 35;

for i = 1:length(events.onsets)
    start = single(events.onsets(i)*fs);
    stop = single(events.offsets(i)*fs);
    if i == length(events.onsets)
        start2 = length(lfp);
    else
        start2 = single(events.onsets(i+1)*fs);
    end
    lfpin{i} = lfp(start:stop);
    lfpout{i} = lfp(stop+1:start2-1);
end

sizein = [];
sizeout = [];

cont = 0;
for i = 1:length(lfpin)
    if length(lfpin{i})>8
        cont = cont+1;
        [fin{cont},Pin{cont}] = psd(lfpin{i}, fs, 1, 'welch', 0);
        sizein = [sizein size(Pin{cont},1)];
    end    
end
cont = 0;
for i = 1:length(lfpout)
    if length(lfpout{i})>8
        cont = cont+1;
        [fout{cont},Pout{cont}] = psd(lfpout{i}, fs, 1, 'welch', 0);
        sizeout = [sizeout size(Pout{cont},1)];
    end
end

[maxL, indmax] = max(sizeout);
for i = 1:length(Pout)
    if sizeout(i)<maxL
        Poutmat(i,:) = interp1(fout{i}, Pout{i}, fout{indmax});
    else
        Poutmat(i,:) = Pout{i}';
    end
end
psdmean_out = mean(Poutmat);

[maxLin, indmaxin] = max(sizein);
for i = 1:length(Pin)

    if sizein(i)<maxLin
        Pinmat(i,:) = interp1(fin{i}, Pin{i}, fin{indmaxin});
    else
        Pinmat(i,:) = Pin{i}';
    end
end
psdmean_in = mean(Pinmat);

FplotOut = fout{indmax};
figure, hold on
plot(fin{indmaxin},10*log10(psdmean_in), 'r');
plot(fout{indmax}, 10*log10(psdmean_out), 'g');
xlabel('f [Hz]')
ylabel('PSD [db/Hz]')
xlim([0 100])
