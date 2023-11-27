function [f,P1] = mpsd(data,fs,chs,method, plt, nfft)
% calculate power spectrum density of signals contained in data [nchs x
% timestamsp] and recorded at sampling frequency fs. Specify chs to perform only on specific channels
%
% Mattia D'Andola, Apr 2020

if nargin<5
    plt = 1;
end
if nargin<6
    nfft = fs;
end

cont = 0;
switch method
    case 'fourier'
        L = size(data,2);  
        f = fs*(0:(L/2))/L;
        for i = chs
            cont = cont+1;
            Y = fft(data(i,:));
            P2 = abs(Y/L);
            P1(cont,:) = P2(1:L/2+1);
            P1(cont,2:end-1) = 2*P1(2:end-1);
        end
    case 'welch'
        for i = chs
            cont = cont+1;
            [P1(cont,:), f] = pwelch(data(i,:),[],[],[],nfft);
        end
end

if plt
    figure
    plot(f,10*log10(P1));
    xlabel('f [Hz]')
    ylabel('PSD [db/Hz]')
    xlim([0 fs/2])
end