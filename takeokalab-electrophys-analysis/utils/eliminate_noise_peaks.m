function eliminate_noise_peaks(datafile)

nch = 32;
file = dir(datafile); 
samples = file.bytes/2/nch;
D1 = memmapfile(datafile, 'Format',{'int16' [nch samples] 'mapped'});
D1.Writable = true;
for i = 1:nch
    lfp = double(D1.Data(1).mapped(i,:));
    lfpfil = data_filter(lfp,'bandpass',[300 6000],10);
    idx = find(abs(lfpfil)>300);
    sx = sign(lfpfil(idx));
    uu = (ones(1,length(idx))*300).*sx;
    lfpfil(idx) = uu;
    D1.Data(1).mapped(i,:) = lfpfil;
end
%lfpnew = reshape(lfptemp,[1, size(lfptemp,1)*size(lfptemp,2)]);

%disp('Writing new data file for spike sorting...')
%fid = fopen('recording_filt_cut.dat', 'w');
%fwrite(fid, lfpfil, 'int16');
%fclose(fid);


