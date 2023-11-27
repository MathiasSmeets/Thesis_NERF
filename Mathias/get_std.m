peak_area_m = [6,14,4,4,2,101,86,0,2,3,5];
peak_area_m = peak_area_m - round(8.1);
peak_area_m(peak_area_m<0) = 0; 

closest_peaks_m = zeros(1,sum(peak_area_m));
cur_index = 1;
for kk = 1:size(peak_area_m,2)
    %if peaks_to_keep_m(jj) - max_peak_distance < 1
    %    latency = kk-4 - (peaks_to_keep_m(jj) - max_peak_distance-1);
    %else
        latency = kk-5;
    %end
    closest_peaks_m(cur_index:cur_index+peak_area_m(kk)-1) = latency;
    cur_index = cur_index + peak_area_m(kk);
end