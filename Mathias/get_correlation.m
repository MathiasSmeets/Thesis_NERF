function correlation = get_correlation(cc, total_candidates)
% calculate strength correlation based on shao and tsau 1996

% y1j = value in row 1, column j
% y0j = value in row 2, column j
% r1 = sum all y1j
% r0 = sum all y0j
% nj = sum y1j and y0j
% N = sum all n

% in our case, n1 = n2 = ... nj = n, so it will be denoted as n

intermediate_sum = 0;
for i = 1:size(cc,2)
    intermediate_sum = intermediate_sum + cc(1,i)^2 + cc(2,i)^2;
end

n = cc(1,1) + cc(2,1);
r1 = sum(cc(1,:));
r0 = sum(cc(2,:));
N = r1+r0;

correlation = ((1/n) * intermediate_sum - ((r1^2+r0^2) / N)) / (N -((r1^2+r0^2)/N));

%correlation = sqrt(correlation * (n/total_candidates));
end