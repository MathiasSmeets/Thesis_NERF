function N = normrobust(A)
%NORMROBUST robustly normalise A over the whole matrix instead of 1 column
%   A is a 2 dimensional matrix
    temp = reshape(A,1,[]);
    N = (A - median(temp,'omitnan')) ./ median(abs(temp - median(temp,'omitnan')),'omitnan');
end

