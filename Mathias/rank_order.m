% Statistical Testing: For each event, the rank-order correlation between cell number and time was calculated, and a probability was determined. The distribution of correlation values of all events was compared to a distribution of shuffled correlation values using statistical tests (e.g., two-sample Kolmogorov–Smirnov test) to assess significance
%Certainly! Here is an example code snippet in MATLAB that demonstrates how you can calculate the rank-order correlation between cell number and time for a given event and determine its significance:

% Generate example data for cell numbers and corresponding times
cell_numbers = [1, 2, 3, 4, 5]; % Example cell numbers
times = [0.1, 0.2, 0.3, 0.4, 0.5]; % Example corresponding times

% Calculate rank-order correlation
correlation = corr(cell_numbers', times', 'Type', 'Spearman');

% Determine significance using a two-tailed test
alpha = 0.05; % Significance level
n = length(cell_numbers); % Number of data points
critical_value = tinv(1 - alpha/2, n-2); % Critical value for two-tailed test
t_stat = correlation / sqrt((1 - correlation^2) / (n - 2)); % Calculate t-statistic = measure for significanse of difference or strength of the relationship
% if this is high, the observed correlation is significantly different from zero
p_value = 2 * (1 - tcdf(abs(t_stat), n-2)); % Calculate two-tailed p-value, low value (below 0.05 suggests strong evidence to reject the null hypothesis)

% Display results
disp(['Rank-order correlation: ', num2str(correlation)]);
disp(['Two-tailed p-value: ', num2str(p_value)]);
if p_value < alpha
    disp('The event is statistically significant.');
else
    disp('The event is not statistically significant.');
end

% Display results with detailed explanations
disp('--- Rank-Order Correlation Analysis ---');
disp(['Rank-order correlation: ', num2str(correlation)]);
disp(['Degrees of freedom: ', num2str(n - 2)]);
disp(['Critical value at alpha = ', num2str(alpha), ': ', num2str(critical_value)]);
disp(['T-statistic: ', num2str(t_stat)]);
disp(['Two-tailed p-value: ', num2str(p_value)]);

% Interpretation based on p-value and critical value
% if abs(t_stat) > critical_value
%     disp('The event is statistically significant at the 5% level.');
%     disp('There is evidence of a significant relationship between cell numbers and times.');
% else
%     disp('The event is not statistically significant at the 5% level.');
%     disp('There is no strong evidence of a relationship between cell numbers and times.');
% end

%You are correct. To compare the calculated correlation value with a distribution of shuffled correlation values, you would need to generate shuffled data and perform the statistical test. Here is an updated version of the MATLAB code snippet that includes the generation of shuffled correlation values and the statistical test:

% Generate example shuffled data for cell numbers and corresponding times
shuffled_cell_numbers = cell_numbers(randperm(length(cell_numbers)));
shuffled_times = times(randperm(length(times)));

% Calculate rank-order correlation for shuffled data
shuffled_correlation = corr(shuffled_cell_numbers', shuffled_times', 'Type', 'Spearman');

% Determine significance using a two-tailed test
shuffled_correlations = zeros(1, 1000); % Number of shuffled correlations to generate
for i = 1:1000
    shuffled_cell_numbers = cell_numbers(randperm(length(cell_numbers)));
    shuffled_times = times(randperm(length(times)));
    shuffled_correlations(i) = corr(shuffled_cell_numbers', shuffled_times', 'Type', 'Spearman');
end

% Calculate p-value based on shuffled correlations
p_value_shuffled = sum(abs(shuffled_correlations) >= abs(correlation)) / length(shuffled_correlations);

% Display results
disp(['Rank-order correlation: ', num2str(correlation)]);
disp(['Shuffled correlation p-value: ', num2str(p_value_shuffled)]);
if p_value_shuffled < alpha
    disp('The event is statistically significant based on shuffled data.');
else
    disp('The event is not statistically significant based on shuffled data.');
end
%his updated code snippet generates shuffled correlation values, calculates the p-value based on the shuffled correlations, and compares the original correlation value with the distribution of shuffled correlation values to determine significance.