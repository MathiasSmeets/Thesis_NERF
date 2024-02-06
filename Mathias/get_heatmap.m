function get_heatmap(data)

% normalize data
for i = 1:size(data,1)
    data(i,:) = data(i,:) / max(data(i,:));
end

heatmap(data)

end