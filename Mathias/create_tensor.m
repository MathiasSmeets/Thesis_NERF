function synthetic_tensor = create_tensor(synthetic_data, interval_length)

synthetic_tensor = zeros(size(synthetic_data,1),interval_length);
for i = 1:size(synthetic_data,2)/interval_length
    synthetic_tensor(:,:,i) = synthetic_data(:,(i-1)*interval_length+1:i*interval_length);
end

end