function switched_data = get_switched_data(data, assembly)

switched_data = data;
for i = 1:length(assembly)
    cur_index = assembly(i);
    if cur_index > length(assembly)
        switched_data([i, cur_index], :) = switched_data([cur_index, i],:);
    end
end

end