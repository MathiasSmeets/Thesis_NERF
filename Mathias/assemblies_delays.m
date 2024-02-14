function assemblies = assemblies_delays(assemblies, max_delay)

new_assembly = cell(1,1);

if iscell(assemblies)
    for i = 1:length(assemblies)
        cur_assembly = unique(ceil(assemblies{i} ./ (max_delay+1)));

        % check if assembly already in new_assembly
        equal = false;
        for j = 1:length(new_assembly)
            if isequal(new_assembly{j}, cur_assembly)
                equal = true;
            end
        end
        if ~equal
            new_assembly = [new_assembly; {cur_assembly}];
        end

    end
    assemblies = new_assembly(2:end);
end

end