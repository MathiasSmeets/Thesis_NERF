function index = get_index(connection, total_nb_neurons)
% get index in matrix where the connection is

previous = 0;
if connection(1) > 1
    for i = 1:connection(1)-1
        previous = previous + total_nb_neurons - i;
    end
end

index = previous + (connection(2) - connection(1));

end