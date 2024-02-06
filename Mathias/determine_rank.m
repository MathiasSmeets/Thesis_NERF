function rank = determine_rank(cur_total_mouse)

max_rank = min(size(cur_total_mouse))-1;
for i = 1:max_rank
    for j = 1:100
        [~,~,D(j,i)] = nnmf(cur_total_mouse, i);
    end
end
[~, rank] = min(mean(D));

end