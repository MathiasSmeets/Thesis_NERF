function [W, H, rank] = best_nmf(cur_total_mouse, iterations)

max_rank = min(size(cur_total_mouse))-1;
smallest_D = Inf;
for i = 1:max_rank
    for j = 1:iterations
        [W,H,D(j,i)] = nnmf(cur_total_mouse, i);
        if D(j,i) < smallest_D
            smallest_D = D(j,i);
            W_best = W;
            H_best = H;
        end
    end
end
[~, rank] = min(mean(D));
W = W_best;
H = H_best;

end