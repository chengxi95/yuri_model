function neuron_index_list = get_target_neuron_idex(firing_mat, base_start, base_end, test_start, test_end)
    if nargin < 5
        base_start = 100000;
        base_end = 200000;
        test_start = 200000;
        test_end = 300000 - 1000;
    end

    baseline_firing_num = sum(firing_mat(:, base_start:base_end), 2);
    test_firing_num = sum(firing_mat(:, test_start:test_end), 2);
    
    firing_ratio = test_firing_num./(baseline_firing_num + 1);

    % find the bin with most number of neurons
    [bin_count, bin_edges] = histcounts(firing_ratio);
    [max_count, max_index] = max(bin_count);

    bin_start = bin_edges(max_index);
    bin_end = bin_edges(max_index+1);

    % get the index for those neurons
    neuron_index_list = find(firing_ratio >= bin_start & firing_ratio <= bin_end);
end

