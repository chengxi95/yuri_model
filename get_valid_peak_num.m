function valid_peak_num = get_valid_peak_num(baseline_firing, test_firing, window_size, minimum_cell_num)
    baseline_sum = movsum(baseline_firing, window_size);
    test_sum = movsum(test_firing, window_size);
    
    baseline_sum(baseline_sum < minimum_cell_num) = 0;
    test_sum(test_sum < minimum_cell_num) = 0;

    
    baseline_peak_num = length(findpeaks(baseline_sum));
    test_peak_num = length(findpeaks(test_sum));
    valid_peak_num = test_peak_num - baseline_peak_num;
end