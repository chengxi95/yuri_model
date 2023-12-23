function success_ratio = get_success_ratio(full_PFC, minimum_cell_num, minimum_peak_count, window_size, trial_num)
    if nargin < 5
        numberofneurons = 50;
        base_start = 100000;
        base_end = 200000;
        test_start = 200000;
        test_end = 300000;
        full_end = 400000;
        window_size = 100;
        minimum_cell_num = 20;
        minimum_peak_count = 4;
        trial_num = 1;
    end
    success_count = 0;
    total_count = 0;

    for i = 1:trial_num
        if ~isnan(full_PFC{i})
            firing_matrix = zeros(numberofneurons, full_end);
            temp = full_PFC{i};
 
            indexs = sub2ind(firing_matrix, temp(:, 1), temp(:, 2));
            firing_matrix(indexs) = firing_matrix(indexs) + 1;
            full_firing_num = sum(firing_matrix, 1);
            valid_peak_count = get_valid_peak_num(full_firing_num(base_start:base_end), full_firing_num(test_start:test_end), window_size, minimum_cell_num);
            
            if valid_peak_count >= minimum_peak_count
                success_count = success_count + 1;
            end
            total_count = total_count + 1;
        end
    end
    success_ratio = success_count/total_count;
end