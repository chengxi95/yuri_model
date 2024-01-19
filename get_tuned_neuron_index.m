function neuron_index_list = get_tuned_neuron_index(firing_index_time, trial_num, numberofneurons, base_start, base_end, test_start, test_end)
% get index of the cells that are tuned, the cells are tuned if they fires
% much more often after stimulation (larger than mean+2*std)
%   Detailed explanation goes here
    if nargin < 7
        base_start = 100000;
        base_end = 200000;
        test_start = 200000;
        test_end = 300000;
    end
    
    base_firing_num = -1 * ones(numberofneurons, trial_num);
    test_firing_num = -1 * ones(numberofneurons, trial_num);
 
    for i = 1:trial_num
        temp = firing_index_time{i};
        if ~isnan(temp)
            for j = 1:numberofneurons
                base_firing_num(j, i) = length(find(temp(:, 1) == j & temp(:, 2) >= base_start & temp(:, 2) < base_end));
                test_firing_num(j, i) = length(find(temp(:, 1) == j & temp(:, 2) >= test_start & temp(:, 2) < test_end));
            end
        end
    end
    
    % remove no data entry
    base_firing_num = base_firing_num(base_firing_num(:, 2) ~= -1, :);
    test_firing_num = test_firing_num(test_firing_num(:, 2) ~= -1, :);

    base_mean = mean(base_firing_num, 2);
    base_std = std(base_firing_num, 1, 2);
    test_mean = mean(test_firing_num, 2);
    
    [neuron_index_list, ] = find(test_mean > (base_mean + 2*base_std));
end