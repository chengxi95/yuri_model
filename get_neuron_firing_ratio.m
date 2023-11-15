function firing_ratio = get_neuron_firing_ratio(firing_list, trial_num, neuron_num, base_start, base_end, test_start, test_end)
    if nargin < 7
        neuron_num = 50;
        base_start = 100000;
        base_end = 200000;
        test_start = 200000;
        test_end = 300000;
    end
    
    firing_ratio = zeros(neuron_num, 1);
    for i = 1:trial_num
        if ~isnan(firing_list{i})
            for j = 1:neuron_num
                base_firing_num = length(find((firing_list{i}(:, 1)==j) & (firing_list{i}(:, 2) > base_start) & (firing_list{i}(:, 2) < base_end)));
                test_firing_num = length(find((firing_list{i}(:, 1)==j) & (firing_list{i}(:, 2) > test_start) & (firing_list{i}(:, 2) < test_end)));
                firing_ratio(j) = firing_ratio(j) + test_firing_num/(base_firing_num+1);
            end
        end
    end

    firing_ratio = firing_ratio ./ trial_num;

