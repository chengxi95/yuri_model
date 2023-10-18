function [outputArg1,outputArg2] = get_target_neuron_id(firing_mat, base_start, base_end, test_start, test_end)
    if nargin < 5
        base_start = 100000;
        base_end = 200000;
        test_start = 200000;
        test_end = 300000 - 1000;
    end

    baseline_firing_num = sum(firing_mat(:, base_start:base_end), 2);
    test_firing_num = sum(firing_mat(:, test_start:test_end), 2);
    
    firing_ratio = test_firing_num./(baseline_firing_num + 1);
    figure(5)
    
    subplot(2,1,1)
    plot(firing_ratio)

    subplot(2,1,2)
    spy(firing_mat)
end

