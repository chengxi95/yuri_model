function [base_si_list, test_si_list] = get_si_tuned_neuron_index(bt_ft, rc_ft, gc_ft, yt_ft, trial_num, numberofneurons, base_start, base_end, test_start, test_end)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    if nargin < 7
        base_start = 100000;
        base_end = 200000;
        test_start = 200000;
        test_end = 300000;
        window_size = 10000;
        full_end = 400000;
    end
    
    % placeholder for one cell's firing across all trials to calculate
    % its spike density function
    bt_firing_matrix = zeros(trial_num, full_end);
    rc_firing_matrix = zeros(trial_num, full_end);
    gc_firing_matrix = zeros(trial_num, full_end);
    yt_firing_matrix = zeros(trial_num, full_end);
    
    % placeholder for the si index for all cells
    base_si_list = zeros(numberofneurons, floor((base_end - base_start)/window_size));
    test_si_list = zeros(numberofneurons, floor((test_end - test_start)/window_size));

    for cell_index = 1:numberofneurons
        for i = 1:trial_num
            temp = bt_ft{i};
            temp_firing_matrix = zeros(numberofneurons, full_end);
            if ~isnan(temp)
                % disp(max(temp(:, 1)));
                % disp(max(temp(:, 2)));
                indexs = sub2ind(size(temp_firing_matrix), temp(:, 1), temp(:, 2));
                temp_firing_matrix(indexs) = 1;
            end
            bt_firing_matrix(i, :) = temp_firing_matrix(cell_index, :);
    
            temp = rc_ft{i};
            temp_firing_matrix = zeros(numberofneurons, full_end);
            if ~isnan(temp)
                indexs = sub2ind(size(temp_firing_matrix), temp(:, 1), temp(:, 2));
                temp_firing_matrix(indexs) = 1;
            end
            rc_firing_matrix(i, :) = temp_firing_matrix(cell_index, :);
    
            temp = gc_ft{i};
            temp_firing_matrix = zeros(numberofneurons, full_end);
            if ~isnan(temp)
                indexs = sub2ind(size(temp_firing_matrix), temp(:, 1), temp(:, 2));
                temp_firing_matrix(indexs) = 1;
            end
            gc_firing_matrix(i, :) = temp_firing_matrix(cell_index, :);
    
            temp = yt_ft{i};
            temp_firing_matrix = zeros(numberofneurons, full_end);
            if ~isnan(temp)
                indexs = sub2ind(size(temp_firing_matrix), temp(:, 1), temp(:, 2));
                temp_firing_matrix(indexs) = 1;
            end
            yt_firing_matrix(i, :) = temp_firing_matrix(cell_index, :);
    
        end
        [base_si, test_si] = get_one_si([max(get_firing_num(bt_firing_matrix, 100), get_firing_num(rc_firing_matrix, 100)); max(get_firing_num(gc_firing_matrix, 100), get_firing_num(yt_firing_matrix, 100))]);
        
        base_si_list(cell_index, :) = base_si;
        test_si_list(cell_index, :) = test_si;
    end
end
    
    
