function filtered_time_series = get_firing_num(firing_mat, window_size)
% current: membrane potiential
% window_size: window size for moveing summation
% thresh: threshold for firing
    firing_time_series = sum(firing_mat, 1);
    filtered_time_series = movsum(firing_time_series, window_size);
end

