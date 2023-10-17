function filtered_time_series = get_firing_num(current, window_size, thresh)
% current: membrane potiential
% window_size: window size for moveing summation
% thresh: threshold for firing
    firing_time_series = sum(current > thresh, 1);
    filtered_time_series = movsum(firing_time_series, window_size);
end

