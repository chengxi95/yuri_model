function [base_si, test_si] = get_one_si(spike_density_list)
    % disp(size(spike_density_list));
    window_size = 10000;
    base_start = 100000;
    base_end = 200000;
    test_start = 200000;
    test_end = 300000;
    
    base_window_num = floor((base_end - base_start)/window_size);
    test_window_num = floor((test_end - test_start)/window_size);
    base_si = zeros(base_window_num, 1);
    test_si = zeros(test_window_num, 1);
    
    n = height(spike_density_list);
    lambda_list = zeros(n, 1);

    for i = 1:base_window_num
        for j = 1:n
            lambda_list(j) = sum(spike_density_list(j, base_start + (i-1)*window_size: base_start+i*window_size));
        end
        si = (n - sum(lambda_list)/max(lambda_list))/(n-1);
        base_si(i) = si;
    end
    
    lambda_list = zeros(n, 1);
    for i = 1:test_window_num
        for j = 1:n
            lambda_list(j) = sum(spike_density_list(j, test_start + (i-1)*window_size: test_start+i*window_size));
        end
        si = (n - sum(lambda_list)/max(lambda_list))/(n-1);
        test_si(i) = si;
    end
end