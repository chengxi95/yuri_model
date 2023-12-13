%   Simulation time
dt = 0.01; %step size ms
t_final = 4000; %simulation time ms
end_time = 3500/dt;
total_trial_num = 200;
numberofneurons = 50;% number of neurons per group

srate=1000;
gauss_width= 100;

full_VA_num =  zeros(total_trial_num, end_time);
full_MD_num = zeros(total_trial_num, end_time);
full_PFC_num = zeros(total_trial_num, end_time);

pfc_s_bt_index = get_tuned_neuron_index(full_PFC_S_BT, total_trial_num, numberofneurons);
pfc_s_rc_index = get_tuned_neuron_index(full_PFC_S_RC, total_trial_num, numberofneurons);
pfc_s_gc_index = get_tuned_neuron_index(full_PFC_S_GC, total_trial_num, numberofneurons);
pfc_s_yt_index = get_tuned_neuron_index(full_PFC_S_YT, total_trial_num, numberofneurons);
pfc_d_bt_index = get_tuned_neuron_index(full_PFC_D_BT, total_trial_num, numberofneurons);
pfc_d_rc_index = get_tuned_neuron_index(full_PFC_D_RC, total_trial_num, numberofneurons);
pfc_d_gc_index = get_tuned_neuron_index(full_PFC_D_GC, total_trial_num, numberofneurons);
pfc_d_yt_index = get_tuned_neuron_index(full_PFC_D_YT, total_trial_num, numberofneurons);
pfc_shape_index = get_tuned_neuron_index(full_PFC_shape_ensemble, total_trial_num, numberofneurons);
pfc_ori_index = get_tuned_neuron_index(full_PFC_ori_ensemble, total_trial_num, numberofneurons);
va_shape_index = get_tuned_neuron_index(full_VA_shape, total_trial_num, numberofneurons);
va_ori_index = get_tuned_neuron_index(full_VA_ori, total_trial_num, numberofneurons);
md_shape_index = get_tuned_neuron_index(full_MD_shape, total_trial_num, numberofneurons);
md_ori_index = get_tuned_neuron_index(full_MD_ori, total_trial_num, numberofneurons);


for i = 1:total_trial_num
    
    one_pfc_s_bt = zeros(numberofneurons, end_time);
    one_pfc_s_rc = zeros(numberofneurons, end_time);
    one_pfc_s_gc = zeros(numberofneurons, end_time);
    one_pfc_s_yt = zeros(numberofneurons, end_time);
    one_pfc_d_bt = zeros(numberofneurons, end_time);
    one_pfc_d_rc = zeros(numberofneurons, end_time);
    one_pfc_d_gc = zeros(numberofneurons, end_time);
    one_pfc_d_yt = zeros(numberofneurons, end_time);
    one_pfc_shape = zeros(numberofneurons, end_time);
    one_pfc_ori = zeros(numberofneurons, end_time);
    one_va_shape = zeros(numberofneurons, end_time);
    one_va_ori = zeros(numberofneurons, end_time);
    one_md_shape = zeros(numberofneurons, end_time);
    one_md_ori = zeros(numberofneurons, end_time);

    temp = full_PFC_S_BT{i};
    if ~isnan(temp)
        % increase the neuron index for neurons in different sections of pfc
        indexs = (temp(:, 2) - 1) * numberofneurons + temp(:, 1);
        one_pfc_s_bt(indexs) = 1;
        one_pfc_s_bt = one_pfc_s_bt(pfc_s_bt_index, :);
    end
    
    temp = full_PFC_S_RC{i};
    if ~isnan(temp)
        indexs = (temp(:, 2) - 1) * numberofneurons + temp(:, 1);
        one_pfc_s_rc(indexs) = 1;
        one_pfc_s_rc = one_pfc_s_rc(pfc_s_rc_index, :);
    end

    temp = full_PFC_S_GC{i};
    if ~isnan(temp)
        temp(:, 1) = temp(:, 1) + numberofneurons;
        indexs = (temp(:, 2) - 1) * numberofneurons + temp(:, 1);
        one_pfc_s_gc(indexs) = 1;
        one_pfc_s_gc = one_pfc_s_gc(pfc_s_gc_index, :);
    end

    temp = full_PFC_S_YT{i};
    if ~isnan(temp)
        indexs = (temp(:, 2) - 1) * numberofneurons + temp(:, 1);
        one_pfc_s_yt(indexs) = 1;
        one_pfc_s_yt = one_pfc_s_yt(pfc_s_yt_index, :);
    end

    temp = full_PFC_D_BT{i};
    if ~isnan(temp)
        indexs = (temp(:, 2) - 1) * numberofneurons + temp(:, 1);
        one_pfc_d_bt(indexs) = 1;
        one_pfc_d_bt = one_pfc_d_bt(pfc_d_bt_index, :);
    end

    temp = full_PFC_D_RC{i};
    if ~isnan(temp)
        indexs = (temp(:, 2) - 1) * numberofneurons + temp(:, 1);
        one_pfc_d_rc(indexs) = 1;
        one_pfc_d_rc = one_pfc_d_rc(pfc_d_rc_index, :);
    end

    temp = full_PFC_D_GC{i};
    if ~isnan(temp)  
        indexs = (temp(:, 2) - 1) * numberofneurons + temp(:, 1);
        one_pfc_d_gc(indexs) = 1;
        one_pfc_d_gc = one_pfc_d_gc(pfc_d_gc_index, :);
    end

    temp = full_PFC_D_YT{i};
    if ~isnan(temp)       
        indexs = (temp(:, 2) - 1) * numberofneurons + temp(:, 1);
        one_pfc_d_yt(indexs) = 1;
        one_pfc_d_yt = one_pfc_d_yt(pfc_d_yt_index, :);
    end 

    temp = full_PFC_shape_ensemble{i};
    if ~isnan(temp)      
        indexs = (temp(:, 2) - 1) * numberofneurons + temp(:, 1);
        one_pfc_shape(indexs) = 1;
        one_pfc_shape = one_pfc_shape(pfc_shape_index, :);
    end
    
    temp = full_PFC_ori_ensemble{i};
    if ~isnan(temp)    
        indexs = (temp(:, 2) - 1) * numberofneurons + temp(:, 1);
        one_pfc_ori(indexs) = 1;
        one_pfc_ori = one_pfc_ori(pfc_ori_index, :);
    end   

    temp = full_VA_shape{i};
    if ~isnan(temp) 
        indexs = (temp(:, 2) - 1) * numberofneurons + temp(:, 1);
        one_va_shape(indexs) = 1;
        one_va_shape = one_va_shape(va_shape_index, :);
    end  

    temp = full_VA_ori{i};
    if ~isnan(temp)    
        indexs = (temp(:, 2) - 1) * numberofneurons + temp(:, 1);
        one_va_ori(indexs) = 1;
        one_va_ori = one_va_ori(va_ori_index, :);
    end

    temp = full_MD_shape{i};
    if ~isnan(temp)        
        indexs = (temp(:, 2) - 1) * numberofneurons + temp(:, 1);
        one_md_shape(indexs) = 1;
        one_md_shpae = one_md_shape(md_shape_index, :);
    end   

    temp = full_MD_ori{i};
    if ~isnan(temp)       
        indexs = (temp(:, 2) - 1) * numberofneurons + temp(:, 1);
        one_md_ori(indexs) = 1;
        one_md_ori = one_md_ori(md_ori_index, :);
    end
    
    full_VA_num(i, :) =  get_firing_num(cat(1, one_va_shape, one_va_ori), gauss_width);
    full_MD_num(i, :) = get_firing_num(cat(1, one_md_shape, one_md_ori), gauss_width);
    full_PFC_num(i, :) = get_firing_num(cat(1, one_pfc_s_bt, one_pfc_s_rc, one_pfc_s_gc, one_pfc_s_yt, one_pfc_d_bt, one_pfc_d_rc, one_pfc_d_gc, one_pfc_d_yt, one_pfc_shape, one_pfc_ori), gauss_width);
end

average_VA_num = sum(full_VA_num, 1);
VA_firing_rate_timeseries=  conv_gaussian(average_VA_num, srate,gauss_width);

average_MD_num = sum(full_MD_num, 1);
MD_firing_rate_timeseries=  conv_gaussian(average_MD_num, srate,gauss_width);

average_PFC_num = sum(full_PFC_num, 1);
PFC_firing_rate_timeseries=  conv_gaussian(average_PFC_num, srate,gauss_width);


figure(2)
plot(VA_firing_rate_timeseries, 'b'), xlim([0, end_time])
hold on
plot(MD_firing_rate_timeseries, 'g'), xlim([0, end_time])
hold on
plot(PFC_firing_rate_timeseries, 'r'), xlim([0, end_time])
legend('VA', 'MD', 'PFC')