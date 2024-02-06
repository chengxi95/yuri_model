% bt_file = 'spike_time_healthy_BT_stim.mat';
% gc_file = 'spike_time_healthy_GC_stim.mat';
% rc_file = 'spike_time_healthy_RC_stim.mat';
% yt_file = 'spike_time_healthy_YT_stim.mat';

% bt_file = 'spike_time_va_off_BT_stim.mat';
% gc_file = 'spike_time_va_off_GC_stim.mat';
% rc_file = 'spike_time_va_off_RC_stim.mat';
% yt_file = 'spike_time_va_off_YT_stim.mat';

bt_file = 'spike_time_md_off_BT_stim.mat';
gc_file = 'spike_time_md_off_GC_stim.mat';
rc_file = 'spike_time_md_off_RC_stim.mat';
yt_file = 'spike_time_md_off_YT_stim.mat';

bt = load(bt_file);
rc = load(rc_file);
gc = load(gc_file);
yt = load(yt_file);

[PFC_shape_base, PFC_shape_test] = get_si_tuned_neuron_index(bt.full_PFC_shape_ensemble, rc.full_PFC_shape_ensemble, gc.full_PFC_shape_ensemble, yt.full_PFC_shape_ensemble, 500, 50);
PFC_shape_si = cat(2, PFC_shape_base, PFC_shape_test);

[PFC_ori_base, PFC_ori_test] = get_si_tuned_neuron_index(bt.full_PFC_ori_ensemble, rc.full_PFC_ori_ensemble, gc.full_PFC_ori_ensemble, yt.full_PFC_ori_ensemble, 500, 50);
PFC_ori_si = cat(2, PFC_ori_base, PFC_ori_test);

[VA_shape_base, VA_shape_test] = get_si_tuned_neuron_index(bt.full_VA_shape, rc.full_VA_shape, gc.full_VA_shape, yt.full_VA_shape, 500, 50);
VA_shape_si = cat(2, VA_shape_base, VA_shape_test);

[VA_ori_base, VA_ori_test] = get_si_tuned_neuron_index(bt.full_VA_ori, rc.full_VA_ori, gc.full_VA_ori, yt.full_VA_ori, 500, 50);
VA_ori_si = cat(2, VA_ori_base, VA_ori_test);

[MD_shape_base, MD_shape_test] = get_si_tuned_neuron_index(bt.full_MD_shape, rc.full_MD_shape, gc.full_MD_shape, yt.full_MD_shape, 500, 50);
MD_shape_si = cat(2, MD_shape_base, MD_shape_test);

[MD_ori_base, MD_ori_test] = get_si_tuned_neuron_index(bt.full_MD_ori, rc.full_MD_ori, gc.full_MD_ori, yt.full_MD_ori, 500, 50);
MD_ori_si = cat(2, MD_ori_base, MD_ori_test);

[PFC_S_BT_base, PFC_S_BT_test] = get_si_tuned_neuron_index(bt.full_PFC_S_BT, rc.full_PFC_S_BT, gc.full_PFC_S_BT, yt.full_PFC_S_BT, 500, 50);
PFC_S_BT_si = cat(2, PFC_S_BT_base, PFC_S_BT_test);

[PFC_S_RC_base, PFC_S_RC_test] = get_si_tuned_neuron_index(bt.full_PFC_S_RC, rc.full_PFC_S_RC, gc.full_PFC_S_RC, yt.full_PFC_S_RC, 500, 50);
PFC_S_RC_si = cat(2, PFC_S_RC_base, PFC_S_RC_test);

[PFC_S_GC_base, PFC_S_GC_test] = get_si_tuned_neuron_index(bt.full_PFC_S_GC, rc.full_PFC_S_GC, gc.full_PFC_S_GC, yt.full_PFC_S_GC, 500, 50);
PFC_S_GC_si = cat(2, PFC_S_GC_base, PFC_S_GC_test);

[PFC_S_YT_base, PFC_S_YT_test] = get_si_tuned_neuron_index(bt.full_PFC_S_YT, rc.full_PFC_S_YT, gc.full_PFC_S_YT, yt.full_PFC_S_YT, 500, 50);
PFC_S_YT_si = cat(2, PFC_S_YT_base, PFC_S_YT_test);

[PFC_D_BT_base, PFC_D_BT_test] = get_si_tuned_neuron_index(bt.full_PFC_D_BT, rc.full_PFC_D_BT, gc.full_PFC_D_BT, yt.full_PFC_D_BT, 500, 50);
PFC_D_BT_si = cat(2, PFC_D_BT_base, PFC_D_BT_test);

[PFC_D_RC_base, PFC_D_RC_test] = get_si_tuned_neuron_index(bt.full_PFC_D_RC, rc.full_PFC_D_RC, gc.full_PFC_D_RC, yt.full_PFC_D_RC, 500, 50);
PFC_D_RC_si = cat(2, PFC_D_RC_base, PFC_D_RC_test);

[PFC_D_GC_base, PFC_D_GC_test] = get_si_tuned_neuron_index(bt.full_PFC_D_GC, rc.full_PFC_D_GC, gc.full_PFC_D_GC, yt.full_PFC_D_GC, 500, 50);
PFC_D_GC_si = cat(2, PFC_D_GC_base, PFC_D_GC_test);

[PFC_D_YT_base, PFC_D_YT_test] = get_si_tuned_neuron_index(bt.full_PFC_D_YT, rc.full_PFC_D_YT, gc.full_PFC_D_YT, yt.full_PFC_D_YT, 500, 50);
PFC_D_YT_si = cat(2, PFC_D_YT_base, PFC_D_YT_test); 