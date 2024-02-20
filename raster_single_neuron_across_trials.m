%% pltting set up
total_trial_num = 500;
%   Simulation time
dt = 0.01; %step size ms
t_final = 4000; %simulation time ms
end_time = 3500/dt;
T = 0:dt:t_final;
gauss_width= 100;

neuron_id = 9;
PFC_S_BT = zeros(total_trial_num, length(T));
PFC_S_RC = zeros(total_trial_num, length(T));
PFC_S_GC = zeros(total_trial_num, length(T));
PFC_S_YT = zeros(total_trial_num, length(T));

PFC_D_BT = zeros(total_trial_num, length(T));
PFC_D_RC = zeros(total_trial_num, length(T));
PFC_D_GC = zeros(total_trial_num, length(T));
PFC_D_YT = zeros(total_trial_num, length(T));

VA_shape = zeros(total_trial_num, length(T));
VA_ori = zeros(total_trial_num, length(T));

PFC_remote_shape = zeros(total_trial_num, length(T));
PFC_remote_ori = zeros(total_trial_num, length(T));

MD_shape = zeros(total_trial_num, length(T));
MD_ori = zeros(total_trial_num, length(T));

for i = 1:total_trial_num
    temp = full_PFC_S_BT{i};
    if ~isnan(temp)
        PFC_S_BT(i, temp(temp(:, 1) == neuron_id, 2)) = 1;
    end
    
    temp = full_PFC_S_RC{i};
    if ~isnan(temp)
        PFC_S_RC(i, temp(temp(:, 1) == neuron_id, 2)) = 1;
    end

    temp = full_PFC_S_GC{i};
    if ~isnan(temp)
        PFC_S_GC(i, temp(temp(:, 1) == neuron_id, 2)) = 1;
    end

    temp = full_PFC_S_YT{i};
    if ~isnan(temp)
        PFC_S_YT(i, temp(temp(:, 1) == neuron_id, 2)) = 1;
    end

    temp = full_PFC_D_BT{i};
    if ~isnan(temp)
        PFC_D_BT(i, temp(temp(:, 1) == neuron_id, 2)) = 1;
    end

    temp = full_PFC_D_RC{i};
    if ~isnan(temp)
        PFC_D_RC(i, temp(temp(:, 1) == neuron_id, 2)) = 1;
    end

    temp = full_PFC_D_GC{i};
    if ~isnan(temp)
        PFC_D_GC(i, temp(temp(:, 1) == neuron_id, 2)) = 1;
    end

    temp = full_PFC_D_YT{i};
    if ~isnan(temp)
        PFC_D_YT(i, temp(temp(:, 1) == neuron_id, 2)) = 1;
    end

    temp = full_VA_shape{i};
    if ~isnan(temp)
        VA_shape(i, temp(temp(:, 1) == neuron_id, 2)) = 1;
    end

    temp = full_VA_ori{i};
    if ~isnan(temp)
        VA_ori(i, temp(temp(:, 1) == neuron_id, 2)) = 1;
    end

    temp = full_PFC_shape_ensemble{i};
    if ~isnan(temp)
        PFC_remote_shape(i, temp(temp(:, 1) == neuron_id, 2)) = 1;
    end

    temp = full_PFC_ori_ensemble{i};
    if ~isnan(temp)
        PFC_remote_ori(i, temp(temp(:, 1) == neuron_id, 2)) = 1;
    end

    temp = full_MD_shape{i};
    if ~isnan(temp)
        MD_shape(i, temp(temp(:, 1) == neuron_id, 2)) = 1;
    end

    temp = full_MD_ori{i};
    if ~isnan(temp)
        MD_ori(i, temp(temp(:, 1) == neuron_id, 2)) = 1;
    end
end

%%

spy( PFC_S_BT(:, t_start:t_end), 8,'b'),title('PFC Superficial blue triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[100 1 1]),ylabel('Trial number', 'FontSize',16),xlabel('time(ms)');
xticks([0 50000, 140000])
xticklabels([-500, 0, 900])
saveas(gcf,'PFC_S_BT.pdf')

spy( PFC_S_RC(:, t_start:t_end), 8,'b'),title('PFC Superficial red circle ', 'FontSize', 16)
set(gca,'DataAspectRatio',[100 1 1], 'XTick', -500:900:100),ylabel('Trial number', 'FontSize',16),xlabel('time(ms)');
xticks([0 50000, 140000])
xticklabels([-500, 0, 900])
saveas(gcf,'PFC_S_RC.pdf')

spy( PFC_S_GC(:, t_start:t_end), 8,'b'),title('PFC Superficial green circle ', 'FontSize', 16)
set(gca,'DataAspectRatio',[100 1 1], 'XTick', -500:900:100),ylabel('Trial number', 'FontSize',16),xlabel('time(ms)');
xticks([0 50000, 140000])
xticklabels([-500, 0, 900])
saveas(gcf,'PFC_S_GC.pdf')

spy( PFC_S_YT(:, t_start:t_end),8,'b'),title('PFC Superficial yellow triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[100 1 1], 'XTick', -500:900:100),ylabel('Trial number', 'FontSize',16),xlabel('time(ms)');
xticks([0 50000, 140000])
xticklabels([-500, 0, 900])
saveas(gcf,'PFC_S_YT.pdf')

spy( PFC_D_BT(:, t_start:t_end),8,'b'),title('PFC Deep blue triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[100 1 1], 'XTick', -500:900:100),ylabel('Trial number', 'FontSize',16),xlabel('time(ms)');
xticks([0 50000, 140000])
xticklabels([-500, 0, 900])
saveas(gcf,'PFC_D_BT.pdf')

spy( PFC_D_RC(:, t_start:t_end),8,'b'),title('PFC Deep red circle', 'FontSize', 16)
set(gca,'DataAspectRatio',[100 1 1]),ylabel('Trial number', 'FontSize',16),xlabel('time(ms)');
xticks([0 50000, 140000])
xticklabels([-500, 0, 900])
saveas(gcf,'PFC_D_RC.pdf')

spy( PFC_D_GC(:, t_start:t_end),8,'b'),title('PFC Deep green circle', 'FontSize', 16)
set(gca,'DataAspectRatio',[100 1 1]),ylabel('Trial number', 'FontSize',16),xlabel('time(ms)');
xticks([0 50000, 140000])
xticklabels([-500, 0, 900])
saveas(gcf,'PFC_D_GC.pdf')

spy( PFC_D_YT(:, t_start:t_end),8,'b'),title('PFC Deep yellow triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[100 1 1]),ylabel('Trial number', 'FontSize',16),xlabel('time(ms)');
xticks([0 50000, 140000])
xticklabels([-500, 0, 900])
saveas(gcf,'PFC_D_YT.pdf')


spy( VA_shape(:, t_start:t_end),8,'b'),title('VA Thalamus Shape', 'FontSize', 16)
set(gca,'DataAspectRatio',[100 1 1]),ylabel('Trial number', 'FontSize',16),xlabel('time(ms)');
xticks([0 50000, 140000])
xticklabels([-500, 0, 900])
saveas(gcf,'VA_shape.pdf')

spy( VA_ori(:, t_start:t_end),8,'b'),title('VA Thalamus Orientation', 'FontSize', 16)
set(gca,'DataAspectRatio',[100 1 1]),ylabel('Trial number', 'FontSize',16),xlabel('time(ms)');
xticks([0 50000, 140000])
xticklabels([-500, 0, 900])
saveas(gcf,'VA_orientation.pdf')

spy( PFC_remote_shape(:, t_start:t_end),8,'b'),title('PFC Shape Ensemble', 'FontSize', 16)
set(gca,'DataAspectRatio',[100 1 1]),ylabel('Trial number', 'FontSize',16),xlabel('time(ms)');
xticks([0 50000, 140000])
xticklabels([-500, 0, 900])
saveas(gcf,'PFC_shape.pdf')

spy( PFC_remote_ori(:, t_start:t_end),8,'b'),title('PFC Ori Ensemble', 'FontSize', 16)
set(gca,'DataAspectRatio',[100 1 1]),ylabel('Trial number', 'FontSize',16),xlabel('time(ms)');
xticks([0 50000, 140000])
xticklabels([-500, 0, 900])
saveas(gcf,'PFC_orientation.pdf')

spy( MD_shape(:, t_start:t_end),8,'b'),title('MD Shape', 'FontSize', 16)
set(gca,'DataAspectRatio',[100 1 1]),ylabel('Trial number', 'FontSize',16),xlabel('time(ms)');
xticks([0 50000, 140000])
xticklabels([-500, 0, 900])
saveas(gcf,'MD_shape.pdf')

spy( MD_ori(:, t_start:t_end),8,'b'),title('MD ori', 'FontSize', 16)
set(gca,'DataAspectRatio',[100 1 1]),ylabel('Trial number', 'FontSize',16),xlabel('time(ms)');
xticks([0 50000, 140000])
xticklabels([-500, 0, 900])
saveas(gcf,'MD_orientation.pdf')