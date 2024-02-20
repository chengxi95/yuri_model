%% pltting set up
trial_num = 1;
%   Simulation time
dt = 0.01; %step size ms
t_final = 4000; %simulation time ms
end_time = 3500/dt;
T = 0:dt:t_final;
gauss_width= 100;

t_PFC_S_BT = zeros(numberofneurons, length(T));
t_PFC_S_RC = zeros(numberofneurons, length(T));
t_PFC_S_GC = zeros(numberofneurons, length(T));
t_PFC_S_YT = zeros(numberofneurons, length(T));

t_PFC_D_BT = zeros(numberofneurons, length(T));
t_PFC_D_RC = zeros(numberofneurons, length(T));
t_PFC_D_GC = zeros(numberofneurons, length(T));
t_PFC_D_YT = zeros(numberofneurons, length(T));

t_VA_shape = zeros(numberofneurons, length(T));
t_VA_ori = zeros(numberofneurons, length(T));

t_PFC_remote_shape = zeros(numberofneurons, length(T));
t_PFC_remote_ori = zeros(numberofneurons, length(T));

t_MD_shape = zeros(numberofneurons, length(T));
t_MD_ori = zeros(numberofneurons, length(T));

[r, ~] = size(full_PFC_S_BT{trial_num});
for i = 1:r
    if ~isnan(full_PFC_S_BT{trial_num})
        t_PFC_S_BT(full_PFC_S_BT{trial_num}(i, 1), full_PFC_S_BT{trial_num}(i, 2)) = 1;
    end
end

[r, ~] = size(full_PFC_S_RC{trial_num});
for i = 1:r
    if ~isnan(full_PFC_S_RC{trial_num})
        t_PFC_S_RC(full_PFC_S_RC{trial_num}(i, 1), full_PFC_S_RC{trial_num}(i, 2)) = 1;
    end
end

[r, ~] = size(full_PFC_S_GC{trial_num});
for i = 1:r
    if ~isnan(full_PFC_S_GC{trial_num})
        t_PFC_S_GC(full_PFC_S_GC{trial_num}(i, 1), full_PFC_S_GC{trial_num}(i, 2)) = 1;
    end
end

[r, ~] = size(full_PFC_S_YT{trial_num});
for i = 1:r
    if ~isnan(full_PFC_S_YT{trial_num})
        t_PFC_S_YT(full_PFC_S_YT{trial_num}(i, 1), full_PFC_S_YT{trial_num}(i, 2)) = 1;
    end
end

[r, ~] = size(full_PFC_D_BT{trial_num});
for i = 1:r
    if ~isnan(full_PFC_D_BT{trial_num})
        t_PFC_D_BT(full_PFC_D_BT{trial_num}(i, 1), full_PFC_D_BT{trial_num}(i, 2)) = 1;
    end
end

[r, ~] = size(full_PFC_D_RC{trial_num});
for i = 1:r
    if ~isnan(full_PFC_D_RC{trial_num})
        t_PFC_D_RC(full_PFC_D_RC{trial_num}(i, 1), full_PFC_D_RC{trial_num}(i, 2)) = 1;
    end
end

[r, ~] = size(full_PFC_D_GC{trial_num});
for i = 1:r
    if ~isnan(full_PFC_D_GC{trial_num})
        t_PFC_D_GC(full_PFC_D_GC{trial_num}(i, 1), full_PFC_D_GC{trial_num}(i, 2)) = 1;
    end
end

[r, ~] = size(full_PFC_D_YT{trial_num});
for i = 1:r
    if ~isnan(full_PFC_D_YT{trial_num})
        t_PFC_D_YT(full_PFC_D_YT{trial_num}(i, 1), full_PFC_D_YT{trial_num}(i, 2)) = 1;
    end
end

[r, ~] = size(full_VA_shape{trial_num});
for i = 1:r
    if ~isnan(full_VA_shape{trial_num})
        t_VA_shape(full_VA_shape{trial_num}(i, 1), full_VA_shape{trial_num}(i, 2)) = 1;
    end
end

[r, ~] = size(full_VA_ori{trial_num});
for i = 1:r
    if ~isnan(full_VA_ori{trial_num})
        t_VA_ori(full_VA_ori{trial_num}(i, 1), full_VA_ori{trial_num}(i, 2)) = 1;
    end
end

[r, ~] = size(full_PFC_shape_ensemble{trial_num});
for i = 1:r
    if ~isnan(full_PFC_shape_ensemble{trial_num})
        t_PFC_remote_shape(full_PFC_shape_ensemble{trial_num}(i, 1), full_PFC_shape_ensemble{trial_num}(i, 2)) = 1;
    end
end

[r, ~] = size(full_PFC_ori_ensemble{trial_num});
for i = 1:r
    if ~isnan(full_PFC_ori_ensemble{trial_num})
        t_PFC_remote_ori(full_PFC_ori_ensemble{trial_num}(i, 1), full_PFC_ori_ensemble{trial_num}(i, 2)) = 1;
    end
end

[r, ~] = size(full_MD_shape{trial_num});
for i = 1:r
    if ~isnan(full_MD_shape{trial_num})
        t_MD_shape(full_MD_shape{trial_num}(i, 1), full_MD_shape{trial_num}(i, 2)) = 1;
    end
end

[r, ~] = size(full_MD_ori{trial_num});
for i = 1:r
    if ~isnan(full_MD_ori{trial_num})
        t_MD_ori(full_MD_ori{trial_num}(i, 1), full_MD_ori{trial_num}(i, 2)) = 1;
    end
end

figure(1)

subplot(4,1,1)
spy( t_PFC_S_BT,8,'b'),title('PFC Superficial blue triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 end_time]);%,xlabel('time')


subplot(4,1,2)
spy( t_PFC_S_RC,8,'b'),title('PFC Superficial red circle ', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 end_time]);%,xlabel('time')

subplot(4,1,3)
spy( t_PFC_S_GC,8,'b'),title('PFC Superficial green circle ', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 end_time]);%,xlabel('time')

subplot(4,1,4)
spy(t_PFC_S_YT,8,'b'),title('PFC Superficial yellow triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 end_time]);%,xlabel('time')
%
figure (2)
subplot(4,1,1)
spy( t_PFC_D_BT,8,'k'),title('PFC Deep blue triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 end_time]);%,xlabel('time')
subplot(4,1,2)
spy( t_PFC_D_RC,8,'k'),title('PFC Deep red circle', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 end_time]);%,xlabel('time')

subplot(4,1,3)
spy( t_PFC_D_GC,8,'k'),title('PFC Deep green circle', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 end_time]);%,xlabel('time')
subplot(4,1,4)
spy( t_PFC_D_YT,8,'k'),title('PFC Deep yellow triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 end_time]);%,xlabel('time')


figure(3)
subplot(6,1,1)
spy( t_VA_shape,8,'k'),title('VA Thalamus Shape', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 end_time]);%,xlabel('time')

subplot(6,1,2)
spy( t_VA_ori,8,'k'),title('VA Thalamus Ori', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 end_time]);%,xlabel('time')

subplot(6,1,3)
spy( t_PFC_remote_shape,8,'b'),title('PFC Shape Ensemble', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 end_time]);%,xlabel('time')

subplot(6,1,4)
spy( t_PFC_remote_ori,8,'b'),title('PFC Ori Ensemble', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 end_time]);%,xlabel('time')

subplot(6,1,5)
spy( t_MD_shape,8,'k'),title('MD shape', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 end_time]);%,xlabel('time')

subplot(6,1,6)
spy( t_MD_ori,8,'k'),title('MD ori', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 end_time]);%,xlabel('time')