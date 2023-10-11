
close all;
clear all;
% this code simulate two populations of cortical cells (aPFC and pPFC superficial and deep layers)  and three population of thalamic cells (aPFC core, pPFC core and Matrix cells)during a working memory task.
% our assumption is that thalamus increases the
% cortico-cortical connectivity gain by a coefiecint an amplification
% factor for local and distant connections.
% which is bigger than one in the awake state and is equal to 1 in the coma state.
% in the awake state thalamic cells have shorter timeconstant so higher firing frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   constants
numberofneurons = 50;% number of neurons per group

%   Time constants
tha = 20; %time constant

%   Simulation time
dt = 0.01; %step size ms
t_final = 3000; %simulation time ms
T = 0:dt:t_final;

%   Intrinsic property of neuron
delay = 5/dt; % 3ms
v_th = -50; %mv
E_L = -65; %mv
RM = 10;
leaky_coef = 1; %
tref = 1/dt; % 1 ms refratory time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neurons talk to each other by sending short pulses, define duration and
% amplitude of these pulses
abs_pfc_width = 1.5;% cortical input liftime
spikewidth = abs_pfc_width/dt;
spikewidth_inh = 100/dt; % assumed that inhibition effect is prolonged due to the synch exc input SOM neurons
spikewidth_MD = 100/dt;

%   Connectivity Map: local PFC and PRL connectivity = all to all connectivity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Stimulation showing the abstract cue
no_of_trl = 12;
column_length = 10;
I_stim = zeros(numberofneurons,length(T));
t_start_stim_1_abs = 1000; % cue time
cue_amp = 0.5;
abs_stim_duration = 450; %0.1;% stim duration is one ms
stim_duration = abs_stim_duration/dt;
abs_stim_intervals = 0;
stim_intervals = abs_stim_intervals/dt;
t_start_stim_1 = t_start_stim_1_abs/dt;
for i=1:no_of_trl
    I_stim(1:column_length, t_start_stim_1: t_start_stim_1+ stim_duration) = cue_amp; % first stimulus delivered to visual starter neuron
    t_start_stim_1 = t_start_stim_1+ 1000/dt;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   local cortical connectivity

W_local = 0.002; % synaptic weights
W_M= 0.02;
matrix_local = zeros(numberofneurons*2,numberofneurons*2);
matrix_M= zeros(numberofneurons*2,numberofneurons*2);
for ii = 1:column_length: 2*numberofneurons-column_length

    matrix_local(ii:ii+ column_length - 1,ii + column_length:ii+ 2*column_length-1) = W_local;
    matrix_M(ii:ii+ column_length - 1, ii :ii+ column_length-1) = W_M;
end
matrix_local = matrix_local - diag(diag(matrix_local));
matrix_local(1:10, 11:20)= 2* matrix_local(1:10, 11:20);
matrix_M(41:50,1:10 )= W_M;
matrix_local(41:50,1:10 )= W_local;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   superficial to deep layers
W_S_D = 0.004;% 0.0075;
matrix_local_S_D = zeros(numberofneurons,numberofneurons);
for ii = 1:column_length: numberofneurons-column_length+1

    matrix_local_S_D(ii:ii+ column_length - 1,ii: ii + column_length - 1) = W_S_D;

end
matrix_local_S_D(1:10, 1:10) = 2*W_S_D;

%%%%%%%%%%%%%%%%%
%connections from superficial layer to the deep layer
W11=2; W12=0.8; W13=0.2; W14= 0.2;
W21=0.8; W22=2; W23=0.2; W24= 0.2;
W31=0.2; W32=0.2; W33=2; W34= 0.8;
W41=0.2; W42=0.2; W43=0.8; W44= 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   Superficial-layers connectivity
% W_D_toremote = 0.002;
% matrix_local_D_toremote = zeros(numberofneurons,numberofneurons);
% matrix_local_D_toremote (20:40, 1:10)= W_D_toremote;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PFC Deep to Thalamus
W_PFC_TH =  0.009;% 0.0175;% 0.05;
W_PFC_MD= 0.002;
PFC_VA_matrix = zeros(numberofneurons,numberofneurons);
PFC_VA_matrix(1: 20, 1:numberofneurons) = W_PFC_TH;
PFC_MD_matrix = zeros(numberofneurons,numberofneurons);
for  ii = 1:column_length: numberofneurons-column_length
    PFC_MD_matrix(ii:ii+ column_length - 1,ii:ii+ column_length - 1) = W_PFC_TH;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PFC Deep to STR
W_PFC_to_str= 0.0145;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %md to inh
% W_MD_to_PV= 0.0015;
% Matrix_MD_to_PV = zeros(numberofneurons,numberofneurons);
% for  ii = 1:column_length: numberofneurons-column_length
%     Matrix_MD_to_PV(ii:ii+ column_length - 1,ii) = W_MD_to_PV;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MD to superficial PFC layers
matrix_MD_to_Crtex = zeros(numberofneurons,numberofneurons);
matrix_VA_to_rule= zeros(numberofneurons,numberofneurons);
matrix_VA_to_rule(1:50, 1:10)= 1;
matrix_MD_to_Crtex(1:10, 1:50) = 1;                                          %????????????????????????????????????????????
matrix_MD_to_shape= matrix_local;
matrix_MD_to_orientation= matrix_local;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Connectivity Map: Excitatory and inhibitory cell connectivity = localy, no long distance inhibitory connection, inhibition only local within population
W_IPL_Inh_to_exc = 0.001;%25;%synaptic weight
W_IPL_Inh_to_exc_pv = 0.01;%25;%synaptic weight
W_IPL_Inh_to_exc_fs= 0.005;%25;%synaptic weight

W_VA_exi = 0.005;
W_MD_inh = 0.003;
W_MD_PFC = 0.012;

for rr= 1


    % Initialize membrane  potential---------------------------------------
    y_PFC_S_BT= zeros(numberofneurons,length(T));
    y_PFC_S_BT(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;%initializing the excitatory PFC neurons Voltages

    y_PFC_S_RC = zeros(numberofneurons,length(T));
    y_PFC_S_RC(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;%initializing the excitatory PFC neurons Voltages

    y_PFC_S_GC = zeros(numberofneurons,length(T));
    y_PFC_S_GC(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;

    y_PFC_S_YT = zeros(numberofneurons,length(T));
    y_PFC_S_YT(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;

    y_PFC_M = zeros(numberofneurons,length(T));
    y_PFC_M(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;

    y_aPFC_D_BT = zeros(numberofneurons,length(T));
    y_aPFC_D_BT(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;%initializing the excitatory PFC neurons Voltages

    y_aPFC_D_RC= zeros(numberofneurons,length(T));
    y_aPFC_D_RC(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;%initializing the excitatory PFC neurons Voltages

    y_aPFC_D_GC= zeros(numberofneurons,length(T));
    y_aPFC_D_GC(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;%initializing the excitatory PFC neurons Voltages

    y_aPFC_D_YT= zeros(numberofneurons,length(T));
    y_aPFC_D_YT(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;%initializing the excitatory PFC neurons Voltages

    y_aPFC_D_6 = zeros(numberofneurons,length(T));
    y_aPFC_D_6(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;%initializing the excitatory PFC neurons Voltages

    y_MD_common = zeros(numberofneurons,length(T));
    y_MD_common(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;%initializing the excitatory PFC neurons Voltages

    y_pPFC_D = zeros(numberofneurons,length(T));
    y_pPFC_D(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;%initializing the excitatory PFC neurons Voltages
    y_aPFC_inh = zeros(numberofneurons,length(T));
    y_aPFC_inh(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5; %initializing inhibitory PFC neurons Voltages
    y_response_Ori_right = zeros(numberofneurons,length(T));
    y_response_Ori_right(1:numberofneurons,1)=-50+rand(numberofneurons,1)*1;
    y_response_Ori_left = zeros(numberofneurons,length(T));
    y_response_Ori_left(1:numberofneurons,1)=-50+rand(numberofneurons,1)*1;
    y_pPFC_ORI = zeros(numberofneurons,length(T));
    y_pPFC_ORI(1:numberofneurons,1)= -50+ rand(numberofneurons,1)*5; %initializing the excitatory IPL neurons Voltages
    y_pPFC_ORI_LEFT = zeros(numberofneurons,length(T));
    y_pPFC_ORI_LEFT(1:numberofneurons,1)= -50 + rand(numberofneurons,1)*5; %initializing the excitatory IPL neurons Voltages
    y_pPFC_Shape = zeros(numberofneurons,length(T));
    y_pPFC_Shape(1:numberofneurons,1)=-50+rand(numberofneurons,1)*5;

    y_pPFC_Orientation = zeros(numberofneurons,length(T));
    y_pPFC_Orientation(1:numberofneurons,1)=-50+rand(numberofneurons,1)*5;

    y_VA_matrix_shape = zeros(numberofneurons,length(T));
    y_VA_matrix_shape(1:numberofneurons,1)=-50+rand(numberofneurons,1)*1; %initializing the thalamic neurons Voltages

    y_VA_matrix_Orientation  = zeros(numberofneurons,length(T));
    y_VA_matrix_Orientation  =-50+rand(numberofneurons,1)*1; %initializing the thalamic neurons Voltages


    y_MD_core_ori = zeros(numberofneurons,length(T));
    y_MD_core_ori(1:numberofneurons,1)=-50+rand(numberofneurons,1)*1; %initializing the thalamic neurons Voltages
    y_MD_core_shape = zeros(numberofneurons,length(T));
    y_MD_core_shape(1:numberofneurons,1)=-50+rand(numberofneurons,1)*1; %initializing the thalamic neurons Voltages
    y_SNpr_inh=zeros(numberofneurons,length(T));
    y_SNpr_inh(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;% initializing the PFC neurons Voltages
    y_ST_inh=zeros(numberofneurons,length(T));
    y_ST_inh(1:numberofneurons,1)=-50+rand(numberofneurons,1)*0;% initializing the PFC neurons Voltages
    y_PV_inh= zeros(numberofneurons,length(T));
    y_PV_inh (1:numberofneurons,1)=-50+rand(numberofneurons,1)*1;
    y_FS_inh= zeros(numberofneurons,length(T));
    y_FS_inh(1:numberofneurons,1)=-50+rand(numberofneurons,1)*1;

    % driving currents-----------------------------------------------------
    Iext_IPL = 1.5 - 0.01 *(rand(numberofneurons,1)); % similar driving current for PFC and IPL cells small variation across cells
    Iext_PFC = 1.5 - 0.01 *(rand(numberofneurons,1)); % similar driving current for PFC and IPL cells small variation across cells
    Iext_PFC_M = 1.5 - 0.01 *(rand(numberofneurons,1));
    Iext_PFC_remote_ORI = 1.5 - 0.01 *(rand(numberofneurons,1)); % similar driving current for PFC and IPL cells small variation across cells
    Iext_PFC_remote_Shape = 1.5 - 0.01 *(rand(numberofneurons,1)); % similar driving current for PFC and IPL cells small variation across cells
    Iext_PFC_D = 1.5 - 0.01 *(rand(numberofneurons,1)); % similar driving current for PFC and IPL cells small variation across cells
    Iext_PFC_D_2 = 1.5 - 0.01 *(rand(numberofneurons,1)); % similar driving current for PFC and IPL cells small variation across cells
    Iext_PFC_D_3 = 1.5 - 0.01 *(rand(numberofneurons,1)); % similar driving current for PFC and IPL cells small variation across cells
    Iext_PFC_2 = 1.5 - 0.01 *(rand(numberofneurons,1)); % similar driving current for PFC and IPL cells small variation across cells
    Iext_PFC_3 = 1.5 - 0.01 *(rand(numberofneurons,1)); % similar driving current for PFC and IPL cells small variation across cells
    Iext_md = 1.5 - 0.01 *rand(numberofneurons,1); % same driving current for all thalamic cells
    Iext_VA = 1.5 - 0.01 *rand(numberofneurons,1); % same driving current for all thalamic cells
    Iext_Inh = 1.5 -  0.01 *(rand(numberofneurons,1)-0.5); % driving current for inhibitory cells similar in PFC and IPL
    Iext_Inh_pPFC = 1.5 - 0.01 *(rand(numberofneurons,1)-0.5); %driving current for inhibitory cells similar in PFC and IPL
    Iext_Inh_PV = 1.49 - 0.0 *(rand(numberofneurons,1));
    Iext_response = 1.5 - 0.01 *(rand(numberofneurons,1));

    % initialize noise input-----------------------------------------------
    I_noise_md = zeros(numberofneurons,length(T));
    I_noise_PFC = zeros(numberofneurons,length(T));

    % initialize spike counting--------------------------------------------
    last_spike_IPL=10^10*ones(numberofneurons,1);
    last_spike_response_ori_right = 10^10*ones(numberofneurons,1);
    last_spike_PFC_S_BT=10^10*ones(numberofneurons,1);
    last_spike_PFC_S_GC=10^10*ones(numberofneurons,1);
    last_spike_PFC_S_YT=10^10*ones(numberofneurons,1);

    last_spike_PFC_S_RC =10^10*ones(numberofneurons,1);
    last_spike_aPFC_2 = 10^10*ones(numberofneurons,1);
    last_spike_aPFC_3 = 10^10*ones(numberofneurons,1);
    last_spike_PFC_D_BT = 10^10*ones(numberofneurons,1);
    last_spike_aPFC_D_YT= 10^10*ones(numberofneurons,1);
    last_spike_aPFC_D_GC= 10^10*ones(numberofneurons,1);
    last_spike_aPFC_D_RC = 10^10*ones(numberofneurons,1);
    last_spike_aPFC_D6 = 10^10*ones(numberofneurons,1);
    last_spike_MD_ori = 10^10*ones(numberofneurons,1);
    last_spike_pPFC_D = 10^10*ones(numberofneurons,1);
    last_spike_VA_matrix_shape = 10^10*ones(numberofneurons,1);
    last_spike_VA_matrix_Orientation= 10^10*ones(numberofneurons,1);
    last_spike_VA_matrix_Ori= 10^10*ones(numberofneurons,1);


    last_spike_aPFC_D_2 = 10^10*ones(numberofneurons,1);
    last_spike_aPFC_D_3 = 10^10*ones(numberofneurons,1);
    last_spike_aPFC_Inh = 10^10*ones(numberofneurons,1);
    last_spike_SNpr_Inh = 10^10*ones(numberofneurons,1);
    last_spike_IPL_remote = 10^10*ones(numberofneurons,1);
    last_spike_MD_shape = 10^10*ones(numberofneurons,1);
    last_spike_MD2_core = 10^10*ones(numberofneurons,1);
    last_spike_CMT = 10^10*ones(numberofneurons,1);
    last_spike_pPFC_remote_ORI = 10^10*ones(numberofneurons,1);
    last_spike_pPFC_remote_ORI_LEFT = 10^10*ones(numberofneurons,1);
    last_spike_pPFC_remote_shape= 10^10*ones(numberofneurons,1);
    last_spike_pPFC_remote_Orientation= 10^10*ones(numberofneurons,1);
    last_spike_PV= 10^10*ones(numberofneurons,1);
    last_spike_FS= 10^10*ones(numberofneurons,1);
    last_spike_ST = 10^10*ones(numberofneurons,1);
    last_spike_PFC_M= 10^10*ones(numberofneurons,1);
    last_spike_MD_common= 10^10*ones(numberofneurons,1);

    % spike holders--------------------------------------------------------
    spiketimes_S = [];
    spiketimes_S_2 = [];
    spiketimes_S_3 = [];
    spiketimes_M =[];
    spiketimes_D_BT= [];
    spiketimes_D_RC = [];
    spiketimes_MD = [];
    spiketimes_MD2 = [];
    spiketimes_MD_common=[];
    spiketimes_VA = [];
    spiketimes_VA_shape =[];
    spiketimes_VA_ori =[];
    spiketimes_D_2 = [];
    spiketimes_D_3 = [];
    spiketimes_pPFC_D = [];
    spiketimes_CMT = [];
    spiketimes_ipl_inh =[];
    spiketimes_IPL =[];
    spiketimes_S_remote_ORI = [];
    spiketimes_S_remote_ORI_LEFT = [];
    spiketimes_S_remote_Shape = [];
    spiketimes_S_remote_Orientation=[];
    spiketimes_PV_inh = [];
    spiketimes_FS_inh=[];
    spiketimes_response = [];
    spiketimes_D6 = [];
    spiketimes_S_blue = [];
    spiketimes_S_red = [];
    spiketimes_S_YT = [];
    spiketimes_D_YT = [];
    spiketimes_D_GC= [];
    spiketimes_D_RC =[];

    % initialize synaptic input current------------------------------------
    Isyn_aPFC_local = zeros(numberofneurons,length(T));
    Isyn_aPFC_local_shared = zeros(numberofneurons,length(T));
    Isyn_aPFC_local_BT= zeros(numberofneurons,length(T));
    Isyn_aPFC_local_RC= zeros(numberofneurons,length(T));
    Isyn_aPFC_local_GC= zeros(numberofneurons,length(T));
    Isyn_aPFC_local_YT= zeros(numberofneurons,length(T));

    Isyn_PFC_M_S= zeros(numberofneurons,length(T));
    Isyn_aPFC_to_D = zeros(numberofneurons,length(T));
    Isyn_aPFC_D = zeros(numberofneurons,length(T));
    Isyn_aPFC_D_shared = zeros(numberofneurons,length(T));

    Isyn_pPFC_D = zeros(numberofneurons,length(T));
    Isyn_pPFC_S_to_D = zeros(numberofneurons,length(T));
    Isyn_pPFC_local = zeros(numberofneurons,length(T));
    Isyn_pPFC_local_left = zeros(numberofneurons,length(T));
    Isyn_pPFCD_to_response = zeros(numberofneurons,length(T));
    Isyn_aPFC_IPL = zeros(numberofneurons,length(T));
    Isyn_aPFC_Exc_to_Inh = zeros(numberofneurons,length(T));
    Isyn_aPFC_Inh_to_exc = zeros(numberofneurons,length(T));
    Isyn_SNpr_to_VA = zeros(numberofneurons,length(T));
    Isyn_IPL = zeros(numberofneurons,length(T));
    Isyn_IPL_PFC = zeros(numberofneurons,length(T));
    Isyn_IPL_exc_inh = zeros(numberofneurons,length(T));
    Isyn_IPL_loc_remote = zeros(numberofneurons,length(T));
    Isyn_IPL_md = zeros(numberofneurons,length(T));
    Isyn_aPFC_cmt = zeros(numberofneurons,length(T));
    Isyn_PFC_D_MD_ori = zeros(numberofneurons,length(T));
    Isyn_PFC_D_MD_shape_common= zeros(numberofneurons,length(T));
    Isyn_PFC_D_MD_shape = zeros(numberofneurons,length(T));
    Isyn_MD_to_inh = zeros(numberofneurons,length(T));
    Isyn_MD_to_inh_PV = zeros(numberofneurons,length(T));
    Isyn_MD_shape_to_PFC= zeros(numberofneurons,length(T));
    Isyn_MD2_to_pPFC = zeros(numberofneurons,length(T));
    Isyn_VA_Matrix_to_PFC = ones(numberofneurons,length(T));
    Isyn_VA_Matrix_to_PFC_exc_shape= zeros(numberofneurons,length(T));
    Isyn_VA_Matrix_to_PFC_exc_orientation= zeros(numberofneurons,length(T));
    Isyn_ST = zeros(numberofneurons,length(T));
    Isyn_PV = zeros(numberofneurons,length(T));
    Isyn_MD_to_Inh = zeros(numberofneurons,length(T));
    Isyn_FS= zeros(numberofneurons,length(T));
    Isyn_aPFC_local_effective = zeros(numberofneurons,length(T));
    Isyn_aPFC_D_effective = zeros(numberofneurons,length(T));
    Isyn_aPFC_D_effective_2 = zeros(numberofneurons,length(T));
    Isyn_aPFC_to_D_effective = zeros(numberofneurons,length(T));
    Isyn_pPFC_local_effective = zeros(numberofneurons,length(T));
    Isyn_pPFC_local_effective_2 = zeros(numberofneurons,length(T));
    Isyn_pPFC_local_effective_3 = zeros(numberofneurons,length(T));
    Isyn_pPFC_S_to_D_effective = zeros(numberofneurons,length(T));
    Isyn_pPFC_D_effective = zeros(numberofneurons,length(T));
    Isyn_pPFCD_to_response_effective = zeros(numberofneurons,length(T));
    Isyn_aPFC_D5 = zeros(numberofneurons,length(T));
    Isyn_aPFC_D5_to_ST = zeros(numberofneurons,length(T));
    Isyn_VA_Matrix_to_Inh= zeros(numberofneurons,length(T));
    Isyn_aPFC_D5_toremote = zeros(numberofneurons,length(T));
    Isyn_aPFC_D5_orientation = zeros(numberofneurons,length(T));
    Isyn_aPFC_to_D_RC= zeros(numberofneurons,length(T));
    Isyn_aPFC_to_D_GC= zeros(numberofneurons,length(T));
    Isyn_aPFC_to_D_YT= zeros(numberofneurons,length(T));
    Isyn_aPFC_to_D_BT= zeros(numberofneurons,length(T));
    Isyn_aPFC_local_shape= zeros(numberofneurons,length(T));

    % amplification = 1;
    Matrix_amplification = 2;


    % simulation starts here
    for i= 2:length(T)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Midle layer (l = 1) to superficial layer (l = 2)

        for jj = 1:numberofneurons
            for kk = 1:numberofneurons
                if (i-last_spike_PFC_M(kk))<(spikewidth) && (i-last_spike_PFC_M(kk))> 0
                    if matrix_M(kk,jj) ~= 0
                        Isyn_PFC_M_S(jj,i+ delay) = Isyn_PFC_M_S(jj,i+ delay)+ matrix_M(kk,jj);
                    end
                end
            end
        end
        %%%%%%%%% Superficial layer (l = 2) to superficial layer (l = 2)
        %%%%%%%%% excitatory to excitatory  within chains

        w_within_chains = 0.6;
        for jj = 1:numberofneurons
            for kk = 1:numberofneurons
                if (i-last_spike_PFC_S_BT(kk))<(spikewidth) && (i-last_spike_PFC_S_BT(kk))> 0
                    if matrix_local(kk,jj) ~= 0
                        Isyn_aPFC_local_BT(jj,i+ delay) = Isyn_aPFC_local_BT(jj,i+ delay)+ w_within_chains* matrix_local(kk,jj);
                    end
                end
            end
        end

        for jj = 1:numberofneurons
            for kk = 1:numberofneurons
                if (i-last_spike_PFC_S_RC(kk))<(spikewidth) && (i-last_spike_PFC_S_RC(kk))> 0
                    if matrix_local(kk,jj) ~= 0
                        Isyn_aPFC_local_RC(jj,i+ delay) = Isyn_aPFC_local_RC(jj,i+ delay)+ w_within_chains* matrix_local(kk,jj);
                    end
                end
            end
        end


        for jj = 1:numberofneurons
            for kk = 1:numberofneurons
                if (i-last_spike_PFC_S_YT(kk))<(spikewidth) && (i-last_spike_PFC_S_YT(kk))> 0
                    if matrix_local(kk,jj) ~= 0
                        Isyn_aPFC_local_YT(jj,i+ delay) = Isyn_aPFC_local_YT(jj,i+ delay)+ w_within_chains* matrix_local(kk,jj);
                    end
                end
            end
        end

        for jj = 1:numberofneurons
            for kk = 1:numberofneurons
                if (i-last_spike_PFC_S_GC(kk))<(spikewidth) && (i-last_spike_PFC_S_GC(kk))> 0
                    if matrix_local(kk,jj) ~= 0
                        Isyn_aPFC_local_GC(jj,i+ delay) = Isyn_aPFC_local_GC(jj,i+ delay)+ w_within_chains* matrix_local(kk,jj);
                    end
                end
            end
        end

        %%%%%%%%%%%Superficial layer (l = 2) to superficial layer (l = 2)
        %%%%%%%%%%% excitatory to excitatory across chains

        w_across_chains= 0.45;
        for jj = 1:numberofneurons
            for kk = 1:numberofneurons
                if (i-last_spike_PFC_S_BT(kk))<(spikewidth) && (i-last_spike_PFC_S_BT(kk))> 0
                    if matrix_local(kk,jj) ~= 0
                        Isyn_aPFC_local_shared(jj,i+ delay) = Isyn_aPFC_local_shared(jj,i+ delay)+ w_across_chains* matrix_local(kk,jj);
                    end
                end
            end
        end

        for jj = 1:numberofneurons
            for kk = 1:numberofneurons
                if (i-last_spike_PFC_S_RC(kk))<(spikewidth) && (i-last_spike_PFC_S_RC(kk))> 0
                    if matrix_local(kk,jj) ~= 0
                        Isyn_aPFC_local_shared(jj,i+ delay) = Isyn_aPFC_local_shared(jj,i+ delay)+ w_across_chains* matrix_local(kk,jj);
                    end
                end
            end
        end


        for jj = 1:numberofneurons
            for kk = 1:numberofneurons
                if (i-last_spike_PFC_S_YT(kk))<(spikewidth) && (i-last_spike_PFC_S_YT(kk))> 0
                    if matrix_local(kk,jj) ~= 0
                        Isyn_aPFC_local_shared(jj,i+ delay) = Isyn_aPFC_local_shared(jj,i+ delay)+ w_across_chains* matrix_local(kk,jj);
                    end
                end
            end
        end

        for jj = 1:numberofneurons
            for kk = 1:numberofneurons
                if (i-last_spike_PFC_S_GC(kk))<(spikewidth) && (i-last_spike_PFC_S_GC(kk))> 0
                    if matrix_local(kk,jj) ~= 0
                        Isyn_aPFC_local_shared(jj,i+ delay) = Isyn_aPFC_local_shared(jj,i+ delay)+ w_across_chains* matrix_local(kk,jj);
                    end
                end
            end
        end

        %%%%%%%%%%%Superficial layer (l = 2) to deep layer (l = 3)

        for jj = 1:numberofneurons
            for kk = 1:numberofneurons
                if matrix_local_S_D(kk,jj) ~= 0
                    if (i-last_spike_PFC_S_BT(kk))<(spikewidth) && (i-last_spike_PFC_S_BT(kk))> 0
                        Isyn_aPFC_to_D_BT(jj,i+ delay) = Isyn_aPFC_to_D_BT(jj,i+ delay) + matrix_local_S_D(kk,jj);
                    end
                end
            end
        end

        for jj = 1:numberofneurons
            for kk = 1:numberofneurons
                if matrix_local_S_D(kk,jj) ~= 0
                    if (i-last_spike_PFC_S_RC(kk))<(spikewidth) && (i-last_spike_PFC_S_RC(kk))> 0
                        Isyn_aPFC_to_D_RC(jj,i+ delay) = Isyn_aPFC_to_D_RC(jj,i+ delay) + matrix_local_S_D(kk,jj);
                    end
                end
            end
        end

        for jj = 1:numberofneurons
            for kk = 1:numberofneurons
                if matrix_local_S_D(kk,jj) ~= 0
                    if (i-last_spike_PFC_S_GC(kk))<(spikewidth) && (i-last_spike_PFC_S_GC(kk))> 0
                        Isyn_aPFC_to_D_GC(jj,i+ delay) = Isyn_aPFC_to_D_GC(jj,i+ delay) + matrix_local_S_D(kk,jj);
                    end
                end
            end
        end

        for jj = 1:numberofneurons
            for kk = 1:numberofneurons
                if matrix_local_S_D(kk,jj) ~= 0
                    if (i-last_spike_PFC_S_YT(kk))<(spikewidth) && (i-last_spike_PFC_S_YT(kk))> 0
                        Isyn_aPFC_to_D_YT(jj,i+ delay) = Isyn_aPFC_to_D_YT(jj,i+ delay) + matrix_local_S_D(kk,jj);
                    end
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%% Deep layer (l = 3) to deep layer (l = 3)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% excitatory to excitatory across chains
        for jj = 1:numberofneurons
            for kk = 1:numberofneurons
                if (i-last_spike_PFC_D_BT(kk))<(spikewidth) && (i-last_spike_PFC_D_BT(kk))> 0
                    if matrix_local(kk,jj) ~= 0
                        Isyn_aPFC_D_shared(jj,i+ delay) = Isyn_aPFC_D_shared(jj,i+ delay)+ w_across_chains* matrix_local(kk,jj);
                    end
                end
            end
        end

        for jj = 1:numberofneurons
            for kk = 1:numberofneurons
                if (i-last_spike_aPFC_D_GC(kk))<(spikewidth) && (i-last_spike_aPFC_D_GC(kk))> 0
                    if matrix_local(kk,jj) ~= 0
                        Isyn_aPFC_D_shared(jj,i+ delay) = Isyn_aPFC_D_shared(jj,i+ delay)+ w_across_chains* matrix_local(kk,jj);
                    end
                end
            end
        end

        for jj = 1:numberofneurons
            for kk = 1:numberofneurons
                if (i-last_spike_aPFC_D_RC(kk))<(spikewidth) && (i-last_spike_aPFC_D_RC(kk))> 0
                    if matrix_local(kk,jj) ~= 0
                        Isyn_aPFC_D_shared(jj,i+ delay) = Isyn_aPFC_D_shared(jj,i+ delay)+ w_across_chains* matrix_local(kk,jj);
                    end
                end
            end
        end

        for jj = 1:numberofneurons
            for kk = 1:numberofneurons
                if (i-last_spike_aPFC_D_YT(kk))<(spikewidth) && (i-last_spike_aPFC_D_YT(kk))> 0
                    if matrix_local(kk,jj) ~= 0
                        Isyn_aPFC_D_shared(jj,i+ delay) = Isyn_aPFC_D_shared(jj,i+ delay)+ w_across_chains* matrix_local(kk,jj);
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Deep layer (l=3) to VA thalamus (l=4) 
        %%%%  SHAPE
        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_PFC_D_BT(kk))<(spikewidth) && (i-last_spike_PFC_D_BT(kk))> 0  %%%%%%%%%%%%%%%%%%%%%%%%%%delay ==0
                    if PFC_VA_matrix(kk,jj)~=0
                        Isyn_aPFC_D5(jj,i+ delay)= Isyn_aPFC_D5(jj,i+ delay) + W_PFC_TH;
                    end
                end
            end
        end
        
        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_aPFC_D_RC(kk))<(spikewidth) && (i-last_spike_aPFC_D_RC(kk))> 0  %%%%%%%%%%%%%%%%%%%%%%%%%%delay ==0
                    if PFC_VA_matrix(kk,jj)~=0
                        Isyn_aPFC_D5(jj,i+ delay)= Isyn_aPFC_D5(jj,i+ delay) + W_PFC_TH;
                    end
                end
            end
        end


       %%%%  ORIENTATION
        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_aPFC_D_GC(kk))<(spikewidth) && (i-last_spike_aPFC_D_GC(kk))> 0  
                    if PFC_VA_matrix(kk,jj)~=0
                        Isyn_aPFC_D5_orientation(jj,i+ delay)= Isyn_aPFC_D5_orientation(jj,i+ delay) + W_PFC_TH;
                    end
                end
            end
        end



        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_aPFC_D_YT(kk))<(spikewidth) && (i-last_spike_aPFC_D_YT(kk))> 0 
                    if PFC_VA_matrix(kk,jj)~=0
                        Isyn_aPFC_D5_orientation(jj,i+ delay)= Isyn_aPFC_D5_orientation(jj,i+ delay) + W_PFC_TH;
                    end
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PFC deep  to MD 
% SHAPE 
        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_PFC_D_BT(kk))<(spikewidth) && (i-last_spike_PFC_D_BT(kk))> 0  
                    if matrix_M(kk,jj)~=0
                        Isyn_PFC_D_MD_shape(jj,i+ delay)= Isyn_PFC_D_MD_shape(jj,i+ delay) + W_PFC_MD;
                    end
                end
            end
        end


        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_aPFC_D_RC(kk))<(spikewidth) && (i-last_spike_aPFC_D_RC(kk))> 0 
                    if matrix_M(kk,jj)~=0
                        Isyn_PFC_D_MD_shape(jj,i+ delay)= Isyn_PFC_D_MD_shape(jj,i+ delay) + W_PFC_MD;
                    end
                end
            end
        end

        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_pPFC_remote_shape(kk))<(spikewidth) && (i-last_spike_pPFC_remote_shape(kk))> 0  
                    if matrix_M(kk,jj)~=0
                        Isyn_PFC_D_MD_shape(jj,i+ delay) = Isyn_PFC_D_MD_shape(jj,i+ delay)+ 3* W_PFC_MD;%%%%%%%%%%%%%%%%%% hard coding 
                    end
                end
            end
        end

% ORIENTATION 


        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_aPFC_D_YT(kk))<(spikewidth) && (i-last_spike_aPFC_D_YT(kk))> 0 
                    if PFC_VA_matrix(kk,jj)~=0
                        Isyn_PFC_D_MD_ori(jj,i+ delay)= Isyn_PFC_D_MD_ori(jj,i+ delay) + W_PFC_MD;
                    end
                end
            end
        end
        
        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_aPFC_D_GC(kk))<(spikewidth) && (i-last_spike_aPFC_D_GC(kk))> 0 
                    if PFC_VA_matrix(kk,jj)~=0
                        Isyn_PFC_D_MD_ori(jj,i+ delay)= Isyn_PFC_D_MD_ori(jj,i+ delay) + W_PFC_MD;
                    end
                end
            end
        end
        
        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_pPFC_remote_Orientation(kk))<(spikewidth) && (i-last_spike_pPFC_remote_Orientation(kk))> 0
                    if matrix_M(kk,jj)~=0
                        Isyn_PFC_D_MD_ori(jj,i+ delay) = Isyn_PFC_D_MD_ori(jj,i+ delay)+ W_PFC_MD;
                    end
                end
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VA thalamus to cortex 



        for jj=1:numberofneurons
            q=0;
            for kk=1:numberofneurons
                if (i-last_spike_VA_matrix_shape(kk))<(spikewidth_MD) && (i-last_spike_VA_matrix_shape(kk))> 0 
                    if matrix_MD_to_Crtex(kk,jj)~=0
                        q=q+1;
                        spike_count(jj, i)=q;
                        Isyn_VA_Matrix_to_PFC(jj,i+ delay)= Matrix_amplification; % md_coefficient* spike_count(jj, i);
                    end
                end
            end
        end


        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_VA_matrix_shape(kk))<(spikewidth_MD) && (i-last_spike_VA_matrix_shape(kk))> 0 %%%%%%%%%%%%%%%%%%%%%%%%%%delay ==0
                    if matrix_VA_to_rule(kk,jj)~=0
                        Isyn_VA_Matrix_to_PFC_exc_shape(jj,i+ delay)= Isyn_VA_Matrix_to_PFC_exc_shape(jj,i+ delay)+ W_VA_exi;% md_coefficient* spike_count(jj, i);
                    end
                end
            end
        end

        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_VA_matrix_Orientation(kk))<(spikewidth_MD) && (i-last_spike_VA_matrix_Orientation(kk))> 0 %%%%%%%%%%%%%%%%%%%%%%%%%%delay ==0
                    if matrix_VA_to_rule(kk,jj)~=0
                        Isyn_VA_Matrix_to_PFC_exc_orientation(jj,i+ delay)= Isyn_VA_Matrix_to_PFC_exc_orientation(jj,i+ delay)+ W_VA_exi;% md_coefficient* spike_count(jj, i);
                    end
                end
            end
        end


        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_VA_matrix_Orientation(kk))<(spikewidth_MD) && (i-last_spike_VA_matrix_Orientation(kk))> 0 
                    if matrix_MD_to_Crtex(kk,jj)~=0
                        Isyn_VA_Matrix_to_PFC(jj,i+ delay)= Matrix_amplification;% md_coefficient* spike_count(jj, i);
                    end
                end
            end
        end

        % VA to Inh 
        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_VA_matrix_shape(kk))<(spikewidth) && (i-last_spike_VA_matrix_shape(kk))> 0  %%%%%%%%%%%%%%%%%%%%%%%%%%delay ==0

                    Isyn_VA_Matrix_to_Inh(jj,i+ delay) = Isyn_VA_Matrix_to_Inh(jj,i+ delay)+ W_PFC_to_str;

                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MD to PFC 

        %SHAPE
        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_MD_shape(kk))<(spikewidth) && (i-last_spike_MD_shape(kk))> 0 
                    if matrix_local(kk,jj)~=0
                        Isyn_MD_shape_to_PFC(jj,i+ delay)= Isyn_MD_shape_to_PFC(jj,i+ delay)+ W_MD_PFC;% md_coefficient* spike_count(jj, i);
                    end
                end
            end
        end


        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PFC deep  to ST 
        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_PFC_D_BT(kk))<(spikewidth) && (i-last_spike_PFC_D_BT(kk))> 0  %%%%%%%%%%%%%%%%%%%%%%%%%%delay ==0

                    Isyn_aPFC_D5_to_ST(jj,i+ delay) = Isyn_aPFC_D5_to_ST(jj,i+ delay)+ W_PFC_to_str;

                end
            end
        end

        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_aPFC_D_YT(kk))<(spikewidth) && (i-last_spike_aPFC_D_YT(kk))> 0  %%%%%%%%%%%%%%%%%%%%%%%%%%delay ==0

                    Isyn_aPFC_D5_to_ST(jj,i+ delay) = Isyn_aPFC_D5_to_ST(jj,i+ delay)+ W_PFC_to_str;

                end
            end
        end

        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_aPFC_D_RC(kk))<(spikewidth) && (i-last_spike_aPFC_D_RC(kk))> 0  %%%%%%%%%%%%%%%%%%%%%%%%%%delay ==0

                    Isyn_aPFC_D5_to_ST(jj,i+ delay) = Isyn_aPFC_D5_to_ST(jj,i+ delay)+ W_PFC_to_str;

                end
            end
        end

        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_aPFC_D_GC(kk))<(spikewidth) && (i-last_spike_aPFC_D_GC(kk))> 0  %%%%%%%%%%%%%%%%%%%%%%%%%%delay ==0

                    Isyn_aPFC_D5_to_ST(jj,i+ delay) = Isyn_aPFC_D5_to_ST(jj,i+ delay)+ W_PFC_to_str;

                end
            end
        end
       
        %%%%%%%%%%%%%%%%%%%%% ST to SNpr
        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_ST(kk))<(spikewidth_inh) && (i-last_spike_ST(kk))> 0  % inserting delay && (i-last_spike_PFC(kk))>4
                    Isyn_ST(jj,i+delay)= Isyn_ST(jj,i+delay)+ W_IPL_Inh_to_exc*exp((1-(i-last_spike_ST(kk)-delay))*2./spikewidth_inh);%$$$$$$$$$$$$$&&&&&&&&&&&&&&&&&FIX the exponent
                end
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SNpr tO vA

        %         IPL Inh neurons to excitatory
        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_SNpr_Inh(kk))<(spikewidth_inh) && (i-last_spike_SNpr_Inh(kk))> 0  
                    Isyn_SNpr_to_VA(jj,i+delay)= Isyn_SNpr_to_VA(jj,i+delay)+ W_IPL_Inh_to_exc*exp((1-(i-last_spike_SNpr_Inh(kk)-delay))*2./spikewidth_inh);%$$$$$$$$$$$$$&&&&&&&&&&&&&&&&&FIX the exponent
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cortical inhibitory cells 
        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_PV(kk))<(spikewidth_inh) && (i-last_spike_PV(kk))> 0  
                    Isyn_PV(jj,i+delay)= Isyn_PV(jj,i+delay)+ W_IPL_Inh_to_exc_pv*exp((1-(i-last_spike_PV(kk)-delay))*2./spikewidth_inh);%$$$$$$$$$$$$$&&&&&&&&&&&&&&&&&FIX the exponent
                end
            end
        end
        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_MD_shape(kk))<(spikewidth_MD) && (i-last_spike_MD_shape(kk))> 0  
                    Isyn_MD_to_Inh(jj,i+delay)= Isyn_MD_to_Inh(jj,i+delay)+ W_MD_inh;%$$$$$$$$$$$$$&&&&&&&&&&&&&&&&&FIX the exponent
                end
            end
        end
        for jj=1:numberofneurons
            for kk=1:numberofneurons
                if (i-last_spike_FS(kk))<(spikewidth_inh) && (i-last_spike_FS(kk))> 0  
                    Isyn_FS(jj,i+delay)= Isyn_FS(jj,i+delay)+ W_IPL_Inh_to_exc_fs*exp((1-(i-last_spike_FS(kk)-delay))*2./spikewidth_inh);%$$$$$$$$$$$$$&&&&&&&&&&&&&&&&&FIX the exponent
                end
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% local excitatory input to
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% shape cells 

        for jj = 1:numberofneurons
            for kk = 1:numberofneurons
                if (i-last_spike_pPFC_remote_shape(kk))<(spikewidth) && (i-last_spike_pPFC_remote_shape(kk))> 0
                    if matrix_local(kk,jj) ~= 0
                        Isyn_aPFC_local_shape(jj,i+ delay) = Isyn_aPFC_local_shape(jj,i+ delay)+  matrix_local(kk,jj);
                    end
                end
            end
        end






        % for jj = 1:numberofneurons
        %     for kk = 1:numberofneurons
        %         if (i-last_spike_pPFC_remote_ORI(kk))<(spikewidth) && (i-last_spike_pPFC_remote_ORI(kk))> 0
        %             if matrix_local(kk,jj) ~= 0
        %                 Isyn_pPFC_local(jj,i+ delay) = Isyn_pPFC_local(jj,i+ delay)+ matrix_local(kk,jj);
        %             end
        %         end
        %     end
        % end

        % for jj = 1:numberofneurons
        %     for kk = 1:numberofneurons
        %         if (i-last_spike_pPFC_D(kk))<(spikewidth) && (i-last_spike_pPFC_D(kk))> 0
        %             if matrix_local_S_D(kk,jj) ~= 0
        %                 Isyn_pPFCD_to_response(jj,i+ delay) = Isyn_pPFCD_to_response(jj,i+ delay)+ matrix_local_S_D(kk,jj);
        %             end
        %         end
        %     end
        % end


        % for jj = 1:numberofneurons
        %     for kk = 1:numberofneurons
        %         if (i-last_spike_pPFC_remote_ORI_LEFT(kk))<(spikewidth) && (i-last_spike_pPFC_remote_ORI_LEFT(kk))> 0
        %             if matrix_local(kk,jj) ~= 0
        %                 Isyn_pPFC_local_left(jj,i+ delay) = Isyn_pPFC_local_left(jj,i+ delay)+ matrix_local(kk,jj);
        %             end
        %         end
        %     end
        % end
        % %



        %
        % for jj = 1:numberofneurons
        %     for kk = 1:numberofneurons
        %         if (i-last_spike_aPFC_D(kk))<(spikewidth) && (i-last_spike_aPFC_D(kk))> 0
        %             if matrix_local_D_toremote(kk,jj) ~= 0
        %                 Isyn_aPFC_D_to_Remote(jj,i+ delay) = Isyn_aPFC_D_to_Remote(jj,i+ delay)+ matrix_local_D_toremote(kk,jj);
        %             end
        %         end
        %     end
        % end



        % for jj=1:numberofneurons
        %     q=0;
        %     for kk=1:numberofneurons
        %         if (i-last_spike_MD2_core(kk))<(spikewidth_MD) && (i-last_spike_MD2_core(kk))> 0 %%%%%%%%%%%%%%%%%%%%%%%%%%delay ==0
        %             if matrix_MD_to_Crtex(kk,jj)~=0
        %                 q=q+1;
        %                 spike_count(jj, i)=q;
        %                 Isyn_MD2_to_pPFC(jj,i+ delay)= amplification;% md_coefficient* spike_count(jj, i);
        %             end
        %         end
        %     end
        % end


        % for jj=1:numberofneurons
        %     for kk=1:numberofneurons
        %         if (i-last_spike_PFC_D_BT(kk))<(spikewidth) && (i-last_spike_PFC_D_BT(kk))> 0  
        % 
        %             Isyn_aPFC_D5_toremote(jj,i+ delay)=Isyn_aPFC_D5_toremote(jj,i+ delay) + 0.3*W_PFC_TH;%%%%%%%%%%%%%%%%%%%%%                      hard coding
        % 
        %         end
        %     end
        % end
        
        %         %
        %current to PV cell
        % for jj=1:numberofneurons
        %     for kk=1:numberofneurons
        %         if (i-last_spike_MD2_core(kk))<(spikewidth_MD) && (i-last_spike_MD2_core(kk))> 0  % inserting delay && (i-last_spike_PFC(kk))>4
        %             Isyn_MD_to_inh_PV(jj,i+delay)= Isyn_MD_to_inh_PV(jj,i+delay)+ Matrix_MD_to_PV(kk, jj);
        %         end
        %     end
        % end
        %
       



     

        


        for j = 1:numberofneurons
            % create noise
            n1=rand;
            noise_amp = 10;
            n1=rand;
            noise_prob = 0.000005;
            noise_prob_PFC=0.000005;
            if n1<noise_prob
                I_noise_PFC_1(j,i)= noise_amp;
            else
                I_noise_PFC_1(j,i)=0;
            end





            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 1 PFC middle_ first  layer that gets the stim

            if (last_spike_PFC_M(j)~=10^10 && (i-last_spike_PFC_M(j))>tref)

                y_PFC_M(j,i) = y_PFC_M(j,i-1)+leaky_coef*((E_L-y_PFC_M(j,i-1))/tha)*dt+(Iext_PFC_M(j)+ I_stim(j,i)+ I_noise_PFC_1(j, i) )*dt*(RM/tha);

            elseif last_spike_PFC_M(j)==10^10

                y_PFC_M(j,i) = y_PFC_M(j,i-1)+leaky_coef*((E_L-y_PFC_M(j,i-1))/tha)*dt+(Iext_PFC_M(j)+ I_stim(j,i)+ I_noise_PFC_1(j, i) )*dt*(RM/tha);

            else

                y_PFC_M(j,i)=E_L;

            end

            if y_PFC_M(j,i)>=v_th
                last_spike_PFC_M(j)=i;
                y_PFC_M(j,i)=0;
                spiketimes_M=[spiketimes_M;i,j];

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Layer 2, where the layer 1 send signal for blue triangle to
            % this cell

            % membrane potential across PFC cells

            if (last_spike_PFC_S_BT(j)~=10^10 && (i-last_spike_PFC_S_BT(j))>tref)

                y_PFC_S_BT(j,i) = y_PFC_S_BT(j,i-1)+leaky_coef*((E_L-y_PFC_S_BT(j,i-1))/tha)*dt+(Iext_PFC(j)+ Isyn_PFC_M_S (j, i)+ (Isyn_aPFC_local_BT(j, i) + Isyn_aPFC_local_shared (j,i))* Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1(j, i)- Isyn_PV(j, i) )*dt*(RM/tha);%- Isyn_PV(j, i)
            elseif last_spike_PFC_S_BT(j)==10^10

                y_PFC_S_BT(j,i) = y_PFC_S_BT(j,i-1)+leaky_coef*((E_L-y_PFC_S_BT(j,i-1))/tha)*dt+(Iext_PFC(j)+ Isyn_PFC_M_S (j, i)+ (Isyn_aPFC_local_BT(j, i) + Isyn_aPFC_local_shared (j,i))* Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1(j, i)- Isyn_PV(j, i))*dt*(RM/tha);%- Isyn_PV(j, i)

            else

                y_PFC_S_BT(j,i)=E_L;

            end

            if y_PFC_S_BT(j,i)>=v_th
                last_spike_PFC_S_BT(j)=i;
                y_PFC_S_BT(j,i)=0;
                spiketimes_S=[spiketimes_S;i,j];

            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 2 PFC superficial  cells in superficial layers

            % membrane potential across PFC cells

            if (last_spike_PFC_S_RC(j)~=10^10 && (i-last_spike_PFC_S_RC(j))>tref)

                y_PFC_S_RC(j,i) = y_PFC_S_RC(j,i-1)+ leaky_coef*((E_L-y_PFC_S_RC(j,i-1))/tha)*dt + (Iext_PFC(j)+ (Isyn_aPFC_local_RC(j, i) + Isyn_aPFC_local_shared (j,i))* Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1(j, i)- Isyn_PV(j, i) )*dt*(RM/tha);%- Isyn_PV(j, i)
            elseif last_spike_PFC_S_RC(j)==10^10

                y_PFC_S_RC(j,i) = y_PFC_S_RC(j,i-1)+ leaky_coef*((E_L-y_PFC_S_RC(j,i-1))/tha)*dt + (Iext_PFC(j)+ (Isyn_aPFC_local_RC(j, i) + Isyn_aPFC_local_shared (j,i))* Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1(j, i)- Isyn_PV(j, i) )*dt*(RM/tha);%- Isyn_PV(j, i)

            else

                y_PFC_S_RC(j,i) = E_L;

            end

            if y_PFC_S_RC(j,i)>= v_th
                last_spike_PFC_S_RC(j)= i;
                y_PFC_S_RC(j,i)= 0;
                spiketimes_S_blue= [spiketimes_S_blue;i,j];

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            if (last_spike_PFC_S_GC(j)~=10^10 && (i-last_spike_PFC_S_GC(j))>tref)

                y_PFC_S_GC(j,i) = y_PFC_S_GC(j,i-1)+ leaky_coef*((E_L-y_PFC_S_GC(j,i-1))/tha)*dt + (Iext_PFC(j) + (Isyn_aPFC_local_GC(j, i) + Isyn_aPFC_local_shared (j,i)) * Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1(j, i)- Isyn_PV(j, i) )*dt*(RM/tha);%- Isyn_PV(j, i)
            elseif last_spike_PFC_S_GC(j)==10^10

                y_PFC_S_GC(j,i) = y_PFC_S_GC(j,i-1)+ leaky_coef*((E_L-y_PFC_S_GC(j,i-1))/tha)*dt + (Iext_PFC(j) + (Isyn_aPFC_local_GC(j, i) + Isyn_aPFC_local_shared (j,i)) * Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1(j, i)- Isyn_PV(j, i) )*dt*(RM/tha);%- Isyn_PV(j, i)

            else

                y_PFC_S_GC(j,i) = E_L;

            end

            if y_PFC_S_GC(j,i)>= v_th
                last_spike_PFC_S_GC(j)= i;
                y_PFC_S_GC(j,i)= 0;
                spiketimes_S_red= [spiketimes_S_red;i,j];

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            if (last_spike_PFC_S_YT(j)~=10^10 && (i-last_spike_PFC_S_YT(j))>tref)

                y_PFC_S_YT(j,i) = y_PFC_S_YT(j,i-1)+ leaky_coef*((E_L-y_PFC_S_YT(j,i-1))/tha)*dt + (Iext_PFC(j) + (Isyn_aPFC_local_YT(j, i) + Isyn_aPFC_local_shared (j,i)) * Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1(j, i) - Isyn_PV(j, i))*dt*(RM/tha);%- Isyn_PV(j, i)
            elseif last_spike_PFC_S_YT(j)==10^10

                y_PFC_S_YT(j,i) = y_PFC_S_YT(j,i-1)+ leaky_coef*((E_L-y_PFC_S_YT(j,i-1))/tha)*dt + (Iext_PFC(j) + (Isyn_aPFC_local_YT(j, i) + Isyn_aPFC_local_shared (j,i)) * Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1(j, i) - Isyn_PV(j, i))*dt*(RM/tha);%- Isyn_PV(j, i)

            else

                y_PFC_S_YT(j,i) = E_L;

            end

            if y_PFC_S_YT(j,i)>= v_th
                last_spike_PFC_S_YT(j)= i;
                y_PFC_S_YT(j,i)= 0;
                spiketimes_S_YT= [spiketimes_S_YT;i,j];

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 3 PFC Deep layer 5



            n1=rand;
            if n1<noise_prob_PFC
                I_noise_PFC_2(j,i)= noise_amp;
            else
                I_noise_PFC_2(j,i)= 0;
            end


            % membrane potential across deep layers of PFC cells

            if (last_spike_PFC_D_BT(j)~=10^10 && (i-last_spike_PFC_D_BT(j))>tref)

                y_aPFC_D_BT(j,i) = y_aPFC_D_BT(j,i-1)+leaky_coef*((E_L-y_aPFC_D_BT(j,i-1))/tha)*dt+(Iext_PFC_D(j)+ (W11 * Isyn_aPFC_to_D_BT(j,i)+ W12 * Isyn_aPFC_to_D_RC(j,i)+ W13 *Isyn_aPFC_to_D_GC(j,i)+ W14* Isyn_aPFC_to_D_YT(j,i)) + Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2(j, i))*dt*(RM/tha); %+Isyn_IPL_PFC(j,i)-Isyn_aPFC_Inh_to_exc(j,i)

            elseif last_spike_PFC_D_BT(j)==10^10

                y_aPFC_D_BT(j,i) = y_aPFC_D_BT(j,i-1)+leaky_coef*((E_L-y_aPFC_D_BT(j,i-1))/tha)*dt+(Iext_PFC_D(j)+ (W11 * Isyn_aPFC_to_D_BT(j,i)+ W12 * Isyn_aPFC_to_D_RC(j,i)+ W13 *Isyn_aPFC_to_D_GC(j,i)+ W14* Isyn_aPFC_to_D_YT(j,i)) + Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2(j, i))*dt*(RM/tha); %+Isyn_IPL_PFC(j,i)-Isyn_aPFC_Inh_to_exc(j,i)

            else

                y_aPFC_D_BT(j,i)= E_L;

            end

            if y_aPFC_D_BT(j,i)>= v_th
                last_spike_PFC_D_BT(j)= i;
                y_aPFC_D_BT(j,i)= 0;
                spiketimes_D_BT= [spiketimes_D_BT;i,j];

            end


            % membrane potential across deep layers of PFC cells

            if (last_spike_aPFC_D_RC(j)~=10^10 && (i-last_spike_aPFC_D_RC(j))>tref)

                y_aPFC_D_RC(j,i) = y_aPFC_D_RC(j,i-1)+leaky_coef*((E_L-y_aPFC_D_RC(j,i-1))/tha)*dt+(Iext_PFC_D(j)+ (W21 * Isyn_aPFC_to_D_BT(j,i)+ W22 * Isyn_aPFC_to_D_RC(j,i)+ W23 *Isyn_aPFC_to_D_GC(j,i)+ W24* Isyn_aPFC_to_D_YT(j,i)) +Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2(j, i))*dt*(RM/tha); %+Isyn_IPL_PFC(j,i)-Isyn_aPFC_Inh_to_exc(j,i)

            elseif last_spike_aPFC_D_RC(j)==10^10

                y_aPFC_D_RC(j,i) = y_aPFC_D_RC(j,i-1)+leaky_coef*((E_L-y_aPFC_D_RC(j,i-1))/tha)*dt+(Iext_PFC_D(j)+ (W21 * Isyn_aPFC_to_D_BT(j,i)+ W22 * Isyn_aPFC_to_D_RC(j,i)+ W23 *Isyn_aPFC_to_D_GC(j,i)+ W24* Isyn_aPFC_to_D_YT(j,i)) +Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2(j, i))*dt*(RM/tha); %+Isyn_IPL_PFC(j,i)-Isyn_aPFC_Inh_to_exc(j,i)

            else

                y_aPFC_D_RC(j,i)= E_L;

            end

            if y_aPFC_D_RC(j,i)>= v_th
                last_spike_aPFC_D_RC(j)= i;
                y_aPFC_D_RC(j,i)= 0;
                spiketimes_D_RC = [spiketimes_D_RC;i,j];

            end

            %%%%%%%%%%%% membrane potential across deep layers of PFC cells

            if (last_spike_aPFC_D_GC(j)~=10^10 && (i-last_spike_aPFC_D_GC(j))>tref)

                y_aPFC_D_GC(j,i) = y_aPFC_D_GC(j,i-1)+leaky_coef*((E_L-y_aPFC_D_GC(j,i-1))/tha)*dt+(Iext_PFC_D(j)+ (W31 * Isyn_aPFC_to_D_BT(j,i)+ W32 * Isyn_aPFC_to_D_RC(j,i)+ W33 *Isyn_aPFC_to_D_GC(j,i)+ W34* Isyn_aPFC_to_D_YT(j,i))+ Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2(j, i))*dt*(RM/tha); %+Isyn_IPL_PFC(j,i)-Isyn_aPFC_Inh_to_exc(j,i)

            elseif last_spike_aPFC_D_GC(j)==10^10

                y_aPFC_D_GC(j,i) = y_aPFC_D_GC(j,i-1)+leaky_coef*((E_L-y_aPFC_D_GC(j,i-1))/tha)*dt+(Iext_PFC_D(j)+ (W31 * Isyn_aPFC_to_D_BT(j,i)+ W32 * Isyn_aPFC_to_D_RC(j,i)+ W33 *Isyn_aPFC_to_D_GC(j,i)+ W34* Isyn_aPFC_to_D_YT(j,i)) + Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2(j, i))*dt*(RM/tha); %+Isyn_IPL_PFC(j,i)-Isyn_aPFC_Inh_to_exc(j,i)

            else

                y_aPFC_D_GC(j,i)= E_L;

            end

            if y_aPFC_D_GC(j,i)>= v_th
                last_spike_aPFC_D_GC(j)= i;
                y_aPFC_D_GC(j,i)= 0;
                spiketimes_D_GC = [spiketimes_D_GC;i,j];

            end

            %%%%%%%%%%%% membrane potential across deep layers of PFC cells

            if (last_spike_aPFC_D_YT(j)~=10^10 && (i-last_spike_aPFC_D_YT(j))>tref)

                y_aPFC_D_YT(j,i) = y_aPFC_D_YT(j,i-1)+leaky_coef*((E_L-y_aPFC_D_YT(j,i-1))/tha)*dt+(Iext_PFC_D(j)+ (W41 * Isyn_aPFC_to_D_BT(j,i)+ W42 * Isyn_aPFC_to_D_RC(j,i)+ W43 *Isyn_aPFC_to_D_GC(j,i)+ W44* Isyn_aPFC_to_D_YT(j,i))+ Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2(j, i))*dt*(RM/tha); %+Isyn_IPL_PFC(j,i)-Isyn_aPFC_Inh_to_exc(j,i)

            elseif last_spike_aPFC_D_YT(j)==10^10

                y_aPFC_D_YT(j,i) = y_aPFC_D_YT(j,i-1)+leaky_coef*((E_L-y_aPFC_D_YT(j,i-1))/tha)*dt+(Iext_PFC_D(j)+ (W41 * Isyn_aPFC_to_D_BT(j,i)+ W42 * Isyn_aPFC_to_D_RC(j,i)+ W43 *Isyn_aPFC_to_D_GC(j,i)+ W44* Isyn_aPFC_to_D_YT(j,i))+ Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2(j, i))*dt*(RM/tha); %+Isyn_IPL_PFC(j,i)-Isyn_aPFC_Inh_to_exc(j,i)

            else

                y_aPFC_D_YT(j,i)= E_L;

            end

            if y_aPFC_D_YT(j,i)>= v_th
                last_spike_aPFC_D_YT(j)= i;
                y_aPFC_D_YT(j,i)= 0;
                spiketimes_D_YT = [spiketimes_D_YT;i,j];

            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VA matrix cells 
            n1 = rand;
            noise_prob_md = noise_prob_PFC;
            if n1<noise_prob_md
                I_noise_VA(j,i) = noise_amp;
            else
                I_noise_VA(j,i) = 0;
            end

            %matrix VA cells
            if (last_spike_VA_matrix_shape(j)~=10^10 && (i-last_spike_VA_matrix_shape(j))>tref)

                y_VA_matrix_shape(j,i) = y_VA_matrix_shape(j,i-1)+((E_L-y_VA_matrix_shape(j,i-1))/tha)*dt+(Iext_VA(j)+ I_noise_VA(j,i)+ Isyn_aPFC_D5(j, i) - Isyn_SNpr_to_VA(j,i))*dt*(RM/tha);%

            elseif last_spike_VA_matrix_shape(j)==10^10

                y_VA_matrix_shape(j,i) = y_VA_matrix_shape(j,i-1)+((E_L-y_VA_matrix_shape(j,i-1))/tha)*dt+(Iext_VA(j)+ I_noise_VA(j,i)+ Isyn_aPFC_D5(j, i) - Isyn_SNpr_to_VA(j,i))*dt*(RM/tha);%

            else

                y_VA_matrix_shape(j,i)=E_L;

            end

            if y_VA_matrix_shape(j,i)>=v_th
                last_spike_VA_matrix_shape(j)=i;
                y_VA_matrix_shape(j,i)=0;
                spiketimes_VA=[spiketimes_VA;i,j];

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VA matrix cells connecting aPFC to pPFC
            n1 = rand;
            noise_prob_md = noise_prob_PFC;
            if n1<noise_prob_md
                I_noise_VA(j,i) = noise_amp;
            else
                I_noise_VA(j,i) = 0;
            end

            %matrix MD  cells
            if (last_spike_VA_matrix_Orientation(j)~=10^10 && (i-last_spike_VA_matrix_Orientation(j))>tref)

                y_VA_matrix_Orientation (j,i) = y_VA_matrix_Orientation (j,i-1)+((E_L-y_VA_matrix_Orientation (j,i-1))/tha)*dt+(Iext_VA(j)+ I_noise_VA(j,i)+ Isyn_aPFC_D5_orientation(j, i) - Isyn_SNpr_to_VA(j,i))*dt*(RM/tha);%

            elseif last_spike_VA_matrix_Orientation(j)==10^10

                y_VA_matrix_Orientation (j,i) = y_VA_matrix_Orientation (j,i-1)+((E_L-y_VA_matrix_Orientation (j,i-1))/tha)*dt+(Iext_VA(j)+ I_noise_VA(j,i)+ Isyn_aPFC_D5_orientation(j, i) - Isyn_SNpr_to_VA(j,i))*dt*(RM/tha);%

            else

                y_VA_matrix_Orientation (j,i)=E_L;

            end

            if y_VA_matrix_Orientation (j,i)>=v_th
                last_spike_VA_matrix_Orientation(j)=i;
                y_VA_matrix_Orientation (j,i)=0;
                spiketimes_VA_shape=[spiketimes_VA_shape;i,j];

            end

              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Core MD cells amplifying pPFC connections
            n1 = rand;
            noise_prob_md = noise_prob_PFC;
            if n1<noise_prob_md
                I_noise_md_2(j,i) = 0;%noise_amp;
            else
                I_noise_md_2(j,i) = 0;
            end


            if (last_spike_MD_shape(j)~=10^10 && (i-last_spike_MD_shape(j))>tref)

                y_MD_core_shape(j,i) = y_MD_core_shape(j,i-1)+((E_L-y_MD_core_shape(j,i-1))/tha)*dt+(Iext_VA(j)+ I_noise_md_2(j,i)+ Isyn_PFC_D_MD_shape(j,i) )*dt*(RM/tha);%

            elseif last_spike_MD_shape(j)==10^10

                y_MD_core_shape(j,i) = y_MD_core_shape(j,i-1)+((E_L-y_MD_core_shape(j,i-1))/tha)*dt+(Iext_VA(j)+ I_noise_md_2(j,i)+ Isyn_PFC_D_MD_shape(j,i) )*dt*(RM/tha);%

            else

                y_MD_core_shape(j,i) = E_L;

            end

            if y_MD_core_shape(j,i) >= v_th
                last_spike_MD_shape(j)=i;
                y_MD_core_shape(j,i) = 0;
                spiketimes_MD2=[spiketimes_MD2;i,j];

            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Core MD cells amplifying aPFC connections
            n1 = rand;
            noise_prob_md = noise_prob_PFC;
            if n1<noise_prob_md
                I_noise_md_1(j,i) =0;% noise_amp;
            else
                I_noise_md_1(j,i) = 0;
            end

            % md cells to IPL;        DEEP LAYER aPFC drives these

            if (last_spike_MD_ori(j)~=10^10 && (i-last_spike_MD_ori(j))>tref)

                y_MD_core_ori(j,i) = y_MD_core_ori(j,i-1)+((E_L-y_MD_core_ori(j,i-1))/tha)*dt+(Iext_VA(j)+I_noise_md_1(j,i)+ Isyn_PFC_D_MD_ori(j, i) )*dt*(RM/tha);%

            elseif last_spike_MD_ori(j)==10^10

                y_MD_core_ori(j,i) = y_MD_core_ori(j,i-1)+((E_L-y_MD_core_ori(j,i-1))/tha)*dt+(Iext_VA(j)+I_noise_md_1(j,i)+ Isyn_PFC_D_MD_ori(j, i) )*dt*(RM/tha);%

            else

                y_MD_core_ori(j,i)=E_L;

            end

            if y_MD_core_ori(j,i)>=v_th
                last_spike_MD_ori(j)=i;
                y_MD_core_ori(j,i)=0;
                spiketimes_MD=[spiketimes_MD;i,j];

            end






            %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % controls shape or orientation
            if (last_spike_ST(j)~=10^10 && (i-last_spike_ST(j))>tref)

                y_ST_inh(j,i) = y_ST_inh(j,i-1)+((E_L-y_ST_inh(j,i-1))/tha)*dt+(Iext_Inh_PV(j) + Isyn_aPFC_D5_to_ST(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            elseif last_spike_ST(j)==10^10

                y_ST_inh(j,i) = y_ST_inh(j,i-1)+((E_L-y_ST_inh(j,i-1))/tha)*dt+(Iext_Inh_PV(j) + Isyn_aPFC_D5_to_ST(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            else

                y_ST_inh(j,i)=E_L;

            end

            if y_ST_inh(j,i)>=v_th
                last_spike_ST(j)=i;
                y_ST_inh(j,i)=0;
                spiketimes_PV_inh=[spiketimes_PV_inh;i,j];


            end

            %%%%%%%%%%%%%%%%%%%%%%%
            %  %som inhibitory controls left or right
            if (last_spike_SNpr_Inh(j)~=10^10 && (i-last_spike_SNpr_Inh(j))>tref)

                y_SNpr_inh(j,i) = y_SNpr_inh(j,i-1)+((E_L-y_SNpr_inh(j,i-1))/tha)*dt+(Iext_Inh_pPFC(j) - Isyn_ST(j, i) )*dt*(RM/tha);%+ Isyn_aPFC_D5(j,i)+ Isyn_pPFC_local(j,i)

            elseif last_spike_SNpr_Inh(j)==10^10

                y_SNpr_inh(j,i) = y_SNpr_inh(j,i-1)+((E_L-y_SNpr_inh(j,i-1))/tha)*dt+(Iext_Inh_pPFC(j) - Isyn_ST(j, i) )*dt*(RM/tha);%+ Isyn_aPFC_D5(j,i)+ Isyn_pPFC_local(j,i)

            else

                y_SNpr_inh(j,i)=E_L;

            end

            if y_SNpr_inh(j,i)>=v_th
                last_spike_SNpr_Inh(j)=i;
                y_SNpr_inh(j,i)=0;
                spiketimes_ipl_inh=[spiketimes_ipl_inh;i,j];


            end
            %%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Shape cells in pPFC
            % Isyn_aPFC_D_to_Remote_effective_3(j,i) = Isyn_aPFC_D_to_Remote(j,i)*Isyn_MD_Matrix_to_PFC(j, i);
            % Isyn_pPFC_local_effective_3(j,i) = Isyn_pPFC_local(j,i)*Isyn_MD2_to_pPFC(j, i);
            n1=rand;
            if n1<noise_prob_PFC
                I_noise_PFC_6(j,i)= 0*noise_amp;
            else
                I_noise_PFC_6(j,i)= 0;
            end
            %             %
            if (last_spike_pPFC_remote_shape(j)~=10^10 && (i-last_spike_pPFC_remote_shape(j))>tref)

                y_pPFC_Shape(j,i) = y_pPFC_Shape(j,i-1)+leaky_coef*((E_L-y_pPFC_Shape(j,i-1))/tha)*dt+(Iext_PFC_remote_Shape(j)+  Isyn_aPFC_local_shared (j,i) * Isyn_VA_Matrix_to_PFC(j, i) + I_noise_PFC_6(j,i)+ Isyn_VA_Matrix_to_PFC_exc_shape(j, i)  + Isyn_MD_shape_to_PFC(j,i) + Isyn_aPFC_local_shape(j,i)* Isyn_VA_Matrix_to_PFC(j,i)  - Isyn_PV(j, i))*dt*(RM/tha);%

            elseif last_spike_pPFC_remote_shape(j)==10^10

                y_pPFC_Shape(j,i) = y_pPFC_Shape(j,i-1)+leaky_coef*((E_L-y_pPFC_Shape(j,i-1))/tha)*dt+(Iext_PFC_remote_Shape(j)+  Isyn_aPFC_local_shared (j,i) * Isyn_VA_Matrix_to_PFC(j, i) + I_noise_PFC_6(j,i)+ Isyn_VA_Matrix_to_PFC_exc_shape(j, i)  + Isyn_MD_shape_to_PFC(j,i) + Isyn_aPFC_local_shape(j,i)* Isyn_VA_Matrix_to_PFC(j,i)  - Isyn_PV(j, i))*dt*(RM/tha);%

            else

                y_pPFC_Shape(j,i)=E_L;

            end

            if y_pPFC_Shape(j,i)>=v_th
                last_spike_pPFC_remote_shape(j)=i;
                y_pPFC_Shape(j,i)=0;
                spiketimes_S_remote_Shape=[spiketimes_S_remote_Shape;i,j];

            end


            %%%%%%%%%%%%%%%%%%%%%%%%
            n1=rand;
            if n1<noise_prob_PFC
                I_noise_PFC_6(j,i)= noise_amp;
            else
                I_noise_PFC_6(j,i)= 0;
            end
            %             %
            if (last_spike_pPFC_remote_Orientation(j)~=10^10 && (i-last_spike_pPFC_remote_Orientation(j))>tref)

                y_pPFC_Orientation(j,i) = y_pPFC_Orientation(j,i-1)+leaky_coef*((E_L-y_pPFC_Orientation(j,i-1))/tha)*dt+(Iext_PFC_remote_Shape(j)+ I_noise_PFC_6(j,i) + Isyn_VA_Matrix_to_PFC_exc_orientation(j, i) - Isyn_PV(j, i) - Isyn_FS(j, i))*dt*(RM/tha);

            elseif last_spike_pPFC_remote_Orientation(j)==10^10

                y_pPFC_Orientation(j,i) = y_pPFC_Orientation(j,i-1)+leaky_coef*((E_L-y_pPFC_Orientation(j,i-1))/tha)*dt+(Iext_PFC_remote_Shape(j)+ I_noise_PFC_6(j,i) + Isyn_VA_Matrix_to_PFC_exc_orientation(j, i) - Isyn_PV(j, i) - Isyn_FS(j, i))*dt*(RM/tha);

            else

                y_pPFC_Orientation(j,i)=E_L;

            end

            if y_pPFC_Orientation(j,i)>=v_th
                last_spike_pPFC_remote_Orientation(j)=i;
                y_pPFC_Orientation(j,i)=0;
                spiketimes_S_remote_Orientation=[spiketimes_S_remote_Orientation;i,j];

            end





            % controls shape or orientation
            if (last_spike_PV(j)~=10^10 && (i-last_spike_PV(j))>tref)

                y_PV_inh(j,i) = y_PV_inh(j,i-1)+((E_L-y_PV_inh(j,i-1))/tha)*dt+(Iext_Inh_PV(j)+ Isyn_VA_Matrix_to_Inh(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            elseif last_spike_PV(j)==10^10

                y_PV_inh(j,i) = y_PV_inh(j,i-1)+((E_L-y_PV_inh(j,i-1))/tha)*dt+(Iext_Inh_PV(j)+ Isyn_VA_Matrix_to_Inh(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            else

                y_PV_inh(j,i)=E_L;

            end

            if y_PV_inh(j,i)>=v_th
                last_spike_PV(j)=i;
                y_PV_inh(j,i)=0;
                spiketimes_PV_inh=[spiketimes_PV_inh;i,j];
            end



            % controls shape or orientation
            if (last_spike_FS(j)~=10^10 && (i-last_spike_FS(j))>tref)

                y_FS_inh(j,i) = y_FS_inh(j,i-1)+((E_L-y_FS_inh(j,i-1))/tha)*dt+(Iext_Inh_PV(j)+ Isyn_MD_to_Inh(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            elseif last_spike_FS(j)==10^10

                y_FS_inh(j,i) = y_FS_inh(j,i-1)+((E_L-y_FS_inh(j,i-1))/tha)*dt+(Iext_Inh_PV(j)+ Isyn_MD_to_Inh(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            else

                y_FS_inh(j,i)=E_L;

            end

            if y_FS_inh(j,i)>=v_th
                last_spike_FS(j)=i;
                y_FS_inh(j,i)=0;
                spiketimes_FS_inh=[spiketimes_FS_inh;i,j];
            end




        end
    end

    % LFP_aPFC_S(rr,:)= smooth(mean(Isyn_aPFC_local_effective)+ mean(Isyn_pPFC_Remote_to_aPFC_effective)+ mean(I_noise_PFC_1));% %%% !!!!!!!!!!!!!!!! add noise
    % LFP_aPFC_D(rr,:)= smooth(mean(Isyn_aPFC_to_D_effective)+ mean(Isyn_aPFC_D_effective_2)+ mean(I_noise_PFC_2));% post-synaptic current in pfc
    % LFP_pPFC_S(rr,:)= smooth(mean(Isyn_pPFC_local_effective)+ mean(Isyn_pPFC_local_effective_2)+mean(Isyn_pPFC_local_effective_3)+ mean(Isyn_aPFC_D_to_Remote_effective)+ mean(Isyn_aPFC_D_to_Remote_effective_2)+ mean(Isyn_aPFC_D_to_Remote_effective_3)+ mean(I_noise_PFC_5)+ mean(I_noise_PFC_4)+ mean(I_noise_PFC_6));
    % LFP_pPFC_D(rr,:)= smooth(mean(Isyn_pPFC_S_to_D_effective)+ mean(Isyn_pPFC_D_effective)+ mean(I_noise_PFC_3));% post-synaptic current in pfc
    % LFP_MD_matrix(rr,:)= smooth(mean(I_noise_md_1)+ mean(I_noise_md_2)+ mean(I_noise_VA)+ mean(Isyn_PFC_D_MD_ori)+mean(Isyn_PFC_D_MD_shape)+ mean(Isyn_PFC_D_MD_ori));% post-synaptic current in pfc
    % LFP_response(rr,:)= smooth(mean(Isyn_pPFCD_to_response_effective)+ mean(I_noise_PFC));

end


%%
figure(1)

subplot(4,1,1)
% spy( y_PFC_M>-50,8,'b'),title('Middle layer sensory', 'FontSize', 16)
% set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
%
%
% subplot(5,1,2)
spy( y_PFC_S_BT>-50,8,'b'),title('PFC Superficial blue triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')


subplot(4,1,2)
spy( y_PFC_S_RC >-50,8,'b'),title('PFC Superficial red circle ', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

subplot(4,1,3)
spy( y_PFC_S_GC>-50,8,'b'),title('PFC Superficial green circle ', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

subplot(4,1,4)
spy(y_PFC_S_YT >-50,8,'b'),title('PFC Superficial yellow triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
%
figure (3)
subplot(4,1,1)
spy( y_aPFC_D_BT>-50,8,'k'),title('PFC Deep blue triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(4,1,2)
spy( y_aPFC_D_RC>-50,8,'k'),title('PFC Deep red circle', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

subplot(4,1,3)
spy( y_aPFC_D_GC>-50,8,'k'),title('PFC Deep green circle', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(4,1,4)
spy( y_aPFC_D_YT>-50,8,'k'),title('PFC Deep yellow triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')


figure( 2)
subplot(5,1,1)
spy( y_ST_inh>-50,8,'r'),title('Str', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(5,1,2)
spy( y_SNpr_inh>-50,8,'r'),title('SNpr', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(5,1,3)
spy( y_VA_matrix_shape>-50,8,'k'),title('VA Thalamus Shape', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(5,1,4)
spy( y_VA_matrix_Orientation >-50,8,'k'),title('VA Thalamus Orientation', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(5,1,5)
spy( y_PV_inh>-50,8,'r'),title('inhibitory PFC cells', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')


figure (4)
subplot(2,1,1)
spy( y_pPFC_Shape>-50,8,'b'),title('remote PFC Shape', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(2,1,2)
spy( y_pPFC_Orientation>-50,8,'b'),title('remote PFC Ori', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')


%
figure(5)

subplot(3,1,1)
spy( y_MD_core_shape>-50,8,'k'),title('MD', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(3,1,2)


spy( y_MD_core_ori>-50,8,'k'),title('MD', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')


subplot(3,1,3)
spy( y_FS_inh>-50,8,'r'),title('MD driven FS cells', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')


%%

subplot(8,1,3)
plot( y_aPFC_D(1, :),'b'),title('layer 5', 'FontSize', 16)
ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

subplot(8,1,4)
plot( y_aPFC_D_6(1, :),'b'),title('layer 6', 'FontSize', 16)
ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

subplot(8,1,5)
plot( y_ST_inh(1, :),'r'),title('Str', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(8,1,6)
plot( y_SNpr_inh(1, :),'r'),title('SNpr', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(8,1,7)
plot( y_VA_matrix_shape(1, :),'k'),title('VA Thalamus', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

subplot(8,1,8)
plot( y_VA_matrix_Orientation (1, :),'k'),title('VA shape Thalamus', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

%% ploting spikes

figure(5)
hold on
srate=1000;
gauss_width=100;
example_neuron=[1 3 5 10 20 ]
hold on
ax1=subplot(2,1,1),
box(ax1,'on');
no_of_trials=7
spiketimes_VA;
spiketimes_S;
spiketimes_D_YT;
spiketimes_S_remote_ORI;
spiketimes_MD2;
spiketimes_response;
spike_times = spiketimes_D_YT;
spiketimes_pPFC_D;
for ii = 1: length(example_neuron)
    rr= example_neuron(ii);
    trl= floor(spike_times(find(spike_times(:,2) ==rr))./100000);
    y = mod (spike_times(find(spike_times(:,2) ==rr), 1),100000 );
    spikes= [trl,y];
    for jj=1:max(trl)+1
        temp = find(spikes(:,1)==jj);
        if isempty(temp) == 0
            x=[];
            x = spikes( temp, 2);
            nspikes = numel(x);
            for k = 1:nspikes
                ax1=subplot(1,5,ii),

                line([x(k) x(k)], [jj-0.5 jj+0.5],'Color','k');  xlim([1 100000]),ylim([0.5  no_of_trials+0.5 ]),ylabel(' # of trials','FontSize',13),title('D PFC ','FontSize',13);

            end
        end

    end
end

%%
a=axes;
set(a,'Units','normalized')
P=get(a,'Position');
set(a,'Position',P+[0 -0.03 0 0])
new_timeseries=zeros(numberofneurons, 100000,no_of_trials );
for hhh=1:numberofneurons
    y=[];
    spikes=[];
    trl= floor(spike_times(find(spike_times(:,2) == hhh), 1)./100000)
    y = mod (spike_times(find(spike_times(:,2) == hhh), 1),100000 )
    spikes= [trl,y];
    for kk=1:max(trl)+1
        temp = find(spikes(:,1)== kk);
        if isempty(temp) == 0
            x =spikes(find(spikes(:,1)== kk),2)';
            new_timeseries(hhh,x+1,kk)=1;
        end
    end
end

srate=1000;
gauss_width= 5000;

for j= 1 :  length(example_neuron)
    temp=[];xxx = [];
    rr= example_neuron(j);

    temp= squeeze(new_timeseries(rr, :, :));
    temp = temp';
    xxx= mean(temp);
    subplot(2,5,j+5)
    % hold on
    %
    % plot(xxx)

    new_matrix_pfc = conv_gaussian(xxx, srate,gauss_width);

    plot((new_matrix_pfc), 'k' ), xlim([1 100000]); ylim([0 1]),xlabel('time','FontSize',13),ylabel('Cumulative activity','FontSize',13)

    new_matrix_pfc=[];

end



%%


