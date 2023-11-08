
close all;
clear all;
% rng('default')
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

total_trial_num = 1;
gauss_width= 100;
 
%   Simulation time
dt = 0.01; %step size ms
t_final = 4000; %simulation time ms
T = 0:dt:t_final;

% legion test
VA_off = 0;
MD_off = 0;
pPFC_off = 0;
PV_off = 0;

%   Intrinsic property of neuron
delay = 5/dt; % 3ms
MD_delay = 10/dt; % 3ms
v_th = -50; %mv
E_L = -65; %mv
RM = 10;
leaky_coef = 1; %
tref = 1/dt; % 1 ms refratory time

% placeholder for multi trials results (first neuron of the section are saved)
% full_PFC_S_BT = []; 
% full_PFC_S_RC = [];
% full_PFC_S_GC = []; 
% full_PFC_S_YT = [];
% 
% full_PFC_D_BT = []; 
% full_PFC_D_RC = [];
% full_PFC_D_GC = []; 
% full_PFC_D_YT = []; 
% 
% full_VA_shape = [];
% full_VA_ori = []; 
% 
% full_MD_shape = [];
% full_MD_ori = []; 
% 
% full_PFC_remote_shape = [];
% full_PFC_remote_ori = [];

full_PFC_S_BT = {}; 
full_PFC_S_RC = {};
full_PFC_S_GC = {}; 
full_PFC_S_YT = {};

full_PFC_D_BT = {}; 
full_PFC_D_RC = {};
full_PFC_D_GC = {}; 
full_PFC_D_YT = {}; 

full_VA_shape = {};
full_VA_ori = {}; 

full_MD_shape = {};
full_MD_ori = {}; 

full_PFC_remote_shape = {};
full_PFC_remote_ori = {};


% save the total firing number for VA, MD and PFC remote shape to compare
% the order
full_VA_shape_num =  zeros(total_trial_num, length(T));
full_MD_shape_num = zeros(total_trial_num, length(T));
full_PFC_remote_shape_num = zeros(total_trial_num, length(T));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neurons talk to each other by sending short pulses, define duration and
% amplitude of these pulses
abs_pfc_width = 1.5;% cortical input liftime
spikewidth = abs_pfc_width/dt;
spikewidth_inh = 100/dt; % assumed that inhibition effect is prolonged due to the synch exc input SOM neurons
spikewidth_VA = 100/dt;
spikewidth_MD = spikewidth;
%   Connectivity Map: local PFC and PRL connectivity = all to all connectivity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Stimulation showing the abstract cue
no_of_trl = 1;
column_length = 10;
I_stim = zeros(numberofneurons,length(T));
t_start_stim_1_abs = 2000; % cue time
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
matrix_local = zeros(numberofneurons,numberofneurons);
matrix_M= zeros(numberofneurons,numberofneurons);
for ii = 1:column_length: numberofneurons-column_length

    matrix_local(ii:ii+ column_length - 1,ii + column_length:ii+ column_length-1) = W_local;
    matrix_M(ii:ii+ column_length - 1, ii :ii+ column_length-1) = W_M;
end
matrix_local = matrix_local - diag(diag(matrix_local));
matrix_local(1:10, 11:20)= 2* matrix_local(1:10, 11:20);
matrix_M(41:50,1:10 )= W_M;
% matrix_M(1:10, 41:50 )= W_M;
Matrix_PFC_to_MD = zeros(numberofneurons,numberofneurons);%matrix_M;
Matrix_PFC_to_MD (:, :)=1;%matrix_M;

Matrix_DPFC_to_MD = zeros(numberofneurons,numberofneurons);
% for i=1:20
%     random_commections= randi(numberofneurons, 1, 1);
% Matrix_DPFC_to_MD (random_commections,:)=1;
% end 
Matrix_DPFC_to_MD(:, 1:15)=1 ; 
% Matrix_PFC_to_MD(1:25, 26:50)=1;
% Matrix_PFC_to_MD(26:50, 1:25)=1;
Matrix_MD_to_PFC= zeros(numberofneurons,numberofneurons);%matrix_M;

% random_commections= randi(numberofneurons, 1, 20);
Matrix_MD_to_PFC(:,1:25)=1;
% Matrix_MD_to_PFC(26:50, 1:25)=1;

Matrix_random_MD_PFC_BT= zeros(numberofneurons,numberofneurons);
Matrix_random_MD_PFC_RC= zeros(numberofneurons,numberofneurons);
Matrix_random_MD_PFC_GC= zeros(numberofneurons,numberofneurons);
Matrix_random_MD_PFC_YT= zeros(numberofneurons,numberofneurons);

for i=1:5
random_commections= randi(numberofneurons, 1, 1);
Matrix_random_MD_PFC_BT (random_commections,random_commections)=1;
random_commections= randi(numberofneurons, 1, 1);
Matrix_random_MD_PFC_RC (random_commections,random_commections)=1;
random_commections= randi(numberofneurons, 1, 1);
Matrix_random_MD_PFC_GC (random_commections,random_commections)=1;
random_commections= randi(numberofneurons, 1, 1);
Matrix_random_MD_PFC_YT (random_commections,random_commections)=1;
end 
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
W_PFC_TH =  0.015;% 0.0175;% 0.05;
W_PFC_MD= 0.0006;
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
matrix_VA_to_rule(1:50, 1:20)= 1;
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


% driving currents-----------------------------------------------------
Iext_SPFC_YT = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
Iext_SPFC_BT = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
Iext_SPFC_RC = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
Iext_SPFC_GC = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
Iext_DPFC_YT = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
Iext_DPFC_BT = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
Iext_DPFC_RC = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
Iext_DPFC_GC = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
Iext_PFC_M = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5);
Iext_PFC_remote_ORI = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
Iext_PFC_remote_Shape = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
Iext_PFC_D = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
Iext_PFC_D_2 = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
Iext_PFC_D_3 = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
Iext_MD = 1.5 - 0.001 *rand(numberofneurons,1); % same driving current for all thalamic cells
Iext_VA = 1.5 - 0.001 *rand(numberofneurons,1); % same driving current for all thalamic cells
Iext_Inh = 1.5 -  0.01 *(rand(numberofneurons,1)-0.5); % driving current for inhibitory cells similar in PFC and IPL
Iext_Inh_pPFC = 1.5 - 0.01 *(rand(numberofneurons,1)-0.5); %driving current for inhibitory cells similar in PFC and IPL
Iext_Inh_PV = 1.49 - 0.0 *(rand(numberofneurons,1));


for rr= 1:total_trial_num
    tic

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
    y_response_Ori_right(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1;
    y_response_Ori_left = zeros(numberofneurons,length(T));
    y_response_Ori_left(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1;
    y_pPFC_ORI = zeros(numberofneurons,length(T));
    y_pPFC_ORI(1:numberofneurons,1)= -55+ rand(numberofneurons,1)*5; %initializing the excitatory IPL neurons Voltages
    y_pPFC_ORI_LEFT = zeros(numberofneurons,length(T));
    y_pPFC_ORI_LEFT(1:numberofneurons,1)= -55 + rand(numberofneurons,1)*5; %initializing the excitatory IPL neurons Voltages
    y_pPFC_Shape = zeros(numberofneurons,length(T));
    y_pPFC_Shape(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;

    y_pPFC_Orientation = zeros(numberofneurons,length(T));
    y_pPFC_Orientation(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;

    y_VA_matrix_shape = zeros(numberofneurons,length(T));
    y_VA_matrix_shape(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1; %initializing the thalamic neurons Voltages

    y_VA_matrix_Orientation  = zeros(numberofneurons,length(T));
    y_VA_matrix_Orientation  =-55+rand(numberofneurons,1)*1; %initializing the thalamic neurons Voltages


    y_MD_core_ori = zeros(numberofneurons,length(T));
    y_MD_core_ori(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1; %initializing the thalamic neurons Voltages
    y_MD_core_shape = zeros(numberofneurons,length(T));
    y_MD_core_shape(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1; %initializing the thalamic neurons Voltages
    y_SNpr_inh=zeros(numberofneurons,length(T));
    y_SNpr_inh(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;% initializing the PFC neurons Voltages
    y_ST_inh=zeros(numberofneurons,length(T));
    y_ST_inh(1:numberofneurons,1)=-55+rand(numberofneurons,1)*0;% initializing the PFC neurons Voltages
    y_PV_inh= zeros(numberofneurons,length(T));
    y_PV_inh (1:numberofneurons,1)=-55+rand(numberofneurons,1)*1;
    y_FS_inh= zeros(numberofneurons,length(T));
    y_FS_inh(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1;


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
    last_spike_PFC_D_YT= 10^10*ones(numberofneurons,1);
    last_spike_PFC_D_GC= 10^10*ones(numberofneurons,1);
    last_spike_PFC_D_RC = 10^10*ones(numberofneurons,1);
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
    last_spike_PFC_remote_shape= 10^10*ones(numberofneurons,1);
    last_spike_PFC_remote_Orientation= 10^10*ones(numberofneurons,1);
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
    Isyn_MD_ori_to_PFC= zeros(numberofneurons,length(T));
    Isyn_MD_shape_to_random_PFC=zeros(numberofneurons,length(T));
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
    Isyn_aPFC_local_ori= zeros(numberofneurons,length(T));
    Isyn_MD_shape_to_random_PFC_BT= zeros(numberofneurons,length(T));
    Isyn_MD_shape_to_random_PFC_RC= zeros(numberofneurons,length(T));
    Isyn_MD_shape_to_random_PFC_GC= zeros(numberofneurons,length(T));
    Isyn_MD_shape_to_random_PFC_YT= zeros(numberofneurons,length(T));
    Isyn_MD_orientation_to_random_PFC_BT= zeros(numberofneurons,length(T));
    Isyn_MD_orientation_to_random_PFC_RC= zeros(numberofneurons,length(T));
    Isyn_MD_orientation_to_random_PFC_GC= zeros(numberofneurons,length(T));
    Isyn_MD_orientation_to_random_PFC_YT= zeros(numberofneurons,length(T));
    % amplification = 1;
    Matrix_amplification = 2;


    % simulation starts here
    for i= 2:length(T)- MD_delay*10


       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Midle layer (l = 1) to superficial layer (l = 2)

        Isyn_PFC_M_S(:,i+ delay) = Isyn_PFC_M_S(:,i+ delay) + sum(((i-last_spike_PFC_M < spikewidth) & (i-last_spike_PFC_M > 0)) .* matrix_M, 1).';

        %%%%%%%%% Superficial layer (l = 2) to superficial layer (l = 2)
        %%%%%%%%% excitatory to excitatory  within chains

        w_within_chains = 0.6;
        Isyn_aPFC_local_BT(:,i+ delay) = Isyn_aPFC_local_BT(:,i+ delay) + w_within_chains .* sum(((i-last_spike_PFC_S_BT<spikewidth) & (i-last_spike_PFC_S_BT>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_RC(:,i+ delay) = Isyn_aPFC_local_RC(:,i+ delay) + w_within_chains .* sum(((i-last_spike_PFC_S_RC<spikewidth) & (i-last_spike_PFC_S_RC>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_YT(:,i+ delay) = Isyn_aPFC_local_YT(:,i+ delay) + w_within_chains .* sum(((i-last_spike_PFC_S_YT<spikewidth) & (i-last_spike_PFC_S_YT>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_GC(:,i+ delay) = Isyn_aPFC_local_GC(:,i+ delay) + w_within_chains .* sum(((i-last_spike_PFC_S_GC<spikewidth) & (i-last_spike_PFC_S_GC>0)) .* matrix_local, 1).';

        %%%%%%%%%%%Superficial layer (l = 2) to superficial layer (l = 2)
        %%%%%%%%%%% excitatory to excitatory across chains

        w_across_chains= 0.45;
        
        Isyn_aPFC_local_shared(:,i+ delay) = Isyn_aPFC_local_shared(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_S_BT<spikewidth) & (i-last_spike_PFC_S_BT>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_shared(:,i+ delay) = Isyn_aPFC_local_shared(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_S_RC<spikewidth) & (i-last_spike_PFC_S_RC>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_shared(:,i+ delay) = Isyn_aPFC_local_shared(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_S_YT<spikewidth) & (i-last_spike_PFC_S_YT>0)) .* matrix_local, 1).';
        
        Isyn_aPFC_local_shared(:,i+ delay) = Isyn_aPFC_local_shared(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_S_GC<spikewidth) & (i-last_spike_PFC_S_GC>0)) .* matrix_local, 1).';

        %%%%%%%%%%%Superficial layer (l = 2) to deep layer (l = 3)

        Isyn_aPFC_to_D_BT(:,i+ delay) = Isyn_aPFC_to_D_BT(:,i+ delay) + sum(((i-last_spike_PFC_S_BT<spikewidth) & (i-last_spike_PFC_S_BT>0)) .* matrix_local_S_D, 1).';
        
        Isyn_aPFC_to_D_RC(:,i+ delay) = Isyn_aPFC_to_D_RC(:,i+ delay) + sum(((i-last_spike_PFC_S_RC<spikewidth) & (i-last_spike_PFC_S_RC>0)) .* matrix_local_S_D, 1).';

        Isyn_aPFC_to_D_GC(:,i+ delay) = Isyn_aPFC_to_D_GC(:,i+ delay) + sum(((i-last_spike_PFC_S_GC<spikewidth) & (i-last_spike_PFC_S_GC>0)) .* matrix_local_S_D, 1).';
        
        Isyn_aPFC_to_D_YT(:,i+ delay) = Isyn_aPFC_to_D_YT(:,i+ delay) + sum(((i-last_spike_PFC_S_YT<spikewidth) & (i-last_spike_PFC_S_YT>0)) .* matrix_local_S_D, 1).';

        %%%%%%%%%%%%%%%%%%%%%%%%%% Deep layer (l = 3) to deep layer (l = 3)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% excitatory to excitatory across chains

        Isyn_aPFC_D_shared(:,i+ delay) = Isyn_aPFC_D_shared(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_D_BT<spikewidth) & (i-last_spike_PFC_D_BT>0)) .* matrix_local, 1).';

        Isyn_aPFC_D_shared(:,i+ delay) = Isyn_aPFC_D_shared(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_D_GC<spikewidth) & (i-last_spike_PFC_D_GC>0)) .* matrix_local, 1).';

        Isyn_aPFC_D_shared(:,i+ delay) = Isyn_aPFC_D_shared(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_D_RC<spikewidth) & (i-last_spike_PFC_D_RC>0)) .* matrix_local, 1).';

        Isyn_aPFC_D_shared(:,i+ delay) = Isyn_aPFC_D_shared(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_D_YT<spikewidth) & (i-last_spike_PFC_D_YT>0)) .* matrix_local, 1).';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Deep layer (l=3) to VA thalamus (l=4) 
        %%%%  SHAPE

        Isyn_aPFC_D5(:,i+ delay) = Isyn_aPFC_D5(:,i+ delay) +  W_PFC_TH .* sum(((i-last_spike_PFC_D_BT<spikewidth) & (i-last_spike_PFC_D_BT>0) & (PFC_VA_matrix ~= 0)), 1).';
        
        Isyn_aPFC_D5(:,i+ delay) = Isyn_aPFC_D5(:,i+ delay) +  W_PFC_TH .* sum(((i-last_spike_PFC_D_RC<spikewidth) & (i-last_spike_PFC_D_RC>0) & (PFC_VA_matrix ~= 0)), 1).';

       %%%%  ORIENTATION
        
        Isyn_aPFC_D5_orientation(:,i+ delay) = Isyn_aPFC_D5_orientation(:,i+ delay) +  W_PFC_TH .* sum((i-last_spike_PFC_D_GC<spikewidth) & (i-last_spike_PFC_D_GC>0) & (PFC_VA_matrix ~= 0), 1).';

        Isyn_aPFC_D5_orientation(:,i+ delay) = Isyn_aPFC_D5_orientation(:,i+ delay) +  W_PFC_TH .* sum((i-last_spike_PFC_D_YT<spikewidth) & (i-last_spike_PFC_D_YT>0) & (PFC_VA_matrix ~= 0), 1).';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PFC deep  to MD 
% SHAPE 
        PFC_MD_initial_coef= 5;
         PFC_MD_secondary_coef= 3;
        Isyn_PFC_D_MD_shape(:,i+ MD_delay) = Isyn_PFC_D_MD_shape(:,i+ MD_delay) + PFC_MD_initial_coef * W_PFC_MD .* sum(((i-last_spike_PFC_D_BT<spikewidth) & (i-last_spike_PFC_D_BT>0) & (Matrix_DPFC_to_MD ~= 0)), 1).';
        
        Isyn_PFC_D_MD_shape(:,i+ MD_delay) = Isyn_PFC_D_MD_shape(:,i+ MD_delay) + PFC_MD_initial_coef * W_PFC_MD .* sum(((i-last_spike_PFC_D_RC<spikewidth) & (i-last_spike_PFC_D_RC>0) & (Matrix_DPFC_to_MD ~= 0)), 1).';

        Isyn_PFC_D_MD_shape(:,i+ 10*MD_delay) = Isyn_PFC_D_MD_shape(:,i+ 10*MD_delay) + PFC_MD_secondary_coef * W_PFC_MD .* sum(((i-last_spike_PFC_remote_shape<spikewidth) & (i-last_spike_PFC_remote_shape>0) & (Matrix_PFC_to_MD ~= 0)), 1).';

% ORIENTATION 

        Isyn_PFC_D_MD_ori(:,i+ MD_delay) = Isyn_PFC_D_MD_ori(:,i+ MD_delay) + PFC_MD_initial_coef * W_PFC_MD .* sum(((i-last_spike_PFC_D_YT<spikewidth) & (i-last_spike_PFC_D_YT>0) & (Matrix_DPFC_to_MD ~= 0)), 1).';

        Isyn_PFC_D_MD_ori(:,i+ MD_delay) = Isyn_PFC_D_MD_ori(:,i+ MD_delay) + PFC_MD_initial_coef * W_PFC_MD .* sum(((i-last_spike_PFC_D_GC<spikewidth) & (i-last_spike_PFC_D_GC>0) & (Matrix_DPFC_to_MD ~= 0)), 1).';

        Isyn_PFC_D_MD_ori(:,i+ 10*MD_delay) = Isyn_PFC_D_MD_ori(:,i+ 10*MD_delay) + PFC_MD_secondary_coef *W_PFC_MD .* sum(((i-last_spike_PFC_remote_Orientation<spikewidth) & (i-last_spike_PFC_remote_Orientation>0) & (Matrix_PFC_to_MD ~= 0)), 1).';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VA thalamus to cortex 
        if ~VA_off
    
            Isyn_VA_Matrix_to_PFC(sum(((i-last_spike_VA_matrix_shape<spikewidth_VA) & (i-last_spike_VA_matrix_shape>0)) .* matrix_MD_to_Crtex, 1)~=0,i+ delay) = Matrix_amplification;
    
            Isyn_VA_Matrix_to_PFC_exc_shape(:,i+ delay) = Isyn_VA_Matrix_to_PFC_exc_shape(:,i+ delay) +  2*W_VA_exi .* sum(((i-last_spike_VA_matrix_shape<spikewidth) & (i-last_spike_VA_matrix_shape>0) & (matrix_VA_to_rule ~= 0)), 1).';
            
            Isyn_VA_Matrix_to_PFC_exc_orientation(:,i+ delay) = Isyn_VA_Matrix_to_PFC_exc_orientation(:,i+ delay) +  W_VA_exi .* sum(((i-last_spike_VA_matrix_Orientation<spikewidth) & (i-last_spike_VA_matrix_Orientation>0) & (matrix_VA_to_rule ~= 0)), 1).';
    
            Isyn_VA_Matrix_to_PFC(sum(((i-last_spike_VA_matrix_Orientation<spikewidth_VA) & (i-last_spike_VA_matrix_Orientation>0)) .* matrix_MD_to_Crtex, 1)~=0,i+ delay) = Matrix_amplification;
    
            % VA to Inh 
    
            Isyn_VA_Matrix_to_Inh(:,i+ delay) = Isyn_VA_Matrix_to_Inh(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_VA_matrix_shape<spikewidth) & (i-last_spike_VA_matrix_shape>0), 1).';
            Isyn_VA_Matrix_to_Inh(:,i+ delay) = Isyn_VA_Matrix_to_Inh(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_VA_matrix_Orientation<spikewidth) & (i-last_spike_VA_matrix_Orientation>0), 1).';

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MD to PFC 

        %SHAPE

        if ~MD_off

            Isyn_MD_shape_to_PFC(:,i+ 8*MD_delay) = Isyn_MD_shape_to_PFC(:,i+ 8*MD_delay) +  0.3*W_MD_PFC .* sum(((i-last_spike_MD_shape<spikewidth_MD) & (i-last_spike_MD_shape>0) & (Matrix_MD_to_PFC ~= 0)), 1).'; 
    
            Isyn_MD_shape_to_random_PFC_BT(:,i+ delay) = Isyn_MD_shape_to_random_PFC_BT(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_shape<spikewidth_MD) & (i-last_spike_MD_shape>0) & (Matrix_random_MD_PFC_BT ~= 0)), 1).'; 
    
            Isyn_MD_shape_to_random_PFC_RC(:,i+ delay) = Isyn_MD_shape_to_random_PFC_RC(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_shape<spikewidth_MD) & (i-last_spike_MD_shape>0) & (Matrix_random_MD_PFC_RC ~= 0)), 1).'; 
    
            Isyn_MD_shape_to_random_PFC_GC(:,i+ delay) = Isyn_MD_shape_to_random_PFC_GC(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_shape<spikewidth_MD) & (i-last_spike_MD_shape>0) & (Matrix_random_MD_PFC_GC ~= 0)), 1).'; 
    
            Isyn_MD_shape_to_random_PFC_YT(:,i+ delay) = Isyn_MD_shape_to_random_PFC_YT(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_shape<spikewidth_MD) & (i-last_spike_MD_shape>0) & (Matrix_random_MD_PFC_YT ~= 0)), 1).'; 

            Isyn_MD_ori_to_PFC(:,i+ 8*MD_delay) = Isyn_MD_ori_to_PFC(:,i+ 8*MD_delay) +  0.3*W_MD_PFC .* sum(((i-last_spike_MD_ori<spikewidth_MD) & (i-last_spike_MD_ori>0) & (Matrix_MD_to_PFC ~= 0)), 1).'; 
            
            Isyn_MD_orientation_to_random_PFC_BT(:,i+ delay) = Isyn_MD_orientation_to_random_PFC_BT(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_ori<spikewidth_MD) & (i-last_spike_MD_ori>0) & (Matrix_random_MD_PFC_BT ~= 0)), 1).'; 
    
            Isyn_MD_orientation_to_random_PFC_RC(:,i+ delay) = Isyn_MD_orientation_to_random_PFC_RC(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_ori<spikewidth_MD) & (i-last_spike_MD_ori>0) & (Matrix_random_MD_PFC_RC ~= 0)), 1).'; 
    
            Isyn_MD_orientation_to_random_PFC_GC(:,i+ delay) = Isyn_MD_orientation_to_random_PFC_GC(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_ori<spikewidth_MD) & (i-last_spike_MD_ori>0) & (Matrix_random_MD_PFC_GC ~= 0)), 1).'; 
    
            Isyn_MD_orientation_to_random_PFC_YT(:,i+ delay) = Isyn_MD_orientation_to_random_PFC_YT(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_ori<spikewidth_MD) & (i-last_spike_MD_ori>0) & (Matrix_random_MD_PFC_YT ~= 0)), 1).'; 


        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PFC deep  to ST 

        Isyn_aPFC_D5_to_ST(:,i+ delay) = Isyn_aPFC_D5_to_ST(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_PFC_D_BT<spikewidth) & (i-last_spike_PFC_D_BT>0), 1).';

        Isyn_aPFC_D5_to_ST(:,i+ delay) = Isyn_aPFC_D5_to_ST(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_PFC_D_YT<spikewidth) & (i-last_spike_PFC_D_YT>0), 1).';

        Isyn_aPFC_D5_to_ST(:,i+ delay) = Isyn_aPFC_D5_to_ST(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_PFC_D_RC<spikewidth) & (i-last_spike_PFC_D_RC>0), 1).';

        Isyn_aPFC_D5_to_ST(:,i+ delay) = Isyn_aPFC_D5_to_ST(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_PFC_D_GC<spikewidth) & (i-last_spike_PFC_D_GC>0), 1).';
       
        %%%%%%%%%%%%%%%%%%%%% ST to SNpr
        Isyn_ST(:, i + delay) =  Isyn_ST(:, i + delay) + W_IPL_Inh_to_exc.*sum(((i-last_spike_ST < spikewidth_inh) & (i-last_spike_ST > 0)).*(exp(((i-last_spike_ST < spikewidth_inh) & (i-last_spike_ST > 0)).*(1-(i-last_spike_ST-delay)).*2./spikewidth_inh)), 1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SNpr tO vA

        %         IPL Inh neurons to excitatory
        Isyn_SNpr_to_VA(:, i + delay) =  Isyn_SNpr_to_VA(:, i + delay) + W_IPL_Inh_to_exc.*sum(((i-last_spike_SNpr_Inh < spikewidth_inh) & (i-last_spike_SNpr_Inh > 0)).*(exp(((i-last_spike_SNpr_Inh < spikewidth_inh) & (i-last_spike_SNpr_Inh > 0)).*(1-(i-last_spike_SNpr_Inh-delay)).*2./spikewidth_inh)), 1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cortical inhibitory cells 

        if ~PV_off
            Isyn_PV(:, i + delay) =  Isyn_PV(:, i + delay) + W_IPL_Inh_to_exc_pv.*sum(((i-last_spike_PV < spikewidth_inh) & (i-last_spike_PV > 0)).*(exp(((i-last_spike_PV < spikewidth_inh) & (i-last_spike_PV > 0)).*(1-(i-last_spike_PV-delay)).*2./spikewidth_inh)), 1);
        end

        if ~MD_off
            Isyn_MD_to_Inh(:,i+ delay) = Isyn_MD_to_Inh(:,i+ delay) + W_MD_inh .* sum(((i-last_spike_MD_shape<spikewidth_VA) & (i-last_spike_MD_shape>0)), 1).';

            Isyn_MD_to_Inh(:,i+ delay) = Isyn_MD_to_Inh(:,i+ delay) + W_MD_inh .* sum(((i-last_spike_MD_ori<spikewidth_VA) & (i-last_spike_MD_ori>0)), 1).';
        end

        Isyn_FS(:, i + delay) =  Isyn_FS(:, i + delay) + W_IPL_Inh_to_exc_fs.*sum(((i-last_spike_FS < spikewidth_inh) & (i-last_spike_FS > 0)).*(exp(((i-last_spike_FS < spikewidth_inh) & (i-last_spike_FS > 0)).*(1-(i-last_spike_FS-delay)).*2./spikewidth_inh)), 1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% local excitatory input to
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% shape cells 
        
        Isyn_aPFC_local_shape(:,i+ delay) = Isyn_aPFC_local_shape(:,i+ delay) + sum(((i-last_spike_PFC_remote_shape<spikewidth) & (i-last_spike_PFC_remote_shape>0)) .* matrix_local, 1).';
    
        Isyn_aPFC_local_ori(:,i+ delay) = Isyn_aPFC_local_ori(:,i+ delay) + sum(((i-last_spike_PFC_remote_Orientation<spikewidth) & (i-last_spike_PFC_remote_Orientation>0)) .* matrix_local, 1).';


        


        for j = 1:numberofneurons
            % create noise
            n1=rand;
            noise_amp = 10;
            n1=rand;
            % noise_prob = 0.000005;
            % noise_prob_PFC=0.000005;
            noise_prob = 0.000002;
            noise_prob_PFC=0.000002;
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
%                 spiketimes_M=[spiketimes_M;i,j];

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Layer 2, where the layer 1 send signal for blue triangle to
            % this cell

            % membrane potential across PFC cells

            if (last_spike_PFC_S_BT(j)~=10^10 && (i-last_spike_PFC_S_BT(j))>tref)

                y_PFC_S_BT(j,i) = y_PFC_S_BT(j,i-1)+leaky_coef*((E_L-y_PFC_S_BT(j,i-1))/tha)*dt+(Iext_SPFC_BT(j)+ Isyn_PFC_M_S (j, i)+ (Isyn_aPFC_local_BT(j, i) + Isyn_aPFC_local_shared (j,i))* Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1(j, i)- Isyn_PV(j, i) + Isyn_MD_shape_to_random_PFC_BT(j,i)+Isyn_MD_orientation_to_random_PFC_BT(j,i))*dt*(RM/tha);%- Isyn_PV(j, i)
            elseif last_spike_PFC_S_BT(j)==10^10

                y_PFC_S_BT(j,i) = y_PFC_S_BT(j,i-1)+leaky_coef*((E_L-y_PFC_S_BT(j,i-1))/tha)*dt+(Iext_SPFC_BT(j)+ Isyn_PFC_M_S (j, i)+ (Isyn_aPFC_local_BT(j, i) + Isyn_aPFC_local_shared (j,i))* Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1(j, i)- Isyn_PV(j, i) +Isyn_MD_shape_to_random_PFC_BT(j,i)+Isyn_MD_orientation_to_random_PFC_BT(j,i))*dt*(RM/tha);%- Isyn_PV(j, i)

            else

                y_PFC_S_BT(j,i)=E_L;

            end

            if y_PFC_S_BT(j,i)>=v_th
                last_spike_PFC_S_BT(j)=i;
                y_PFC_S_BT(j,i)=0;
%                 spiketimes_S=[spiketimes_S;i,j];

            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 2 PFC superficial  cells in superficial layers

            % membrane potential across PFC cells

            if (last_spike_PFC_S_RC(j)~=10^10 && (i-last_spike_PFC_S_RC(j))>tref)

                y_PFC_S_RC(j,i) = y_PFC_S_RC(j,i-1)+ leaky_coef*((E_L-y_PFC_S_RC(j,i-1))/tha)*dt + (Iext_SPFC_RC(j)+ (Isyn_aPFC_local_RC(j, i) + Isyn_aPFC_local_shared (j,i))* Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1(j, i)- Isyn_PV(j, i) +Isyn_MD_shape_to_random_PFC_RC(j,i)+Isyn_MD_orientation_to_random_PFC_RC(j,i))*dt*(RM/tha);%- Isyn_PV(j, i)
            elseif last_spike_PFC_S_RC(j)==10^10

                y_PFC_S_RC(j,i) = y_PFC_S_RC(j,i-1)+ leaky_coef*((E_L-y_PFC_S_RC(j,i-1))/tha)*dt + (Iext_SPFC_RC(j)+ (Isyn_aPFC_local_RC(j, i) + Isyn_aPFC_local_shared (j,i))* Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1(j, i)- Isyn_PV(j, i) +Isyn_MD_shape_to_random_PFC_RC(j,i)+Isyn_MD_orientation_to_random_PFC_RC(j,i))*dt*(RM/tha);%- Isyn_PV(j, i)

            else

                y_PFC_S_RC(j,i) = E_L;

            end

            if y_PFC_S_RC(j,i)>= v_th
                last_spike_PFC_S_RC(j)= i;
                y_PFC_S_RC(j,i)= 0;
%                 spiketimes_S_blue= [spiketimes_S_blue;i,j];

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            if (last_spike_PFC_S_GC(j)~=10^10 && (i-last_spike_PFC_S_GC(j))>tref)

                y_PFC_S_GC(j,i) = y_PFC_S_GC(j,i-1)+ leaky_coef*((E_L-y_PFC_S_GC(j,i-1))/tha)*dt + (Iext_SPFC_GC(j) + (Isyn_aPFC_local_GC(j, i) + Isyn_aPFC_local_shared (j,i)) * Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1(j, i)- Isyn_PV(j, i) +Isyn_MD_shape_to_random_PFC_GC(j,i)+Isyn_MD_orientation_to_random_PFC_GC(j,i))*dt*(RM/tha);%- Isyn_PV(j, i)
            elseif last_spike_PFC_S_GC(j)==10^10

                y_PFC_S_GC(j,i) = y_PFC_S_GC(j,i-1)+ leaky_coef*((E_L-y_PFC_S_GC(j,i-1))/tha)*dt + (Iext_SPFC_GC(j) + (Isyn_aPFC_local_GC(j, i) + Isyn_aPFC_local_shared (j,i)) * Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1(j, i)- Isyn_PV(j, i)+Isyn_MD_shape_to_random_PFC_GC(j,i)+Isyn_MD_orientation_to_random_PFC_GC(j,i) )*dt*(RM/tha);%- Isyn_PV(j, i)

            else

                y_PFC_S_GC(j,i) = E_L;

            end

            if y_PFC_S_GC(j,i)>= v_th
                last_spike_PFC_S_GC(j)= i;
                y_PFC_S_GC(j,i)= 0;
%                 spiketimes_S_red= [spiketimes_S_red;i,j];

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            if (last_spike_PFC_S_YT(j)~=10^10 && (i-last_spike_PFC_S_YT(j))>tref)

                y_PFC_S_YT(j,i) = y_PFC_S_YT(j,i-1)+ leaky_coef*((E_L-y_PFC_S_YT(j,i-1))/tha)*dt + (Iext_SPFC_YT(j) + (Isyn_aPFC_local_YT(j, i) + Isyn_aPFC_local_shared (j,i)) * Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1(j, i) - Isyn_PV(j, i)+Isyn_MD_shape_to_random_PFC_YT(j,i)+Isyn_MD_orientation_to_random_PFC_YT(j,i))*dt*(RM/tha);%- Isyn_PV(j, i)
            elseif last_spike_PFC_S_YT(j)==10^10

                y_PFC_S_YT(j,i) = y_PFC_S_YT(j,i-1)+ leaky_coef*((E_L-y_PFC_S_YT(j,i-1))/tha)*dt + (Iext_SPFC_YT(j) + (Isyn_aPFC_local_YT(j, i) + Isyn_aPFC_local_shared (j,i)) * Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1(j, i) - Isyn_PV(j, i)+Isyn_MD_shape_to_random_PFC_YT(j,i)+Isyn_MD_orientation_to_random_PFC_YT(j,i))*dt*(RM/tha);%- Isyn_PV(j, i)

            else

                y_PFC_S_YT(j,i) = E_L;

            end

            if y_PFC_S_YT(j,i)>= v_th
                last_spike_PFC_S_YT(j)= i;
                y_PFC_S_YT(j,i)= 0;
%                 spiketimes_S_YT= [spiketimes_S_YT;i,j];

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

                y_aPFC_D_BT(j,i) = y_aPFC_D_BT(j,i-1)+leaky_coef*((E_L-y_aPFC_D_BT(j,i-1))/tha)*dt+(Iext_DPFC_BT(j)+ (W11 * Isyn_aPFC_to_D_BT(j,i)+ W12 * Isyn_aPFC_to_D_RC(j,i)+ W13 *Isyn_aPFC_to_D_GC(j,i)+ W14* Isyn_aPFC_to_D_YT(j,i)) + Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2(j, i))*dt*(RM/tha); %+Isyn_IPL_PFC(j,i)-Isyn_aPFC_Inh_to_exc(j,i)

            elseif last_spike_PFC_D_BT(j)==10^10

                y_aPFC_D_BT(j,i) = y_aPFC_D_BT(j,i-1)+leaky_coef*((E_L-y_aPFC_D_BT(j,i-1))/tha)*dt+(Iext_DPFC_BT(j)+ (W11 * Isyn_aPFC_to_D_BT(j,i)+ W12 * Isyn_aPFC_to_D_RC(j,i)+ W13 *Isyn_aPFC_to_D_GC(j,i)+ W14* Isyn_aPFC_to_D_YT(j,i)) + Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2(j, i))*dt*(RM/tha); %+Isyn_IPL_PFC(j,i)-Isyn_aPFC_Inh_to_exc(j,i)

            else

                y_aPFC_D_BT(j,i)= E_L;

            end

            if y_aPFC_D_BT(j,i)>= v_th
                last_spike_PFC_D_BT(j)= i;
                y_aPFC_D_BT(j,i)= 0;
%                 spiketimes_D_BT= [spiketimes_D_BT;i,j];

            end


            % membrane potential across deep layers of PFC cells

            if (last_spike_PFC_D_RC(j)~=10^10 && (i-last_spike_PFC_D_RC(j))>tref)

                y_aPFC_D_RC(j,i) = y_aPFC_D_RC(j,i-1)+leaky_coef*((E_L-y_aPFC_D_RC(j,i-1))/tha)*dt+(Iext_DPFC_RC(j)+ (W21 * Isyn_aPFC_to_D_BT(j,i)+ W22 * Isyn_aPFC_to_D_RC(j,i)+ W23 *Isyn_aPFC_to_D_GC(j,i)+ W24* Isyn_aPFC_to_D_YT(j,i)) +Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2(j, i))*dt*(RM/tha); %+Isyn_IPL_PFC(j,i)-Isyn_aPFC_Inh_to_exc(j,i)

            elseif last_spike_PFC_D_RC(j)==10^10

                y_aPFC_D_RC(j,i) = y_aPFC_D_RC(j,i-1)+leaky_coef*((E_L-y_aPFC_D_RC(j,i-1))/tha)*dt+(Iext_DPFC_RC(j)+ (W21 * Isyn_aPFC_to_D_BT(j,i)+ W22 * Isyn_aPFC_to_D_RC(j,i)+ W23 *Isyn_aPFC_to_D_GC(j,i)+ W24* Isyn_aPFC_to_D_YT(j,i)) +Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2(j, i))*dt*(RM/tha); %+Isyn_IPL_PFC(j,i)-Isyn_aPFC_Inh_to_exc(j,i)

            else

                y_aPFC_D_RC(j,i)= E_L;

            end

            if y_aPFC_D_RC(j,i)>= v_th
                last_spike_PFC_D_RC(j)= i;
                y_aPFC_D_RC(j,i)= 0;
%                 spiketimes_D_RC = [spiketimes_D_RC;i,j];

            end

            %%%%%%%%%%%% membrane potential across deep layers of PFC cells

            if (last_spike_PFC_D_GC(j)~=10^10 && (i-last_spike_PFC_D_GC(j))>tref)

                y_aPFC_D_GC(j,i) = y_aPFC_D_GC(j,i-1)+leaky_coef*((E_L-y_aPFC_D_GC(j,i-1))/tha)*dt+(Iext_DPFC_GC(j)+ (W31 * Isyn_aPFC_to_D_BT(j,i)+ W32 * Isyn_aPFC_to_D_RC(j,i)+ W33 *Isyn_aPFC_to_D_GC(j,i)+ W34* Isyn_aPFC_to_D_YT(j,i))+ Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2(j, i))*dt*(RM/tha); %+Isyn_IPL_PFC(j,i)-Isyn_aPFC_Inh_to_exc(j,i)

            elseif last_spike_PFC_D_GC(j)==10^10

                y_aPFC_D_GC(j,i) = y_aPFC_D_GC(j,i-1)+leaky_coef*((E_L-y_aPFC_D_GC(j,i-1))/tha)*dt+(Iext_DPFC_GC(j)+ (W31 * Isyn_aPFC_to_D_BT(j,i)+ W32 * Isyn_aPFC_to_D_RC(j,i)+ W33 *Isyn_aPFC_to_D_GC(j,i)+ W34* Isyn_aPFC_to_D_YT(j,i)) + Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2(j, i))*dt*(RM/tha); %+Isyn_IPL_PFC(j,i)-Isyn_aPFC_Inh_to_exc(j,i)

            else

                y_aPFC_D_GC(j,i)= E_L;

            end

            if y_aPFC_D_GC(j,i)>= v_th
                last_spike_PFC_D_GC(j)= i;
                y_aPFC_D_GC(j,i)= 0;
%                 spiketimes_D_GC = [spiketimes_D_GC;i,j];

            end

            %%%%%%%%%%%% membrane potential across deep layers of PFC cells

            if (last_spike_PFC_D_YT(j)~=10^10 && (i-last_spike_PFC_D_YT(j))>tref)

                y_aPFC_D_YT(j,i) = y_aPFC_D_YT(j,i-1)+leaky_coef*((E_L-y_aPFC_D_YT(j,i-1))/tha)*dt+(Iext_DPFC_YT(j)+ (W41 * Isyn_aPFC_to_D_BT(j,i)+ W42 * Isyn_aPFC_to_D_RC(j,i)+ W43 *Isyn_aPFC_to_D_GC(j,i)+ W44* Isyn_aPFC_to_D_YT(j,i))+ Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2(j, i))*dt*(RM/tha); %+Isyn_IPL_PFC(j,i)-Isyn_aPFC_Inh_to_exc(j,i)

            elseif last_spike_PFC_D_YT(j)==10^10

                y_aPFC_D_YT(j,i) = y_aPFC_D_YT(j,i-1)+leaky_coef*((E_L-y_aPFC_D_YT(j,i-1))/tha)*dt+(Iext_DPFC_YT(j)+ (W41 * Isyn_aPFC_to_D_BT(j,i)+ W42 * Isyn_aPFC_to_D_RC(j,i)+ W43 *Isyn_aPFC_to_D_GC(j,i)+ W44* Isyn_aPFC_to_D_YT(j,i))+ Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2(j, i))*dt*(RM/tha); %+Isyn_IPL_PFC(j,i)-Isyn_aPFC_Inh_to_exc(j,i)

            else

                y_aPFC_D_YT(j,i)= E_L;

            end

            if y_aPFC_D_YT(j,i)>= v_th
                last_spike_PFC_D_YT(j)= i;
                y_aPFC_D_YT(j,i)= 0;
%                 spiketimes_D_YT = [spiketimes_D_YT;i,j];

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
%                 spiketimes_VA=[spiketimes_VA;i,j];

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % VA matrix cells connecting aPFC to pPFC
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

                y_MD_core_shape(j,i) = y_MD_core_shape(j,i-1)+((E_L-y_MD_core_shape(j,i-1))/tha)*dt+(Iext_MD(j)+ I_noise_md_2(j,i)+ Isyn_PFC_D_MD_shape(j,i) )*dt*(RM/tha);%

            elseif last_spike_MD_shape(j)==10^10

                y_MD_core_shape(j,i) = y_MD_core_shape(j,i-1)+((E_L-y_MD_core_shape(j,i-1))/tha)*dt+(Iext_MD(j)+ I_noise_md_2(j,i)+ Isyn_PFC_D_MD_shape(j,i) )*dt*(RM/tha);%

            else

                y_MD_core_shape(j,i) = E_L;

            end

            if y_MD_core_shape(j,i) >= v_th
                last_spike_MD_shape(j)=i;
                y_MD_core_shape(j,i) = 0;
%                 spiketimes_MD2=[spiketimes_MD2;i,j];

            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % % Core MD cells amplifying aPFC connections
            n1 = rand;
            noise_prob_md = noise_prob_PFC;
            if n1<noise_prob_md
                I_noise_md_1(j,i) =0;% noise_amp;
            else
                I_noise_md_1(j,i) = 0;
            end

            % md cells to IPL;        DEEP LAYER aPFC drives these

            if (last_spike_MD_ori(j)~=10^10 && (i-last_spike_MD_ori(j))>tref)

                y_MD_core_ori(j,i) = y_MD_core_ori(j,i-1)+((E_L-y_MD_core_ori(j,i-1))/tha)*dt+(Iext_MD(j)+I_noise_md_1(j,i)+ Isyn_PFC_D_MD_ori(j, i) )*dt*(RM/tha);%

            elseif last_spike_MD_ori(j)==10^10

                y_MD_core_ori(j,i) = y_MD_core_ori(j,i-1)+((E_L-y_MD_core_ori(j,i-1))/tha)*dt+(Iext_MD(j)+I_noise_md_1(j,i)+ Isyn_PFC_D_MD_ori(j, i) )*dt*(RM/tha);%

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
%                 spiketimes_PV_inh=[spiketimes_PV_inh;i,j];


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
%                 spiketimes_ipl_inh=[spiketimes_ipl_inh;i,j];


            end
            %%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Shape cells in pPFC
            % Isyn_aPFC_D_to_Remote_effective_3(j,i) = Isyn_aPFC_D_to_Remote(j,i)*Isyn_MD_Matrix_to_PFC(j, i);
            % Isyn_pPFC_local_effective_3(j,i) = Isyn_pPFC_local(j,i)*Isyn_MD2_to_pPFC(j, i);
            n1=rand;
            if n1<noise_prob_PFC
                I_noise_PFC_6(j,i)= noise_amp;
            else
                I_noise_PFC_6(j,i)= 0;
            end
            %             %
            if (last_spike_PFC_remote_shape(j)~=10^10 && (i-last_spike_PFC_remote_shape(j))>tref)

                y_pPFC_Shape(j,i) = y_pPFC_Shape(j,i-1)+leaky_coef*((E_L-y_pPFC_Shape(j,i-1))/tha)*dt+(Iext_PFC_remote_Shape(j)+ Isyn_aPFC_local_shared (j,i) * Isyn_VA_Matrix_to_PFC(j, i) + I_noise_PFC_6(j,i)+ Isyn_VA_Matrix_to_PFC_exc_shape(j, i)  + Isyn_MD_shape_to_PFC(j,i) + Isyn_aPFC_local_shape(j,i)* Isyn_VA_Matrix_to_PFC(j,i)  - 0* Isyn_FS(j, i) )*dt*(RM/tha);%

            elseif last_spike_PFC_remote_shape(j)==10^10

                y_pPFC_Shape(j,i) = y_pPFC_Shape(j,i-1)+leaky_coef*((E_L-y_pPFC_Shape(j,i-1))/tha)*dt+(Iext_PFC_remote_Shape(j)+ Isyn_aPFC_local_shared (j,i) * Isyn_VA_Matrix_to_PFC(j, i) + I_noise_PFC_6(j,i)+ Isyn_VA_Matrix_to_PFC_exc_shape(j, i)  +  Isyn_MD_shape_to_PFC(j,i) + Isyn_aPFC_local_shape(j,i)* Isyn_VA_Matrix_to_PFC(j,i)  - 0*  Isyn_FS(j, i) )*dt*(RM/tha);%

            else

                y_pPFC_Shape(j,i)=E_L;

            end

            if y_pPFC_Shape(j,i)>=v_th
                last_spike_PFC_remote_shape(j)=i;
                y_pPFC_Shape(j,i)=0;
%                 spiketimes_S_remote_Shape=[spiketimes_S_remote_Shape;i,j];

            end


            %%%%%%%%%%%%%%%%%%%%%%%%
            n1=rand;
            if n1<noise_prob_PFC
                I_noise_PFC_6(j,i)= noise_amp;
            else
                I_noise_PFC_6(j,i)= 0;
            end
            %             %
            if (last_spike_PFC_remote_Orientation(j)~=10^10 && (i-last_spike_PFC_remote_Orientation(j))>tref)

                y_pPFC_Orientation(j,i) = y_pPFC_Orientation(j,i-1)+leaky_coef*((E_L-y_pPFC_Orientation(j,i-1))/tha)*dt+(Iext_PFC_remote_Shape(j)+ I_noise_PFC_6(j,i)+ Isyn_aPFC_local_shared (j,i) * Isyn_VA_Matrix_to_PFC(j, i) + Isyn_VA_Matrix_to_PFC_exc_orientation(j, i) +Isyn_MD_ori_to_PFC(j,i) + Isyn_aPFC_local_ori(j,i)* Isyn_VA_Matrix_to_PFC(j,i)   - Isyn_FS(j, i))*dt*(RM/tha);

            elseif last_spike_PFC_remote_Orientation(j)==10^10

                y_pPFC_Orientation(j,i) = y_pPFC_Orientation(j,i-1)+leaky_coef*((E_L-y_pPFC_Orientation(j,i-1))/tha)*dt+(Iext_PFC_remote_Shape(j)+ I_noise_PFC_6(j,i)+ Isyn_aPFC_local_shared (j,i) * Isyn_VA_Matrix_to_PFC(j, i) + Isyn_VA_Matrix_to_PFC_exc_orientation(j, i) +Isyn_MD_ori_to_PFC(j,i) + Isyn_aPFC_local_ori(j,i)* Isyn_VA_Matrix_to_PFC(j,i)  - Isyn_FS(j, i))*dt*(RM/tha);

            else

                y_pPFC_Orientation(j,i)=E_L;

            end

            if y_pPFC_Orientation(j,i)>=v_th
                last_spike_PFC_remote_Orientation(j)=i;
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
%                 spiketimes_PV_inh=[spiketimes_PV_inh;i,j];
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
%                 spiketimes_FS_inh=[spiketimes_FS_inh;i,j];
            end




        end
    end

    % LFP_aPFC_S(rr,:)= smooth(mean(Isyn_aPFC_local_effective)+ mean(Isyn_pPFC_Remote_to_aPFC_effective)+ mean(I_noise_PFC_1));% %%% !!!!!!!!!!!!!!!! add noise
    % LFP_aPFC_D(rr,:)= smooth(mean(Isyn_aPFC_to_D_effective)+ mean(Isyn_aPFC_D_effective_2)+ mean(I_noise_PFC_2));% post-synaptic current in pfc
    % LFP_pPFC_S(rr,:)= smooth(mean(Isyn_pPFC_local_effective)+ mean(Isyn_pPFC_local_effective_2)+mean(Isyn_pPFC_local_effective_3)+ mean(Isyn_aPFC_D_to_Remote_effective)+ mean(Isyn_aPFC_D_to_Remote_effective_2)+ mean(Isyn_aPFC_D_to_Remote_effective_3)+ mean(I_noise_PFC_5)+ mean(I_noise_PFC_4)+ mean(I_noise_PFC_6));
    % LFP_pPFC_D(rr,:)= smooth(mean(Isyn_pPFC_S_to_D_effective)+ mean(Isyn_pPFC_D_effective)+ mean(I_noise_PFC_3));% post-synaptic current in pfc
    % LFP_response(rr,:)= smooth(mean(Isyn_pPFCD_to_response_effective)+ mean(I_noise_PFC));
        
    [neuo_id, firing_time] = find(get_firing(y_PFC_S_BT, v_th, 1));
    full_PFC_S_BT = [full_PFC_S_BT; cat(2, neuo_id, firing_time)];

    [neuo_id, firing_time] = find(get_firing(y_PFC_S_RC, v_th, 1));
    full_PFC_S_RC = [full_PFC_S_RC; cat(2, neuo_id, firing_time)];

    [neuo_id, firing_time] = find(get_firing(y_PFC_S_GC, v_th, 1));
    full_PFC_S_GC = [full_PFC_S_GC; cat(2, neuo_id, firing_time)];

    [neuo_id, firing_time] = find(get_firing(y_PFC_S_YT, v_th, 1));
    full_PFC_S_YT = [full_PFC_S_YT; cat(2, neuo_id, firing_time)];
    
    [neuo_id, firing_time] = find(get_firing(y_aPFC_D_BT, v_th, 1));
    full_PFC_D_BT = [full_PFC_D_BT; cat(2, neuo_id, firing_time)];

    [neuo_id, firing_time] = find(get_firing(y_aPFC_D_RC, v_th, 1));
    full_PFC_D_RC = [full_PFC_D_RC; cat(2, neuo_id, firing_time)];

    [neuo_id, firing_time] = find(get_firing(y_aPFC_D_GC, v_th, 1));
    full_PFC_D_GC = [full_PFC_D_GC; cat(2, neuo_id, firing_time)];

    [neuo_id, firing_time] = find(get_firing(y_aPFC_D_YT, v_th, 1));
    full_PFC_D_YT = [full_PFC_D_YT; cat(2, neuo_id, firing_time)];
    
    [neuo_id, firing_time] = find(get_firing(y_VA_matrix_shape, v_th, 1));
    full_VA_shape = [full_VA_shape; cat(2, neuo_id, firing_time)];

    [neuo_id, firing_time] = find(get_firing(y_VA_matrix_Orientation, v_th, 1));
    full_VA_ori = [full_VA_ori; cat(2, neuo_id, firing_time)];
    
    [neuo_id, firing_time] = find(get_firing(y_MD_core_shape, v_th, 1));
    full_MD_shape = [full_MD_shape; cat(2, neuo_id, firing_time)];
     
    [neuo_id, firing_time] = find(get_firing(y_MD_core_ori, v_th, 1));
    full_MD_ori = [full_MD_ori; cat(2, neuo_id, firing_time)];
    
    [neuo_id, firing_time] = find(get_firing(y_pPFC_Shape, v_th, 1));
    full_PFC_remote_shape = [full_PFC_remote_shape; cat(2, neuo_id, firing_time)];
    
    [neuo_id, firing_time] = find(get_firing(y_pPFC_Orientation, v_th, 1));
    full_PFC_remote_ori = [full_PFC_remote_ori; cat(2, neuo_id, firing_time)];

    VA_neuron_index = get_target_neuron_idex(get_firing(y_VA_matrix_shape, v_th, 1));
    MD_neuron_index = get_target_neuron_idex(get_firing(y_MD_core_shape, v_th, 1));
    pPFC_neuron_index = get_target_neuron_idex(get_firing(y_pPFC_Shape, v_th, 1));

    % only select target neuron set, movsum in get_firing_num keeps the dimenstion same   
    full_VA_shape_num(rr, :) =  get_firing_num(get_firing(y_VA_matrix_shape(VA_neuron_index, :), v_th, 1), gauss_width);
    full_MD_shape_num(rr, :) = get_firing_num(get_firing(y_MD_core_shape(MD_neuron_index, :), v_th, 1), gauss_width);
    full_PFC_remote_shape_num(rr, :) = get_firing_num(get_firing(y_pPFC_Shape(pPFC_neuron_index, :), v_th, 1), gauss_width);
    
    disp(rr);
    toc
        
end


save("spike_time", ...
    'full_PFC_S_BT', ...
    'full_PFC_S_RC', ...
    'full_PFC_S_GC', ...
    'full_PFC_S_YT', ...
    'full_PFC_D_BT', ...
    'full_PFC_D_RC', ...
    'full_PFC_D_GC', ...
    'full_PFC_D_YT', ...
    'full_VA_shape', ...
    'full_VA_ori', ...
    'full_MD_shape', ...
    'full_MD_ori', ...
    'full_PFC_remote_shape', ...
    'full_PFC_remote_ori', ...
    'full_VA_shape_num', ...
    'full_MD_shape_num', ...
    'full_PFC_remote_shape_num', ...
    '-v7.3');

% trial_num = 1;
% 
% t_PFC_S_BT = zeros(numberofneurons, length(T));
% t_PFC_S_RC = zeros(numberofneurons, length(T));
% t_PFC_S_GC = zeros(numberofneurons, length(T));
% t_PFC_S_YT = zeros(numberofneurons, length(T));
% 
% t_PFC_D_BT = zeros(numberofneurons, length(T));
% t_PFC_D_RC = zeros(numberofneurons, length(T));
% t_PFC_D_GC = zeros(numberofneurons, length(T));
% t_PFC_D_YT = zeros(numberofneurons, length(T));
% 
% t_VA_shape = zeros(numberofneurons, length(T));
% t_VA_ori = zeros(numberofneurons, length(T));
% 
% t_PFC_remote_shape = zeros(numberofneurons, length(T));
% t_PFC_remote_ori = zeros(numberofneurons, length(T));
% 
% t_MD_shape = zeros(numberofneurons, length(T));
% t_MD_ori = zeros(numberofneurons, length(T));
% 
% for i = 1:length(full_PFC_S_BT{trial_num})
%     t_PFC_S_BT(full_PFC_S_BT{trial_num}(i, 1), full_PFC_S_BT{trial_num}(i, 2)) = 1;
% end
% 
% for i = 1:length(full_PFC_S_RC{trial_num})
%     t_PFC_S_RC(full_PFC_S_RC{trial_num}(i, 1), full_PFC_S_RC{trial_num}(i, 2)) = 1;
% end
% 
% for i = 1:length(full_PFC_S_GC{trial_num})
%     t_PFC_S_GC(full_PFC_S_GC{trial_num}(i, 1), full_PFC_S_GC{trial_num}(i, 2)) = 1;
% end
% 
% for i = 1:length(full_PFC_S_YT{trial_num})
%     t_PFC_S_YT(full_PFC_S_YT{trial_num}(i, 1), full_PFC_S_YT{trial_num}(i, 2)) = 1;
% end
% 
% for i = 1:length(full_PFC_D_BT{trial_num})
%     t_PFC_D_BT(full_PFC_D_BT{trial_num}(i, 1), full_PFC_D_BT{trial_num}(i, 2)) = 1;
% end
% 
% for i = 1:length(full_PFC_D_RC{trial_num})
%     t_PFC_D_RC(full_PFC_D_RC{trial_num}(i, 1), full_PFC_D_RC{trial_num}(i, 2)) = 1;
% end
% 
% for i = 1:length(full_PFC_D_GC{trial_num})
%     t_PFC_D_GC(full_PFC_D_GC{trial_num}(i, 1), full_PFC_D_GC{trial_num}(i, 2)) = 1;
% end
% 
% for i = 1:length(full_PFC_D_YT{trial_num})
%     t_PFC_D_YT(full_PFC_D_YT{trial_num}(i, 1), full_PFC_D_YT{trial_num}(i, 2)) = 1;
% end
% 
% for i = 1:length(full_VA_shape{trial_num})
%     t_VA_shape(full_VA_shape{trial_num}(i, 1), full_VA_shape{trial_num}(i, 2)) = 1;
% end
% 
% for i = 1:length(full_PFC_remote_shape{trial_num})
%     t_PFC_remote_shape(full_PFC_remote_shape{trial_num}(i, 1), full_PFC_remote_shape{trial_num}(i, 2)) = 1;
% end
% 
% for i = 1:length(full_MD_shape{trial_num})
%     t_MD_shape(full_MD_shape{trial_num}(i, 1), full_MD_shape{trial_num}(i, 2)) = 1;
% end
% 
% figure(1)
% 
% subplot(4,1,1)
% spy( t_PFC_S_BT,8,'b'),title('PFC Superficial blue triangle', 'FontSize', 16)
% set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
% 
% 
% subplot(4,1,2)
% spy( t_PFC_S_RC,8,'b'),title('PFC Superficial red circle ', 'FontSize', 16)
% set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
% 
% subplot(4,1,3)
% spy( t_PFC_S_GC,8,'b'),title('PFC Superficial green circle ', 'FontSize', 16)
% set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
% 
% subplot(4,1,4)
% spy(t_PFC_S_YT,8,'b'),title('PFC Superficial yellow triangle', 'FontSize', 16)
% set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
% %
% figure (2)
% subplot(4,1,1)
% spy( t_PFC_D_BT,8,'k'),title('PFC Deep blue triangle', 'FontSize', 16)
% set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
% subplot(4,1,2)
% spy( t_PFC_D_RC,8,'k'),title('PFC Deep red circle', 'FontSize', 16)
% set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
% 
% subplot(4,1,3)
% spy( t_PFC_D_GC,8,'k'),title('PFC Deep green circle', 'FontSize', 16)
% set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
% subplot(4,1,4)
% spy( t_PFC_D_YT,8,'k'),title('PFC Deep yellow triangle', 'FontSize', 16)
% set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
% 
% 
% figure(3)
% subplot(3,1,1)
% spy( t_VA_shape,8,'k'),title('VA Thalamus Shape', 'FontSize', 16)
% set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
% subplot(3,1,2)
% spy( t_PFC_remote_shape,8,'b'),title('remote PFC Shape', 'FontSize', 16)
% set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
% subplot(3,1,3)
% spy( t_MD_shape,8,'k'),title('MD', 'FontSize', 16)
% set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')


%%
neuron_id = 1;
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
    PFC_S_BT(i, temp(temp(:, 1) == neuron_id, 2)) = 1;

    temp = full_PFC_S_RC{i};
    PFC_S_RC(i, temp(temp(:, 1) == neuron_id, 2)) = 1;

    temp = full_PFC_S_GC{i};
    PFC_S_GC(i, temp(temp(:, 1) == neuron_id, 2)) = 1;

    temp = full_PFC_S_YT{i};
    PFC_S_YT(i, temp(temp(:, 1) == neuron_id, 2)) = 1;

    temp = full_PFC_D_BT{i};
    PFC_D_BT(i, temp(temp(:, 1) == neuron_id, 2)) = 1;

    temp = full_PFC_D_RC{i};
    PFC_D_RC(i, temp(temp(:, 1) == neuron_id, 2)) = 1;

    temp = full_PFC_D_GC{i};
    PFC_D_GC(i, temp(temp(:, 1) == neuron_id, 2)) = 1;

    temp = full_PFC_D_YT{i};
    PFC_D_YT(i, temp(temp(:, 1) == neuron_id, 2)) = 1;

    temp = full_VA_shape{i};
    VA_shape(i, temp(temp(:, 1) == neuron_id, 2)) = 1;

    % temp = full_VA_ori{i};
    % VA_ori(i, temp(temp(:, 1) == neuron_id, 2)) = 1;

    temp = full_PFC_remote_shape{i};
    PFC_remote_shape(i, temp(temp(:, 1) == neuron_id, 2)) = 1;

    % temp = full_PFC_remote_ori{i};
    % PFC_remote_ori(i, temp(temp(:, 1) == neuron_id, 2)) = 1;

    temp = full_MD_shape{i};
    MD_shape(i, temp(temp(:, 1) == neuron_id, 2)) = 1;

    % temp = full_MD_ori{i};
    % MD_ori(i, temp(temp(:, 1) == neuron_id, 2)) = 1;
end

figure(1)

subplot(4,1,1)
spy( PFC_S_BT, 8,'b'),title('PFC Superficial blue triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[10000 1 1]),ylabel('Trial number', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')


subplot(4,1,2)
spy( PFC_S_RC, 8,'b'),title('PFC Superficial red circle ', 'FontSize', 16)
set(gca,'DataAspectRatio',[10000 1 1]),ylabel('Trial number', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

subplot(4,1,3)
spy( PFC_S_GC, 8,'b'),title('PFC Superficial green circle ', 'FontSize', 16)
set(gca,'DataAspectRatio',[10000 1 1]),ylabel('Trial number', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

subplot(4,1,4)
spy( PFC_S_YT,8,'b'),title('PFC Superficial yellow triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[10000 1 1]),ylabel('Trial number', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

%
figure (2)
subplot(4,1,1)
spy( PFC_D_BT,8,'k'),title('PFC Deep blue triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[10000 1 1]),ylabel('Trial number', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(4,1,2)
spy( PFC_D_RC,8,'k'),title('PFC Deep red circle', 'FontSize', 16)
set(gca,'DataAspectRatio',[10000 1 1]),ylabel('Trial number', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

subplot(4,1,3)
spy( PFC_D_GC,8,'k'),title('PFC Deep green circle', 'FontSize', 16)
set(gca,'DataAspectRatio',[10000 1 1]),ylabel('Trial number', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(4,1,4)
spy( PFC_D_YT,8,'k'),title('PFC Deep yellow triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[10000 1 1]),ylabel('Trial number', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')


figure( 3)
subplot(6,1,1)
spy( VA_shape,8,'k'),title('VA Thalamus Shape', 'FontSize', 16)
set(gca,'DataAspectRatio',[10000 1 1]),ylabel('Trial number', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(6,1,2)
spy( VA_ori,8,'k'),title('VA Thalamus Orientation', 'FontSize', 16)
set(gca,'DataAspectRatio',[10000 1 1]),ylabel('Trial number', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

subplot(6,1,3)
spy( PFC_remote_shape,8,'b'),title('remote PFC Shape', 'FontSize', 16)
set(gca,'DataAspectRatio',[10000 1 1]),ylabel('Trial number', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(6,1,4)
spy( PFC_remote_ori,8,'b'),title('remote PFC Ori', 'FontSize', 16)
set(gca,'DataAspectRatio',[10000 1 1]),ylabel('Trial number', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

subplot(6,1,5)
spy( MD_shape,8,'k'),title('MD Shape', 'FontSize', 16)
set(gca,'DataAspectRatio',[10000 1 1]),ylabel('Trial number', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

subplot(6,1,6)
spy( MD_ori,8,'k'),title('MD ori', 'FontSize', 16)
set(gca,'DataAspectRatio',[10000 1 1]),ylabel('Trial number', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

figure(4)
srate=1000;
st= 10000;
final=250000;

average_VA_shape_num = sum(full_VA_shape_num, 1);
VA_firing_rate_timeseries=  conv_gaussian(average_VA_shape_num, srate,gauss_width);

average_MD_shape_num = sum(full_MD_shape_num, 1);
MD_firing_rate_timeseries=  conv_gaussian(average_MD_shape_num, srate,gauss_width);

average_pPFC_shape_num = sum(full_PFC_remote_shape_num, 1);
pPFC_firing_rate_timeseries=  conv_gaussian(average_pPFC_shape_num, srate,gauss_width);

plot(VA_firing_rate_timeseries, 'b'), xlim([st, final])
hold on
plot(MD_firing_rate_timeseries, 'g'), xlim([st, final])
hold on
plot(pPFC_firing_rate_timeseries, 'r'), xlim([st, final])
legend('VA', 'MD', 'PFC remote')


%%
trial_num = 10;

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

for i = 1:length(full_PFC_S_BT{trial_num})
    t_PFC_S_BT(full_PFC_S_BT{trial_num}(i, 1), full_PFC_S_BT{trial_num}(i, 2)) = 1;
end

for i = 1:length(full_PFC_S_RC{trial_num})
    t_PFC_S_RC(full_PFC_S_RC{trial_num}(i, 1), full_PFC_S_RC{trial_num}(i, 2)) = 1;
end

for i = 1:length(full_PFC_S_GC{trial_num})
    t_PFC_S_GC(full_PFC_S_GC{trial_num}(i, 1), full_PFC_S_GC{trial_num}(i, 2)) = 1;
end

for i = 1:length(full_PFC_S_YT{trial_num})
    t_PFC_S_YT(full_PFC_S_YT{trial_num}(i, 1), full_PFC_S_YT{trial_num}(i, 2)) = 1;
end

for i = 1:length(full_PFC_D_BT{trial_num})
    t_PFC_D_BT(full_PFC_D_BT{trial_num}(i, 1), full_PFC_D_BT{trial_num}(i, 2)) = 1;
end

for i = 1:length(full_PFC_D_RC{trial_num})
    t_PFC_D_RC(full_PFC_D_RC{trial_num}(i, 1), full_PFC_D_RC{trial_num}(i, 2)) = 1;
end

for i = 1:length(full_PFC_D_GC{trial_num})
    t_PFC_D_GC(full_PFC_D_GC{trial_num}(i, 1), full_PFC_D_GC{trial_num}(i, 2)) = 1;
end

for i = 1:length(full_PFC_D_YT{trial_num})
    t_PFC_D_YT(full_PFC_D_YT{trial_num}(i, 1), full_PFC_D_YT{trial_num}(i, 2)) = 1;
end

for i = 1:length(full_VA_shape{trial_num})
    t_VA_shape(full_VA_shape{trial_num}(i, 1), full_VA_shape{trial_num}(i, 2)) = 1;
end

for i = 1:length(full_PFC_remote_shape{trial_num})
    t_PFC_remote_shape(full_PFC_remote_shape{trial_num}(i, 1), full_PFC_remote_shape{trial_num}(i, 2)) = 1;
end

for i = 1:length(full_MD_shape{trial_num})
    t_MD_shape(full_MD_shape{trial_num}(i, 1), full_MD_shape{trial_num}(i, 2)) = 1;
end

figure(1)

subplot(4,1,1)
spy( t_PFC_S_BT,8,'b'),title('PFC Superficial blue triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')


subplot(4,1,2)
spy( t_PFC_S_RC,8,'b'),title('PFC Superficial red circle ', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

subplot(4,1,3)
spy( t_PFC_S_GC,8,'b'),title('PFC Superficial green circle ', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

subplot(4,1,4)
spy(t_PFC_S_YT,8,'b'),title('PFC Superficial yellow triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
%
figure (2)
subplot(4,1,1)
spy( t_PFC_D_BT,8,'k'),title('PFC Deep blue triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(4,1,2)
spy( t_PFC_D_RC,8,'k'),title('PFC Deep red circle', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

subplot(4,1,3)
spy( t_PFC_D_GC,8,'k'),title('PFC Deep green circle', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(4,1,4)
spy( t_PFC_D_YT,8,'k'),title('PFC Deep yellow triangle', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')


figure(3)
subplot(3,1,1)
spy( t_VA_shape,8,'k'),title('VA Thalamus Shape', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(3,1,2)
spy( t_PFC_remote_shape,8,'b'),title('remote PFC Shape', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(3,1,3)
spy( t_MD_shape,8,'k'),title('MD', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')


%%
% save("good_initial", ...
%     'Iext_SPFC_YT', ...
%     'Iext_SPFC_BT', ...
%     'Iext_SPFC_RC', ...
%     'Iext_SPFC_GC', ...
%     'Iext_DPFC_YT', ...
%     'Iext_DPFC_BT', ...
%     'Iext_DPFC_RC', ...
%     'Iext_DPFC_GC', ...
%     'Iext_PFC_M', ...
%     'Iext_PFC_remote_ORI', ...
%     'Iext_PFC_remote_Shape', ...
%     'Iext_PFC_D', ...
%     'Iext_PFC_D_2', ...
%     'Iext_PFC_D_3', ...
%     'Iext_MD', ...
%     'Iext_VA', ...
%     'Iext_Inh', ...
%     'Iext_Inh_pPFC', ...
%     'Iext_Inh_PV', ...
%     'Iext_response');
% Iext_SPFC_YT = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
% Iext_SPFC_BT = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
% Iext_SPFC_RC = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
% Iext_SPFC_GC = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
% Iext_DPFC_YT = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
% Iext_DPFC_BT = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
% Iext_DPFC_RC = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
% Iext_DPFC_GC = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
% Iext_PFC_M = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5);
% Iext_PFC_remote_ORI = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
% Iext_PFC_remote_Shape = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
% Iext_PFC_D = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
% Iext_PFC_D_2 = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
% Iext_PFC_D_3 = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); % similar driving current for PFC and IPL cells small variation across cells
% Iext_MD = 1.5 - 0.005 *rand(numberofneurons,1); % same driving current for all thalamic cells
% Iext_VA = 1.5 - 0.01 *rand(numberofneurons,1); % same driving current for all thalamic cells
% Iext_Inh = 1.5 -  0.01 *(rand(numberofneurons,1)-0.5); % driving current for inhibitory cells similar in PFC and IPL
% Iext_Inh_pPFC = 1.5 - 0.01 *(rand(numberofneurons,1)-0.5); %driving current for inhibitory cells similar in PFC and IPL
% Iext_Inh_PV = 1.49 - 0.0 *(rand(numberofneurons,1));
% Iext_response = 1.5 - 0.01 *(rand(numberofneurons,1));

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
spy( y_VA_matrix_Orientation>-50,8,'b'),title('VA Thalamus Orientation', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(5,1,5)
spy( y_PV_inh>-50,8,'r'),title('inhibitory PFC cells', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')



%
figure(5)

subplot(5,1,1)
spy( y_MD_core_shape>-50,8,'k'),title('MD shape', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(5,1,2)

spy( y_MD_core_ori>-50,8,'k'),title('MD Ori', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

subplot(5,1,3)
spy( y_pPFC_Shape>-50,8,'b'),title('remote PFC Shape', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
subplot(5,1,4)
spy( y_pPFC_Orientation>-50,8,'b'),title('remote PFC Ori', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

subplot(5,1,5)
spy( y_FS_inh>-50,8,'r'),title('MD driven FS cells', 'FontSize', 16)
set(gca,'DataAspectRatio',[1000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
%%
close all

figure(11)
spy( y_VA_matrix_shape>-50,8,'k'),title('VA Thalamus Shape', 'FontSize', 16)
set(gca,'DataAspectRatio',[10000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')
hold on 
spy( y_MD_core_shape>-50,8,'r'),title('MD', 'FontSize', 16)
set(gca,'DataAspectRatio',[10000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')

hold on
spy( y_pPFC_Shape>-50,8,'b'),title('remote PFC Shape', 'FontSize', 16)
set(gca,'DataAspectRatio',[10000 1 1]),ylabel('Neuron ID', 'FontSize',16), xlim([0 t_final/dt]);%,xlabel('time')


%% ploting spikes
%%


%%
close all
figure(13)
srate=1000;
gauss_width= 100;
st= 10000;
final=800000;
firing_rate_timeseries= get_firing_num(y_pPFC_Shape, gauss_width, v_th);
firing_rate_timeseries=  conv_gaussian(firing_rate_timeseries, srate,gauss_width);


plot(firing_rate_timeseries, 'b'), xlim([st, final])

firing_rate_timeseries=[];
firing_rate_timeseries= get_firing_num(y_VA_matrix_shape, gauss_width, v_th);
firing_rate_timeseries=  conv_gaussian(firing_rate_timeseries, srate,gauss_width);

hold on

plot(firing_rate_timeseries, 'k'), xlim([st, final])

firing_rate_timeseries=[];
firing_rate_timeseries= get_firing_num(y_MD_core_shape, gauss_width, v_th);
firing_rate_timeseries=  conv_gaussian(firing_rate_timeseries, srate,gauss_width);

hold on

plot(firing_rate_timeseries, 'r'), xlim([st, final])


function firing = get_firing(current, threshold, nanvalue)
% current: membrane potiential
% window_size: threshold for firing
% thresh: predefined nan value
    firing = (current > threshold) & (current ~= nanvalue);
end
