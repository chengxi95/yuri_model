% use fixed initialization to generate certain number of trials data 
close all;

%%%%%%%%%%%%%%%%%% Simulation Setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numberofneurons = 50;% number of neurons per group
initialization = 0; % whether randomize the initialization (must load an initialization file if set to 0)
total_trial_num = 10; % total number of trials 
stimulation_type = 1; % Stimulation type 1:BT, 2:RC, 3:GC, 4:YT
con_stimulation_type = 2; 

% lesion test
VA_off = 0;
MD_off = 0;
pPFC_off = 0;
PV_off = 0;  
ST_off = 0;
SNPR_off = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tha = 20; % time constant 
 
% simulation time
dt = 0.01; %step size ms
t_final = 4500; % simulation time ms
T = 0:dt:t_final;
end_time = 4000/dt; % there is a delay between spike and corresponding current, only saves the results till end_time

%   Intrinsic property of neuron
delay = 5/dt; % 5ms
MD_delay = 10/dt; % 10ms
v_th = -50; %mv
E_L = -65; %mv
RM = 10;
leaky_coef = 1; %
tref = 1/dt; % 1 ms refratory time

% neurons talk to each other by sending short pulses, define duration and
% amplitude of these pulses
abs_pfc_width = 1.5; % cortical input liftime
spikewidth = abs_pfc_width/dt;
spikewidth_inh = 100/dt; % assumed that inhibition effect is prolonged due to the synch exc input SOM neurons
spikewidth_inh_FS= 1000/dt; 
spikewidth_VA = 100/dt;
spikewidth_MD = spikewidth;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Stimulation showing the abstract cue
no_of_trl = 1;
column_length = 10;
I_stim = zeros(numberofneurons,length(T));
t_start_stim_1_abs = 2000; % cue time
cue_amp = 0.5;
abs_stim_duration = 450; % stim duration is 450 ms
stim_duration = abs_stim_duration/dt;
t_start_stim_1 = t_start_stim_1_abs/dt;
for i=1:no_of_trl
    I_stim(1:column_length, t_start_stim_1: t_start_stim_1+ stim_duration) = cue_amp; % first stimulus delivered to visual starter neuron
    t_start_stim_1 = t_start_stim_1+ 1000/dt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Stimulation showing the concrete cue
no_of_trl = 1;
column_length = 10;
I_stim_con = zeros(numberofneurons,length(T));
t_start_stim_con = 2950; % cue time
cue_amp = 0.5;
abs_stim_duration = 450; %0.1;% stim duration is one ms
stim_duration = abs_stim_duration/dt;
t_start_stim_con = t_start_stim_con/dt;
for i=1:no_of_trl
    I_stim_con(1:column_length, t_start_stim_con: t_start_stim_con+ stim_duration) = cue_amp; % first stimulus delivered to visual starter neuron
    t_start_stim_con = t_start_stim_con + 1000/dt;
end

%%%%%%%%%%%%%%%%%% Connectivity Maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matrix_M= zeros(numberofneurons,numberofneurons); % connectivity from PFC middle to superficial
matrix_local = zeros(numberofneurons,numberofneurons); % connectivity within PFC superficial and deep
W_local = 0.002; % synaptic weights
W_M= 0.02;
for ii = 1:column_length: numberofneurons-column_length

    matrix_local(ii:ii+ column_length - 1,ii :ii+ column_length-1) = W_local;
    matrix_M(ii:ii+ column_length - 1, ii :ii+ column_length-1) = W_M;
end
matrix_local = matrix_local - diag(diag(matrix_local));
matrix_local(1:10, 11:20)= 2* matrix_local(1:10, 11:20);
matrix_local(41:50,1:10 )= W_local;
matrix_M(41:50,1:10 )= W_M;

% PFC superficial to PFC deep cells
W_S_D = 0.004;
matrix_local_S_D = zeros(numberofneurons,numberofneurons);
for ii = 1:column_length: numberofneurons-column_length+1
    matrix_local_S_D(ii:ii+ column_length - 1,ii: ii + column_length - 1) = W_S_D;
end
matrix_local_S_D(1:10, 1:10) = 2*W_S_D;

% weights from superficial layer to the deep layer across difference chains
W11=2; W12=0.8; W13=0.2; W14= 0.2;
W21=0.8; W22=2; W23=0.2; W24= 0.2;
W31=0.2; W32=0.2; W33=2; W34= 0.8;
W41=0.2; W42=0.2; W43=0.8; W44= 2;

% PFC deep cells to MD cells
W_PFC_MD= 0.0006;
Matrix_DPFC_to_MD = zeros(numberofneurons,numberofneurons);
Matrix_DPFC_to_MD(:, 1:15)=1 ; 

% PFC Deep to VA cells
W_PFC_TH =  0.015;
PFC_VA_matrix = zeros(numberofneurons,numberofneurons);
PFC_VA_matrix(1: 20, 1:numberofneurons) = 1;

% PFC Deep to STR
W_PFC_to_str= 0.0145;

Matrix_amplification = 2;
matrix_MD_to_Crtex = zeros(numberofneurons,numberofneurons); % modulatory effect of VA cells to PFC superficial cells and rule cells
matrix_VA_to_rule= zeros(numberofneurons,numberofneurons); % VA cells to PFC rule cells
matrix_VA_to_rule(1:50, 1:20)= 1;
matrix_MD_to_Crtex(1:10, 1:50) = 1;     

% MD cells to PFC superficial cells
W_MD_PFC = 0.012;
Matrix_random_MD_PFC_BT= zeros(numberofneurons,numberofneurons);
Matrix_random_MD_PFC_RC= zeros(numberofneurons,numberofneurons);
Matrix_random_MD_PFC_GC= zeros(numberofneurons,numberofneurons);
Matrix_random_MD_PFC_YT= zeros(numberofneurons,numberofneurons);

Matrix_random_MD_PFC_Bow = zeros(numberofneurons,numberofneurons);
Matrix_random_MD_PFC_Rec = zeros(numberofneurons,numberofneurons);
Matrix_random_MD_PFC_Oval = zeros(numberofneurons,numberofneurons);
Matrix_random_MD_PFC_SW = zeros(numberofneurons,numberofneurons);
Matrix_random_MD_PFC_N = zeros(numberofneurons,numberofneurons);
Matrix_random_MD_PFC_SE = zeros(numberofneurons,numberofneurons);

for i=1:5
    random_commections= randi(numberofneurons, 1, 1);
    Matrix_random_MD_PFC_BT (random_commections,random_commections)=1;
    random_commections= randi(numberofneurons, 1, 1);
    Matrix_random_MD_PFC_RC (random_commections,random_commections)=1;
    random_commections= randi(numberofneurons, 1, 1);
    Matrix_random_MD_PFC_GC (random_commections,random_commections)=1;
    random_commections= randi(numberofneurons, 1, 1);
    Matrix_random_MD_PFC_YT (random_commections,random_commections)=1;

    random_commections= randi(numberofneurons, 1, 1);
    Matrix_random_MD_PFC_Bow (random_commections,random_commections)=1;
    random_commections= randi(numberofneurons, 1, 1);
    Matrix_random_MD_PFC_Rec (random_commections,random_commections)=1;
    random_commections= randi(numberofneurons, 1, 1);
    Matrix_random_MD_PFC_Oval (random_commections,random_commections)=1;
    random_commections= randi(numberofneurons, 1, 1);
    Matrix_random_MD_PFC_SW (random_commections,random_commections)=1;
    random_commections= randi(numberofneurons, 1, 1);
    Matrix_random_MD_PFC_N (random_commections,random_commections)=1;
    random_commections= randi(numberofneurons, 1, 1);
    Matrix_random_MD_PFC_SE (random_commections,random_commections)=1;
end 

% MD cells to PFC rule cells
Matrix_MD_to_PFC= zeros(numberofneurons,numberofneurons);
Matrix_MD_to_PFC(:,1:25)=1;

% PFC rules cells to MD cells
Matrix_PFC_to_MD = zeros(numberofneurons,numberofneurons);
Matrix_PFC_to_MD (:, :)=1;                                      

% Excitatory and inhibitory cell connectivity = localy, no long distance inhibitory connection, inhibition only local within population
W_IPL_Inh_to_exc = 0.001;
W_IPL_Inh_to_exc_pv = 0.01;
W_IPL_Inh_to_exc_fs= 0.0025;

W_VA_exi = 0.01;
W_MD_inh = 0.001;

% placeholder for multi trials results (first neuron of the section are saved)
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

full_PFC_shape_ensemble = {};
full_PFC_ori_ensemble = {};

full_con_inh_shape = {};
full_con_inh_ori = {};
full_PFC_bow = {};
full_PFC_rec = {};
full_PFC_oval = {};
full_PFC_SW = {};
full_PFC_N = {};
full_PFC_SE = {};
full_MD_left = {};
full_MD_mid = {};
full_MD_right = {};
full_PFC_left = {};
full_PFC_mid = {};
full_PFC_right = {};
full_output_left = {};
full_output_mid = {};
full_output_right = {};

full_Isyn_FS_left= {};
full_Isyn_FS_mid= {};
full_Isyn_FS_right= {};
full_Isyn_MD_left_to_PFC = {};
full_Isyn_MD_mid_to_PFC = {};
full_Isyn_MD_right_to_PFC = {};


if ~initialization
    load("initialization_con_output.mat")
end

% % concrete rule
% Iext_PFC_M_con = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5);
% Iext_SPFC_Bow = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5);
% Iext_SPFC_Rec = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5);
% Iext_SPFC_Oval = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5);
% Iext_SPFC_SW = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5);
% Iext_SPFC_N = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5);
% Iext_SPFC_SE = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5);
% Iext_Inh_pPFC_con = 1.5 - 0.01 *(rand(numberofneurons,1)-0.5); 
% Iext_Inh_PV_con = 1.49 - 0.0 *(rand(numberofneurons,1));
% Iext_PFC_remote_left = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); 
% Iext_PFC_remote_right = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); 
% Iext_PFC_remote_mid = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5);
% Iext_PFC_output_left = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); 
% Iext_PFC_output_right = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); 
% Iext_PFC_output_mid = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); 

for rr= 1:total_trial_num
    disp(rr)
    if initialization
    
        % driving currents-----------------------------------------------------
        % similar driving current for PFC and IPL cells small variation across cells
        Iext_SPFC_YT = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); 
        Iext_SPFC_BT = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); 
        Iext_SPFC_RC = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); 
        Iext_SPFC_GC = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); 
        Iext_DPFC_YT = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); 
        Iext_DPFC_BT = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); 
        Iext_DPFC_RC = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); 
        Iext_DPFC_GC = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); 
        Iext_PFC_M = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5);
        Iext_PFC_remote = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5);               

        % same driving current for all thalamic cells
        Iext_MD = 1.5 - 0.001 *rand(numberofneurons,1); 
        Iext_VA = 1.5 - 0.001 *rand(numberofneurons,1); 
        Iext_Inh_pPFC = 1.5 - 0.01 *(rand(numberofneurons,1)-0.5); 
        Iext_Inh_PV = 1.49 - 0.0 *(rand(numberofneurons,1));

        % concrete rule
        Iext_PFC_M_con = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5);
        Iext_SPFC_Bow = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5);
        Iext_SPFC_Rec = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5);
        Iext_SPFC_Oval = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5);
        Iext_SPFC_SW = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5);
        Iext_SPFC_N = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5);
        Iext_SPFC_SE = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5);
        Iext_Inh_pPFC_con = 1.5 - 0.01 *(rand(numberofneurons,1)-0.5); 
        Iext_Inh_PV_con = 1.49 - 0.0 *(rand(numberofneurons,1));
        Iext_PFC_remote_left = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); 
        Iext_PFC_remote_right = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); 
        Iext_PFC_remote_mid = 1.4955 +  0.001 *(rand(numberofneurons,1)-0.5); 
    end

    % Initialize membrane potential---------------------------------------
    y_PFC_S_BT= zeros(numberofneurons,length(T));
    y_PFC_S_BT(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;
    y_PFC_S_RC = zeros(numberofneurons,length(T));
    y_PFC_S_RC(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;
    y_PFC_S_GC = zeros(numberofneurons,length(T));
    y_PFC_S_GC(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;
    y_PFC_S_YT = zeros(numberofneurons,length(T));
    y_PFC_S_YT(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;

    y_PFC_M = zeros(numberofneurons,length(T));
    y_PFC_M(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;

    y_aPFC_D_BT = zeros(numberofneurons,length(T));
    y_aPFC_D_BT(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;
    y_aPFC_D_RC= zeros(numberofneurons,length(T));
    y_aPFC_D_RC(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;
    y_aPFC_D_GC= zeros(numberofneurons,length(T));
    y_aPFC_D_GC(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;
    y_aPFC_D_YT= zeros(numberofneurons,length(T));
    y_aPFC_D_YT(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;

    y_PFC_Shape = zeros(numberofneurons,length(T));
    y_PFC_Shape(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;
    y_PFC_Orientation = zeros(numberofneurons,length(T));
    y_PFC_Orientation(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;

    y_VA_matrix_shape = zeros(numberofneurons,length(T));
    y_VA_matrix_shape(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1; 
    y_VA_matrix_Orientation  = zeros(numberofneurons,length(T));
    y_VA_matrix_Orientation(1:numberofneurons,1)  =-55+rand(numberofneurons,1)*1; 

    y_MD_core_ori = zeros(numberofneurons,length(T));
    y_MD_core_ori(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1;
    y_MD_core_shape = zeros(numberofneurons,length(T));
    y_MD_core_shape(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1; 

    y_SNpr_inh=zeros(numberofneurons,length(T));
    y_SNpr_inh(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;
    y_ST_inh=zeros(numberofneurons,length(T));
    y_ST_inh(1:numberofneurons,1)=-55+rand(numberofneurons,1)*0;
    y_PV_inh= zeros(numberofneurons,length(T));
    y_PV_inh (1:numberofneurons,1)=-55+rand(numberofneurons,1)*1;
    y_FS_inh_shape = zeros(numberofneurons,length(T));
    y_FS_inh_shape(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1;
    y_FS_inh_ori= zeros(numberofneurons,length(T));
    y_FS_inh_ori(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1;

    % concrete rule
    y_PFC_M_con = zeros(numberofneurons,length(T));
    y_PFC_M_con(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;
    y_PFC_S_Bow= zeros(numberofneurons,length(T));
    y_PFC_S_Bow(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;
    y_PFC_S_Rec= zeros(numberofneurons,length(T));
    y_PFC_S_Rec(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;
    y_PFC_S_Oval= zeros(numberofneurons,length(T));
    y_PFC_S_Oval(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;
    y_PFC_S_SW= zeros(numberofneurons,length(T));
    y_PFC_S_SW(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;
    y_PFC_S_N= zeros(numberofneurons,length(T));
    y_PFC_S_N(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;
    y_PFC_S_SE= zeros(numberofneurons,length(T));
    y_PFC_S_SE(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;
    
    y_VA_matrix_left = zeros(numberofneurons,length(T));
    y_VA_matrix_left(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1; 
    y_VA_matrix_mid = zeros(numberofneurons,length(T));
    y_VA_matrix_mid(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1; 
    y_VA_matrix_right = zeros(numberofneurons,length(T));
    y_VA_matrix_right(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1; 

    y_MD_core_left = zeros(numberofneurons,length(T));
    y_MD_core_left(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1; 
    y_MD_core_mid = zeros(numberofneurons,length(T));
    y_MD_core_mid(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1; 
    y_MD_core_right = zeros(numberofneurons,length(T));
    y_MD_core_right(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1; 

    y_SNpr_inh_con =zeros(numberofneurons,length(T));
    y_SNpr_inh_con(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;
    y_ST_inh_con =zeros(numberofneurons,length(T));
    y_ST_inh_con(1:numberofneurons,1)=-55+rand(numberofneurons,1)*0;% 

    y_pPFC_left = zeros(numberofneurons,length(T));
    y_pPFC_left(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;
    y_pPFC_mid = zeros(numberofneurons,length(T));
    y_pPFC_mid(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;
    y_pPFC_right = zeros(numberofneurons,length(T));
    y_pPFC_right(1:numberofneurons,1)=-55+rand(numberofneurons,1)*5;

    y_PV_inh_con = zeros(numberofneurons,length(T));
    y_PV_inh_con(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1;
    y_FS_inh_left= zeros(numberofneurons,length(T));
    y_FS_inh_left(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1;
    y_FS_inh_mid= zeros(numberofneurons,length(T));
    y_FS_inh_mid(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1;
    y_FS_inh_right= zeros(numberofneurons,length(T));
    y_FS_inh_right(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1;

    y_con_inh_shape= zeros(numberofneurons,length(T));
    y_con_inh_shape(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1;
    y_con_inh_ori= zeros(numberofneurons,length(T));
    y_con_inh_ori(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1;

    y_PFC_output_left = zeros(numberofneurons,length(T));
    y_PFC_output_left(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1;
    y_PFC_output_mid = zeros(numberofneurons,length(T));
    y_PFC_output_mid(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1;
    y_PFC_output_right = zeros(numberofneurons,length(T));
    y_PFC_output_right(1:numberofneurons,1)=-55+rand(numberofneurons,1)*1;

    % initialize spike counting--------------------------------------------
    last_spike_PFC_M = 10^10*ones(numberofneurons,1);

    last_spike_PFC_S_BT=10^10*ones(numberofneurons,1);
    last_spike_PFC_S_GC=10^10*ones(numberofneurons,1);
    last_spike_PFC_S_YT=10^10*ones(numberofneurons,1);
    last_spike_PFC_S_RC =10^10*ones(numberofneurons,1);

    last_spike_PFC_D_BT = 10^10*ones(numberofneurons,1);
    last_spike_PFC_D_YT= 10^10*ones(numberofneurons,1);
    last_spike_PFC_D_GC= 10^10*ones(numberofneurons,1);
    last_spike_PFC_D_RC = 10^10*ones(numberofneurons,1);
    
    last_spike_MD_shape = 10^10*ones(numberofneurons,1);
    last_spike_MD_ori = 10^10*ones(numberofneurons,1);

    last_spike_VA_matrix_shape = 10^10*ones(numberofneurons,1);
    last_spike_VA_matrix_Orientation= 10^10*ones(numberofneurons,1);

    last_spike_SNpr_Inh = 10^10*ones(numberofneurons,1);

    last_spike_PFC_remote_shape= 10^10*ones(numberofneurons,1);
    last_spike_PFC_remote_Orientation= 10^10*ones(numberofneurons,1);

    last_spike_PV= 10^10*ones(numberofneurons,1);
    last_spike_FS_shape = 10^10*ones(numberofneurons,1);
    last_spike_FS_ori = 10^10*ones(numberofneurons,1); 
    last_spike_ST = 10^10*ones(numberofneurons,1);

    % concrete rule 
    last_spike_PFC_S_Bow=10^10*ones(numberofneurons,1);
    last_spike_PFC_S_Rec=10^10*ones(numberofneurons,1);
    last_spike_PFC_S_Oval=10^10*ones(numberofneurons,1);
    last_spike_PFC_S_SW =10^10*ones(numberofneurons,1);
    last_spike_PFC_S_N =10^10*ones(numberofneurons,1);
    last_spike_PFC_S_SE =10^10*ones(numberofneurons,1);
    last_spike_PFC_M_con= 10^10*ones(numberofneurons,1);

    last_spike_PFC_remote_left= 10^10*ones(numberofneurons,1);
    last_spike_PFC_remote_mid= 10^10*ones(numberofneurons,1);
    last_spike_PFC_remote_right= 10^10*ones(numberofneurons,1);

    last_spike_ST_con = 10^10*ones(numberofneurons,1);
    last_spike_SNpr_Inh_con = 10^10*ones(numberofneurons,1);

    last_spike_VA_matrix_left = 10^10*ones(numberofneurons,1);
    last_spike_VA_matrix_mid= 10^10*ones(numberofneurons,1);
    last_spike_VA_matrix_right = 10^10*ones(numberofneurons,1);

    last_spike_MD_left = 10^10*ones(numberofneurons,1);
    last_spike_MD_mid = 10^10*ones(numberofneurons,1);
    last_spike_MD_right = 10^10*ones(numberofneurons,1);

    last_spike_FS_left = 10^10*ones(numberofneurons,1);
    last_spike_FS_mid = 10^10*ones(numberofneurons,1);
    last_spike_FS_right = 10^10*ones(numberofneurons,1);

    last_spike_PV_con = 10^10*ones(numberofneurons,1);

    last_spike_abs_inh_shape = 10^10*ones(numberofneurons,1);
    last_spike_abs_inh_ori = 10^10*ones(numberofneurons,1);

    last_spike_PFC_output_left = 10^10*ones(numberofneurons,1);
    last_spike_PFC_output_mid = 10^10*ones(numberofneurons,1);
    last_spike_PFC_output_right = 10^10*ones(numberofneurons,1);

    % initialize synaptic input current------------------------------------
    Isyn_aPFC_local_shared = zeros(numberofneurons,length(T));
    Isyn_aPFC_local_BT= zeros(numberofneurons,length(T));
    Isyn_aPFC_local_RC= zeros(numberofneurons,length(T));
    Isyn_aPFC_local_GC= zeros(numberofneurons,length(T));
    Isyn_aPFC_local_YT= zeros(numberofneurons,length(T));

    Isyn_PFC_M_S= zeros(numberofneurons,length(T));
    Isyn_PFC_M_BT= zeros(numberofneurons,length(T));
    Isyn_PFC_M_RC= zeros(numberofneurons,length(T));
    Isyn_PFC_M_YT= zeros(numberofneurons,length(T));
    Isyn_PFC_M_GC= zeros(numberofneurons,length(T));

    Isyn_aPFC_D_shared = zeros(numberofneurons,length(T));
    Isyn_aPFC_to_D_RC= zeros(numberofneurons,length(T));
    Isyn_aPFC_to_D_GC= zeros(numberofneurons,length(T));
    Isyn_aPFC_to_D_YT= zeros(numberofneurons,length(T));
    Isyn_aPFC_to_D_BT= zeros(numberofneurons,length(T));
    Isyn_aPFC_D5_shape = zeros(numberofneurons,length(T));
    Isyn_aPFC_D5_orientation = zeros(numberofneurons,length(T));

    Isyn_PFC_D_MD_ori = zeros(numberofneurons,length(T));
    Isyn_PFC_D_MD_shape = zeros(numberofneurons,length(T));

    Isyn_SNpr_to_VA = zeros(numberofneurons,length(T));
   
    Isyn_MD_shape_to_PFC= zeros(numberofneurons,length(T));
    Isyn_MD_ori_to_PFC= zeros(numberofneurons,length(T));

    Isyn_VA_Matrix_to_PFC = ones(numberofneurons,length(T));
    Isyn_VA_Matrix_to_PFC_exc_shape= zeros(numberofneurons,length(T));
    Isyn_VA_Matrix_to_PFC_exc_orientation= zeros(numberofneurons,length(T));

    Isyn_ST = zeros(numberofneurons,length(T));
    Isyn_PV = zeros(numberofneurons,length(T));

    Isyn_MD_shape_to_Inh = zeros(numberofneurons,length(T));
    Isyn_MD_ori_to_Inh = zeros(numberofneurons,length(T));

    Isyn_FS_shape= zeros(numberofneurons,length(T));
    Isyn_FS_ori= zeros(numberofneurons,length(T));

    Isyn_aPFC_D5_to_ST = zeros(numberofneurons,length(T));
    Isyn_VA_Matrix_to_Inh= zeros(numberofneurons,length(T));

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

    % Concrete Rule
    Isyn_PFC_M_S_con= zeros(numberofneurons,length(T));
    Isyn_PFC_M_Bow = zeros(numberofneurons,length(T));
    Isyn_PFC_M_Rec = zeros(numberofneurons,length(T));
    Isyn_PFC_M_Oval = zeros(numberofneurons,length(T));
    Isyn_PFC_M_SW = zeros(numberofneurons,length(T));
    Isyn_PFC_M_N = zeros(numberofneurons,length(T));
    Isyn_PFC_M_SE = zeros(numberofneurons,length(T));

    Isyn_aPFC_D_left = zeros(numberofneurons,length(T));
    Isyn_aPFC_D_mid = zeros(numberofneurons,length(T));
    Isyn_aPFC_D_right = zeros(numberofneurons,length(T));

    Isyn_PFC_D_MD_left = zeros(numberofneurons,length(T));
    Isyn_PFC_D_MD_mid = zeros(numberofneurons,length(T));
    Isyn_PFC_D_MD_right = zeros(numberofneurons,length(T));
    
    Isyn_aPFC_D_con_to_ST = zeros(numberofneurons,length(T));
    Isyn_ST_con = zeros(numberofneurons,length(T));
    Isyn_SNpr_to_VA_con = zeros(numberofneurons,length(T));

    Isyn_VA_Matrix_to_PFC_con = ones(numberofneurons,length(T));
    Isyn_VA_Matrix_to_PFC_exc_left = zeros(numberofneurons,length(T));
    Isyn_VA_Matrix_to_PFC_exc_mid = zeros(numberofneurons,length(T));
    Isyn_VA_Matrix_to_PFC_exc_right = zeros(numberofneurons,length(T));
    Isyn_VA_Matrix_to_Inh_con = zeros(numberofneurons,length(T));

    Isyn_MD_left_to_PFC= zeros(numberofneurons,length(T));
    Isyn_MD_mid_to_PFC= zeros(numberofneurons,length(T));
    Isyn_MD_right_to_PFC= zeros(numberofneurons,length(T));
    
    Isyn_MD_left_to_random_PFC_Bow= zeros(numberofneurons,length(T));
    Isyn_MD_left_to_random_PFC_Rec= zeros(numberofneurons,length(T));
    Isyn_MD_left_to_random_PFC_Oval= zeros(numberofneurons,length(T));
    Isyn_MD_left_to_random_PFC_SW= zeros(numberofneurons,length(T));
    Isyn_MD_left_to_random_PFC_N= zeros(numberofneurons,length(T));
    Isyn_MD_left_to_random_PFC_SE= zeros(numberofneurons,length(T));

    Isyn_MD_mid_to_random_PFC_Bow= zeros(numberofneurons,length(T));
    Isyn_MD_mid_to_random_PFC_Rec= zeros(numberofneurons,length(T));
    Isyn_MD_mid_to_random_PFC_Oval= zeros(numberofneurons,length(T));
    Isyn_MD_mid_to_random_PFC_SW= zeros(numberofneurons,length(T));
    Isyn_MD_mid_to_random_PFC_N= zeros(numberofneurons,length(T));
    Isyn_MD_mid_to_random_PFC_SE= zeros(numberofneurons,length(T));

    Isyn_MD_right_to_random_PFC_Bow= zeros(numberofneurons,length(T));
    Isyn_MD_right_to_random_PFC_Rec= zeros(numberofneurons,length(T));
    Isyn_MD_right_to_random_PFC_Oval= zeros(numberofneurons,length(T));
    Isyn_MD_right_to_random_PFC_SW= zeros(numberofneurons,length(T));
    Isyn_MD_right_to_random_PFC_N= zeros(numberofneurons,length(T));
    Isyn_MD_right_to_random_PFC_SE= zeros(numberofneurons,length(T));

    Isyn_MD_left_to_Inh = zeros(numberofneurons,length(T));
    Isyn_MD_mid_to_Inh = zeros(numberofneurons,length(T));
    Isyn_MD_right_to_Inh = zeros(numberofneurons,length(T));
    Isyn_FS_left = zeros(numberofneurons,length(T));
    Isyn_FS_mid = zeros(numberofneurons,length(T));
    Isyn_FS_right = zeros(numberofneurons,length(T));

    Isyn_aPFC_local_left = zeros(numberofneurons,length(T));
    Isyn_aPFC_local_mid = zeros(numberofneurons,length(T));
    Isyn_aPFC_local_right = zeros(numberofneurons,length(T));

    Isyn_aPFC_local_shared_con = zeros(numberofneurons,length(T));
    Isyn_aPFC_local_Bow= zeros(numberofneurons,length(T));
    Isyn_aPFC_local_Rec= zeros(numberofneurons,length(T));
    Isyn_aPFC_local_Oval= zeros(numberofneurons,length(T));
    Isyn_aPFC_local_SW= zeros(numberofneurons,length(T));
    Isyn_aPFC_local_N= zeros(numberofneurons,length(T));
    Isyn_aPFC_local_SE= zeros(numberofneurons,length(T));

    Isyn_PV_con= zeros(numberofneurons,length(T));

    Isyn_abs_inh_cell_shape = zeros(numberofneurons,length(T)); % driving current from pfc rule ensemble to inhibitory cell
    Isyn_abs_inh_cell_ori = zeros(numberofneurons,length(T)); % driving current from pfc rule ensemble to inhibitory cell
    Isyn_abs_cell_shape_inh_con = zeros(numberofneurons,length(T)); % inhibitory current from inhibitory cell to pfc concrete cell
    Isyn_abs_cell_ori_inh_con = zeros(numberofneurons,length(T)); % inhibitory current from inhibitory cell to pfc concrete cell

    Isyn_PFC_output_left = zeros(numberofneurons,length(T));
    Isyn_PFC_output_mid = zeros(numberofneurons,length(T));
    Isyn_PFC_output_right = zeros(numberofneurons,length(T));

    % simulation starts here
    for i= 2:length(T)- MD_delay*10

       %%%%%% Midle layer (l = 1) to superficial layer (l = 2) %%%%%%%%%%%%

        Isyn_PFC_M_S(:,i+ delay) = Isyn_PFC_M_S(:,i+ delay) + sum(((i-last_spike_PFC_M < spikewidth) & (i-last_spike_PFC_M > 0)) .* matrix_M, 1).';
            
        if stimulation_type == 1
            Isyn_PFC_M_BT(:, i+delay) = Isyn_PFC_M_S(:,i+ delay);
        elseif stimulation_type == 2
            Isyn_PFC_M_RC(:, i+delay) = Isyn_PFC_M_S(:,i+ delay);
        elseif stimulation_type == 3
            Isyn_PFC_M_GC(:, i+delay) = Isyn_PFC_M_S(:,i+ delay);
        elseif stimulation_type == 4
            Isyn_PFC_M_YT(:, i+delay) = Isyn_PFC_M_S(:,i+ delay);
        else
            disp("invalid stimulation type")
            return
        end

        %%%%%%%%% Superficial layer (l = 2) to superficial layer (l = 2)
        %%%%%%%%% excitatory to excitatory  within chains

        w_within_chains = 0.6;
        Isyn_aPFC_local_BT(:,i+ delay) = Isyn_aPFC_local_BT(:,i+ delay) + w_within_chains .* sum(((i-last_spike_PFC_S_BT<spikewidth) & (i-last_spike_PFC_S_BT>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_RC(:,i+ delay) = Isyn_aPFC_local_RC(:,i+ delay) + w_within_chains .* sum(((i-last_spike_PFC_S_RC<spikewidth) & (i-last_spike_PFC_S_RC>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_YT(:,i+ delay) = Isyn_aPFC_local_YT(:,i+ delay) + w_within_chains .* sum(((i-last_spike_PFC_S_YT<spikewidth) & (i-last_spike_PFC_S_YT>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_GC(:,i+ delay) = Isyn_aPFC_local_GC(:,i+ delay) + w_within_chains .* sum(((i-last_spike_PFC_S_GC<spikewidth) & (i-last_spike_PFC_S_GC>0)) .* matrix_local, 1).';

        %%%%%%%%%% Superficial layer (l = 2) to superficial layer (l = 2)
        %%%%%%%%%%% excitatory to excitatory across chains

        w_across_chains= 0.45;
        
        Isyn_aPFC_local_shared(:,i+ delay) = Isyn_aPFC_local_shared(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_S_BT<spikewidth) & (i-last_spike_PFC_S_BT>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_shared(:,i+ delay) = Isyn_aPFC_local_shared(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_S_RC<spikewidth) & (i-last_spike_PFC_S_RC>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_shared(:,i+ delay) = Isyn_aPFC_local_shared(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_S_YT<spikewidth) & (i-last_spike_PFC_S_YT>0)) .* matrix_local, 1).';
        
        Isyn_aPFC_local_shared(:,i+ delay) = Isyn_aPFC_local_shared(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_S_GC<spikewidth) & (i-last_spike_PFC_S_GC>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_shared(:,i+ delay) = Isyn_aPFC_local_shared(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_remote_shape<spikewidth) & (i-last_spike_PFC_remote_shape>0)) .* matrix_local, 1).';
        
        Isyn_aPFC_local_shared(:,i+ delay) = Isyn_aPFC_local_shared(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_remote_Orientation<spikewidth) & (i-last_spike_PFC_remote_Orientation>0)) .* matrix_local, 1).';


        %%%%%%%% Superficial layer (l = 2) to deep layer (l = 3) %%%%%%%%%%

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
        
        %%%%%%%%%%%%%%%% Deep layer (l=3) to VA thalamus (l=4) %%%%%%%%%%%%
        %%%%  SHAPE

        Isyn_aPFC_D5_shape(:,i+ delay) = Isyn_aPFC_D5_shape(:,i+ delay) +  W_PFC_TH .* sum(((i-last_spike_PFC_D_BT<spikewidth) & (i-last_spike_PFC_D_BT>0) & (PFC_VA_matrix ~= 0)), 1).';
        
        Isyn_aPFC_D5_shape(:,i+ delay) = Isyn_aPFC_D5_shape(:,i+ delay) +  W_PFC_TH .* sum(((i-last_spike_PFC_D_RC<spikewidth) & (i-last_spike_PFC_D_RC>0) & (PFC_VA_matrix ~= 0)), 1).';

        %%%%  ORIENTATION
        
        Isyn_aPFC_D5_orientation(:,i+ delay) = Isyn_aPFC_D5_orientation(:,i+ delay) +  W_PFC_TH .* sum((i-last_spike_PFC_D_GC<spikewidth) & (i-last_spike_PFC_D_GC>0) & (PFC_VA_matrix ~= 0), 1).';

        Isyn_aPFC_D5_orientation(:,i+ delay) = Isyn_aPFC_D5_orientation(:,i+ delay) +  W_PFC_TH .* sum((i-last_spike_PFC_D_YT<spikewidth) & (i-last_spike_PFC_D_YT>0) & (PFC_VA_matrix ~= 0), 1).';

        %%%%%%%%%%%%%%%%%%%%%% PFC deep  to MD %%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        %%%%%%%%%%%%%%%%%%% VA thalamus to cortex %%%%%%%%%%%%%%%%%%%%%%%%%
    
        % VA to excitatory both modulatory and driving 
        
        Isyn_VA_Matrix_to_PFC(sum(((i-last_spike_VA_matrix_shape<spikewidth_VA) & (i-last_spike_VA_matrix_shape>0)) .* matrix_MD_to_Crtex, 1)~=0,i+ delay) = Matrix_amplification; % modulatory efect of VA that goes to superficial + rule selective ensembles shape 

        Isyn_VA_Matrix_to_PFC(sum(((i-last_spike_VA_matrix_Orientation<spikewidth_VA) & (i-last_spike_VA_matrix_Orientation>0)) .* matrix_MD_to_Crtex, 1)~=0,i+ delay) = Matrix_amplification;% modulatory efect of VA that goes to superficial + rule selective ensembles orientation 
        
        Isyn_VA_Matrix_to_PFC_exc_shape(:,i+ delay) = Isyn_VA_Matrix_to_PFC_exc_shape(:,i+ delay) +  W_VA_exi .* sum(((i-last_spike_VA_matrix_shape<spikewidth) & (i-last_spike_VA_matrix_shape>0) & (matrix_VA_to_rule ~= 0)), 1).';% driving current from VA to rule selective cells (shape) 
        
        Isyn_VA_Matrix_to_PFC_exc_orientation(:,i+ delay) = Isyn_VA_Matrix_to_PFC_exc_orientation(:,i+ delay) +  W_VA_exi .* sum(((i-last_spike_VA_matrix_Orientation<spikewidth) & (i-last_spike_VA_matrix_Orientation>0) & (matrix_VA_to_rule ~= 0)), 1).';% driving current from VA to rule selective cells (orientation)

        % VA to Inh 

        Isyn_VA_Matrix_to_Inh(:,i+ delay) = Isyn_VA_Matrix_to_Inh(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_VA_matrix_shape<spikewidth) & (i-last_spike_VA_matrix_shape>0), 1).';
        Isyn_VA_Matrix_to_Inh(:,i+ delay) = Isyn_VA_Matrix_to_Inh(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_VA_matrix_Orientation<spikewidth) & (i-last_spike_VA_matrix_Orientation>0), 1).';

        %%%%%%%%%%%%%%%%%%%%%%% MD to PFC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %SHAPE

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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PFC deep  to ST %%%%%%%%%%%%%%%%%%%%

        Isyn_aPFC_D5_to_ST(:,i+ delay) = Isyn_aPFC_D5_to_ST(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_PFC_D_BT<spikewidth) & (i-last_spike_PFC_D_BT>0), 1).';

        Isyn_aPFC_D5_to_ST(:,i+ delay) = Isyn_aPFC_D5_to_ST(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_PFC_D_YT<spikewidth) & (i-last_spike_PFC_D_YT>0), 1).';

        Isyn_aPFC_D5_to_ST(:,i+ delay) = Isyn_aPFC_D5_to_ST(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_PFC_D_RC<spikewidth) & (i-last_spike_PFC_D_RC>0), 1).';

        Isyn_aPFC_D5_to_ST(:,i+ delay) = Isyn_aPFC_D5_to_ST(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_PFC_D_GC<spikewidth) & (i-last_spike_PFC_D_GC>0), 1).';
       
        %%%%%%%%%%%%%%%%%%%%% ST to SNpr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Isyn_ST(:, i + delay) =  Isyn_ST(:, i + delay) + W_IPL_Inh_to_exc.*sum(((i-last_spike_ST < spikewidth_inh) & (i-last_spike_ST > 0)).*(exp(((i-last_spike_ST < spikewidth_inh) & (i-last_spike_ST > 0)).*(1-(i-last_spike_ST-delay)).*2./spikewidth_inh)), 1);

        %%%%%%%%%%%%%%%%%%%%% SNpr to VA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %   IPL Inh neurons to excitatory

        Isyn_SNpr_to_VA(:, i + delay) =  Isyn_SNpr_to_VA(:, i + delay) + W_IPL_Inh_to_exc.*sum(((i-last_spike_SNpr_Inh < spikewidth_inh) & (i-last_spike_SNpr_Inh > 0)).*(exp(((i-last_spike_SNpr_Inh < spikewidth_inh) & (i-last_spike_SNpr_Inh > 0)).*(1-(i-last_spike_SNpr_Inh-delay)).*2./spikewidth_inh)), 1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cortical inhibitory cells 
        %   Va cells drive general inhibition of PV cells 
        if ~PV_off
            Isyn_PV(:, i + delay) =  Isyn_PV(:, i + delay) + W_IPL_Inh_to_exc_pv.*sum(((i-last_spike_PV < spikewidth_inh) & (i-last_spike_PV > 0)).*(exp(((i-last_spike_PV < spikewidth_inh) & (i-last_spike_PV > 0)).*(1-(i-last_spike_PV-delay)).*2./spikewidth_inh)), 1);
        end

        Isyn_MD_shape_to_Inh(:,i+ delay) = Isyn_MD_shape_to_Inh(:,i+ delay) + W_MD_inh .* sum(((i-last_spike_MD_shape<spikewidth_VA) & (i-last_spike_MD_shape>0)), 1).';

        Isyn_MD_ori_to_Inh(:,i+ delay) = Isyn_MD_ori_to_Inh(:,i+ delay) + W_MD_inh .* sum(((i-last_spike_MD_ori<spikewidth_VA) & (i-last_spike_MD_ori>0)), 1).';

        Isyn_FS_shape(:, i + delay) =  Isyn_FS_shape(:, i + delay) + W_IPL_Inh_to_exc_fs.*sum(((i-last_spike_FS_shape < spikewidth_inh_FS) & (i-last_spike_FS_shape > 0)).*(exp(((i-last_spike_FS_shape < spikewidth_inh_FS) & (i-last_spike_FS_shape > 0)).*(1-(i-last_spike_FS_shape-delay)).*2./spikewidth_inh_FS)), 1);

        Isyn_FS_ori(:, i + delay) =  Isyn_FS_ori(:, i + delay) + W_IPL_Inh_to_exc_fs.*sum(((i-last_spike_FS_ori < spikewidth_inh_FS) & (i-last_spike_FS_ori > 0)).*(exp(((i-last_spike_FS_ori < spikewidth_inh_FS) & (i-last_spike_FS_ori > 0)).*(1-(i-last_spike_FS_ori-delay)).*2./spikewidth_inh_FS)), 1);
        
        %%%%%%%%%%% local excitatory input to shape cells %%%%%%%%%%%%%%%%%
                
        Isyn_aPFC_local_shape(:,i+ delay) = Isyn_aPFC_local_shape(:,i+ delay) + sum(((i-last_spike_PFC_remote_shape<spikewidth) & (i-last_spike_PFC_remote_shape>0)) .* matrix_local, 1).';
    
        Isyn_aPFC_local_ori(:,i+ delay) = Isyn_aPFC_local_ori(:,i+ delay) + sum(((i-last_spike_PFC_remote_Orientation<spikewidth) & (i-last_spike_PFC_remote_Orientation>0)) .* matrix_local, 1).';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Concrete Rule Midle layer (l = 1) to superficial layer (l = 2)

        W_abstract_inh = 0.005;
        abs_strength = 0.5;
        spikewidth_inh_abs = 1000/dt;

        Isyn_abs_inh_cell_shape(:, i + delay) = Isyn_abs_inh_cell_shape(:, i + delay) +  W_abstract_inh .* sum(((i-last_spike_PFC_remote_shape<spikewidth_VA) & (i-last_spike_PFC_remote_shape>0)), 1).';
        
        Isyn_abs_inh_cell_ori(:, i + delay) = Isyn_abs_inh_cell_ori(:, i + delay) +  W_abstract_inh .* sum(((i-last_spike_PFC_remote_Orientation<spikewidth_VA) & (i-last_spike_PFC_remote_Orientation>0)), 1).';
        
        Isyn_abs_cell_shape_inh_con(:, i + delay) =  Isyn_abs_cell_shape_inh_con(:, i + delay) + abs_strength * W_abstract_inh.*sum(((i-last_spike_abs_inh_shape < spikewidth_inh_abs) & (i-last_spike_abs_inh_shape > 0)).*(exp(((i-last_spike_abs_inh_shape < spikewidth_inh_abs) & (i-last_spike_abs_inh_shape > 0)).*(1-(i-last_spike_abs_inh_shape-delay)).*2./spikewidth_inh_abs)), 1);

        Isyn_abs_cell_ori_inh_con(:, i + delay) =  Isyn_abs_cell_ori_inh_con(:, i + delay) + abs_strength * W_abstract_inh.*sum(((i-last_spike_abs_inh_ori < spikewidth_inh_abs) & (i-last_spike_abs_inh_ori > 0)).*(exp(((i-last_spike_abs_inh_ori < spikewidth_inh_abs) & (i-last_spike_abs_inh_ori > 0)).*(1-(i-last_spike_abs_inh_ori-delay)).*2./spikewidth_inh_abs)), 1);        

        Isyn_PFC_M_S_con(:,i+ delay) = Isyn_PFC_M_S_con(:,i+ delay) + sum(((i-last_spike_PFC_M_con < spikewidth) & (i-last_spike_PFC_M_con > 0)) .* matrix_M, 1).';
        
        % SW Bow
        if con_stimulation_type == 1
            Isyn_PFC_M_SW(:, i+delay) = Isyn_PFC_M_S_con(:,i+ delay);
            Isyn_PFC_M_Bow(:, i+delay) = Isyn_PFC_M_S_con(:,i+ delay);
        % N Bow
        elseif con_stimulation_type == 2
            Isyn_PFC_M_N(:, i+delay) = Isyn_PFC_M_S_con(:,i+ delay);
            Isyn_PFC_M_Bow(:, i+delay) = Isyn_PFC_M_S_con(:,i+ delay);
        % SE Bow
        elseif con_stimulation_type == 3
            Isyn_PFC_M_SE(:, i+delay) = Isyn_PFC_M_S_con(:,i+ delay);
            Isyn_PFC_M_Bow(:, i+delay) = Isyn_PFC_M_S_con(:,i+ delay);
        % SW Rec
        elseif con_stimulation_type == 4
            Isyn_PFC_M_SW(:, i+delay) = Isyn_PFC_M_con(:,i+ delay);
            Isyn_PFC_M_Rec(:, i+delay) = Isyn_PFC_M_con(:,i+ delay);
        % N Rec
        elseif con_stimulation_type == 5
            Isyn_PFC_M_N(:, i+delay) = Isyn_PFC_M_con(:,i+ delay);
            Isyn_PFC_M_Rec(:, i+delay) = Isyn_PFC_M_con(:,i+ delay);
        % SE Rec
        elseif con_stimulation_type == 6
            Isyn_PFC_M_SE(:, i+delay) = Isyn_PFC_M_con(:,i+ delay);
            Isyn_PFC_M_Rec(:, i+delay) = Isyn_PFC_M_con(:,i+ delay);
        % SW Oval
        elseif con_stimulation_type == 7
            Isyn_PFC_M_SW(:, i+delay) = Isyn_PFC_M_con(:,i+ delay);
            Isyn_PFC_M_Oval(:, i+delay) = Isyn_PFC_M_con(:,i+ delay);
        % N Oval
        elseif con_stimulation_type == 8
            Isyn_PFC_M_N(:, i+delay) = Isyn_PFC_M_con(:,i+ delay);
            Isyn_PFC_M_Oval(:, i+delay) = Isyn_PFC_M_con(:,i+ delay);
        % SE Oval
        elseif con_stimulation_type == 9
            Isyn_PFC_M_SE(:, i+delay) = Isyn_PFC_M_con(:,i+ delay);
            Isyn_PFC_M_Oval(:, i+delay) = Isyn_PFC_M_con(:,i+ delay);
        else
            disp("invalid stimulation type")
            return
        end
        

        %%%%%%%%% Concrete rule Superficial layer (l = 2) to superficial layer (l = 2)
        %%%%%%%%% excitatory to excitatory  within chains

        w_within_chains = 0.6;
        Isyn_aPFC_local_Bow(:,i+ delay) = Isyn_aPFC_local_Bow(:,i+ delay) + w_within_chains .* sum(((i-last_spike_PFC_S_Bow<spikewidth) & (i-last_spike_PFC_S_Bow>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_Rec(:,i+ delay) = Isyn_aPFC_local_Rec(:,i+ delay) + w_within_chains .* sum(((i-last_spike_PFC_S_Rec<spikewidth) & (i-last_spike_PFC_S_Rec>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_Oval(:,i+ delay) = Isyn_aPFC_local_Oval(:,i+ delay) + w_within_chains .* sum(((i-last_spike_PFC_S_Oval<spikewidth) & (i-last_spike_PFC_S_Oval>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_SW(:,i+ delay) = Isyn_aPFC_local_SW(:,i+ delay) + w_within_chains .* sum(((i-last_spike_PFC_S_SW<spikewidth) & (i-last_spike_PFC_S_SW>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_N(:,i+ delay) = Isyn_aPFC_local_N(:,i+ delay) + w_within_chains .* sum(((i-last_spike_PFC_S_N<spikewidth) & (i-last_spike_PFC_S_N>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_SE(:,i+ delay) = Isyn_aPFC_local_SE(:,i+ delay) + w_within_chains .* sum(((i-last_spike_PFC_S_SE<spikewidth) & (i-last_spike_PFC_S_SE>0)) .* matrix_local, 1).';

        %%%%%%%%%%%Superficial layer (l = 2) to superficial layer (l = 2)
        %%%%%%%%%%% excitatory to excitatory across chains

        w_across_chains= 0.45;
        
        Isyn_aPFC_local_shared_con(:,i+ delay) = Isyn_aPFC_local_shared_con(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_S_Bow<spikewidth) & (i-last_spike_PFC_S_Bow>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_shared_con(:,i+ delay) = Isyn_aPFC_local_shared_con(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_S_Rec<spikewidth) & (i-last_spike_PFC_S_Rec>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_shared_con(:,i+ delay) = Isyn_aPFC_local_shared_con(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_S_Oval<spikewidth) & (i-last_spike_PFC_S_Oval>0)) .* matrix_local, 1).';
        
        Isyn_aPFC_local_shared_con(:,i+ delay) = Isyn_aPFC_local_shared_con(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_S_SW<spikewidth) & (i-last_spike_PFC_S_SW>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_shared_con(:,i+ delay) = Isyn_aPFC_local_shared_con(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_S_N<spikewidth) & (i-last_spike_PFC_S_N>0)) .* matrix_local, 1).';
        
        Isyn_aPFC_local_shared_con(:,i+ delay) = Isyn_aPFC_local_shared_con(:,i+ delay) + w_across_chains .* sum(((i-last_spike_PFC_S_SE<spikewidth) & (i-last_spike_PFC_S_SE>0)) .* matrix_local, 1).';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Concrete Rule Superficial layer (l=2) to VA thalamus (l=4) 
        %%%%  Left

        Isyn_aPFC_D_left(:,i+ delay) = Isyn_aPFC_D_left(:,i+ delay) +  W_PFC_TH .* sum(((i-last_spike_PFC_S_Bow<spikewidth) & (i-last_spike_PFC_S_Bow>0) & (PFC_VA_matrix ~= 0)), 1).';
        
        Isyn_aPFC_D_left(:,i+ delay) = Isyn_aPFC_D_left(:,i+ delay) +  W_PFC_TH .* sum(((i-last_spike_PFC_S_SW<spikewidth) & (i-last_spike_PFC_S_SW>0) & (PFC_VA_matrix ~= 0)), 1).';

        %%%%  Middle

        Isyn_aPFC_D_mid(:,i+ delay) = Isyn_aPFC_D_mid(:,i+ delay) +  W_PFC_TH .* sum(((i-last_spike_PFC_S_Rec<spikewidth) & (i-last_spike_PFC_S_Rec>0) & (PFC_VA_matrix ~= 0)), 1).';
        
        Isyn_aPFC_D_mid(:,i+ delay) = Isyn_aPFC_D_mid(:,i+ delay) +  W_PFC_TH .* sum(((i-last_spike_PFC_S_N<spikewidth) & (i-last_spike_PFC_S_N>0) & (PFC_VA_matrix ~= 0)), 1).';

        %%%%  Right

        Isyn_aPFC_D_right(:,i+ delay) = Isyn_aPFC_D_right(:,i+ delay) +  W_PFC_TH .* sum(((i-last_spike_PFC_S_Oval<spikewidth) & (i-last_spike_PFC_S_Oval>0) & (PFC_VA_matrix ~= 0)), 1).';
        
        Isyn_aPFC_D_right(:,i+ delay) = Isyn_aPFC_D_right(:,i+ delay) +  W_PFC_TH .* sum(((i-last_spike_PFC_S_SE<spikewidth) & (i-last_spike_PFC_S_SE>0) & (PFC_VA_matrix ~= 0)), 1).';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Concrete Rule PFC deep/superficial to MD 
        
        % Left
        PFC_MD_initial_coef= 5;
        PFC_MD_secondary_coef= 1.5;
        Isyn_PFC_D_MD_left(:,i+ MD_delay) = Isyn_PFC_D_MD_left(:,i+MD_delay) + PFC_MD_secondary_coef * W_PFC_MD .* sum(((i-last_spike_PFC_S_Bow<spikewidth) & (i-last_spike_PFC_S_Bow>0) & (Matrix_PFC_to_MD ~= 0)), 1).';
        
        Isyn_PFC_D_MD_left(:,i+ MD_delay) = Isyn_PFC_D_MD_left(:,i+ MD_delay) + PFC_MD_secondary_coef * W_PFC_MD .* sum(((i-last_spike_PFC_S_SW<spikewidth) & (i-last_spike_PFC_S_SW>0) & (Matrix_PFC_to_MD ~= 0)), 1).';

        Isyn_PFC_D_MD_left(:,i+ 10*MD_delay) = Isyn_PFC_D_MD_left(:,i+ 10*MD_delay) + PFC_MD_secondary_coef * W_PFC_MD .* sum(((i-last_spike_PFC_remote_left<spikewidth) & (i-last_spike_PFC_remote_left>0) & (Matrix_PFC_to_MD ~= 0)), 1).';

        % Mid

        Isyn_PFC_D_MD_mid(:,i+ MD_delay) = Isyn_PFC_D_MD_mid(:,i+ MD_delay) + PFC_MD_secondary_coef * W_PFC_MD .* sum(((i-last_spike_PFC_S_Rec<spikewidth) & (i-last_spike_PFC_S_Rec>0) & (Matrix_PFC_to_MD ~= 0)), 1).';
        
        Isyn_PFC_D_MD_mid(:,i+ MD_delay) = Isyn_PFC_D_MD_mid(:,i+ MD_delay) + PFC_MD_secondary_coef * W_PFC_MD .* sum(((i-last_spike_PFC_S_N<spikewidth) & (i-last_spike_PFC_S_N>0) & (Matrix_PFC_to_MD ~= 0)), 1).';

        Isyn_PFC_D_MD_mid(:,i+ 10*MD_delay) = Isyn_PFC_D_MD_mid(:,i+ 10*MD_delay) + PFC_MD_secondary_coef * W_PFC_MD .* sum(((i-last_spike_PFC_remote_mid<spikewidth) & (i-last_spike_PFC_remote_mid>0) & (Matrix_PFC_to_MD ~= 0)), 1).';

        % Right

        Isyn_PFC_D_MD_right(:,i+ MD_delay) = Isyn_PFC_D_MD_right(:,i+ MD_delay) + PFC_MD_secondary_coef * W_PFC_MD .* sum(((i-last_spike_PFC_S_Oval<spikewidth) & (i-last_spike_PFC_S_Oval>0) & (Matrix_PFC_to_MD ~= 0)), 1).';

        Isyn_PFC_D_MD_right(:,i+ MD_delay) = Isyn_PFC_D_MD_right(:,i+ MD_delay) + PFC_MD_secondary_coef * W_PFC_MD .* sum(((i-last_spike_PFC_S_SE<spikewidth) & (i-last_spike_PFC_S_SE>0) & (Matrix_PFC_to_MD ~= 0)), 1).';

        Isyn_PFC_D_MD_right(:,i+ 10*MD_delay) = Isyn_PFC_D_MD_right(:,i+ 10*MD_delay) + PFC_MD_secondary_coef *W_PFC_MD .* sum(((i-last_spike_PFC_remote_right<spikewidth) & (i-last_spike_PFC_remote_right>0) & (Matrix_PFC_to_MD ~= 0)), 1).';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Concrete Rule VA thalamus to cortex 
    
        Isyn_VA_Matrix_to_PFC_con(sum(((i-last_spike_VA_matrix_left<spikewidth_VA) & (i-last_spike_VA_matrix_left>0)) .* matrix_MD_to_Crtex, 1)~=0,i+ delay) = Matrix_amplification;

        Isyn_VA_Matrix_to_PFC_con(sum(((i-last_spike_VA_matrix_mid<spikewidth_VA) & (i-last_spike_VA_matrix_mid>0)) .* matrix_MD_to_Crtex, 1)~=0,i+ delay) = Matrix_amplification;

        Isyn_VA_Matrix_to_PFC_con(sum(((i-last_spike_VA_matrix_right<spikewidth_VA) & (i-last_spike_VA_matrix_right>0)) .* matrix_MD_to_Crtex, 1)~=0,i+ delay) = Matrix_amplification;

        Isyn_VA_Matrix_to_PFC_exc_left(:,i+ delay) = Isyn_VA_Matrix_to_PFC_exc_left(:,i+ delay) +  2*W_VA_exi .* sum(((i-last_spike_VA_matrix_left<spikewidth) & (i-last_spike_VA_matrix_left>0) & (matrix_VA_to_rule ~= 0)), 1).';
        
        Isyn_VA_Matrix_to_PFC_exc_mid(:,i+ delay) = Isyn_VA_Matrix_to_PFC_exc_mid(:,i+ delay) +  2*W_VA_exi .* sum(((i-last_spike_VA_matrix_mid<spikewidth) & (i-last_spike_VA_matrix_mid>0) & (matrix_VA_to_rule ~= 0)), 1).';

        Isyn_VA_Matrix_to_PFC_exc_right(:,i+ delay) = Isyn_VA_Matrix_to_PFC_exc_right(:,i+ delay) +  2*W_VA_exi .* sum(((i-last_spike_VA_matrix_right<spikewidth) & (i-last_spike_VA_matrix_right>0) & (matrix_VA_to_rule ~= 0)), 1).';

        % Concrete Rule VA to Inh 

        Isyn_VA_Matrix_to_Inh_con(:,i+ delay) = Isyn_VA_Matrix_to_Inh_con(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_VA_matrix_left<spikewidth) & (i-last_spike_VA_matrix_left>0), 1).';
        Isyn_VA_Matrix_to_Inh_con(:,i+ delay) = Isyn_VA_Matrix_to_Inh_con(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_VA_matrix_mid<spikewidth) & (i-last_spike_VA_matrix_mid>0), 1).';
        Isyn_VA_Matrix_to_Inh_con(:,i+ delay) = Isyn_VA_Matrix_to_Inh_con(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_VA_matrix_right<spikewidth) & (i-last_spike_VA_matrix_right>0), 1).';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Concrete Rule MD to PFC 

        Isyn_MD_left_to_PFC(:,i+ 8*MD_delay) = Isyn_MD_left_to_PFC(:,i+ 8*MD_delay) +  0.3*W_MD_PFC .* sum(((i-last_spike_MD_left<spikewidth_MD) & (i-last_spike_MD_left>0) & (Matrix_MD_to_PFC ~= 0)), 1).'; 

        Isyn_MD_mid_to_PFC(:,i+ 8*MD_delay) = Isyn_MD_mid_to_PFC(:,i+ 8*MD_delay) +  0.3*W_MD_PFC .* sum(((i-last_spike_MD_mid<spikewidth_MD) & (i-last_spike_MD_mid>0) & (Matrix_MD_to_PFC ~= 0)), 1).'; 

        Isyn_MD_right_to_PFC(:,i+ 8*MD_delay) = Isyn_MD_right_to_PFC(:,i+ 8*MD_delay) +  0.3*W_MD_PFC .* sum(((i-last_spike_MD_right<spikewidth_MD) & (i-last_spike_MD_right>0) & (Matrix_MD_to_PFC ~= 0)), 1).'; 

        Isyn_MD_left_to_random_PFC_Bow(:,i+ delay) = Isyn_MD_left_to_random_PFC_Bow(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_left<spikewidth_MD) & (i-last_spike_MD_left>0) & (Matrix_random_MD_PFC_Bow ~= 0)), 1).'; 

        Isyn_MD_left_to_random_PFC_Rec(:,i+ delay) = Isyn_MD_left_to_random_PFC_Rec(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_left<spikewidth_MD) & (i-last_spike_MD_left>0) & (Matrix_random_MD_PFC_Rec ~= 0)), 1).'; 

        Isyn_MD_left_to_random_PFC_Oval(:,i+ delay) = Isyn_MD_left_to_random_PFC_Oval(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_left<spikewidth_MD) & (i-last_spike_MD_left>0) & (Matrix_random_MD_PFC_Oval ~= 0)), 1).'; 

        Isyn_MD_left_to_random_PFC_SW(:,i+ delay) = Isyn_MD_left_to_random_PFC_SW(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_left<spikewidth_MD) & (i-last_spike_MD_left>0) & (Matrix_random_MD_PFC_SW ~= 0)), 1).'; 

        Isyn_MD_left_to_random_PFC_N(:,i+ delay) = Isyn_MD_left_to_random_PFC_N(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_left<spikewidth_MD) & (i-last_spike_MD_left>0) & (Matrix_random_MD_PFC_N ~= 0)), 1).'; 

        Isyn_MD_left_to_random_PFC_SE(:,i+ delay) = Isyn_MD_left_to_random_PFC_SE(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_left<spikewidth_MD) & (i-last_spike_MD_left>0) & (Matrix_random_MD_PFC_SE ~= 0)), 1).'; 

        Isyn_MD_mid_to_random_PFC_Bow(:,i+ delay) = Isyn_MD_mid_to_random_PFC_Bow(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_mid<spikewidth_MD) & (i-last_spike_MD_mid>0) & (Matrix_random_MD_PFC_Bow ~= 0)), 1).'; 

        Isyn_MD_mid_to_random_PFC_Rec(:,i+ delay) = Isyn_MD_mid_to_random_PFC_Rec(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_mid<spikewidth_MD) & (i-last_spike_MD_mid>0) & (Matrix_random_MD_PFC_Rec ~= 0)), 1).'; 

        Isyn_MD_mid_to_random_PFC_Oval(:,i+ delay) = Isyn_MD_mid_to_random_PFC_Oval(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_mid<spikewidth_MD) & (i-last_spike_MD_mid>0) & (Matrix_random_MD_PFC_Oval ~= 0)), 1).'; 

        Isyn_MD_mid_to_random_PFC_SW(:,i+ delay) = Isyn_MD_mid_to_random_PFC_SW(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_mid<spikewidth_MD) & (i-last_spike_MD_mid>0) & (Matrix_random_MD_PFC_SW ~= 0)), 1).'; 

        Isyn_MD_mid_to_random_PFC_N(:,i+ delay) = Isyn_MD_mid_to_random_PFC_N(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_mid<spikewidth_MD) & (i-last_spike_MD_mid>0) & (Matrix_random_MD_PFC_N ~= 0)), 1).'; 

        Isyn_MD_mid_to_random_PFC_SE(:,i+ delay) = Isyn_MD_mid_to_random_PFC_SE(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_mid<spikewidth_MD) & (i-last_spike_MD_mid>0) & (Matrix_random_MD_PFC_SE ~= 0)), 1).'; 
        
        Isyn_MD_right_to_random_PFC_Bow(:,i+ delay) = Isyn_MD_right_to_random_PFC_Bow(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_right<spikewidth_MD) & (i-last_spike_MD_right>0) & (Matrix_random_MD_PFC_Bow ~= 0)), 1).'; 

        Isyn_MD_right_to_random_PFC_Rec(:,i+ delay) = Isyn_MD_right_to_random_PFC_Rec(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_right<spikewidth_MD) & (i-last_spike_MD_right>0) & (Matrix_random_MD_PFC_Rec ~= 0)), 1).'; 

        Isyn_MD_right_to_random_PFC_Oval(:,i+ delay) = Isyn_MD_right_to_random_PFC_Oval(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_right<spikewidth_MD) & (i-last_spike_MD_right>0) & (Matrix_random_MD_PFC_Oval ~= 0)), 1).'; 

        Isyn_MD_right_to_random_PFC_SW(:,i+ delay) = Isyn_MD_right_to_random_PFC_SW(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_right<spikewidth_MD) & (i-last_spike_MD_right>0) & (Matrix_random_MD_PFC_SW ~= 0)), 1).'; 

        Isyn_MD_right_to_random_PFC_N(:,i+ delay) = Isyn_MD_right_to_random_PFC_N(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_right<spikewidth_MD) & (i-last_spike_MD_right>0) & (Matrix_random_MD_PFC_N ~= 0)), 1).'; 

        Isyn_MD_right_to_random_PFC_SE(:,i+ delay) = Isyn_MD_right_to_random_PFC_SE(:,i+ delay) +  W_MD_PFC .* sum(((i-last_spike_MD_right<spikewidth_MD) & (i-last_spike_MD_right>0) & (Matrix_random_MD_PFC_SE ~= 0)), 1).'; 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Concrete Rule PFC deep/superficial  to ST 

        Isyn_aPFC_D_con_to_ST(:,i+ delay) = Isyn_aPFC_D_con_to_ST(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_PFC_S_Bow<spikewidth) & (i-last_spike_PFC_S_Bow>0), 1).';

        Isyn_aPFC_D_con_to_ST(:,i+ delay) = Isyn_aPFC_D_con_to_ST(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_PFC_S_Rec<spikewidth) & (i-last_spike_PFC_S_Rec>0), 1).';

        Isyn_aPFC_D_con_to_ST(:,i+ delay) = Isyn_aPFC_D_con_to_ST(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_PFC_S_Oval<spikewidth) & (i-last_spike_PFC_S_Oval>0), 1).';

        Isyn_aPFC_D_con_to_ST(:,i+ delay) = Isyn_aPFC_D_con_to_ST(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_PFC_S_SW<spikewidth) & (i-last_spike_PFC_S_SW>0), 1).';
       
        Isyn_aPFC_D_con_to_ST(:,i+ delay) = Isyn_aPFC_D_con_to_ST(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_PFC_S_N<spikewidth) & (i-last_spike_PFC_S_N>0), 1).';

        Isyn_aPFC_D_con_to_ST(:,i+ delay) = Isyn_aPFC_D_con_to_ST(:,i+ delay) +  W_PFC_to_str .* sum((i-last_spike_PFC_S_SE<spikewidth) & (i-last_spike_PFC_S_SE>0), 1).';

        %%%%%%%%%%%%%%%%%%%%% ST to SNpr
        Isyn_ST_con(:, i + delay) =  Isyn_ST_con(:, i + delay) + W_IPL_Inh_to_exc.*sum(((i-last_spike_ST_con < spikewidth_inh) & (i-last_spike_ST_con > 0)).*(exp(((i-last_spike_ST_con < spikewidth_inh) & (i-last_spike_ST_con > 0)).*(1-(i-last_spike_ST_con-delay)).*2./spikewidth_inh)), 1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SNpr tO vA
        %         IPL Inh neurons to excitatory
        Isyn_SNpr_to_VA_con(:, i + delay) =  Isyn_SNpr_to_VA_con(:, i + delay) + W_IPL_Inh_to_exc.*sum(((i-last_spike_SNpr_Inh_con < spikewidth_inh) & (i-last_spike_SNpr_Inh_con > 0)).*(exp(((i-last_spike_SNpr_Inh_con < spikewidth_inh) & (i-last_spike_SNpr_Inh_con > 0)).*(1-(i-last_spike_SNpr_Inh_con-delay)).*2./spikewidth_inh)), 1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Concrete Rule cortical inhibitory cells 

        if ~PV_off
            Isyn_PV_con(:, i + delay) =  Isyn_PV_con(:, i + delay) + W_IPL_Inh_to_exc_pv.*sum(((i-last_spike_PV_con < spikewidth_inh) & (i-last_spike_PV_con > 0)).*(exp(((i-last_spike_PV_con < spikewidth_inh) & (i-last_spike_PV_con > 0)).*(1-(i-last_spike_PV_con-delay)).*2./spikewidth_inh)), 1);
        end

        Isyn_MD_left_to_Inh(:,i+ delay) = Isyn_MD_left_to_Inh(:,i+ delay) + W_MD_inh .* sum(((i-last_spike_MD_left<spikewidth_VA) & (i-last_spike_MD_left>0)), 1).';
        
        Isyn_MD_mid_to_Inh(:,i+ delay) = Isyn_MD_mid_to_Inh(:,i+ delay) + W_MD_inh .* sum(((i-last_spike_MD_mid<spikewidth_VA) & (i-last_spike_MD_mid>0)), 1).';

        Isyn_MD_right_to_Inh(:,i+ delay) = Isyn_MD_right_to_Inh(:,i+ delay) + W_MD_inh .* sum(((i-last_spike_MD_right<spikewidth_VA) & (i-last_spike_MD_right>0)), 1).';

        Isyn_FS_left(:, i + delay) =  Isyn_FS_left(:, i + delay) + W_IPL_Inh_to_exc_fs * W_IPL_Inh_to_exc_fs.*sum(((i-last_spike_FS_left < spikewidth_inh_FS) & (i-last_spike_FS_left > 0)).*(exp(((i-last_spike_FS_left < spikewidth_inh_FS) & (i-last_spike_FS_left > 0)).*(1-(i-last_spike_FS_left-delay)).*2./spikewidth_inh_FS)), 1);

        Isyn_FS_mid(:, i + delay) =  Isyn_FS_mid(:, i + delay) + W_IPL_Inh_to_exc_fs * W_IPL_Inh_to_exc_fs.*sum(((i-last_spike_FS_mid < spikewidth_inh_FS) & (i-last_spike_FS_mid > 0)).*(exp(((i-last_spike_FS_mid < spikewidth_inh_FS) & (i-last_spike_FS_mid > 0)).*(1-(i-last_spike_FS_mid-delay)).*2./spikewidth_inh_FS)), 1);

        Isyn_FS_right(:, i + delay) =  Isyn_FS_right(:, i + delay) + W_IPL_Inh_to_exc_fs * W_IPL_Inh_to_exc_fs.*sum(((i-last_spike_FS_right < spikewidth_inh_FS) & (i-last_spike_FS_right > 0)).*(exp(((i-last_spike_FS_right < spikewidth_inh_FS) & (i-last_spike_FS_right > 0)).*(1-(i-last_spike_FS_right-delay)).*2./spikewidth_inh_FS)), 1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% local excitatory input to
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% shape cells 
        
        Isyn_aPFC_local_left(:,i+ delay) = Isyn_aPFC_local_left(:,i+ delay) + sum(((i-last_spike_PFC_remote_left<spikewidth) & (i-last_spike_PFC_remote_left>0)) .* matrix_local, 1).';

        Isyn_aPFC_local_mid(:,i+ delay) = Isyn_aPFC_local_mid(:,i+ delay) + sum(((i-last_spike_PFC_remote_mid<spikewidth) & (i-last_spike_PFC_remote_mid>0)) .* matrix_local, 1).';
    
        Isyn_aPFC_local_right(:,i+ delay) = Isyn_aPFC_local_right(:,i+ delay) + sum(((i-last_spike_PFC_remote_right<spikewidth) & (i-last_spike_PFC_remote_right>0)) .* matrix_local, 1).';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Concrete PFC rule cells to output cells, use the connectivity as PFC_M to PFC_superfical 
        
        Isyn_PFC_output_left(:,i+ delay) = Isyn_PFC_output_left(:,i+ delay) + sum(((i-last_spike_PFC_remote_left<spikewidth) & (i-last_spike_PFC_remote_left>0)) .* matrix_M, 1).';

        Isyn_PFC_output_mid(:,i+ delay) = Isyn_PFC_output_mid(:,i+ delay) + sum(((i-last_spike_PFC_remote_mid<spikewidth) & (i-last_spike_PFC_remote_mid>0)) .* matrix_M, 1).';
    
        Isyn_PFC_output_right(:,i+ delay) = Isyn_PFC_output_right(:,i+ delay) + sum(((i-last_spike_PFC_remote_right<spikewidth) & (i-last_spike_PFC_remote_right>0)) .* matrix_M, 1).';

        for j = 1:numberofneurons
            % create noise
            noise_amp = 10;
            noise_prob = 0.000002;
            noise_prob_PFC=0.000002;
            n1=rand;
            if n1<noise_prob
                I_noise_PFC_1= noise_amp;
            else
                I_noise_PFC_1=0;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 1 PFC middle layer first that gets the stim

            if (last_spike_PFC_M(j)~=10^10 && (i-last_spike_PFC_M(j))>tref)

                y_PFC_M(j,i) = y_PFC_M(j,i-1)+leaky_coef*((E_L-y_PFC_M(j,i-1))/tha)*dt+(Iext_PFC_M(j)+ I_stim(j,i)+ I_noise_PFC_1 )*dt*(RM/tha);

            elseif last_spike_PFC_M(j)==10^10

                y_PFC_M(j,i) = y_PFC_M(j,i-1)+leaky_coef*((E_L-y_PFC_M(j,i-1))/tha)*dt+(Iext_PFC_M(j)+ I_stim(j,i)+ I_noise_PFC_1 )*dt*(RM/tha);

            else

                y_PFC_M(j,i)=E_L;

            end

            if y_PFC_M(j,i)>=v_th
                last_spike_PFC_M(j)=i;
                y_PFC_M(j,i)=0;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Layer 2, where the layer 1 send signal to this cell
            % membrane potential across PFC cells
            % BT
            n1=rand;
            if n1<noise_prob
                I_noise_PFC_1= noise_amp;
            else
                I_noise_PFC_1=0;
            end

            if (last_spike_PFC_S_BT(j)~=10^10 && (i-last_spike_PFC_S_BT(j))>tref)

                y_PFC_S_BT(j,i) = y_PFC_S_BT(j,i-1)+leaky_coef*((E_L-y_PFC_S_BT(j,i-1))/tha)*dt+(Iext_SPFC_BT(j)+ Isyn_PFC_M_BT (j, i)+ (Isyn_aPFC_local_BT(j, i) + Isyn_aPFC_local_shared (j,i))* Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1- Isyn_PV(j, i) + Isyn_MD_shape_to_random_PFC_BT(j,i)+Isyn_MD_orientation_to_random_PFC_BT(j,i))*dt*(RM/tha);
            elseif last_spike_PFC_S_BT(j)==10^10

                y_PFC_S_BT(j,i) = y_PFC_S_BT(j,i-1)+leaky_coef*((E_L-y_PFC_S_BT(j,i-1))/tha)*dt+(Iext_SPFC_BT(j)+ Isyn_PFC_M_BT (j, i)+ (Isyn_aPFC_local_BT(j, i) + Isyn_aPFC_local_shared (j,i))* Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1- Isyn_PV(j, i) +Isyn_MD_shape_to_random_PFC_BT(j,i)+Isyn_MD_orientation_to_random_PFC_BT(j,i))*dt*(RM/tha);

            else

                y_PFC_S_BT(j,i)=E_L;

            end

            if y_PFC_S_BT(j,i)>=v_th
                last_spike_PFC_S_BT(j)=i;
                y_PFC_S_BT(j,i)=0;
            end

            % RC

            n1=rand;
            if n1<noise_prob
                I_noise_PFC_1= noise_amp;
            else
                I_noise_PFC_1=0;
            end

            if (last_spike_PFC_S_RC(j)~=10^10 && (i-last_spike_PFC_S_RC(j))>tref)

                y_PFC_S_RC(j,i) = y_PFC_S_RC(j,i-1)+ leaky_coef*((E_L-y_PFC_S_RC(j,i-1))/tha)*dt + (Iext_SPFC_RC(j)+ Isyn_PFC_M_RC (j, i)+ (Isyn_aPFC_local_RC(j, i) + Isyn_aPFC_local_shared (j,i))* Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1- Isyn_PV(j, i) +Isyn_MD_shape_to_random_PFC_RC(j,i)+Isyn_MD_orientation_to_random_PFC_RC(j,i))*dt*(RM/tha);
            elseif last_spike_PFC_S_RC(j)==10^10

                y_PFC_S_RC(j,i) = y_PFC_S_RC(j,i-1)+ leaky_coef*((E_L-y_PFC_S_RC(j,i-1))/tha)*dt + (Iext_SPFC_RC(j)+ Isyn_PFC_M_RC (j, i)+ (Isyn_aPFC_local_RC(j, i) + Isyn_aPFC_local_shared (j,i))* Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1- Isyn_PV(j, i) +Isyn_MD_shape_to_random_PFC_RC(j,i)+Isyn_MD_orientation_to_random_PFC_RC(j,i))*dt*(RM/tha);

            else

                y_PFC_S_RC(j,i) = E_L;

            end

            if y_PFC_S_RC(j,i)>= v_th
                last_spike_PFC_S_RC(j)= i;
                y_PFC_S_RC(j,i)= 0;
            end
            
            % GC

            n1=rand;
            if n1<noise_prob
                I_noise_PFC_1= noise_amp;
            else
                I_noise_PFC_1=0;
            end

            if (last_spike_PFC_S_GC(j)~=10^10 && (i-last_spike_PFC_S_GC(j))>tref)

                y_PFC_S_GC(j,i) = y_PFC_S_GC(j,i-1)+ leaky_coef*((E_L-y_PFC_S_GC(j,i-1))/tha)*dt + (Iext_SPFC_GC(j)+ Isyn_PFC_M_GC (j, i) + (Isyn_aPFC_local_GC(j, i) + Isyn_aPFC_local_shared (j,i)) * Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1- Isyn_PV(j, i) +Isyn_MD_shape_to_random_PFC_GC(j,i)+Isyn_MD_orientation_to_random_PFC_GC(j,i))*dt*(RM/tha);
            elseif last_spike_PFC_S_GC(j)==10^10

                y_PFC_S_GC(j,i) = y_PFC_S_GC(j,i-1)+ leaky_coef*((E_L-y_PFC_S_GC(j,i-1))/tha)*dt + (Iext_SPFC_GC(j)+ Isyn_PFC_M_GC (j, i) + (Isyn_aPFC_local_GC(j, i) + Isyn_aPFC_local_shared (j,i)) * Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1- Isyn_PV(j, i)+Isyn_MD_shape_to_random_PFC_GC(j,i)+Isyn_MD_orientation_to_random_PFC_GC(j,i) )*dt*(RM/tha);

            else

                y_PFC_S_GC(j,i) = E_L;

            end

            if y_PFC_S_GC(j,i)>= v_th
                last_spike_PFC_S_GC(j)= i;
                y_PFC_S_GC(j,i)= 0;
            end
            
            % YT

            n1=rand;
            if n1<noise_prob
                I_noise_PFC_1= noise_amp;
            else
                I_noise_PFC_1=0;
            end

            if (last_spike_PFC_S_YT(j)~=10^10 && (i-last_spike_PFC_S_YT(j))>tref)

                y_PFC_S_YT(j,i) = y_PFC_S_YT(j,i-1)+ leaky_coef*((E_L-y_PFC_S_YT(j,i-1))/tha)*dt + (Iext_SPFC_YT(j)+ Isyn_PFC_M_YT (j, i) + (Isyn_aPFC_local_YT(j, i) + Isyn_aPFC_local_shared (j,i)) * Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1 - Isyn_PV(j, i)+Isyn_MD_shape_to_random_PFC_YT(j,i)+Isyn_MD_orientation_to_random_PFC_YT(j,i))*dt*(RM/tha);
            elseif last_spike_PFC_S_YT(j)==10^10

                y_PFC_S_YT(j,i) = y_PFC_S_YT(j,i-1)+ leaky_coef*((E_L-y_PFC_S_YT(j,i-1))/tha)*dt + (Iext_SPFC_YT(j)+ Isyn_PFC_M_YT (j, i) + (Isyn_aPFC_local_YT(j, i) + Isyn_aPFC_local_shared (j,i)) * Isyn_VA_Matrix_to_PFC(j, i)+ I_noise_PFC_1 - Isyn_PV(j, i)+Isyn_MD_shape_to_random_PFC_YT(j,i)+Isyn_MD_orientation_to_random_PFC_YT(j,i))*dt*(RM/tha);

            else

                y_PFC_S_YT(j,i) = E_L;

            end

            if y_PFC_S_YT(j,i)>= v_th
                last_spike_PFC_S_YT(j)= i;
                y_PFC_S_YT(j,i)= 0;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 3 PFC Deep layer
            % membrane potential across deep layers of PFC cells

            % BT

            n1=rand;
            if n1<noise_prob_PFC
                I_noise_PFC_2= noise_amp;
            else
                I_noise_PFC_2= 0;
            end

            if (last_spike_PFC_D_BT(j)~=10^10 && (i-last_spike_PFC_D_BT(j))>tref)

                y_aPFC_D_BT(j,i) = y_aPFC_D_BT(j,i-1)+leaky_coef*((E_L-y_aPFC_D_BT(j,i-1))/tha)*dt+(Iext_DPFC_BT(j)+ (W11 * Isyn_aPFC_to_D_BT(j,i)+ W12 * Isyn_aPFC_to_D_RC(j,i)+ W13 *Isyn_aPFC_to_D_GC(j,i)+ W14* Isyn_aPFC_to_D_YT(j,i)) + Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2)*dt*(RM/tha);

            elseif last_spike_PFC_D_BT(j)==10^10

                y_aPFC_D_BT(j,i) = y_aPFC_D_BT(j,i-1)+leaky_coef*((E_L-y_aPFC_D_BT(j,i-1))/tha)*dt+(Iext_DPFC_BT(j)+ (W11 * Isyn_aPFC_to_D_BT(j,i)+ W12 * Isyn_aPFC_to_D_RC(j,i)+ W13 *Isyn_aPFC_to_D_GC(j,i)+ W14* Isyn_aPFC_to_D_YT(j,i)) + Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2)*dt*(RM/tha); 

            else

                y_aPFC_D_BT(j,i)= E_L;

            end

            if y_aPFC_D_BT(j,i)>= v_th
                last_spike_PFC_D_BT(j)= i;
                y_aPFC_D_BT(j,i)= 0;
            end


            % RC

            n1=rand;
            if n1<noise_prob_PFC
                I_noise_PFC_2= noise_amp;
            else
                I_noise_PFC_2= 0;
            end

            if (last_spike_PFC_D_RC(j)~=10^10 && (i-last_spike_PFC_D_RC(j))>tref)

                y_aPFC_D_RC(j,i) = y_aPFC_D_RC(j,i-1)+leaky_coef*((E_L-y_aPFC_D_RC(j,i-1))/tha)*dt+(Iext_DPFC_RC(j)+ (W21 * Isyn_aPFC_to_D_BT(j,i)+ W22 * Isyn_aPFC_to_D_RC(j,i)+ W23 *Isyn_aPFC_to_D_GC(j,i)+ W24* Isyn_aPFC_to_D_YT(j,i)) +Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2)*dt*(RM/tha);

            elseif last_spike_PFC_D_RC(j)==10^10

                y_aPFC_D_RC(j,i) = y_aPFC_D_RC(j,i-1)+leaky_coef*((E_L-y_aPFC_D_RC(j,i-1))/tha)*dt+(Iext_DPFC_RC(j)+ (W21 * Isyn_aPFC_to_D_BT(j,i)+ W22 * Isyn_aPFC_to_D_RC(j,i)+ W23 *Isyn_aPFC_to_D_GC(j,i)+ W24* Isyn_aPFC_to_D_YT(j,i)) +Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2)*dt*(RM/tha); 

            else

                y_aPFC_D_RC(j,i)= E_L;

            end

            if y_aPFC_D_RC(j,i)>= v_th
                last_spike_PFC_D_RC(j)= i;
                y_aPFC_D_RC(j,i)= 0;
            end

            % GC

            n1=rand;
            if n1<noise_prob_PFC
                I_noise_PFC_2= noise_amp;
            else
                I_noise_PFC_2= 0;
            end

            if (last_spike_PFC_D_GC(j)~=10^10 && (i-last_spike_PFC_D_GC(j))>tref)

                y_aPFC_D_GC(j,i) = y_aPFC_D_GC(j,i-1)+leaky_coef*((E_L-y_aPFC_D_GC(j,i-1))/tha)*dt+(Iext_DPFC_GC(j)+ (W31 * Isyn_aPFC_to_D_BT(j,i)+ W32 * Isyn_aPFC_to_D_RC(j,i)+ W33 *Isyn_aPFC_to_D_GC(j,i)+ W34* Isyn_aPFC_to_D_YT(j,i))+ Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2)*dt*(RM/tha); 

            elseif last_spike_PFC_D_GC(j)==10^10

                y_aPFC_D_GC(j,i) = y_aPFC_D_GC(j,i-1)+leaky_coef*((E_L-y_aPFC_D_GC(j,i-1))/tha)*dt+(Iext_DPFC_GC(j)+ (W31 * Isyn_aPFC_to_D_BT(j,i)+ W32 * Isyn_aPFC_to_D_RC(j,i)+ W33 *Isyn_aPFC_to_D_GC(j,i)+ W34* Isyn_aPFC_to_D_YT(j,i)) + Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2)*dt*(RM/tha);

            else

                y_aPFC_D_GC(j,i)= E_L;

            end

            if y_aPFC_D_GC(j,i)>= v_th
                last_spike_PFC_D_GC(j)= i;
                y_aPFC_D_GC(j,i)= 0;
            end

            % YT

            n1=rand;
            if n1<noise_prob_PFC
                I_noise_PFC_2= noise_amp;
            else
                I_noise_PFC_2= 0;
            end

            if (last_spike_PFC_D_YT(j)~=10^10 && (i-last_spike_PFC_D_YT(j))>tref)

                y_aPFC_D_YT(j,i) = y_aPFC_D_YT(j,i-1)+leaky_coef*((E_L-y_aPFC_D_YT(j,i-1))/tha)*dt+(Iext_DPFC_YT(j)+ (W41 * Isyn_aPFC_to_D_BT(j,i)+ W42 * Isyn_aPFC_to_D_RC(j,i)+ W43 *Isyn_aPFC_to_D_GC(j,i)+ W44* Isyn_aPFC_to_D_YT(j,i))+ Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2)*dt*(RM/tha);

            elseif last_spike_PFC_D_YT(j)==10^10

                y_aPFC_D_YT(j,i) = y_aPFC_D_YT(j,i-1)+leaky_coef*((E_L-y_aPFC_D_YT(j,i-1))/tha)*dt+(Iext_DPFC_YT(j)+ (W41 * Isyn_aPFC_to_D_BT(j,i)+ W42 * Isyn_aPFC_to_D_RC(j,i)+ W43 *Isyn_aPFC_to_D_GC(j,i)+ W44* Isyn_aPFC_to_D_YT(j,i))+ Isyn_aPFC_D_shared(j,i) +I_noise_PFC_2)*dt*(RM/tha); 

            else

                y_aPFC_D_YT(j,i)= E_L;

            end

            if y_aPFC_D_YT(j,i)>= v_th
                last_spike_PFC_D_YT(j)= i;
                y_aPFC_D_YT(j,i)= 0;
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VA cells 
            if ~VA_off

                % shape
                n1 = rand;
                noise_prob_md = noise_prob_PFC;
                if n1<noise_prob_md
                    I_noise_VA = noise_amp;
                else
                    I_noise_VA = 0;
                end
    
                if (last_spike_VA_matrix_shape(j)~=10^10 && (i-last_spike_VA_matrix_shape(j))>tref)
    
                    y_VA_matrix_shape(j,i) = y_VA_matrix_shape(j,i-1)+((E_L-y_VA_matrix_shape(j,i-1))/tha)*dt+(Iext_VA(j)+ I_noise_VA+ Isyn_aPFC_D5_shape(j, i) - Isyn_SNpr_to_VA(j,i))*dt*(RM/tha);
    
                elseif last_spike_VA_matrix_shape(j)==10^10
    
                    y_VA_matrix_shape(j,i) = y_VA_matrix_shape(j,i-1)+((E_L-y_VA_matrix_shape(j,i-1))/tha)*dt+(Iext_VA(j)+ I_noise_VA+ Isyn_aPFC_D5_shape(j, i) - Isyn_SNpr_to_VA(j,i))*dt*(RM/tha);
    
                else
    
                    y_VA_matrix_shape(j,i)=E_L;
    
                end
    
                if y_VA_matrix_shape(j,i)>=v_th
                    last_spike_VA_matrix_shape(j)=i;
                    y_VA_matrix_shape(j,i)=0;
                end

                % orientation
                n1 = rand;
                noise_prob_md = noise_prob_PFC;
                if n1<noise_prob_md
                    I_noise_VA = noise_amp;
                else
                    I_noise_VA = 0;
                end
    
                if (last_spike_VA_matrix_Orientation(j)~=10^10 && (i-last_spike_VA_matrix_Orientation(j))>tref)
    
                    y_VA_matrix_Orientation (j,i) = y_VA_matrix_Orientation (j,i-1)+((E_L-y_VA_matrix_Orientation (j,i-1))/tha)*dt+(Iext_VA(j)+ I_noise_VA+ Isyn_aPFC_D5_orientation(j, i) - Isyn_SNpr_to_VA(j,i))*dt*(RM/tha);
    
                elseif last_spike_VA_matrix_Orientation(j)==10^10
    
                    y_VA_matrix_Orientation (j,i) = y_VA_matrix_Orientation (j,i-1)+((E_L-y_VA_matrix_Orientation (j,i-1))/tha)*dt+(Iext_VA(j)+ I_noise_VA+ Isyn_aPFC_D5_orientation(j, i) - Isyn_SNpr_to_VA(j,i))*dt*(RM/tha);
    
                else
    
                    y_VA_matrix_Orientation (j,i)=E_L;
    
                end
    
                if y_VA_matrix_Orientation (j,i)>=v_th
                    last_spike_VA_matrix_Orientation(j)=i;
                    y_VA_matrix_Orientation (j,i)=0;
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % MD cells 
            if ~MD_off

                % shape
                n1 = rand;
                noise_prob_md = noise_prob_PFC;
                I_noise_md_2 = 0;
    
    
                if (last_spike_MD_shape(j)~=10^10 && (i-last_spike_MD_shape(j))>tref)
    
                    y_MD_core_shape(j,i) = y_MD_core_shape(j,i-1)+((E_L-y_MD_core_shape(j,i-1))/tha)*dt+(Iext_MD(j)+ I_noise_md_2+ Isyn_PFC_D_MD_shape(j,i) )*dt*(RM/tha);
    
                elseif last_spike_MD_shape(j)==10^10
    
                    y_MD_core_shape(j,i) = y_MD_core_shape(j,i-1)+((E_L-y_MD_core_shape(j,i-1))/tha)*dt+(Iext_MD(j)+ I_noise_md_2+ Isyn_PFC_D_MD_shape(j,i) )*dt*(RM/tha);
    
                else
    
                    y_MD_core_shape(j,i) = E_L;
    
                end
    
                if y_MD_core_shape(j,i) >= v_th
                    last_spike_MD_shape(j)=i;
                    y_MD_core_shape(j,i) = 0;
                end

                % orientation
                n1 = rand;
                noise_prob_md = noise_prob_PFC;
                I_noise_md_1 = 0;
    
                if (last_spike_MD_ori(j)~=10^10 && (i-last_spike_MD_ori(j))>tref)
    
                    y_MD_core_ori(j,i) = y_MD_core_ori(j,i-1)+((E_L-y_MD_core_ori(j,i-1))/tha)*dt+(Iext_MD(j)+I_noise_md_1+ Isyn_PFC_D_MD_ori(j, i) )*dt*(RM/tha);
    
                elseif last_spike_MD_ori(j)==10^10
    
                    y_MD_core_ori(j,i) = y_MD_core_ori(j,i-1)+((E_L-y_MD_core_ori(j,i-1))/tha)*dt+(Iext_MD(j)+I_noise_md_1+ Isyn_PFC_D_MD_ori(j, i) )*dt*(RM/tha);
    
                else
    
                    y_MD_core_ori(j,i)=E_L;
    
                end
    
                if y_MD_core_ori(j,i)>=v_th
                    last_spike_MD_ori(j)=i;
                    y_MD_core_ori(j,i)=0;
                end
            end

            %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % controls shape or orientation

            if ~ST_off
                if (last_spike_ST(j)~=10^10 && (i-last_spike_ST(j))>tref)
    
                    y_ST_inh(j,i) = y_ST_inh(j,i-1)+((E_L-y_ST_inh(j,i-1))/tha)*dt+(Iext_Inh_PV(j) + Isyn_aPFC_D5_to_ST(j,i))*dt*(RM/tha);
    
                elseif last_spike_ST(j)==10^10
    
                    y_ST_inh(j,i) = y_ST_inh(j,i-1)+((E_L-y_ST_inh(j,i-1))/tha)*dt+(Iext_Inh_PV(j) + Isyn_aPFC_D5_to_ST(j,i))*dt*(RM/tha);
    
                else
    
                    y_ST_inh(j,i)=E_L;
    
                end
    
                if y_ST_inh(j,i)>=v_th
                    last_spike_ST(j)=i;
                    y_ST_inh(j,i)=0;
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%

            if ~SNPR_off
                if (last_spike_SNpr_Inh(j)~=10^10 && (i-last_spike_SNpr_Inh(j))>tref)
    
                    y_SNpr_inh(j,i) = y_SNpr_inh(j,i-1)+((E_L-y_SNpr_inh(j,i-1))/tha)*dt+(Iext_Inh_pPFC(j) - Isyn_ST(j, i) )*dt*(RM/tha);
    
                elseif last_spike_SNpr_Inh(j)==10^10
    
                    y_SNpr_inh(j,i) = y_SNpr_inh(j,i-1)+((E_L-y_SNpr_inh(j,i-1))/tha)*dt+(Iext_Inh_pPFC(j) - Isyn_ST(j, i) )*dt*(RM/tha);
    
                else
    
                    y_SNpr_inh(j,i)=E_L;
    
                end
    
                if y_SNpr_inh(j,i)>=v_th
                    last_spike_SNpr_Inh(j)=i;
                    y_SNpr_inh(j,i)=0;
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % PFC rule cells         
            % shape
            n1=rand;
            if n1<noise_prob_PFC
                I_noise_PFC_6= noise_amp;
            else
                I_noise_PFC_6= 0;
            end

            if (last_spike_PFC_remote_shape(j)~=10^10 && (i-last_spike_PFC_remote_shape(j))>tref)

                y_PFC_Shape(j,i) = y_PFC_Shape(j,i-1)+leaky_coef*((E_L-y_PFC_Shape(j,i-1))/tha)*dt+(Iext_PFC_remote(j)+ Isyn_aPFC_local_shared (j,i) * Isyn_VA_Matrix_to_PFC(j, i) + I_noise_PFC_6+ Isyn_VA_Matrix_to_PFC_exc_shape(j, i)  + Isyn_MD_shape_to_PFC(j,i) + Isyn_aPFC_local_shape(j,i)* Isyn_VA_Matrix_to_PFC(j,i) - Isyn_FS_shape(j, i) )*dt*(RM/tha);

            elseif last_spike_PFC_remote_shape(j)==10^10

                y_PFC_Shape(j,i) = y_PFC_Shape(j,i-1)+leaky_coef*((E_L-y_PFC_Shape(j,i-1))/tha)*dt+(Iext_PFC_remote(j)+ Isyn_aPFC_local_shared (j,i) * Isyn_VA_Matrix_to_PFC(j, i) + I_noise_PFC_6+ Isyn_VA_Matrix_to_PFC_exc_shape(j, i)  +  Isyn_MD_shape_to_PFC(j,i) + Isyn_aPFC_local_shape(j,i)* Isyn_VA_Matrix_to_PFC(j,i) - Isyn_FS_shape(j, i) )*dt*(RM/tha);

            else

                y_PFC_Shape(j,i)=E_L;

            end

            if y_PFC_Shape(j,i)>=v_th
                last_spike_PFC_remote_shape(j)=i;
                y_PFC_Shape(j,i)=0;
            end


            % orientation
            n1=rand;
            if n1<noise_prob_PFC
                I_noise_PFC_6= noise_amp;
            else
                I_noise_PFC_6= 0;
            end

            if (last_spike_PFC_remote_Orientation(j)~=10^10 && (i-last_spike_PFC_remote_Orientation(j))>tref)

                y_PFC_Orientation(j,i) = y_PFC_Orientation(j,i-1)+leaky_coef*((E_L-y_PFC_Orientation(j,i-1))/tha)*dt+(Iext_PFC_remote(j)+ I_noise_PFC_6+ Isyn_aPFC_local_shared (j,i) * Isyn_VA_Matrix_to_PFC(j, i) + Isyn_VA_Matrix_to_PFC_exc_orientation(j, i) +Isyn_MD_ori_to_PFC(j,i) + Isyn_aPFC_local_ori(j,i)* Isyn_VA_Matrix_to_PFC(j,i) - Isyn_FS_ori(j, i))*dt*(RM/tha);

            elseif last_spike_PFC_remote_Orientation(j)==10^10

                y_PFC_Orientation(j,i) = y_PFC_Orientation(j,i-1)+leaky_coef*((E_L-y_PFC_Orientation(j,i-1))/tha)*dt+(Iext_PFC_remote(j)+ I_noise_PFC_6+ Isyn_aPFC_local_shared (j,i) * Isyn_VA_Matrix_to_PFC(j, i) + Isyn_VA_Matrix_to_PFC_exc_orientation(j, i) +Isyn_MD_ori_to_PFC(j,i) + Isyn_aPFC_local_ori(j,i)* Isyn_VA_Matrix_to_PFC(j,i) - Isyn_FS_ori(j, i))*dt*(RM/tha);

            else

                y_PFC_Orientation(j,i)=E_L;

            end

            if y_PFC_Orientation(j,i)>=v_th
                last_spike_PFC_remote_Orientation(j)=i;
                y_PFC_Orientation(j,i)=0;
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
            end

            % shape firing inhibit orientation
            if (last_spike_FS_ori(j)~=10^10 && (i-last_spike_FS_ori(j))>tref)

                y_FS_inh_ori(j,i) = y_FS_inh_ori(j,i-1)+((E_L-y_FS_inh_ori(j,i-1))/tha)*dt+(Iext_Inh_PV(j)+ Isyn_MD_shape_to_Inh(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            elseif last_spike_FS_ori(j)==10^10

                y_FS_inh_ori(j,i) = y_FS_inh_ori(j,i-1)+((E_L-y_FS_inh_ori(j,i-1))/tha)*dt+(Iext_Inh_PV(j)+ Isyn_MD_shape_to_Inh(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            else

                y_FS_inh_ori(j,i)=E_L;

            end
            
            if y_FS_inh_ori(j,i)>=v_th
                last_spike_FS_ori(j)=i;
                y_FS_inh_ori(j,i)=0;
            end
            
            % orientation firing inhibit shape
            if (last_spike_FS_shape(j)~=10^10 && (i-last_spike_FS_shape(j))>tref)

                y_FS_inh_shape(j,i) = y_FS_inh_shape(j,i-1)+((E_L-y_FS_inh_shape(j,i-1))/tha)*dt+(Iext_Inh_PV(j)+ Isyn_MD_ori_to_Inh(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            elseif last_spike_FS_shape(j)==10^10

                y_FS_inh_shape(j,i) = y_FS_inh_shape(j,i-1)+((E_L-y_FS_inh_shape(j,i-1))/tha)*dt+(Iext_Inh_PV(j)+ Isyn_MD_ori_to_Inh(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            else

                y_FS_inh_shape(j,i)=E_L;

            end

            if y_FS_inh_shape(j,i)>=v_th
                last_spike_FS_shape(j)=i;
                y_FS_inh_shape(j,i)=0;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % concrete rule
            % orientation firing inhibit shape
            if (last_spike_abs_inh_shape(j)~=10^10 && (i-last_spike_abs_inh_shape(j))>tref)

                y_con_inh_shape(j,i) = y_con_inh_shape(j,i-1)+((E_L-y_con_inh_shape(j,i-1))/tha)*dt+(Iext_Inh_PV(j)+ Isyn_abs_inh_cell_ori(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            elseif last_spike_abs_inh_shape(j)==10^10

                y_con_inh_shape(j,i) = y_con_inh_shape(j,i-1)+((E_L-y_con_inh_shape(j,i-1))/tha)*dt+(Iext_Inh_PV(j)+ Isyn_abs_inh_cell_ori(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            else

                y_con_inh_shape(j,i)=E_L;

            end

            if y_con_inh_shape(j,i)>=v_th
                last_spike_abs_inh_shape(j)=i;
                y_con_inh_shape(j,i)=0;
%                 spiketimes_FS_inh=[spiketimes_FS_inh;i,j];
            end
            
            % shape firing inhibit orientation
            if (last_spike_abs_inh_ori(j)~=10^10 && (i-last_spike_abs_inh_ori(j))>tref)

                y_con_inh_ori(j,i) = y_con_inh_ori(j,i-1)+((E_L-y_con_inh_ori(j,i-1))/tha)*dt+(Iext_Inh_PV(j)+ Isyn_abs_inh_cell_shape(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            elseif last_spike_abs_inh_ori(j)==10^10

                y_con_inh_ori(j,i) = y_con_inh_ori(j,i-1)+((E_L-y_con_inh_ori(j,i-1))/tha)*dt+(Iext_Inh_PV(j)+ Isyn_abs_inh_cell_shape(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            else

                y_con_inh_ori(j,i)=E_L;

            end

            if y_con_inh_ori(j,i)>=v_th
                last_spike_abs_inh_ori(j)=i;
                y_con_inh_ori(j,i)=0;
%                 spiketimes_FS_inh=[spiketimes_FS_inh;i,j];
            end


            % 1 PFC middle_ first  layer that gets the concrete stim

            if (last_spike_PFC_M_con(j)~=10^10 && (i-last_spike_PFC_M_con(j))>tref)

                y_PFC_M_con(j,i) = y_PFC_M_con(j,i-1)+leaky_coef*((E_L-y_PFC_M_con(j,i-1))/tha)*dt+(Iext_PFC_M_con(j)+ I_stim_con(j,i)+ I_noise_PFC_1 )*dt*(RM/tha);

            elseif last_spike_PFC_M_con(j)==10^10

                y_PFC_M_con(j,i) = y_PFC_M_con(j,i-1)+leaky_coef*((E_L-y_PFC_M_con(j,i-1))/tha)*dt+(Iext_PFC_M_con(j)+ I_stim_con(j,i)+ I_noise_PFC_1 )*dt*(RM/tha);

            else

                y_PFC_M_con(j,i)=E_L;

            end

            if y_PFC_M_con(j,i)>=v_th
                last_spike_PFC_M_con(j)=i;
                y_PFC_M_con(j,i)=0;
%                 spiketimes_M=[spiketimes_M;i,j];

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Layer 2, where the layer 1 send signal for correpond concrete rule to
            % this cell

            % membrane potential across PFC cells
            % Bow
            n1=rand;
            if n1<noise_prob
                I_noise_PFC_1= noise_amp;
            else
                I_noise_PFC_1=0;
            end


            if (last_spike_PFC_S_Bow(j)~=10^10 && (i-last_spike_PFC_S_Bow(j))>tref)

                y_PFC_S_Bow(j,i) = y_PFC_S_Bow(j,i-1)+leaky_coef*((E_L-y_PFC_S_Bow(j,i-1))/tha)*dt+(Iext_SPFC_Bow(j)+ Isyn_PFC_M_Bow (j, i)+ (Isyn_aPFC_local_Bow(j, i) + Isyn_aPFC_local_shared_con(j,i))* Isyn_VA_Matrix_to_PFC_con(j, i)+ I_noise_PFC_1- Isyn_PV_con(j, i) + Isyn_MD_left_to_random_PFC_Bow(j,i)+Isyn_MD_mid_to_random_PFC_Bow(j,i)+Isyn_MD_right_to_random_PFC_Bow(j,i) - Isyn_abs_cell_shape_inh_con(j, i))*dt*(RM/tha);%- Isyn_PV(j, i)
            elseif last_spike_PFC_S_Bow(j)==10^10

                y_PFC_S_Bow(j,i) = y_PFC_S_Bow(j,i-1)+leaky_coef*((E_L-y_PFC_S_Bow(j,i-1))/tha)*dt+(Iext_SPFC_Bow(j)+ Isyn_PFC_M_Bow (j, i)+ (Isyn_aPFC_local_Bow(j, i) + Isyn_aPFC_local_shared_con(j,i))* Isyn_VA_Matrix_to_PFC_con(j, i)+ I_noise_PFC_1- Isyn_PV_con(j, i) + Isyn_MD_left_to_random_PFC_Bow(j,i)+Isyn_MD_mid_to_random_PFC_Bow(j,i)+Isyn_MD_right_to_random_PFC_Bow(j,i) - Isyn_abs_cell_shape_inh_con(j, i))*dt*(RM/tha);%- Isyn_PV(j, i)

            else

                y_PFC_S_Bow(j,i)=E_L;

            end

            if y_PFC_S_Bow(j,i)>=v_th
                last_spike_PFC_S_Bow(j)=i;
                y_PFC_S_Bow(j,i)=0;
%                 spiketimes_S=[spiketimes_S;i,j];
            end
            
            % Rec
            n1=rand;
            if n1<noise_prob
                I_noise_PFC_1= noise_amp;
            else
                I_noise_PFC_1=0;
            end


            if (last_spike_PFC_S_Rec(j)~=10^10 && (i-last_spike_PFC_S_Rec(j))>tref)

                y_PFC_S_Rec(j,i) = y_PFC_S_Rec(j,i-1)+leaky_coef*((E_L-y_PFC_S_Rec(j,i-1))/tha)*dt+(Iext_SPFC_Rec(j)+ Isyn_PFC_M_Rec(j, i)+ (Isyn_aPFC_local_Rec(j, i) + Isyn_aPFC_local_shared_con(j,i))* Isyn_VA_Matrix_to_PFC_con(j, i)+ I_noise_PFC_1- Isyn_PV_con(j, i) + Isyn_MD_left_to_random_PFC_Rec(j,i)+Isyn_MD_mid_to_random_PFC_Rec(j,i)+Isyn_MD_right_to_random_PFC_Rec(j,i) - Isyn_abs_cell_shape_inh_con(j, i))*dt*(RM/tha);%- Isyn_PV(j, i)
                
            elseif last_spike_PFC_S_Rec(j)==10^10

                y_PFC_S_Rec(j,i) = y_PFC_S_Rec(j,i-1)+leaky_coef*((E_L-y_PFC_S_Rec(j,i-1))/tha)*dt+(Iext_SPFC_Rec(j)+ Isyn_PFC_M_Rec(j, i)+ (Isyn_aPFC_local_Rec(j, i) + Isyn_aPFC_local_shared_con(j,i))* Isyn_VA_Matrix_to_PFC_con(j, i)+ I_noise_PFC_1- Isyn_PV_con(j, i) + Isyn_MD_left_to_random_PFC_Rec(j,i)+Isyn_MD_mid_to_random_PFC_Rec(j,i)+Isyn_MD_right_to_random_PFC_Rec(j,i) - Isyn_abs_cell_shape_inh_con(j, i))*dt*(RM/tha);%- Isyn_PV(j, i)

            else

                y_PFC_S_Rec(j,i)=E_L;

            end

            if y_PFC_S_Rec(j,i)>=v_th
                last_spike_PFC_S_Rec(j)=i;
                y_PFC_S_Rec(j,i)=0;
%                 spiketimes_S=[spiketimes_S;i,j];
            end

            % Oval
            n1=rand;
            if n1<noise_prob
                I_noise_PFC_1= noise_amp;
            else
                I_noise_PFC_1=0;
            end


            if (last_spike_PFC_S_Oval(j)~=10^10 && (i-last_spike_PFC_S_Oval(j))>tref)

                y_PFC_S_Oval(j,i) = y_PFC_S_Oval(j,i-1)+leaky_coef*((E_L-y_PFC_S_Oval(j,i-1))/tha)*dt+(Iext_SPFC_Oval(j)+ Isyn_PFC_M_Oval(j, i)+ (Isyn_aPFC_local_Oval(j, i) + Isyn_aPFC_local_shared_con(j,i))* Isyn_VA_Matrix_to_PFC_con(j, i)+ I_noise_PFC_1- Isyn_PV_con(j, i) + Isyn_MD_left_to_random_PFC_Oval(j,i)+Isyn_MD_mid_to_random_PFC_Oval(j,i)+Isyn_MD_right_to_random_PFC_Oval(j,i) - Isyn_abs_cell_shape_inh_con(j, i))*dt*(RM/tha);%- Isyn_PV(j, i)
            elseif last_spike_PFC_S_Oval(j)==10^10

                y_PFC_S_Oval(j,i) = y_PFC_S_Oval(j,i-1)+leaky_coef*((E_L-y_PFC_S_Oval(j,i-1))/tha)*dt+(Iext_SPFC_Oval(j)+ Isyn_PFC_M_Oval(j, i)+ (Isyn_aPFC_local_Oval(j, i) + Isyn_aPFC_local_shared_con(j,i))* Isyn_VA_Matrix_to_PFC_con(j, i)+ I_noise_PFC_1- Isyn_PV_con(j, i) + Isyn_MD_left_to_random_PFC_Oval(j,i)+Isyn_MD_mid_to_random_PFC_Oval(j,i)+Isyn_MD_right_to_random_PFC_Oval(j,i) - Isyn_abs_cell_shape_inh_con(j, i))*dt*(RM/tha);%- Isyn_PV(j, i)

            else

                y_PFC_S_Oval(j,i)=E_L;

            end

            if y_PFC_S_Oval(j,i)>=v_th
                last_spike_PFC_S_Oval(j)=i;
                y_PFC_S_Oval(j,i)=0;
%                 spiketimes_S=[spiketimes_S;i,j];
            end
            
            % SW
            n1=rand;
            if n1<noise_prob
                I_noise_PFC_1= noise_amp;
            else
                I_noise_PFC_1=0;
            end


            if (last_spike_PFC_S_SW(j)~=10^10 && (i-last_spike_PFC_S_SW(j))>tref)

                y_PFC_S_SW(j,i) = y_PFC_S_SW(j,i-1)+leaky_coef*((E_L-y_PFC_S_SW(j,i-1))/tha)*dt+(Iext_SPFC_SW(j)+ Isyn_PFC_M_SW(j, i)+ (Isyn_aPFC_local_SW(j, i) + Isyn_aPFC_local_shared_con(j,i))* Isyn_VA_Matrix_to_PFC_con(j, i)+ I_noise_PFC_1- Isyn_PV_con(j, i) + Isyn_MD_left_to_random_PFC_SW(j,i)+Isyn_MD_mid_to_random_PFC_SW(j,i)+Isyn_MD_right_to_random_PFC_SW(j,i) - Isyn_abs_cell_ori_inh_con(j, i))*dt*(RM/tha);%- Isyn_PV(j, i)
            elseif last_spike_PFC_S_SW(j)==10^10

                y_PFC_S_SW(j,i) = y_PFC_S_SW(j,i-1)+leaky_coef*((E_L-y_PFC_S_SW(j,i-1))/tha)*dt+(Iext_SPFC_SW(j)+ Isyn_PFC_M_SW(j, i)+ (Isyn_aPFC_local_SW(j, i) + Isyn_aPFC_local_shared_con(j,i))* Isyn_VA_Matrix_to_PFC_con(j, i)+ I_noise_PFC_1- Isyn_PV_con(j, i) + Isyn_MD_left_to_random_PFC_SW(j,i)+Isyn_MD_mid_to_random_PFC_SW(j,i)+Isyn_MD_right_to_random_PFC_SW(j,i) - Isyn_abs_cell_ori_inh_con(j, i))*dt*(RM/tha);%- Isyn_PV(j, i)

            else

                y_PFC_S_SW(j,i)=E_L;

            end

            if y_PFC_S_SW(j,i)>=v_th
                last_spike_PFC_S_SW(j)=i;
                y_PFC_S_SW(j,i)=0;
%                 spiketimes_S=[spiketimes_S;i,j];
            end

            % N
            n1=rand;
            if n1<noise_prob
                I_noise_PFC_1= noise_amp;
            else
                I_noise_PFC_1=0;
            end


            if (last_spike_PFC_S_N(j)~=10^10 && (i-last_spike_PFC_S_N(j))>tref)

                y_PFC_S_N(j,i) = y_PFC_S_N(j,i-1)+leaky_coef*((E_L-y_PFC_S_N(j,i-1))/tha)*dt+(Iext_SPFC_N(j)+ Isyn_PFC_M_N(j, i)+ (Isyn_aPFC_local_N(j, i) + Isyn_aPFC_local_shared_con(j,i))* Isyn_VA_Matrix_to_PFC_con(j, i)+ I_noise_PFC_1- Isyn_PV_con(j, i) + Isyn_MD_left_to_random_PFC_N(j,i)+Isyn_MD_mid_to_random_PFC_N(j,i)+Isyn_MD_right_to_random_PFC_N(j,i) - Isyn_abs_cell_ori_inh_con(j, i))*dt*(RM/tha);%- Isyn_PV(j, i)
            elseif last_spike_PFC_S_N(j)==10^10

                y_PFC_S_N(j,i) = y_PFC_S_N(j,i-1)+leaky_coef*((E_L-y_PFC_S_N(j,i-1))/tha)*dt+(Iext_SPFC_N(j)+ Isyn_PFC_M_N(j, i)+ (Isyn_aPFC_local_N(j, i) + Isyn_aPFC_local_shared_con(j,i))* Isyn_VA_Matrix_to_PFC_con(j, i)+ I_noise_PFC_1- Isyn_PV_con(j, i) + Isyn_MD_left_to_random_PFC_N(j,i)+Isyn_MD_mid_to_random_PFC_N(j,i)+Isyn_MD_right_to_random_PFC_N(j,i) - Isyn_abs_cell_ori_inh_con(j, i))*dt*(RM/tha);%- Isyn_PV(j, i)

            else

                y_PFC_S_N(j,i)=E_L;

            end

            if y_PFC_S_N(j,i)>=v_th
                last_spike_PFC_S_N(j)=i;
                y_PFC_S_N(j,i)=0;
%                 spiketimes_S=[spiketimes_S;i,j];
            end

            % SE
            n1=rand;
            if n1<noise_prob
                I_noise_PFC_1= noise_amp;
            else
                I_noise_PFC_1=0;
            end


            if (last_spike_PFC_S_SE(j)~=10^10 && (i-last_spike_PFC_S_SE(j))>tref)

                y_PFC_S_SE(j,i) = y_PFC_S_SE(j,i-1)+leaky_coef*((E_L-y_PFC_S_SE(j,i-1))/tha)*dt+(Iext_SPFC_SE(j)+ Isyn_PFC_M_SE(j, i)+ (Isyn_aPFC_local_SE(j, i) + Isyn_aPFC_local_shared_con(j,i))* Isyn_VA_Matrix_to_PFC_con(j, i)+ I_noise_PFC_1- Isyn_PV_con(j, i) + Isyn_MD_left_to_random_PFC_SE(j,i)+Isyn_MD_mid_to_random_PFC_SE(j,i)+Isyn_MD_right_to_random_PFC_SE(j,i) - Isyn_abs_cell_ori_inh_con(j, i))*dt*(RM/tha);%- Isyn_PV(j, i)

            elseif last_spike_PFC_S_SE(j)==10^10

                y_PFC_S_SE(j,i) = y_PFC_S_SE(j,i-1)+leaky_coef*((E_L-y_PFC_S_SE(j,i-1))/tha)*dt+(Iext_SPFC_SE(j)+ Isyn_PFC_M_SE(j, i)+ (Isyn_aPFC_local_SE(j, i) + Isyn_aPFC_local_shared_con(j,i))* Isyn_VA_Matrix_to_PFC_con(j, i)+ I_noise_PFC_1- Isyn_PV_con(j, i) + Isyn_MD_left_to_random_PFC_SE(j,i)+Isyn_MD_mid_to_random_PFC_SE(j,i)+Isyn_MD_right_to_random_PFC_SE(j,i) - Isyn_abs_cell_ori_inh_con(j, i))*dt*(RM/tha);%- Isyn_PV(j, i)

            else

                y_PFC_S_SE(j,i)=E_L;

            end

            if y_PFC_S_SE(j,i)>=v_th
                last_spike_PFC_S_SE(j)=i;
                y_PFC_S_SE(j,i)=0;
%                 spiketimes_S=[spiketimes_S;i,j];
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Concrete rule VA matrix cells (VA is off for concrete rule)
            if 0

                % left
                n1 = rand;
                noise_prob_md = noise_prob_PFC;
                if n1<noise_prob_md
                    I_noise_VA = noise_amp;
                else
                    I_noise_VA = 0;
                end
    
                %matrix VA cells
                if (last_spike_VA_matrix_left(j)~=10^10 && (i-last_spike_VA_matrix_left(j))>tref)
    
                    y_VA_matrix_left(j,i) = y_VA_matrix_left(j,i-1)+((E_L-y_VA_matrix_left(j,i-1))/tha)*dt+(Iext_VA(j)+ I_noise_VA+ Isyn_aPFC_D_left(j, i) - Isyn_SNpr_to_VA_con(j,i))*dt*(RM/tha);%
    
                elseif last_spike_VA_matrix_left(j)==10^10
    
                    y_VA_matrix_left(j,i) = y_VA_matrix_left(j,i-1)+((E_L-y_VA_matrix_left(j,i-1))/tha)*dt+(Iext_VA(j)+ I_noise_VA+ Isyn_aPFC_D_left(j, i) - Isyn_SNpr_to_VA_con(j,i))*dt*(RM/tha);%
    
                else
    
                    y_VA_matrix_left(j,i)=E_L;
    
                end
    
                if y_VA_matrix_left(j,i)>=v_th
                    last_spike_VA_matrix_left(j)=i;
                    y_VA_matrix_left(j,i)=0;
    %                 spiketimes_VA=[spiketimes_VA;i,j];
                end
                
                % middle
                n1 = rand;
                noise_prob_md = noise_prob_PFC;
                if n1<noise_prob_md
                    I_noise_VA = noise_amp;
                else
                    I_noise_VA = 0;
                end
    
                %matrix VA cells
                if (last_spike_VA_matrix_mid(j)~=10^10 && (i-last_spike_VA_matrix_mid(j))>tref)
    
                    y_VA_matrix_mid(j,i) = y_VA_matrix_mid(j,i-1)+((E_L-y_VA_matrix_mid(j,i-1))/tha)*dt+(Iext_VA(j)+ I_noise_VA+ Isyn_aPFC_D_mid(j, i) - Isyn_SNpr_to_VA_con(j,i))*dt*(RM/tha);%
    
                elseif last_spike_VA_matrix_mid(j)==10^10
                        
                    y_VA_matrix_mid(j,i) = y_VA_matrix_mid(j,i-1)+((E_L-y_VA_matrix_mid(j,i-1))/tha)*dt+(Iext_VA(j)+ I_noise_VA+ Isyn_aPFC_D_mid(j, i) - Isyn_SNpr_to_VA_con(j,i))*dt*(RM/tha);%
    
                else
    
                    y_VA_matrix_mid(j,i)=E_L;
    
                end
    
                if y_VA_matrix_mid(j,i)>=v_th
                    last_spike_VA_matrix_mid(j)=i;
                    y_VA_matrix_mid(j,i)=0;
    %                 spiketimes_VA=[spiketimes_VA;i,j];
                end

                % right
                n1 = rand;
                noise_prob_md = noise_prob_PFC;
                if n1<noise_prob_md
                    I_noise_VA = noise_amp;
                else
                    I_noise_VA = 0;
                end
    
                %matrix VA cells
                if (last_spike_VA_matrix_right(j)~=10^10 && (i-last_spike_VA_matrix_right(j))>tref)
    
                    y_VA_matrix_right(j,i) = y_VA_matrix_right(j,i-1)+((E_L-y_VA_matrix_right(j,i-1))/tha)*dt+(Iext_VA(j)+ I_noise_VA+ Isyn_aPFC_D_right(j, i) - Isyn_SNpr_to_VA_con(j,i))*dt*(RM/tha);%
    
                elseif last_spike_VA_matrix_right(j)==10^10
                        
                    y_VA_matrix_right(j,i) = y_VA_matrix_right(j,i-1)+((E_L-y_VA_matrix_right(j,i-1))/tha)*dt+(Iext_VA(j)+ I_noise_VA+ Isyn_aPFC_D_right(j, i) - Isyn_SNpr_to_VA_con(j,i))*dt*(RM/tha);%
    
                else
    
                    y_VA_matrix_right(j,i)=E_L;
    
                end
    
                if y_VA_matrix_right(j,i)>=v_th
                    last_spike_VA_matrix_right(j)=i;
                    y_VA_matrix_right(j,i)=0;
    %                 spiketimes_VA=[spiketimes_VA;i,j];
                end
            end

              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Core MD cells amplifying pPFC connections
            if ~MD_off

                % left
                n1 = rand;
                noise_prob_md = noise_prob_PFC;
                if n1<noise_prob_md
                    I_noise_md_2 = 0;%noise_amp;
                else
                    I_noise_md_2 = 0;
                end
    
    
                if (last_spike_MD_left(j)~=10^10 && (i-last_spike_MD_left(j))>tref)
    
                    y_MD_core_left(j,i) = y_MD_core_left(j,i-1)+((E_L-y_MD_core_left(j,i-1))/tha)*dt+(Iext_MD(j)+ I_noise_md_2+ Isyn_PFC_D_MD_left(j,i) )*dt*(RM/tha);%
    
                elseif last_spike_MD_left(j)==10^10
    
                    y_MD_core_left(j,i) = y_MD_core_left(j,i-1)+((E_L-y_MD_core_left(j,i-1))/tha)*dt+(Iext_MD(j)+ I_noise_md_2+ Isyn_PFC_D_MD_left(j,i) )*dt*(RM/tha);%
    
                else
    
                    y_MD_core_left(j,i) = E_L;
    
                end
    
                if y_MD_core_left(j,i) >= v_th
                    last_spike_MD_left(j)=i;
                    y_MD_core_left(j,i) = 0;
    %                 spiketimes_MD2=[spiketimes_MD2;i,j];
    
                end

                % middle
                n1 = rand;
                noise_prob_md = noise_prob_PFC;
                if n1<noise_prob_md
                    I_noise_md_2 = 0;%noise_amp;
                else
                    I_noise_md_2 = 0;
                end
    
    
                if (last_spike_MD_mid(j)~=10^10 && (i-last_spike_MD_mid(j))>tref)
    
                    y_MD_core_mid(j,i) = y_MD_core_mid(j,i-1)+((E_L-y_MD_core_mid(j,i-1))/tha)*dt+(Iext_MD(j)+ I_noise_md_2+ Isyn_PFC_D_MD_mid(j,i) )*dt*(RM/tha);%
    
                elseif last_spike_MD_mid(j)==10^10
    
                    y_MD_core_mid(j,i) = y_MD_core_mid(j,i-1)+((E_L-y_MD_core_mid(j,i-1))/tha)*dt+(Iext_MD(j)+ I_noise_md_2+ Isyn_PFC_D_MD_mid(j,i) )*dt*(RM/tha);%
    
                else
    
                    y_MD_core_mid(j,i) = E_L;
    
                end
    
                if y_MD_core_mid(j,i) >= v_th
                    last_spike_MD_mid(j)=i;
                    y_MD_core_mid(j,i) = 0;
    %                 spiketimes_MD2=[spiketimes_MD2;i,j];
    
                end

                % right
                n1 = rand;
                noise_prob_md = noise_prob_PFC;
                if n1<noise_prob_md
                    I_noise_md_2 = 0;%noise_amp;
                else
                    I_noise_md_2 = 0;
                end
    
    
                if (last_spike_MD_right(j)~=10^10 && (i-last_spike_MD_right(j))>tref)
    
                    y_MD_core_right(j,i) = y_MD_core_right(j,i-1)+((E_L-y_MD_core_right(j,i-1))/tha)*dt+(Iext_MD(j)+ I_noise_md_2+ Isyn_PFC_D_MD_right(j,i) )*dt*(RM/tha);%
    
                elseif last_spike_MD_right(j)==10^10
    
                    y_MD_core_right(j,i) = y_MD_core_right(j,i-1)+((E_L-y_MD_core_right(j,i-1))/tha)*dt+(Iext_MD(j)+ I_noise_md_2+ Isyn_PFC_D_MD_right(j,i) )*dt*(RM/tha);%
    
                else
    
                    y_MD_core_right(j,i) = E_L;
    
                end
    
                if y_MD_core_right(j,i) >= v_th
                    last_spike_MD_right(j)=i;
                    y_MD_core_right(j,i) = 0;
    %                 spiketimes_MD2=[spiketimes_MD2;i,j];
    
                end                
            end

            %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % controls shape or orientation
            if (last_spike_ST_con(j)~=10^10 && (i-last_spike_ST_con(j))>tref)

                y_ST_inh_con(j,i) = y_ST_inh_con(j,i-1)+((E_L-y_ST_inh_con(j,i-1))/tha)*dt+(Iext_Inh_PV_con(j) + Isyn_aPFC_D_con_to_ST(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            elseif last_spike_ST_con(j)==10^10

                y_ST_inh_con(j,i) = y_ST_inh_con(j,i-1)+((E_L-y_ST_inh_con(j,i-1))/tha)*dt+(Iext_Inh_PV_con(j) + Isyn_aPFC_D_con_to_ST(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            else

                y_ST_inh_con(j,i)=E_L;

            end

            if y_ST_inh_con(j,i)>=v_th
                last_spike_ST_con(j)=i;
                y_ST_inh_con(j,i)=0;
%                 spiketimes_PV_inh=[spiketimes_PV_inh;i,j];


            end

            %%%%%%%%%%%%%%%%%%%%%%%
            %  %som inhibitory controls left or right
            if (last_spike_SNpr_Inh_con(j)~=10^10 && (i-last_spike_SNpr_Inh_con(j))>tref)

                y_SNpr_inh_con(j,i) = y_SNpr_inh_con(j,i-1)+((E_L-y_SNpr_inh_con(j,i-1))/tha)*dt+(Iext_Inh_pPFC_con(j) - Isyn_ST_con(j, i) )*dt*(RM/tha);%+ Isyn_aPFC_D5(j,i)+ Isyn_pPFC_local(j,i)

            elseif last_spike_SNpr_Inh_con(j)==10^10

                y_SNpr_inh_con(j,i) = y_SNpr_inh_con(j,i-1)+((E_L-y_SNpr_inh_con(j,i-1))/tha)*dt+(Iext_Inh_pPFC_con(j) - Isyn_ST_con(j, i) )*dt*(RM/tha);%+ Isyn_aPFC_D5(j,i)+ Isyn_pPFC_local(j,i)

            else

                y_SNpr_inh_con(j,i)=E_L;

            end

            if y_SNpr_inh_con(j,i)>=v_th
                last_spike_SNpr_Inh_con(j)=i;
                y_SNpr_inh_con(j,i)=0;
%                 spiketimes_ipl_inh=[spiketimes_ipl_inh;i,j];


            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Shape cells in pPFC
            % Isyn_aPFC_D_to_Remote_effective_3(j,i) = Isyn_aPFC_D_to_Remote(j,i)*Isyn_MD_Matrix_to_PFC(j, i);
            % Isyn_pPFC_local_effective_3(j,i) = Isyn_pPFC_local(j,i)*Isyn_MD2_to_pPFC(j, i);
            % left
            n1=rand;
            if n1<noise_prob_PFC
                I_noise_PFC_6= noise_amp;
            else
                I_noise_PFC_6= 0;
            end
            %             %
            if (last_spike_PFC_remote_left(j)~=10^10 && (i-last_spike_PFC_remote_left(j))>tref)

                y_pPFC_left(j,i) = y_pPFC_left(j,i-1)+leaky_coef*((E_L-y_pPFC_left(j,i-1))/tha)*dt+(Iext_PFC_remote_left(j)+ Isyn_aPFC_local_shared_con (j,i) * Isyn_VA_Matrix_to_PFC_con(j, i) + I_noise_PFC_6+ Isyn_VA_Matrix_to_PFC_exc_left(j, i)  + Isyn_MD_left_to_PFC(j,i) + Isyn_aPFC_local_left(j,i)* Isyn_VA_Matrix_to_PFC_con(j,i)  - Isyn_FS_left(j, i) )*dt*(RM/tha);%

            elseif last_spike_PFC_remote_left(j)==10^10

                y_pPFC_left(j,i) = y_pPFC_left(j,i-1)+leaky_coef*((E_L-y_pPFC_left(j,i-1))/tha)*dt+(Iext_PFC_remote_left(j)+ Isyn_aPFC_local_shared_con (j,i) * Isyn_VA_Matrix_to_PFC_con(j, i) + I_noise_PFC_6+ Isyn_VA_Matrix_to_PFC_exc_left(j, i)  + Isyn_MD_left_to_PFC(j,i) + Isyn_aPFC_local_left(j,i)* Isyn_VA_Matrix_to_PFC_con(j,i)  - Isyn_FS_left(j, i) )*dt*(RM/tha);%

            else

                y_pPFC_left(j,i)=E_L;

            end

            if y_pPFC_left(j,i)>=v_th
                last_spike_PFC_remote_left(j)=i;
                y_pPFC_left(j,i)=0;
%                 spiketimes_S_remote_Shape=[spiketimes_S_remote_Shape;i,j];
            end


            %%%%%%%%%%%%%%%%%%%%%%%%
            % mid
            n1=rand;
            if n1<noise_prob_PFC
                I_noise_PFC_6= noise_amp;
            else
                I_noise_PFC_6= 0;
            end
        
            if (last_spike_PFC_remote_mid(j)~=10^10 && (i-last_spike_PFC_remote_mid(j))>tref)

                y_pPFC_mid(j,i) = y_pPFC_mid(j,i-1)+leaky_coef*((E_L-y_pPFC_mid(j,i-1))/tha)*dt+(Iext_PFC_remote_mid(j)+ Isyn_aPFC_local_shared_con (j,i) * Isyn_VA_Matrix_to_PFC_con(j, i) + I_noise_PFC_6+ Isyn_VA_Matrix_to_PFC_exc_mid(j, i)  + Isyn_MD_mid_to_PFC(j,i) + Isyn_aPFC_local_mid(j,i)* Isyn_VA_Matrix_to_PFC_con(j,i)  - Isyn_FS_mid(j, i) )*dt*(RM/tha);%

            elseif last_spike_PFC_remote_mid(j)==10^10

                y_pPFC_mid(j,i) = y_pPFC_mid(j,i-1)+leaky_coef*((E_L-y_pPFC_mid(j,i-1))/tha)*dt+(Iext_PFC_remote_mid(j)+ Isyn_aPFC_local_shared_con (j,i) * Isyn_VA_Matrix_to_PFC_con(j, i) + I_noise_PFC_6+ Isyn_VA_Matrix_to_PFC_exc_mid(j, i)  + Isyn_MD_mid_to_PFC(j,i) + Isyn_aPFC_local_mid(j,i)* Isyn_VA_Matrix_to_PFC_con(j,i)  - Isyn_FS_mid(j, i) )*dt*(RM/tha);%

            else

                y_pPFC_mid(j,i)=E_L;

            end

            if y_pPFC_mid(j,i)>=v_th
                last_spike_PFC_remote_mid(j)=i;
                y_pPFC_mid(j,i)=0;
%                 spiketimes_S_remote_Shape=[spiketimes_S_remote_Shape;i,j];
            end

            %%%%%%%%%%%%%%%%%%%%%%%%
            % right
            n1=rand;
            if n1<noise_prob_PFC
                I_noise_PFC_6= noise_amp;
            else
                I_noise_PFC_6= 0;
            end
        
            if (last_spike_PFC_remote_right(j)~=10^10 && (i-last_spike_PFC_remote_right(j))>tref)

                y_pPFC_right(j,i) = y_pPFC_right(j,i-1)+leaky_coef*((E_L-y_pPFC_right(j,i-1))/tha)*dt+(Iext_PFC_remote_right(j)+ Isyn_aPFC_local_shared_con (j,i) * Isyn_VA_Matrix_to_PFC_con(j, i) + I_noise_PFC_6+ Isyn_VA_Matrix_to_PFC_exc_right(j, i)  + Isyn_MD_right_to_PFC(j,i) + Isyn_aPFC_local_right(j,i)* Isyn_VA_Matrix_to_PFC_con(j,i)  - Isyn_FS_right(j, i) )*dt*(RM/tha);%

            elseif last_spike_PFC_remote_right(j)==10^10

                y_pPFC_right(j,i) = y_pPFC_right(j,i-1)+leaky_coef*((E_L-y_pPFC_right(j,i-1))/tha)*dt+(Iext_PFC_remote_right(j)+ Isyn_aPFC_local_shared_con (j,i) * Isyn_VA_Matrix_to_PFC_con(j, i) + I_noise_PFC_6+ Isyn_VA_Matrix_to_PFC_exc_right(j, i)  + Isyn_MD_right_to_PFC(j,i) + Isyn_aPFC_local_right(j,i)* Isyn_VA_Matrix_to_PFC_con(j,i)  - Isyn_FS_right(j, i) )*dt*(RM/tha);%

            else

                y_pPFC_right(j,i)=E_L;

            end

            if y_pPFC_right(j,i)>=v_th
                last_spike_PFC_remote_right(j)=i;
                y_pPFC_right(j,i)=0;
%                 spiketimes_S_remote_Shape=[spiketimes_S_remote_Shape;i,j];
            end

             % Concrete rule 
            if (last_spike_PV_con(j)~=10^10 && (i-last_spike_PV_con(j))>tref)

                y_PV_inh_con(j,i) = y_PV_inh_con(j,i-1)+((E_L-y_PV_inh_con(j,i-1))/tha)*dt+(Iext_Inh_PV_con(j)+ Isyn_VA_Matrix_to_Inh_con(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            elseif last_spike_PV_con(j)==10^10

                y_PV_inh_con(j,i) = y_PV_inh_con(j,i-1)+((E_L-y_PV_inh_con(j,i-1))/tha)*dt+(Iext_Inh_PV_con(j)+ Isyn_VA_Matrix_to_Inh_con(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            else

                y_PV_inh_con(j,i)=E_L;

            end

            if y_PV_inh_con(j,i)>=v_th
                last_spike_PV_con(j)=i;
                y_PV_inh_con(j,i)=0;
%                 spiketimes_PV_inh=[spiketimes_PV_inh;i,j];
            end



            % mid and right firing inhibit left
            if (last_spike_FS_left(j)~=10^10 && (i-last_spike_FS_left(j))>tref)

                y_FS_inh_left(j,i) = y_FS_inh_left(j,i-1)+((E_L-y_FS_inh_left(j,i-1))/tha)*dt+(Iext_Inh_PV_con(j)+ Isyn_MD_mid_to_Inh(j,i) + Isyn_MD_right_to_Inh(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            elseif last_spike_FS_left(j)==10^10

                y_FS_inh_left(j,i) = y_FS_inh_left(j,i-1)+((E_L-y_FS_inh_left(j,i-1))/tha)*dt+(Iext_Inh_PV_con(j)+ Isyn_MD_mid_to_Inh(j,i) + Isyn_MD_right_to_Inh(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            else

                y_FS_inh_left(j,i)=E_L;

            end
            
            if y_FS_inh_left(j,i)>=v_th
                last_spike_FS_left(j)=i;
                y_FS_inh_left(j,i)=0;
%                 spiketimes_FS_inh=[spiketimes_FS_inh;i,j];
            end
            
            % left and right firing inhibit mid
            if (last_spike_FS_mid(j)~=10^10 && (i-last_spike_FS_mid(j))>tref)

                y_FS_inh_mid(j,i) = y_FS_inh_mid(j,i-1)+((E_L-y_FS_inh_mid(j,i-1))/tha)*dt+(Iext_Inh_PV_con(j)+ Isyn_MD_left_to_Inh(j,i) + Isyn_MD_right_to_Inh(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            elseif last_spike_FS_mid(j)==10^10

                y_FS_inh_mid(j,i) = y_FS_inh_mid(j,i-1)+((E_L-y_FS_inh_mid(j,i-1))/tha)*dt+(Iext_Inh_PV_con(j)+ Isyn_MD_left_to_Inh(j,i) + Isyn_MD_right_to_Inh(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            else

                y_FS_inh_mid(j,i)=E_L;

            end
            
            if y_FS_inh_mid(j,i)>=v_th
                last_spike_FS_mid(j)=i;
                y_FS_inh_mid(j,i)=0;
%                 spiketimes_FS_inh=[spiketimes_FS_inh;i,j];
            end

            % mid and right firing inhibit left
            if (last_spike_FS_left(j)~=10^10 && (i-last_spike_FS_left(j))>tref)

                y_FS_inh_left(j,i) = y_FS_inh_left(j,i-1)+((E_L-y_FS_inh_left(j,i-1))/tha)*dt+(Iext_Inh_PV_con(j)+ Isyn_MD_mid_to_Inh(j,i) + Isyn_MD_right_to_Inh(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            elseif last_spike_FS_left(j)==10^10

                y_FS_inh_left(j,i) = y_FS_inh_left(j,i-1)+((E_L-y_FS_inh_left(j,i-1))/tha)*dt+(Iext_Inh_PV_con(j)+ Isyn_MD_mid_to_Inh(j,i) + Isyn_MD_right_to_Inh(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            else

                y_FS_inh_left(j,i)=E_L;

            end
            
            if y_FS_inh_left(j,i)>=v_th
                last_spike_FS_left(j)=i;
                y_FS_inh_left(j,i)=0;
%                 spiketimes_FS_inh=[spiketimes_FS_inh;i,j];
            end
            
            % left and mid firing inhibit right
            if (last_spike_FS_right(j)~=10^10 && (i-last_spike_FS_right(j))>tref)

                y_FS_inh_right(j,i) = y_FS_inh_right(j,i-1)+((E_L-y_FS_inh_right(j,i-1))/tha)*dt+(Iext_Inh_PV_con(j)+ Isyn_MD_left_to_Inh(j,i) + Isyn_MD_mid_to_Inh(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            elseif last_spike_FS_right(j)==10^10

                y_FS_inh_right(j,i) = y_FS_inh_right(j,i-1)+((E_L-y_FS_inh_right(j,i-1))/tha)*dt+(Iext_Inh_PV_con(j)+ Isyn_MD_left_to_Inh(j,i) + Isyn_MD_mid_to_Inh(j,i))*dt*(RM/tha);%+Isyn_IPL_exc_inh(j,i)

            else

                y_FS_inh_right(j,i)=E_L;

            end
            
            if y_FS_inh_right(j,i)>=v_th
                last_spike_FS_right(j)=i;
                y_FS_inh_right(j,i)=0;
%                 spiketimes_FS_inh=[spiketimes_FS_inh;i,j];
            end

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Output layer, where the PFC concrete rule ensemble send signal to this cell
            % left
            n1=rand;
            if n1<noise_prob
                I_noise_PFC_1= noise_amp;
            else
                I_noise_PFC_1=0;
            end

            if (last_spike_PFC_output_left(j)~=10^10 && (i-last_spike_PFC_output_left(j))>tref)

                y_PFC_output_left(j,i) = y_PFC_output_left(j,i-1)+leaky_coef*((E_L-y_PFC_output_left(j,i-1))/tha)*dt+(Iext_PFC_output_left(j)+ Isyn_PFC_output_left (j, i))*dt*(RM/tha);
            elseif last_spike_PFC_output_left(j)==10^10

                y_PFC_output_left(j,i) = y_PFC_output_left(j,i-1)+leaky_coef*((E_L-y_PFC_output_left(j,i-1))/tha)*dt+(Iext_PFC_output_left(j)+ Isyn_PFC_output_left (j, i))*dt*(RM/tha);

            else

                y_PFC_output_left(j,i)=E_L;

            end

            if y_PFC_output_left(j,i)>=v_th
                last_spike_PFC_output_left(j)=i;
                y_PFC_output_left(j,i)=0;
            end
            
            % mid
            n1=rand;
            if n1<noise_prob
                I_noise_PFC_1= noise_amp;
            else
                I_noise_PFC_1=0;
            end

            if (last_spike_PFC_output_mid(j)~=10^10 && (i-last_spike_PFC_output_mid(j))>tref)

                y_PFC_output_mid(j,i) = y_PFC_output_mid(j,i-1)+leaky_coef*((E_L-y_PFC_output_mid(j,i-1))/tha)*dt+(Iext_PFC_output_mid(j)+ Isyn_PFC_output_mid(j, i))*dt*(RM/tha);
            elseif last_spike_PFC_output_mid(j)==10^10

                y_PFC_output_mid(j,i) = y_PFC_output_mid(j,i-1)+leaky_coef*((E_L-y_PFC_output_mid(j,i-1))/tha)*dt+(Iext_PFC_output_mid(j)+ Isyn_PFC_output_mid(j, i))*dt*(RM/tha);

            else

                y_PFC_output_mid(j,i)=E_L;

            end

            if y_PFC_output_mid(j,i)>=v_th
                last_spike_PFC_output_mid(j)=i;
                y_PFC_output_mid(j,i)=0;
            end

            % right
            n1=rand;
            if n1<noise_prob
                I_noise_PFC_1= noise_amp;
            else
                I_noise_PFC_1=0;
            end

            if (last_spike_PFC_output_right(j)~=10^10 && (i-last_spike_PFC_output_right(j))>tref)

                y_PFC_output_right(j,i) = y_PFC_output_right(j,i-1)+leaky_coef*((E_L-y_PFC_output_right(j,i-1))/tha)*dt+(Iext_PFC_output_right(j)+ Isyn_PFC_output_right(j, i))*dt*(RM/tha);
            elseif last_spike_PFC_output_right(j)==10^10

                y_PFC_output_right(j,i) = y_PFC_output_right(j,i-1)+leaky_coef*((E_L-y_PFC_output_right(j,i-1))/tha)*dt+(Iext_PFC_output_right(j)+ Isyn_PFC_output_right(j, i))*dt*(RM/tha);

            else

                y_PFC_output_right(j,i)=E_L;

            end

            if y_PFC_output_right(j,i)>=v_th
                last_spike_PFC_output_right(j)=i;
                y_PFC_output_right(j,i)=0;
            end
        end
    end

    if VA_off
        y_VA_matrix_shape = -10^10 * ones(numberofneurons,length(T));
        y_VA_matrix_Orientation = -10^10 * ones(numberofneurons,length(T));
    end

    if MD_off
        y_MD_core_shape = -10^10 * ones(numberofneurons,length(T));
        y_MD_core_ori = -10^10 * ones(numberofneurons,length(T));
    end

    if ST_off
        y_ST_inh = -10^10 * ones(numberofneurons,length(T));
    end

    if SNPR_off
        Isyn_SNpr_to_VA = -10^10 * ones(numberofneurons,length(T));
    end

    full_PFC_S_BT = save_to_cell(y_PFC_S_BT(:, 1:end_time), full_PFC_S_BT, v_th);

    full_PFC_S_RC = save_to_cell(y_PFC_S_RC(:, 1:end_time), full_PFC_S_RC, v_th);

    full_PFC_S_GC = save_to_cell(y_PFC_S_GC(:, 1:end_time), full_PFC_S_GC, v_th);

    full_PFC_S_YT = save_to_cell(y_PFC_S_YT(:, 1:end_time), full_PFC_S_YT, v_th);
    
    full_PFC_D_BT = save_to_cell(y_aPFC_D_BT(:, 1:end_time), full_PFC_D_BT, v_th);

    full_PFC_D_RC = save_to_cell(y_aPFC_D_RC(:, 1:end_time), full_PFC_D_RC, v_th);

    full_PFC_D_GC = save_to_cell(y_aPFC_D_GC(:, 1:end_time), full_PFC_D_GC, v_th);

    full_PFC_D_YT = save_to_cell(y_aPFC_D_YT(:, 1:end_time), full_PFC_D_YT, v_th);
    
    full_VA_shape = save_to_cell(y_VA_matrix_shape(:, 1:end_time), full_VA_shape, v_th);

    full_VA_ori = save_to_cell(y_VA_matrix_Orientation(:, 1:end_time), full_VA_ori, v_th);
    
    full_MD_shape = save_to_cell(y_MD_core_shape(:, 1:end_time), full_MD_shape, v_th);
     
    full_MD_ori = save_to_cell(y_MD_core_ori(:, 1:end_time), full_MD_ori, v_th);
    
    full_PFC_shape_ensemble = save_to_cell(y_PFC_Shape(:, 1:end_time), full_PFC_shape_ensemble, v_th);
    
    full_PFC_ori_ensemble = save_to_cell(y_PFC_Orientation(:, 1:end_time), full_PFC_ori_ensemble, v_th);    
    
    full_con_inh_shape = save_to_cell(y_con_inh_shape(:, 1:end_time), full_con_inh_shape, v_th);

    full_con_inh_ori = save_to_cell(y_con_inh_ori(:, 1:end_time), full_con_inh_ori, v_th);

    full_PFC_bow = save_to_cell(y_PFC_S_Bow(:, 1:end_time), full_PFC_bow, v_th);

    full_PFC_rec = save_to_cell(y_PFC_S_Rec(:, 1:end_time), full_PFC_rec, v_th);

    full_PFC_oval = save_to_cell(y_PFC_S_Oval(:, 1:end_time), full_PFC_oval, v_th);

    full_PFC_SW = save_to_cell(y_PFC_S_SW(:, 1:end_time), full_PFC_SW, v_th);

    full_PFC_N = save_to_cell(y_PFC_S_N(:, 1:end_time), full_PFC_N, v_th);

    full_PFC_SE = save_to_cell(y_PFC_S_SE(:, 1:end_time), full_PFC_SE, v_th);
    
    full_MD_left = save_to_cell(y_MD_core_left(:, 1:end_time), full_MD_left, v_th);
    full_MD_mid = save_to_cell(y_MD_core_mid(:, 1:end_time), full_MD_mid, v_th);
    full_MD_right = save_to_cell(y_MD_core_right(:, 1:end_time), full_MD_right, v_th);
    
    full_PFC_left = save_to_cell(y_pPFC_left(:, 1:end_time), full_PFC_left, v_th);
    full_PFC_mid = save_to_cell(y_pPFC_mid(:, 1:end_time), full_PFC_mid, v_th);
    full_PFC_right = save_to_cell(y_pPFC_right(:, 1:end_time), full_PFC_right, v_th);
    
    full_output_left = save_to_cell(y_PFC_output_left(:, 1:end_time), full_output_left, v_th);
    full_output_mid = save_to_cell(y_PFC_output_mid(:, 1:end_time), full_output_mid, v_th);
    full_output_right = save_to_cell(y_PFC_output_right(:, 1:end_time), full_output_right, v_th);

    full_Isyn_FS_left= save_to_cell(Isyn_FS_left(:, 1:end_time), full_Isyn_FS_left, v_th);
    full_Isyn_FS_mid= save_to_cell(Isyn_FS_mid(:, 1:end_time), full_Isyn_FS_mid, v_th);
    full_Isyn_FS_right= save_to_cell(Isyn_FS_right(:, 1:end_time), full_Isyn_FS_right, v_th);
    full_Isyn_MD_left_to_PFC = save_to_cell(Isyn_MD_left_to_PFC(:, 1:end_time), full_Isyn_MD_left_to_PFC, v_th);
    full_Isyn_MD_mid_to_PFC = save_to_cell(Isyn_MD_mid_to_PFC(:, 1:end_time), full_Isyn_MD_mid_to_PFC, v_th);
    full_Isyn_MD_right_to_PFC = save_to_cell(Isyn_MD_right_to_PFC(:, 1:end_time), full_Isyn_MD_right_to_PFC, v_th);
end
 
save("spike_time_temp", ...
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
    'full_PFC_shape_ensemble', ...
    'full_PFC_ori_ensemble', ...
    'full_con_inh_shape', ...
    'full_con_inh_ori', ...
    'full_PFC_bow', ...
    'full_PFC_rec', ...
    'full_PFC_oval', ...
    'full_PFC_SW', ...
    'full_PFC_N', ...
    'full_PFC_SE', ...
    'full_MD_left', ...
    'full_MD_mid', ...
    'full_MD_right', ...
    'full_PFC_left', ...
    'full_PFC_mid', ...
    'full_PFC_right', ...
    'full_output_left', ...
    'full_output_mid', ...
    'full_output_right', ...
    'full_Isyn_FS_left', ...
    'full_Isyn_FS_mid', ...
    'full_Isyn_FS_right', ...
    'full_Isyn_MD_left_to_PFC', ...
    'full_Isyn_MD_mid_to_PFC', ...
    'full_Isyn_MD_right_to_PFC', ...
    '-v7.3');

disp('Completed!');

%%
% function firing = get_firing(current, threshold, nanvalue)
% % current: membrane potiential
% % window_size: threshold for firing
% % thresh: predefined nan value
%     firing = (current > threshold) & (current ~= nanvalue);
% end

function output_cell = save_to_cell(current, cell, thresh)
    [neuo_id, firing_time] = find(current > thresh);
%     disp(neuo_id);
%     disp(firing_time);
    if isempty(neuo_id)
        output_cell = [cell; {NaN}];
    else
        output_cell = [cell; cat(2, neuo_id, firing_time)];
    end
end
