% This script creates a dataset from powerflow simulations for SIGA
%%
clear all;
close all;
clc;

%% -----CONFIGS-----
%% Simulation parameters

% excel file of load information
file_name = 'load_data.xlsx';
sheet_name = 'load';

num_periods = 24*1;
case_name = 'case14';
mpc = loadcase(case_name);
[num_branches,~] = size(mpc.branch);

res_base = 60; 
res_target = 15;

opf_opt = true; % run optimal power flow?
load_randomness_opt = true; % change all loads by the following limits
% load_bus_bounds = [0.95,1.05]; % bounds for load changing
load_bus_bounds = [0.99,1.01]; % bounds for load changing
noise_opt = 1; % add noise to "measurements"?

if opf_opt
    pf_name = '_opf';
else
    pf_name = '_pf';
end

% attack options
t_att_start = 7;
t_att_end = 15;
attack_smooth_opt = 'smooth'; % 'step', 'linear', 'smooth'

% index of state variable to attack
idx_target = 10;
% target attack value for selected state variable
c_target_idx = [-0.35]';

% measurements to plot
mea_sel = 1:num_branches;

%% Setting up data

% reading data
file_name = fullfile('data',file_name);

data_table = readtable(file_name,'Sheet',sheet_name);

[num_samples,~] = size(data_table);

data_array = table2array(data_table);
data_array = data_array(:,1);

% normalizing data
load_factor_base = data_array ./ mean(sum(data_array,2));

% saving folder
save_path = fullfile('results');

if not(isfolder(save_path))
    mkdir(save_path);
end

%% Setting ts data

num_samples_base = 60/res_base;
dt_base = 1/num_samples_base;

num_samples_target = 60/res_target;
dt_target = 1/num_samples_target;

t_base = [0:dt_base:num_periods-1]';
load_factor_base = load_factor_base(1:num_periods);

t_target = [0:dt_target:num_periods-dt_target]';
load_factor_target = interp1(t_base,load_factor_base,t_target,'spline');

%% 'Create' more data based on different days

mpc = loadcase(case_name);
mpc.branch(:,9) = 0;

[num_buses,~] = size(mpc.bus);

data_array_norm = data_array ./ mean(sum(data_array,2));

[num_samples_data,~] = size(data_array_norm);
num_samples_day = 24*60/res_base; %(min)(1hr/60min)
num_days = num_samples_data/num_samples_day;
num_days = floor(num_days);

load_days = data_array_norm(1:num_days*num_samples_day);
load_days = reshape(load_days,[num_samples_day,num_days]); % [num_samples_day x num_days]

load_days_sel = load_days;

%% Target ts from different days
load_days_sel_target = zeros(length(t_target),num_buses);

for bus = 1:num_buses
    % load
    load_factor_base_tmp = load_days_sel(:,bus);
    load_days_sel_target(:,bus) = interp1(t_base,load_factor_base_tmp,t_target,'spline');

end

%% Setting timestamps

year_now = year(now);
month_now = month(now);
day_now = day(now);

% Base ts
h_start = 0;
h_end = num_periods-1;

date_start_base = datetime(year_now,month_now,day_now,h_start,0,0);
date_end_base = datetime(year_now,month_now,day_now,h_end,0,0);

ts_base = [date_start_base:hours(dt_base):date_end_base]';

% Target ts
date_start_target = datetime(year_now,month_now,day_now,h_start,0,0);
date_end_target = datetime(year_now,month_now,day_now,h_end+1,0,0) - minutes(res_target);

ts_target = [date_start_target:minutes(res_target):date_end_target]';

num_periods_target = length(ts_target);

%% Creating and Saving timeseries

col_names = ["Timestamp","Load"];
ts_base_table = table(ts_base,load_factor_base,'VariableNames',col_names);

col_names = ["Timestamp","Load"];
ts_target_table = table(ts_target,load_factor_target,'VariableNames',col_names);

% Saving all in a xlsx file
file_name = ['ts_',case_name,'.xlsx'];
file_name = fullfile(save_path, file_name);

if isfile(file_name)
    delete(file_name)
end

writetable(ts_base_table,file_name,'Sheet','ts_base')
writetable(ts_target_table,file_name,'Sheet','ts_target')

%% Setting up attack

mask_attack = (t_att_start <= t_target) & (t_target < t_att_end);
num_periods_att = sum(mask_attack);

alpha_attack_up_tmp = linspace(0,1,num_periods_att/2);
alpha_attack_down_tmp = linspace(1,0,num_periods_att/2);
alpha_attack_tmp = [alpha_attack_up_tmp,alpha_attack_down_tmp];
alpha_attack_linear = zeros(num_periods_target,1);
alpha_attack_linear(mask_attack) = alpha_attack_tmp;

delta_target = zeros(num_periods_target,1);
delta_target(mask_attack,1) = 1;

%% New way to compute 'alpha_attack'

idx_attack = find(mask_attack);
idx_attack_low = idx_attack(1);
idx_attack_high = idx_attack(end);
idx_attack_middle = (idx_attack_high + idx_attack_low)/2;
idx_attack_middle = round(idx_attack_middle);

num_attack_intervals = sum(mask_attack);
t_tmp = 1:length(mask_attack);

alpha_attack_smooth = exp(-((t_tmp-idx_attack_middle)/(0.5*num_attack_intervals/2)).^2);

if strcmp(attack_smooth_opt,'step')
    alpha_attack = mask_attack;

elseif strcmp(attack_smooth_opt,'linear')
    alpha_attack = alpha_attack_linear;

elseif strcmp(attack_smooth_opt,'smooth')
    alpha_attack = alpha_attack_smooth;   

else
        msg = ['Option: ',attack_smooth_opt, 'no available!'];
        error(msg)
end

%% System info data
mpc = loadcase(case_name);
mpc.branch(:,9) = 0;

[num_buses,~] = size(mpc.bus);
[num_branches,~] = size(mpc.branch);

sigma_mea = 0.02;

mpopt = mpoption('verbose', 0, 'out.all', 0); %0 for not printing

%% Compute estimation matrices
H = zeros(num_branches,num_buses);
idx_from = mpc.branch(:,1);
idx_to = mpc.branch(:,2);
x_branches = mpc.branch(:,4);

for branch = 1:num_branches
    from_tmp = idx_from(branch);
    to_tmp = idx_to(branch);
    x_tmp = x_branches(branch);

    H(branch,from_tmp) = 1 / x_tmp;
    H(branch,to_tmp) = -1 / x_tmp;
end

results = rundcopf(mpc, mpopt);

pf_base = results.branch(:,14) / results.baseMVA;
z_base = pf_base;
delta_base = results.bus(:,9)*pi/180; 

sigma_squared = ones(num_branches,1)*sigma_mea^2;
R = diag(sigma_squared);
Rinv = inv(R);

Hred = H(:,2:end);
G = Hred'*Rinv*Hred;
Ginv = inv(G);
K = Hred*Ginv*Hred'*Rinv;

delta_base_hat_red = Ginv*Hred'*Rinv*z_base;
delta_base_hat = [0; delta_base_hat_red];

z_base_hat = H*delta_base_hat;
z_base_hat_K = K*z_base;

% x_target = zeros(num_buses,1);
x_target = delta_base;
x_target(idx_target) = c_target_idx;

%%
pf_array = zeros(num_branches, num_periods_target);
delta_array = zeros(num_buses, num_periods_target);

pf_att_array = zeros(num_branches, num_periods_target);
delta_att_array = zeros(num_buses, num_periods_target);

% getting base pf data

f = waitbar(0,'Please wait...');
for t = 1:num_periods_target
    msg = ['Running simulations [',num2str(t),'/',num2str(num_periods_target),']...'];
    waitbar(t/num_periods_target,f,msg);
    
    % 1 - Obtain "z"
    % active load ts
    mpc_mod = mpc;  

    % changing only active power
    % % using one ts
    % mpc_mod.bus(:,3) = mpc_mod.bus(:,3)*load_factor_target(t);
    % using multiple days of ts
    mpc_mod.bus(:,3) = mpc_mod.bus(:,3).*load_days_sel_target(t,:)';
    
    if load_randomness_opt
        load_random_factor = unifrnd(load_bus_bounds(1),load_bus_bounds(2),[num_buses 1]);
        mpc_mod.bus(:,3) = mpc_mod.bus(:,3) .* load_random_factor;
    end

    % % adding solar generation
    % solar_alpha = 0;
    % mpc_mod.gen(5,2) = solar_days_sel_target(t,1)*solar_alpha;
    % mpc_mod.gen(5,9) = solar_days_sel_target(t,1)*solar_alpha;
    % % mpc_mod.gencost(5,5:end) = 0;
    
    if opf_opt
        results = rundcopf(mpc_mod, mpopt);
    else
        results = rundcpf(mpc_mod, mpopt);
    end
    
    % getting results
    success_tmp = results.success;
    pf_clean_tmp = results.branch(:,14) / results.baseMVA;
    delta_clean_tmp = results.bus(:,9)*pi/180;

    e_pf = sigma_mea * noise_opt * randn(num_branches,1);
    pf_tmp = pf_clean_tmp + e_pf;
    
    delta_tmp = Ginv*Hred'*Rinv*pf_tmp;
    delta_tmp = [0; delta_tmp];

    % FDIA  
    % Computing attack vector
    delta_target_tmp = delta_tmp;
    delta_target_tmp(idx_target) = c_target_idx;
    c_t = (delta_target_tmp - delta_tmp)*alpha_attack(t);                

    delta_att_tmp = delta_tmp + c_t;
    pf_att_tmp = H*delta_att_tmp;

    pf_hat_tmp = H*delta_tmp;

    % Saving results in arrays
    % base ts
    pf_array(:,t) = pf_hat_tmp*success_tmp;
    delta_array(:,t) = delta_tmp*success_tmp;

    % attacked ts
    pf_att_array(:,t) = pf_att_tmp;
    delta_att_array(:,t) = delta_att_tmp;

end
waitbar(1,f,'Finishing');
close(f)

% attack vector
a_pf_array = abs(pf_att_array - pf_array);

%% Saving results [BASE]

% setting data to save

% pf table [num_branches x num_samples]
col_names = "pf_" + [1:num_branches];
col_names = ["Timestamp", col_names];
row_names = "sim_" + [1:num_periods_target];

pf_table = [table(ts_target) array2table(pf_array')];
pf_table.Properties.VariableNames = col_names;

% v_ang table [num_buses x num_samples]
col_names = "v_ang_" + [1:num_buses];
col_names = ["Timestamp", col_names];

v_ang_table = [table(ts_target) array2table(delta_array')];
v_ang_table.Properties.VariableNames = col_names;

% Saving all in a xlsx file
file_name = ['ds_',case_name,pf_name,'_base','.xlsx'];
file_name = fullfile(save_path, file_name);

if isfile(file_name)
    delete(file_name)
end

writetable(pf_table,file_name,'Sheet','pf','WriteRowNames',true)
writetable(v_ang_table,file_name,'Sheet','v_ang')

%% Saving results [ATTACKED]

% setting data to save

% pf table [num_branches x num_samples]
col_names = "pf_" + [1:num_branches];
col_names = ["Timestamp", col_names];
row_names = "sim_" + [1:num_periods_target];

pf_table_att = [table(ts_target) array2table(pf_att_array')];
pf_table_att.Properties.VariableNames = col_names;

% v_ang table [num_buses x num_samples]
col_names = "v_ang_" + [1:num_buses];
col_names = ["Timestamp", col_names];

v_ang_table_att = [table(ts_target) array2table(delta_att_array')];
v_ang_table_att.Properties.VariableNames = col_names;

% a_pf table [num_branches x num_samples]
col_names = "pf_" + [1:num_branches];
col_names = ["Timestamp", col_names];
row_names = "sim_" + [1:num_periods_target];

a_pf_table_att = [table(ts_target) array2table(a_pf_array')];
a_pf_table_att.Properties.VariableNames = col_names;

% Saving all in a xlsx file
file_name = ['ds_',case_name,pf_name,'_attack','.xlsx'];
file_name = fullfile(save_path, file_name);

if isfile(file_name)
    delete(file_name)
end

writetable(pf_table_att,file_name,'Sheet','pf','WriteRowNames',true)
writetable(v_ang_table_att,file_name,'Sheet','v_ang')
writetable(a_pf_table_att,file_name,'Sheet','a_pf','WriteRowNames',true)

%% Plotting measurements [Real vs Fake]
fig = figure;
fig_name = ['measurements_ts_[',case_name,']'];
fig.Name=[fig_name];

hold on;

% Base measurments
mea_data = pf_array;
for mea = mea_sel
    mea_tmp = mea_data(mea,:);
    % line plot
    base_plt = plot(t_target, mea_tmp,...
        'color', 'blue',...
        'LineStyle','-',...
        'DisplayName','Base',...
        'Visible','on');
end

% Attack measurments
mea_data = pf_att_array;
for mea = mea_sel
    mea_tmp = mea_data(mea,:);
    % line plot
    att_plt = plot(t_target, mea_tmp,...
        'color', 'red',...
        'LineStyle','-.',...
        'DisplayName','Attacked');
end

xlabel('time, h')
ylabel('p.u.')

lgd = [base_plt,att_plt];

lgd = legend(lgd);
set(lgd);

%% Plotting state variables [Real vs Fake]
fig = figure;
fig_name = ['states_ts_[',case_name,']'];
fig.Name=[fig_name];

hold on;

% Base
mea_data = delta_array;
for mea = 1:num_buses
    mea_tmp = mea_data(mea,:);
    % line plot
    base_plt = plot(t_target, mea_tmp,...
        'color', 'blue',...
        'LineStyle','-',...
        'DisplayName','Base');
end

% Attack
mea_data = delta_att_array;
for mea = 1:num_buses
    mea_tmp = mea_data(mea,:);
    % line plot
    att_plt = plot(t_target, mea_tmp,...
        'color', 'red',...
        'LineStyle','-.',...
        'DisplayName','Attacked');
end

xlabel('time, h')
ylabel('rad')

lgd = [base_plt,att_plt];

lgd = legend(lgd);
set(lgd);

