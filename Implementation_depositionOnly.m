%% Implementation of the model with deposition-only (only positive sed. rates):
% dG/dt = w * dG/dx - k*G 

% Analytical solution (Berner 1980, page 151) 
% G = G_0 * exp(-k.*x/omega); % Eq 6-21
% --------------------------------------------------------
% 2021-9-18 adapted from previous code

clear

% load sedimentation rate time-series and parameters
load 'SedRate.mat'    % load sedimentation time-series that has only positive sedimentation rates
load 'parameters.mat' % load model parameters

% Define which time-series/set of parameters will be run, between 1-8, 4 is baseline model
model_run = 4; 

% select one of the sed. rate timeseries
% w_ts = omega_exp{model_run};
w_ts = omega{model_run};
T = length (w_ts); % [yr], time of simulation, should be 100,000
k_tr = k_list (model_run);
w_tr = omega_list (model_run);
Zod_tr = depth_O2depletion(model_run);
OC_mean_tr = OC_preserved(model_run);
Da_tr = Da (model_run);

% clear unnessassry variables
clearvars Da depth_O2depletion k_list omega_exp omega_list OC_preserved

% calculate model domain (total depth accumulated)
depth_total = round(sum(w_ts),2);       

% boundary condition and discretization
G_0 = 1;       % OC at sediment water interface 
dx = 0.01;     % grid size, can't be larger to insure stability

% calculate the corresponding initial condition (steady state profile)
domain = round(Zod_tr+depth_total,2);      % define model domain length
x_grid = linspace(0,domain-dx,domain/dx);  % define model domain grid
G_init = G_0 * exp ((-k_tr/w_tr).*x_grid);       % steady-state solution
G_init(Zod_tr/dx+2:domain/dx)= G_init(Zod_tr/dx+1); % constant below Z_od

%% Transient solutions 

G_transient = Sweby(G_init, w_ts(1), k_tr, dx, domain, G_0 , Zod_tr);
for i =2:T
    G_transient = Sweby(G_transient, w_ts(i), k_tr, dx, domain, G_0 , Zod_tr);
end

%% save data to file

filename = append("Exponential_ModelRun",num2str(model_run),".mat");
save(filename) 
    
