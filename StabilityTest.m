%% Stability Test of the numerical scheme
% Boudreau, 1997, page 347; Sweby, 1984, 

% dG/dt = w * dG/dx - k*G 

% Analytical solution (Berner 1980, page 151) 
% G = G_0 * exp(-k.*x/omega); % Eq 6-21
% --------------------------------------------------------
% 2021-9-18 adapted from previous code

clear

% Load the model parameters, choose between 1-8, 4 is baseline model
% run the baseline model for stability test
load 'parameters.mat'
model_run = 4;

% select the corresponding sed. rate timeseries, and associated parameters
T = 100000; % [yr], time of simulation

% define constant sed. rate
w_tr = omega_list (model_run);
w_ts = ones(1,T).* w_tr;

k_tr = k_list (model_run);
Zod_tr = depth_O2depletion(model_run);
OC_mean_tr = OC_preserved(model_run);
Da_tr = Da (model_run);

% clear unnessassry variables
clearvars Da depth_O2depletion k_list omega_list OC_preserved

depth_total = round(sum (w_ts),2);       % total depth accumulated

G_0 = 1;       % OC at sediment water interface
dx = 0.01;      % grid size, can't be larger to insure stability


% calculate the corresponding initial condition (steady state profile)
domain = round(Zod_tr+depth_total,2);      % define model domain length
x_grid = linspace(0,domain-dx,domain/dx);  % define model domain grid
G_init = G_0 * exp ((-k_tr/w_tr).*x_grid);       % steady-state solution
G_init(Zod_tr/dx+2:domain/dx)= G_init(Zod_tr/dx+1); % constant below Z_od


%%
G_transient = Sweby(G_init, w_ts(1), k_tr, dx, domain, G_0 , Zod_tr);
for i =2:T
    G_transient = Sweby(G_transient, w_ts(i), k_tr, dx, domain, G_0 ,Zod_tr);
end 

save 'StabilityTest.mat'


