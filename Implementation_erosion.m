%% Implementation of the model with exponential distributed sed. rates:
% dG/dt = w * dG/dx - k*G 

% Analytical solution (Berner 1980, page 151) 
% G = G_0 * exp(-k.*x/omega); % Eq 6-21
% --------------------------------------------------------
% 2021-9-18 adapted from previous code

clear
% Load sed. rate time series (erosion) and model parameters
load 'SedRate.mat' % load sedimentation rates
load 'parameters.mat'

% select the set of model parameters will be run, between 1-8, 4 is baseline model
model_run = 4; 

% select the corresponding sed. rate timeseries, and associated parameters
sed_run = 11;

w_ts = omega{sed_run};
T = length (w_ts); % [yr], time of simulation, should be 100,000

k_tr = k_list (model_run);
w_tr = omega_list (model_run);

Zod_tr = depth_O2depletion(model_run);
OC_mean_tr = OC_preserved(model_run);
Da_tr = Da (model_run);
% clear unnessassry variables
clearvars Da depth_O2depletion k_list omega_double omega_list OC_preserved

% calculate the transient depth accumulated
depth_trans = zeros(T,1);
for i = 1:T
    depth_trans (i)= sum(w_ts(1:i));
end
depth_total = max(depth_trans);

% boundary condition and discretization
G_0 = 1;       % OC at sediment water interface 
dx = 0.01;     % grid size, can't be larger to insure stability

% calculate the corresponding initial condition (steady state profile)
domain = round(Zod_tr+depth_total,2);      % define model domain length
x_grid = linspace(0,domain-dx,domain/dx);  % define model domain grid
G_init = G_0 * exp ((-k_tr/w_tr).*x_grid);       % steady-state solution
G_init(round(Zod_tr/dx)+2:domain/dx)= G_init(round(Zod_tr/dx)+1); % constant below Z_od

% avoid glitch
n_domain = round(domain/dx);
n_Zod_tr = round(Zod_tr/dx);


%% Transient solutions
if w_ts(1)>=0
    G_transient = Sweby(G_init, w_ts(1), k_tr, dx, domain, G_0 , Zod_tr);
else
    G_transient = G_init;
    n_par =round( w_ts(1)/dx);
    G_transient(1:n_domain+n_par) = G_transient(-n_par+1:n_domain);
        
    G_transient(n_domain+n_par+1:n_domain) = G_transient(n_domain);
        
    G_transient (1: n_Zod_tr) = G_transient (1: n_Zod_tr)-G_transient (1: n_Zod_tr).* k_tr;

end  

for i = 2:T
    if  w_ts(i)>=0
        G_transient = Sweby(G_transient, w_ts(i), k_tr, dx, domain, G_0 , Zod_tr);
    else
        n_par =round( w_ts(i)/dx);
        G_transient(1:n_domain+n_par) = G_transient(-n_par+1:n_domain);
        
        G_transient(n_domain+n_par+1:n_domain) = G_transient(n_domain);
        
        G_transient (1: n_Zod_tr) = G_transient (1: n_Zod_tr)-G_transient (1: n_Zod_tr).* k_tr;
        
    end
end

%% save data to file

filename = append("DoublePareto_ModelRun",num2str(model_run),"_", num2str(sed_run),".mat");
save(filename) 


