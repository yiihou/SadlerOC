%% Implementation of the model with exponential distributed sed. rates, 
%% and time-dependent reaction rate constant:
% dG/dt = w * dG/dx - k*G 

% Analytical solution (Berner 1980, page 151) 
% G = G_0 * exp(-k.*x/omega); % Eq 6-21

% k =  0.125*(2 + t)^(-1); (the original should be (0.0003+t),
% according to Boudreau et al. comments on Rothman&Forney 2007 Science
% but for model stability, we modified the "age" parameter
% --------------------------------------------------------
% 2021-9-18 adapted from previous code

clear

load 'SedRate.mat'
load 'parameters_20210917.mat'

% 4 is the baseline model
model_run = 4;

sed_run = 23;
w_ts = omega{sed_run};

% select the corresponding sed. rate timeseries, and associated parameters
T = length (w_ts); % [yr], time of simulation, should be 100,000
w_tr = omega_list (model_run);
w_ts = w_ts .* (w_tr./mean(w_ts));

% clear unnessassry variables
clearvars Da depth_O2depletion k_list omega_exp omega_list OC_preserved...
            double*
      
% boundary condition and discretization
G_0 = 1;       % OC at sediment water interface 
dx = 0.01;     % grid size, can't be larger to insure stability


%%

% calculate the transient depth accumulated
depth_trans = zeros(T,1);
for i = 1:T
    depth_trans (i)= sum(w_ts(1:i));
end
depth_total = max(depth_trans);

domain = depth_total;  % model domain 
x_grid = linspace(0,domain-dx,domain/dx);

[x,y] = ode45(@(x,y) ssprofile(x,y,w_tr), x_grid, G_0, w_tr);

% initial conditions and k array 
k_array = 0.125*(2 + x_grid ./ w_tr).^(-1);

G_init = y;
k_array = k_array';

n_domain = round(domain/dx);

%% Transient Solution

if w_ts(1)>=0
    G_transient = Sweby_tdk(G_init, w_ts(1), dx, domain, G_0, k_array);
else
    G_transient = G_init;
    n_par =round(w_ts(1)/dx);
    G_transient(1:n_domain+n_par) = G_transient(-n_par+1:n_domain);
        
    G_transient(n_domain+n_par+1:n_domain) = G_transient(n_domain);
        
    % G_transient (1: n_Zod_tr) = G_transient (1: n_Zod_tr)-G_transient (1: n_Zod_tr).* k_tr;
    G_transient = G_transient - G_transient .* k_array;
end

for i =2:T
    if w_ts(i)>=0        
        G_transient = Sweby_tdk(G_transient, w_ts(i), dx, domain, G_0, k_array);
    else
        n_par =round(w_ts(i)/dx);
        G_transient(1:n_domain+n_par) = G_transient(-n_par+1:n_domain);
        
        G_transient(n_domain+n_par+1:n_domain) = G_transient(n_domain);
        G_transient = G_transient - G_transient .* k_array;
        % G_transient (1: n_Zod_tr) = G_transient (1: n_Zod_tr)-G_transient (1: n_Zod_tr).* k_tr;
    end   
end 

%% Save data to file

filename = append("TDK_erosion_ModelRun_",num2str(sed_run),".mat");
save(filename)

%%
function dydx = ssprofile(x,y,w)
t = x/w;
k =  0.125*(2 + t)^(-1);
dydx = (-k*y)/w;
end
