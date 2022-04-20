% ---------------------------------------
% Calculating:  
% (1) the depth of oxygen depletion with the analytical solutions
% given in Berner, 1980 (page 151)
% (2) Da number
% (3) fraction of OC buried
% ----------------------------------------
% 2021-09-17, adapted from previous codes


% C: conc. of oxidant 
% G: conc. of OC
% F: (1-phi)/phi * density (conversion between solid and dissolved phases
% L: stoichiometric coefficient of the oxidation reaction
% Ds: difussivity
% k: reaction rate constant
% omega: sedimentation rate


%%% C = (omega^2 * F * L * G_0)/(omega^2 + k*Ds) * exp(-k.*x/omega) + C_inf; % Eq 6-20
%%% G = G_0 * exp(-k.*x/omega); % Eq 6-21
%%% C_inf = C_0 - (omega^2 * F * L * G_0)/(omega^2 + k*Ds); % Eq 6-22

clear

phi = 0.8; % dimensionless, porosity
L = 1; % stoichiometric coefficient
density = 2.5; % g/cm3
F = (1-phi)/phi * density; % g/cm3
Ds = 0.955; % cm2/d Soetaert et al., 1996 Molecular diffusion coefficient for O2 at 0 degreeC
Ds = Ds * 365; % cm2/yr 

% NOTE: the choice of OC and oxidant concentration is sort of arbitrary -
% they are fixed for almost all our model simulations (except when
% experimenting with oscillating boundary conditions and when simulating
% the dilution effect by sedimentation), so they don't impact our 
% interpretations on model results. They are just one of the many many
% possible combinations of OC and oxidant concentration at the
% sediment-water interface. 

G_0 = 1; % moles per unit mass of total solids = wt%/molar mass
C_0 = 0.0001*7; % molar

x = 0:0.1:100; % cm, spatial domain

% Define functions to calculate:
% C_inf: asymptotic concentration of oxidant at depth, can be negative,
% see Berner 1980 for explanation
% C_profile: depth profile of oxidant
% G_profile: depth profile of OC

Calc_C_inf = @(omega, F, L, G_0, k, Ds, C_0) C_0 - (omega^2 * F * L * G_0)/(omega^2 + k*Ds);
Calc_G_profile = @(omega, k, x, G_0) G_0 * exp(-k.*x/omega);
Calc_C_profile = @(omega, F, L, G_0, k, Ds, C_inf)...
    (omega^2 * F * L * G_0)/(omega^2 + k*Ds) * exp(-k.*x/omega) + C_inf;


%% Pt 1: Calculate the z_OD, Da, and fraction of OC buried for any given k, omega

k = 0.002; % yr-1, reaction rate constant
omega = 0.03;  % cm yr-1, sedimentation rate 


% Calculate C_inf, C_profile, and G_profile
C_inf = Calc_C_inf(omega, F, L, G_0, k, Ds, C_0);
G_profile = Calc_G_profile(omega, k, x, G_0);
C_profile = Calc_C_profile(omega, F, L, G_0, k, Ds, C_inf);

% find the depth where all oxidant has been consumed
O2depletion = find(C_profile < 0, 1);
depth_O2depletion = x(O2depletion);

% find the asymptotic OC fraction and calculate Da
OC_preserved = G_profile(O2depletion);
Da = k * depth_O2depletion /omega;


%% Pt 2: Calculate from a list of selected k, omega values

% This list covers the range of parameters in Figure 3a, 3b where mean OC burial
% efficiency and variance of the simulated record are plotted against Da of
% the system, respectively. In this list, min(Da) = 0.25, max(Da) = 2.65
 

k_list = [0.0005 0.001 0.002 0.0015 0.002 0.002 0.002 0.002];
omega_list = [0.03 0.03 0.04 0.03 0.034 0.031 0.03 0.029];

for i = 1: length(k_list)
    % Calculate C_inf, C_profile, and G_profile
    C_inf = Calc_C_inf(omega_list(i), F, L, G_0, k_list(i), Ds, C_0);
    G_profile = Calc_G_profile(omega_list(i), k_list(i), x, G_0);
    C_profile = Calc_C_profile(omega_list(i), F, L, G_0, k_list(i), Ds, C_inf);

    % find the depth where all oxidant has been consumed
    O2depletion(i) = find(C_profile < 0, 1);
    depth_O2depletion(i) = x(O2depletion(i));

    % find the asymptotic G and calculate Da
    OC_preserved(i) = G_profile(O2depletion(i));
    Da(i) = k_list(i) * depth_O2depletion(i) /omega_list(i);

end

%% save lists of k, omega and calculated Da, z_OD, OC burial to file

clearvars -except k_list omega_list Da depth_O2depletion OC_preserved
save 'parameters.mat'

