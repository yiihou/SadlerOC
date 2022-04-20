# SadlerOC

codes and data of a reactive-transport model incorporating stochastic sedimentation rates

Model Parameters (Supp Tables 2 and 3) are given/calculated in: 

    'Calc_O2Depletion.m', saved as 'parameters.mat'
    
Stochastic sedimentation rate generated from:

    Exponential/Pareto: matlab built in functions
    Double Truncated Pareto: VGAM package - R, modified from code in Trampush & Hajek, 2017

All model runs conducted listed in Supp Table 4,

    sed rate time-sereis given in 'SedRate.mat'
    Raw Output for each simulations in 'RawResults' folder
    Compiled model output in Comp21.mat

Reactive-transport model implementation code:

    'Implementation_*.m':
    with numerical schemes: 'Sweby.m' and 'Sweby_tdk.m' (for reactive continuum kinetics)
  
Spectral analysis:

    linkedearth.github.io/pyleoclim_util/
