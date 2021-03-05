% Solver for 1D isothermal and isobaric dispersion PFR
% 
% Solves the governing eqautions for any kinetic mechanism
% with any number of species.
% Can be used to observe the temporal evolution of the system 
% or to analyze the Steady-State solution.
% 
% Numerical detailes:
% - Finite Volumes domain discretization
% - Central Differencing Scheme for spatial discretization
% - Backward Differencing Scheme for time discretization
%
% Enrico Bussetti 210305
%
% -------------------------------------------------------------------------

function C = dispersion_pfr(opt, tspan)

    if ~exist('opt', 'var')
        opt = struct;       % Define the options structure
    end
    
    if ~exist('dt', 'var')
        dt = min();        % Define the integration time step
    end
    
    if ~exist('tspan', 'var')
        tspan = [0 10000*dt];       % Define the time range to be covered
    end




end

