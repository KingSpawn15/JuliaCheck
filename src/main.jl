using LazyGrids, Trapz, Interpolations, FFTW, JuliaInterpreter, Infiltrator
using Base.Threads
using TimerOutputs
using ThreadsX
include("constants.jl")
include("IndiumArsenide.jl")
include("init.jl")
include("Discretization.jl")
include("UTEMelectron.jl")
include("ChargeDynamics.jl")
include("EELS.jl")
include("index_calculator.jl")
# using Plots


function  default_parameters_2()
    # %DEFAULT_PARAMETERS Summary of this function goes here
    # %   Detailed explanation goes here

    laser = laser_params()
    laser.pulse_energy_experiment = 10 * 1e-9
    laser.pulse_energy_gain_factor = 0.014
    laser.laser_spot_fwhm = 40e-6
    laser.theta_pol = 45*pi/180
    laser.laser_pulse_time_fwhm = 50e-15
    laser.pulse_energy = laser.pulse_energy_gain_factor * laser.pulse_energy_experiment;
    laser.laser_spot_sigma = calculate_sigma(laser.laser_spot_fwhm);
    laser.laser_pulse_time_sigma = calculate_sigma(laser.laser_pulse_time_fwhm);#[s]

    discretization = discretization_params()
    discretization.x0 = 0.
    discretization.y0 = -1e-6
    discretization.ddt = 10e-15   
    discretization.delay_max = 2.5e-12
    discretization.fs = 2.4e15    
    # discretization.l = 1.0e4
    discretization.l = 2.4e4
    discretization.t0 = -0.5e-12
    discretization.xprime_max = 3 * calculate_sigma(laser.laser_spot_fwhm)
    discretization.d_xprime = 4e-2 * 3 * calculate_sigma(laser.laser_spot_fwhm)
    discretization.yprime_max = 1e-6
    discretization.d_yprime = 4e-2 * 1e-6
    discretization.zprime_max = discretization.xprime_max
    discretization.d_zprime = discretization.d_xprime
    discretization.ddz = 2e-6     
    discretization.zmax = 1e-4
    discretization.z_max = 30e-6

    utem = utem_params()
    utem.electron_total_energy = 1.1
    utem.electron_total_time_fs = 360
    utem.electron_time_coherent_fwhm_fs = 50
    utem.electron_theta = -7*pi/180
    utem.electron_velocity_c = 0.7

    numericalp = numerical_params()
    numericalp.tc_subsampling = 30
    numericalp.subsampling_factor = 60

    default_parameters = default_params()
    default_parameters.laser_parameters = laser;
    default_parameters.discretization_parameters = discretization;
    default_parameters.utem_parameters = utem;
    default_parameters.numerical_parameters = numericalp;

    return default_parameters
end

function calculate_sigma(fwhm::Float64)
    return fwhm./(2*sqrt(2*log(2)))
end

function excited_carriers(laser::laser_params, alpha::Float64, hnew::Float64)
            
    consts = constants_fundamental();
    Q_E = consts.Q_E;
    excited_volume = (pi/4)*laser.laser_spot_fwhm^2*alpha^(-1);#[m^3]
    n_exc = laser.pulse_energy/Q_E/hnew/excited_volume;#[m^-3]
    return n_exc

end

function fn()
    parameters = default_parameters_2()
    material = IndiumArsenide()
    dis_params = parameters.discretization_parameters
    dis = Discretization(dis_params)

    elec= ElectronSetup(parameters.utem_parameters)

    interaction_potential_photodember()

    @time t_c_subsampled, t_c, v_pd = interaction_potential_photodember(dis, material, parameters.laser_parameters ,
        parameters.numerical_parameters)


    return v_pd, t_c, t_c_subsampled, elec
end

# v_pd, t_c, t_c_subsampled, elec = fn();

# # # v_pd, t_c, t_c_subsampled, elec = fn()
parameters = default_parameters_2()
material = IndiumArsenide()
dis_params = parameters.discretization_parameters
dis = Discretization(dis_params)
elec= ElectronSetup(parameters.utem_parameters)
# # # # %Time domain Fourier transform
# @time f_t = calculate_ft(dis,1e-1 .* v_pd,elec);

# @time f_t_fast = calculate_ft_fast(dis,1e-1 .* v_pd,elec);

# w, e_w, t_w = energy_time_grid(elec,parameters.numerical_parameters.subsampling_factor,dis.energy, dis.deltat);

# @time psi = calculate_psi_coherent(dis,elec, f_t);
# psi_sub = psi_subsampled(parameters.numerical_parameters.subsampling_factor,psi, e_w);
# @time psi_incoherent = incoherent_convolution(psi_sub, w, t_w, e_w);
# @time psi_incoherent_fast = incoherent_convolution_fast(psi_sub, w, t_w, e_w);

# using Plots
# figure = (; resolution=(600, 400), font="CMU Serif")
# axis = (; xlabel="x", ylabel="y", aspect=1)
# heatmap(vec(e_w), vec(t_w), reverse(psi_incoherent_fast,dims = 1), c =:jet, aspect = 1)

