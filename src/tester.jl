include("Super.jl")
using .SuperModule.LaserMod
using .SuperModule.ElectronMod

las = set_laser!()
set_laser!(;laser=las,
    pulse_energy_experiment = 10 * 1e-9,
    pulse_energy_gain_factor = 0.014,
    laser_spot_fwhm = 40e-6,
    theta_pol = 45*pi/180,
    laser_pulse_time_fwhm = 50e-15)

elec = set_electron!()
set_electron!(;electron=elec,
    electron_total_energy = 1.1,
    electron_total_time_fs = 360.0,
    electron_time_coherent_fwhm_fs = 50.0,
    electron_theta = -7*pi/180,
    electron_velocity_c = 0.7)