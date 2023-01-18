include("CDEM.jl")

using .SuperModule.LaserMod
using .SuperModule.ElectronMod
using .SuperModule.DiscretizationEnergyTime


las = set_laser!()
set_laser!(;laser=las,
    pulse_energy_experiment = 10 * 1e-9,
    pulse_energy_gain_factor = 0.014,
    laser_spot_fwhm = 40e-6,
    theta_pol = 45*pi/180,
    laser_pulse_time_fwhm = 50e-15)

dis_sp = space_config(;x0 = 0.,
    y0 = -1e-6,
    ddt = 10e-15   ,
    delay_max = 2.5e-12,
    xprime_max = 3 * calculate_sigma(las.laser_spot_fwhm),
    d_xprime = 4e-2 * 3 * calculate_sigma(las.laser_spot_fwhm),
    yprime_max = 1e-6,
    d_yprime = 4e-2 * 1e-6,
    zprime_max = xprime_max,
    d_zprime = d_xprime,
    ddz = 2e-6,
    zmax = 1e-4,
    z_max = 30e-6
)

dis_et = energytime_config(;ddt = 10e-15,
    delay_max = 2.5e-12,
    fs = 2.4e15,
    l = 2.4e4,
    t0 = -0.5e-12)



elec = set_electron!()
set_electron!(;electron=elec,
    electron_total_energy = 1.1,
    electron_total_time_fs = 360.0,
    electron_time_coherent_fwhm_fs = 50.0,
    electron_theta = -7*pi/180,
    electron_velocity_c = 0.7)