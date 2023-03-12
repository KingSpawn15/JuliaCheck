using Plots
using JLD2, FileIO, CodecZlib
using .cdem_julia: mod_laser, mod_physicalconstants, mod_utils, mod_customtypes, mod_discrete,
mod_cdem, mod_eels, mod_electron, mod_material

using ..mod_customtypes:Discretization, NumericalParameters
using ..mod_electron: Electron, set_electron!
using ..mod_material: Material, IndiumArsenide
using ..mod_laser: Laser

las = mod_laser.set_laser!()
mod_laser.set_laser!(;laser=las,
    pulse_energy_experiment = 1 * 1e-9,
    # pulse_energy_experiment = 10 * 1e-9,
    pulse_energy_gain_factor = 0.014,
    laser_spot_fwhm = 40e-6,
    theta_pol = 90*pi/180,
    laser_pulse_time_fwhm = 650e-15);
    # laser_pulse_time_fwhm = 50e-15);

dis_sp = mod_discrete.discretization_setup(;x0 = 0.,
    y0 = -1e-6,
    d_xprime = 4e-2 * 3 * mod_utils.calculate_sigma(las.laser_spot_fwhm),
    d_yprime = 4e-2 * 1e-6,
    d_zprime = 4e-2 * 3 * mod_utils.calculate_sigma(las.laser_spot_fwhm),
    xprime_max = 3 * mod_utils.calculate_sigma(las.laser_spot_fwhm),
    yprime_max = 1e-6,
    zprime_max = 3 * mod_utils.calculate_sigma(las.laser_spot_fwhm),
    d_z = 2e-6,
    zmax = 1e-4,
    z_max = 30e-6,
    t0 = -0.5e-12,
    ddt = 10e-15,
    delay_max = 3e-12,
    # delay_max = 2.5e-12,
    fs = 2.4e15,
    # l = 2.4e4)
    l = 1.08e4)

# mat = IndiumArsenide()
elec = set_electron!()
set_electron!(;electron=elec,
    electron_total_energy = 0.94,
    # electron_total_energy = 1.1,
    electron_total_time_fs = 360.0,
    electron_time_coherent_fwhm_fs = 50.0,
    electron_theta = -7*pi/180,
    electron_velocity_c = 0.7)

numericalp = NumericalParameters(;tc_subsampling = 30,subsampling_factor = 60)