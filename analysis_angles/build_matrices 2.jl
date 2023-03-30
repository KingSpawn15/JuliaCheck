using Plots
using JLD2, FileIO, CodecZlib
using .cdem_julia: mod_laser, mod_physicalconstants, mod_utils, mod_customtypes, mod_discrete,
mod_cdem, mod_eels, mod_electron, mod_material

using ..mod_customtypes:Discretization, NumericalParameters
using ..mod_electron: Electron, set_electron!
using ..mod_material: Material, IndiumArsenide
using ..mod_laser: Laser

include("../src/utils.jl");

directory="analysis_angles/saved-matrices/"
simulation_name = "testing"
database_name = "data_base"

las = mod_laser.set_laser!()
mod_laser.set_laser!(;laser=las,
    pulse_energy_experiment = 1 * 1e-9,
    pulse_energy_gain_factor = 0.014,
    laser_spot_fwhm = 100e-6,
    theta_pol = 90*pi/180,
    laser_pulse_time_fwhm = 650e-15,
    pulse_type = true);
   
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
    z_max = 80e-6,
    t0 = -0.5e-12,
    ddt = 10e-15,
    delay_max = 3e-12,
    fs = 2.4e15,
    l = 1.08e4)

    
mat = mod_material.IndiumArsenide()
elec = mod_electron.set_electron!()
mod_electron.set_electron!(;electron=elec,
    electron_total_energy = 0.94,
    electron_total_time_fs = 360.0,
    electron_time_coherent_fwhm_fs = 50.0,
    electron_theta = -7*pi/180,
    electron_velocity_c = 0.7)

numericalp = mod_customtypes.NumericalParameters(;tc_subsampling = 30,subsampling_factor = 60)
# numericalp = mod_customtypes.NumericalParameters(;tc_subsampling = 10, subsampling_factor = 60)

save_to_database(directory, database_name, simulation_name, dis_sp, las, numericalp)

base = directory*simulation_name
# angle_array = vcat([0:10:180;],[45 , 135])
angle_array = [90]

function rectification_builder(base::String, 
    angle_array::Array{Int64,1}, 
    las::Laser,
    dis_sp::Discretization,
    elec::Electron,
    numericalp::NumericalParameters)
 
    for angle in angle_array
        mod_laser.set_laser!(;laser = las,
        theta_pol = angle*pi/180)
        @time _, _, interaction_v=mod_cdem.interaction_potential_rectification(dis_sp, mat, las , elec, numericalp)

        jldopen(base*".jld2", "a+"; compress = true) do f
        f["angle_"*string(angle)] = interaction_v
        end

    end

end

@time t_c_subsampled, t_c, interaction_v_pd=mod_cdem.interaction_potential_photodember(dis_sp, mat, las, numericalp)


jldopen(base*".jld2", "a+"; compress = true) do f

    fieldname = "photodember";
    filename = base*".jld2";
    if haskey(f, fieldname)
        println("The field $fieldname already exists in $filename.")
        println("Do you want to skip (s) or overwrite (o)?")
        choice = String(readline())
        println(choice)
        if choice == "s"
            println("skipped")
            return
        elseif choice == "o"
            delete!(f, fieldname)
        else choice != "o"
            error("Invalid choice"*choice)
        end       
    end
    f[fieldname] = interaction_v_pd
    # f[fieldname] = interaction_v_pd

end




rectification_builder(base,angle_array[1:1],las,dis_sp,elec,numericalp)

rectification_builder(base,angle_array[2:end],las,dis_sp,elec,numericalp)
