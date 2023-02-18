include("./init.jl");

#original values
mat = IndiumArsenide()

pulse_energy_experiment_o = 1 * 1e-9;

base="analysis_photodember/saved-matrices/analysis_"
base = base*"laser_pulse_energy"

# array_pulse_energy_experiment = pulse_energy_experiment_o.*10.0.^(range(-3,stop=-2.5,length=3))

array_pulse_energy_experiment = [0.5 * pulse_energy_experiment_o, 1 * pulse_energy_experiment_o, 2 * pulse_energy_experiment_o]

function build_matrices_laser(array_pulse_energy_experiment::Array{Float64, 1}, base::String, las::Laser, 
    dis_sp::Discretization, mat::Material, numericalp::NumericalParameters, elec::Electron; startind::Int64 = 0)

    for (ind, pulse_energy_experiment) in enumerate(array_pulse_energy_experiment) 
        
        mod_laser.set_laser!(;laser=las,
        pulse_energy_experiment = pulse_energy_experiment,
        pulse_energy_gain_factor = 0.014,
        theta_pol = 90.0*pi/180);

        println(las.pulse_energy)
        @time _, _, interaction_v_pd=mod_cdem.interaction_potential_photodember(dis_sp, mat, las, numericalp)
        @time _, _, interaction_v_or=mod_cdem.interaction_potential_rectification(dis_sp, mat, las , elec, numericalp)

        jldopen(base*".jld2", "a+"; compress = true) do f
            f["photodember_"*string(ind + startind)] = interaction_v_pd
            f["rectification_"*string(ind  + startind)] = interaction_v_or
        end



    end

end

build_matrices_laser(array_pulse_energy_experiment[1:1],base,las,dis_sp,mat,numericalp,elec;startind=0)
build_matrices_laser(array_pulse_energy_experiment[2:3],base,las,dis_sp,mat,numericalp,elec;startind=1)



