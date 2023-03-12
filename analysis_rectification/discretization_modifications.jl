include("../src/cdem.jl")
include("./init.jl");

#original values
mat = IndiumArsenide()


base="analysis_rectification/saved-matrices/zmax_analysis_laser_spot_40_laser_time_50"



mod_laser.set_laser!(;laser=las,
laser_pulse_time_fwhm  = 50 * 1e-15,
laser_spot_fwhm = 40e-6,
theta_pol = 90.0*pi/180);

set_electron!(;electron=elec,
    electron_total_energy = 0.94,
    # electron_total_energy = 1.1,
    electron_total_time_fs = 360.0,
    electron_time_coherent_fwhm_fs = 50.0,
    electron_theta = -7*pi/180,
    electron_velocity_c = 0.7)

function discretization_mod(;x0::Float64 = 0.0, y0::Float64 = -1e-6, z_max::Float64 = 30e-6)

    return mod_discrete.discretization_setup(;x0 = x0,
    y0 = y0,
    d_xprime = 4e-2 * 3 * mod_utils.calculate_sigma(las.laser_spot_fwhm),
    d_yprime = 4e-2 * 1e-6,
    d_zprime = 4e-2 * 3 * mod_utils.calculate_sigma(las.laser_spot_fwhm),
    xprime_max = 3 * mod_utils.calculate_sigma(las.laser_spot_fwhm),
    yprime_max = 1e-6,
    zprime_max = 3 * mod_utils.calculate_sigma(las.laser_spot_fwhm),
    d_z = 2e-6,
    zmax = 1e-4,
    z_max = z_max,
    t0 = -0.5e-12,
    ddt = 10e-15,
    delay_max = 3e-12,
    fs = 2.4e15,
    l = 1.08e4)

end



# y0_array = [-1e-7, -1e-6, -1e-5]
# x0_array = [-1e-7, 0.0, 1e-7]
# z_max_array = [20e-6, 30e-6, 50e-6]

# array_discretization_x0 = []
# array_discretization_y0 = []
# array_discretization_z_max = []

# dict = Dict("y0" => y0_array,
# "x0" => x0_array,
# "z_max" => z_max_array
# )

# y0_array = [-1e-5, -1e-6, -1e-7]
# x0_array = [-1e-7, 0.0, 1e-7]
z_max_array = [20e-6, 30e-6, 50e-6]


# dict = Dict("y0" => y0_array,
# )

dict = Dict("z_max" => z_max_array
)

function get_discretization(array_discretization::Vector{Any}, ss::Int64)::Discretization
    return array_discretization[ss]
end

function generate_matrices(dict:: Dict{String, Vector{Float64}}, 
    mat::Material, las::Laser, numericalp::NumericalParameters, elec::Electron)

    for field in keys(dict)

        println(field)
        base_i = base*"_"*field
    
        for (ind,val) in enumerate(dict[field])
            
     
            ff = Symbol(field)
            dis_sp = eval(:(discretization_mod(;$ff = $val)))
            @time _, _, interaction_v_pd=mod_cdem.interaction_potential_photodember(dis_sp, mat, las, numericalp; cutOff = 0.0)
            @time _, _, interaction_v_or=mod_cdem.interaction_potential_rectification(dis_sp, mat, las , elec, numericalp)

            jldopen(base_i*".jld2", "a+"; compress = true) do f
                f["photodember_"*string(ind)] = interaction_v_pd
                f["rectification_"*string(ind)] = interaction_v_or
            end

        end

        
    end

end

generate_matrices(dict,  mat, las, numericalp, elec)
