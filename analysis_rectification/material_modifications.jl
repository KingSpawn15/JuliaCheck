include("./init.jl");

#original values
mat = IndiumArsenide()
orig_alpha = mat.alpha
orig_kappa = mat.kappa
orig_gamma = mat.gamma
orig_me = mat.me
orig_mh = mat.mh
orig_hnew = mat.hnew
orig_alpha_gamma = mat.alpha_gamma
orig_n_eq = mat.n_eq
orig_n_eq = mat.n_eq

base="analysis_rectification/saved-matrices/materialanalysis_withcutoff_laserspot_40"

array_alpha =[ 1.2 * orig_alpha, 1 * orig_alpha, 0.8 * orig_alpha]
# array_kappa =[ 0.5 * orig_kappa , 1 * orig_kappa, 1.02 * orig_kappa]
# array_gamma =[ 0.5 * orig_gamma, 1 * orig_gamma, 1.01 * orig_gamma]
# array_me = [0.8 * orig_me, 1 * orig_me, 1.2 * orig_me]
# array_mh = [0.8 * orig_mh, 1 * orig_mh, 1.02 * orig_mh]
# array_hnew = [0.8 * orig_hnew, 1 * orig_hnew, 1.2 * orig_hnew]
# array_alpha_gamma = [0.8 * orig_alpha_gamma, 1 * orig_alpha_gamma, 1.2 * orig_alpha_gamma]
# array_n_eq = [0.98 * orig_n_eq, 1 * orig_n_eq, 1.2 * orig_n_eq]



dict = Dict("alpha"=>array_alpha)


# dict2 = Dict("kappa"=>array_kappa ,
# "me"=>array_me ,
# "mh"=>array_mh ,
# "hnew"=>array_hnew ,
# "alpha_gamma"=>array_alpha_gamma ,
# "n_eq"=>array_n_eq)

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

function generate_matrices(dict:: Dict{String, Vector{Float64}}, 
    dis_sp::Discretization, las::Laser, numericalp::NumericalParameters)

    for field in keys(dict)

        println(field)
        base_i = base*"_"*field
        mat = IndiumArsenide()
        for (ind,val) in enumerate(dict[field])
            
            setproperty!(mat,Symbol(field),val)
            mat.m_eq = mod_material.mass_equilibrium_carrier(mat);
            mat.epsilon_e = mat.epsilon_e = mod_material.photoexcited_electron_energy(mat);
            mat.v_t = mod_material.velocity_t(mat)
            mat.me0tilda = mod_material.photoelectron_mass(mat)

            @time _, _, interaction_v_pd=mod_cdem.interaction_potential_photodember(dis_sp, mat, las, numericalp; cutOff = -0.0)
            @time _, _, interaction_v_or=mod_cdem.interaction_potential_rectification(dis_sp, mat, las , elec, numericalp)

            jldopen(base_i*".jld2", "a+"; compress = true) do f
                f["photodember_"*string(ind)] = interaction_v_pd
                f["rectification_"*string(ind)] = interaction_v_or
            end

        end

        
    end

end

generate_matrices(dict, dis_sp, las, numericalp)
generate_matrices(dict2, dis_sp, las, numericalp)
# base_i = base*"_"*strarray[1]
# println(strarray[1])
# for (ind, alpha) in enumerate(array_alpha) 
    
#     mat = IndiumArsenide()
#     mat.alpha = alpha
#     mat.m_eq = mod_material.mass_equilibrium_carrier(mat);
#     mat.epsilon_e = mat.epsilon_e = mod_material.photoexcited_electron_energy(mat);
#     mat.v_t = mod_material.velocity_t(mat)
#     mat.me0tilda = mod_material.photoelectron_mass(mat)

#     @time _, _, interaction_v_pd=mod_cdem.interaction_potential_photodember(dis_sp, mat, las, numericalp)
#     @time _, _, interaction_v_or=mod_cdem.interaction_potential_rectification(dis_sp, mat, las , elec, numericalp)

#     jldopen(base_i*".jld2", "a+"; compress = true) do f
#         f["photodember_"*string(ind)] = interaction_v_pd
#         f["rectification_"*string(ind)] = interaction_v_or
#     end

# end
