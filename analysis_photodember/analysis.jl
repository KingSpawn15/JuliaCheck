include("./init_2.jl");

w, e_w, t_w = mod_electron.energy_time_grid(elec,numericalp.subsampling_factor,dis_sp.energy, dis_sp.deltat);


las_pulse_array = [50.0, 100.0, 150.0, 250.0, 450.0, 650.0]
# las_pulse_array = [50.0]


for las_pulse in las_pulse_array
    base="analysis_photodember/saved-matrices/interact_v_low_power_single_pulse" * string("_",Int(las_pulse),"fs")
    angle_array = vcat([0:10:180;],[45 , 135])
    for angle in angle_array
        mod_laser.set_laser!(;laser = las,
        theta_pol = angle*pi/180,
        laser_pulse_time_fwhm = las_pulse * 1e-15)
        @time _, _, interaction_v=mod_cdem.interaction_potential_rectification(dis_sp, mat, las , elec, numericalp)

        jldopen(base*".jld2", "a+"; compress = true) do f
        f["angle_"*string(angle)] = interaction_v
        end

    end


end


# @time t_c_subsampled, t_c, interaction_v_pd=mod_cdem.interaction_potential_photodember(dis_sp, mat, las, numericalp)
# jldopen(base*".jld2", "a+"; compress = true) do f
#     f["photodember"] = interaction_v_pd
# end


# @time _, _, interact_v=mod_cdem.interaction_potential_rectification(dis_sp, mat, las , elec, numericalp)

# f_t = mod_eels.calculate_ft(discretization, interact_v , elec);
# psi = mod_eels.calculate_psi_coherent(dis_sp, elec, f_t);
# psi_sub = mod_eels.psi_subsampled(numericalp.subsampling_factor,psi, e_w);
# psi_incoherent = mod_eels.incoherent_convolution_fast(psi_sub, w, t_w, e_w);

