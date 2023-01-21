using Plots
using JLD2, FileIO, CodecZlib
using .cdem_julia: mod_laser, mod_physicalconstants, mod_utils, mod_customtypes, mod_discrete,
mod_cdem, mod_eels, mod_electron, mod_material

las = mod_laser.set_laser!()
mod_laser.set_laser!(;laser=las,
    pulse_energy_experiment = 1 * 1e-9,
    pulse_energy_gain_factor = 0.014,
    laser_spot_fwhm = 40e-6,
    theta_pol = 90*pi/180,
    laser_pulse_time_fwhm = 650e-15);

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

base="saved_matrices/interact_v"


# for angle in angle_array
#     mod_laser.set_laser!(;laser = las,
#     theta_pol = angle*pi/180)
#     @time _, _, interaction_v=mod_cdem.interaction_potential_rectification(dis_sp, mat, las , elec, numericalp)

#     jldopen(base*".jld2", "a+"; compress = true) do f
#     f["angle_"*string(angle)] = interaction_v
#     end

# end

# @time t_c_subsampled, t_c, interaction_v_pd=mod_cdem.interaction_potential_photodember(dis_sp, mat, las, numericalp)
# jldopen(base*".jld2", "a+"; compress = true) do f
#     f["photodember"] = interaction_v_pd
# end

v_struct = load(base*".jld2");

t_c = dis_sp.t[dis_sp.t .> dis_sp.t0];
t_c_subsampled = t_c[1:numericalp.tc_subsampling:end];
w, e_w, t_w = mod_electron.energy_time_grid(elec,numericalp.subsampling_factor,dis_sp.energy, dis_sp.deltat);

interaction_v_pd = v_struct["photodember"]

alpha_pd_0 =  .07;
alpha_or_0 = 20 * 1.3 * 1.2;

interact_v_pd = circshift(interaction_v_pd, (18,0));

angle_array = sort(vcat([0:20:180;],[45, 135]))
plot_agg = []

bb=gr(size=(850,300))
for angle in angle_array
    interaction_v_or = v_struct["angle_"*string(angle)]
    interact_v_or = circshift(interaction_v_or, (-15,0));

    interact_v_combined = alpha_pd_0 .* interact_v_pd .+ alpha_or_0 .*interact_v_or

    @time f_t_fast = mod_eels.calculate_ft(dis_sp, interact_v_combined , elec);

    psi = mod_eels.calculate_psi_coherent(dis_sp, elec, f_t_fast);

    psi_sub = mod_eels.psi_subsampled(numericalp.subsampling_factor,psi, e_w);
    psi_incoherent_fast = mod_eels.incoherent_convolution_fast(psi_sub, w, t_w, e_w);

    # t_w_store = t_w;
    # t_w .-= 0.2;

    # figure = (; resolution=(400,400), font="CMU Serif")
    
    
    p = heatmap(e_w, t_w .- 0.2, psi_incoherent_fast, c =:jet, 
    aspect_ratio = 9/2.5,
    xlims=[-4.5,4.5],
    ylims=[-1,1.5],
    colorbar=false,
    xaxis=nothing,
    yaxis=nothing,
    right_margin = -0Plots.mm)
    yflip!(true)
    push!(plot_agg,p)
    
end
plot(plot_agg...,layout=(2,6))

# b = gr(size=(900,300))
# p1 = heatmap(e_w, t_w .- 0.2, rand(length(e_w),length(t_w)), c =:jet, 
#     aspect_ratio = 9/2.5,
#     xlims=[-4.5,4.5],
#     ylims=[-1,1.5],
#     colorbar=false,
#     xaxis=nothing,
#     yaxis=nothing,
#     right_margin = -0Plots.mm)
# yflip!(true)
# # push!(plot_agg,p1)

# plot(p1,p1,p1,p1,p1,p1,p1,p1,p1,p1,p1,p1,layout=(2,6))

# heatmap(e_w, t_w .- 0.2, rand(length(e_w),length(t_w)), c =:jet, 
#     aspect_ratio = 10/2.5,
#     xlims=[-4.5,4.5],
#     ylims=[-1,1.5],
#     colorbar=false,
#     ylabel=L"\Delta t",
#     xlabel=L"\rm{Energy Loss} (eV)")
