using Plots
using JLD2, FileIO, CodecZlib
using .cdem_julia: mod_laser, mod_physicalconstants, mod_utils, mod_customtypes, mod_discrete,
mod_cdem, mod_eels, mod_electron, mod_material
using ..mod_customtypes:Discretization
using ..mod_electron:Electron

include("../src/utils.jl");

#original values
# las = mod_laser.set_laser!()
# mod_laser.set_laser!(;laser=las,
#     pulse_energy_experiment = 1 * 1e-9,
#     pulse_energy_gain_factor = 0.014,
#     laser_spot_fwhm = 40e-6,
#     theta_pol = 90*pi/180,
#     laser_pulse_time_fwhm = 300e-15);


# dis_sp = mod_discrete.discretization_setup(;x0 = 0.,
#     y0 = -1e-6,
#     d_xprime = 4e-2 * 3 * mod_utils.calculate_sigma(las.laser_spot_fwhm),
#     d_yprime = 4e-2 * 1e-6,
#     d_zprime = 4e-2 * 3 * mod_utils.calculate_sigma(las.laser_spot_fwhm),
#     xprime_max = 3 * mod_utils.calculate_sigma(las.laser_spot_fwhm),
#     yprime_max = 1e-6,
#     zprime_max = 3 * mod_utils.calculate_sigma(las.laser_spot_fwhm),
#     d_z = 2e-6,
#     zmax = 1e-4,
#     z_max = 50e-6,
#     t0 = -0.5e-12,
#     ddt = 10e-15,
#     delay_max = 3e-12,
#     fs = 2.4e15,
#     l = 1.08e4)

# mat = mod_material.IndiumArsenide()
# elec = mod_electron.set_electron!()
# mod_electron.set_electron!(;electron=elec,
#     electron_total_energy = 0.94,
#     electron_total_time_fs = 360.0,
#     electron_time_coherent_fwhm_fs = 50.0,
#     electron_theta = -7*pi/180,
#     electron_velocity_c = 0.7)

# numericalp = mod_customtypes.NumericalParameters(;tc_subsampling = 30,subsampling_factor = 60)

base="./analysis_angles/saved-matrices/CLEO"


v_struct = load(base*".jld2");

shift_pd = 18
shift_or = -15
alpha_pd_factor = 0.07
alpha_or_factor = 31.2

array_shift_pd = [14, 18, 22]

ang = 90
# interact_v_pd = circshift(v_struct["photodember"], (shift_pd,0));
# interact_v = alpha_or_0 .* circshift(v_struct["angle_"*string(ang)], (0,0))

w, e_w, t_w = mod_electron.energy_time_grid(elec,numericalp.subsampling_factor,dis_sp.energy, dis_sp.deltat);

psi_sub_array = zeros(length(t_w),length(e_w),length(array_shift_pd));
psi_incoherent_array =similar(psi_sub_array)
psi_sub_array_rect =similar(psi_sub_array)
psi_incoherent_array_rect =similar(psi_sub_array)
psi_sub_array_pd =similar(psi_sub_array)
psi_incoherent_array_pd =similar(psi_sub_array)



for (ind,shift_pd) in enumerate(array_shift_pd)
    
        
    v_pd = alpha_pd_factor * circshift(v_struct["photodember"],(shift_pd,0))
    v_or = alpha_or_factor * circshift(v_struct["angle_"*string(ang)],(shift_or,0))
    v_comb = v_pd .+ v_or
    
    @time loss_spectrum!(psi_sub_array,
    psi_incoherent_array,ind,v_comb, dis_sp, elec, w, t_w, e_w)
    
    @time loss_spectrum!(psi_sub_array_rect,
    psi_incoherent_array_rect,ind,v_or, dis_sp, elec ,w, t_w, e_w)

    @time loss_spectrum!(psi_sub_array_pd,
    psi_incoherent_array_pd,ind,v_pd, dis_sp, elec, w, t_w, e_w)
    

end

using GR 
plot_agg_incoherent = [];
plot_agg_sub = [];
gr(size=(1200,650));


for (ind, _) in enumerate(array_shift_pd)
    
    xl = ind==length(array_shift_pd);
    
    p_comb = my_heatmap_l(e_w, t_w .- 0.2, psi_incoherent_array[:,:,ind]; yl = true, xl=xl, rmargin=-5, lmargin = 10)
    p_or = my_heatmap(e_w, t_w .- 0.2, psi_incoherent_array_rect[:,:,ind]; xl=xl, rmargin=-5)
    p_pd = my_heatmap(e_w, t_w .- 0.2, psi_incoherent_array_pd[:,:,ind]; xl=xl, rmargin=10)
    
    push!(plot_agg_incoherent, p_comb, p_or, p_pd)
    
    p_comb = my_heatmap(e_w, t_w .- 0.2, psi_sub_array[:,:,ind]; xl=xl, rmargin=-5)
    p_or = my_heatmap(e_w, t_w .- 0.2, psi_sub_array_rect[:,:,ind]; xl=xl, rmargin=-5)
    p_pd = my_heatmap_o(e_w, t_w .- 0.2, psi_sub_array_pd[:,:,ind]; xl=xl)
    
    push!(plot_agg_incoherent, p_comb, p_or, p_pd)
end

base_sub = "shift_pd "
titlep = plot(title= base_sub*"Dependence", grid = false, showaxis = false, bottom_margin = -50Plots.px)
l = @layout [A{0.01h}; (length(array_shift_pd),6)]
pp= plot(titlep,plot_agg_incoherent...,layout=l)
savefig("./analysis_angles/saved-plots/CLEO"*string(ang)*".svg")
    
    







# function my_heatmap_o(e_w::Array{Float64,1}, t_w::Array{Float64,1}, psi::Array{Float64,2}; yl::Bool=false, xl::Bool=false)

#     p = heatmap(e_w, t_w, psi, c =:jet, 
#     aspect_ratio = 9/2.5,
#     xlims=[-4.5,4.5],
#     ylims=[-1,1.5],
#     colorbar=false,
#     xaxis=xl,
#     yaxis=yl,
#     yguidefontsize=12,
#     xguidefontsize=8)
#     yflip!(true)
#     if (yl)
#         ylabel!(L"\Delta t (ps)")
#     end
#     if (xl)
#         xlabel!("Energy Loss (ev)")
#     end

#     return p
# end

# function my_heatmap_l(e_w::Array{Float64,1}, t_w::Array{Float64,1}, psi::Array{Float64,2}; yl::Bool=false, xl::Bool=false, 
#     rmargin::Int64=-25, lmargin::Int64=5)

#     p = heatmap(e_w, t_w, psi, c =:jet, 
#     aspect_ratio = 9/2.5,
#     xlims=[-4.5,4.5],
#     ylims=[-1,1.5],
#     colorbar=false,
#     xaxis=xl,
#     yaxis=yl,
#     right_margin = eval(rmargin)Plots.mm,
#     left_margin = eval(lmargin)Plots.mm,
#     yguidefontsize=12,
#     xguidefontsize=8)
#     yflip!(true)
#     if (yl)
#         ylabel!(L"\Delta t (ps)")
#     end
#     if (xl)
#         xlabel!("Energy Loss (ev)")
#     end

#     return p
# end

# function my_heatmap(e_w::Array{Float64,1}, t_w::Array{Float64,1}, psi::Array{Float64,2}; yl::Bool=false, xl::Bool=false, rmargin::Int64=-25)

#     p = heatmap(e_w, t_w, psi, c =:jet, 
#     aspect_ratio = 9/2.5,
#     xlims=[-4.5,4.5],
#     ylims=[-1,1.5],
#     colorbar=false,
#     xaxis=xl,
#     yaxis=yl,
#     right_margin = eval(rmargin)Plots.mm,
#     yguidefontsize=12,
#     xguidefontsize=8)
#     yflip!(true)
#     if (yl)
#         ylabel!(L"\Delta t (ps)")
#     end
#     if (xl)
#         xlabel!("Energy Loss (ev)")
#     end

#     return p
# end