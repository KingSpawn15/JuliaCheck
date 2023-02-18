include("./init_2.jl")

base="analysis_photodember/saved-matrices/interact_v_low_power_triple_pulse"

v_struct = load(base*".jld2");

t_c = dis_sp.t[dis_sp.t .> dis_sp.t0];
t_c_subsampled = t_c[1:numericalp.tc_subsampling:end];
w, e_w, t_w = mod_electron.energy_time_grid(elec,numericalp.subsampling_factor,dis_sp.energy, dis_sp.deltat);

interaction_v_pd = v_struct["photodember"]

# alpha_pd_0 =  .07;
# alpha_or_0 = 20 * 1.3 * 1.2;

alpha_pd_0 =  .5;
alpha_or_0 = 20 * 1.3 * 1.2 / 0.0017;
shift_or = -15;
shift_pd = 30;

interact_v_pd = circshift(interaction_v_pd, (shift_pd,0));

angle_array = sort(vcat([0:20:180;],[45, 135]))

function loss_spectrum!(psi_sub_array::Array{Float64,3},
    psi_incoherent_array::Array{Float64,3},
    ind::Int64,
    interact_v::Array{Float64,2}, 
    discretization::Discretization,
    elec::Electron)
        
    f_t = mod_eels.calculate_ft(discretization, interact_v , elec);
    
    psi = mod_eels.calculate_psi_coherent(dis_sp, elec, f_t);
    psi_sub = mod_eels.psi_subsampled(numericalp.subsampling_factor,psi, e_w);
    @views @inbounds psi_sub_array[:,:,ind] .= psi_sub
    @views @inbounds psi_incoherent_array[:,:,ind] .= mod_eels.incoherent_convolution_fast(psi_sub, w, t_w, e_w);
end

psi_sub_array = zeros(length(t_w),length(e_w),length(angle_array));
psi_incoherent_array = zeros(length(t_w),length(e_w),length(angle_array));

for (ind, angle) in enumerate(angle_array)
    
    
    interact_v = alpha_or_0 .* circshift(v_struct["angle_"*string(angle)], (shift_or,0)) .+
    alpha_pd_0 .* interact_v_pd;
    
    @time loss_spectrum!(psi_sub_array,
    psi_incoherent_array,ind,interact_v, dis_sp, elec)
    
end
# @code_warntype loss_spectrum!(psi_sub_array, psi_incoherent_array,5,interact_v_pd, dis_sp, elec);
plot_agg = [];

gr(size=(850,300));

for (ind, angle) in enumerate(angle_array)
    
    
    @views p = heatmap(e_w, t_w .- 0.2, psi_incoherent_array[:,:,ind], c =:jet, 
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
