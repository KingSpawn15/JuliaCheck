include("./init_2.jl");

las_pulse = 450.0
rectification_base = "analysis_photodember/saved-matrices/interact_v_low_power_triple_pulse" * string("_",Int(las_pulse),"fs");
# rectification_base ="saved_matrices/interact_v"
photodember_base = "analysis_photodember/saved-matrices/laseranalysis_withoutcutoff"

alpha_pd_0 =  .07;
alpha_or_0 = 20 * 1.3 * 1.2;

shift_pd = 20;
shift_or = -20;

#=
=#


base="saved_matrices/interact_v"
v_struct_or = load(rectification_base*".jld2");
v_struct_pd = load(photodember_base*".jld2");

v_pd = v_struct_pd["photodember_pulse_energy_experiment"*string(1)]
#=
=#

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

function loss_spectrum!(psi_sub_array::Array{Float64,2},
    psi_incoherent_array::Array{Float64,2},
    interact_v::Array{Float64,2}, 
    discretization::Discretization,
    elec::Electron)
        
    f_t = mod_eels.calculate_ft(discretization, interact_v , elec);
    
    psi = mod_eels.calculate_psi_coherent(dis_sp, elec, f_t);
    psi_sub = mod_eels.psi_subsampled(numericalp.subsampling_factor,psi, e_w);
    psi_sub_array .= psi_sub
    psi_incoherent_array .= mod_eels.incoherent_convolution_fast(psi_sub, w, t_w, e_w);
end

#=
=#

w, e_w, t_w = mod_electron.energy_time_grid(elec,numericalp.subsampling_factor,dis_sp.energy, dis_sp.deltat);
# angle_array = sort(vcat([0:20:180;],[45, 135]))
angle_array = sort([90])

psi_sub_array = zeros(length(t_w),length(e_w),length(angle_array));
psi_incoherent_array = zeros(length(t_w),length(e_w),length(angle_array));

psi_sub_array_rect = zeros(length(t_w),length(e_w),length(angle_array));
psi_incoherent_array_rect = zeros(length(t_w),length(e_w),length(angle_array));

psi_sub_pd = zeros(length(t_w),length(e_w));
psi_incoherent_pd = zeros(length(t_w),length(e_w));

v_pd = circshift(v_pd, (shift_pd,0));
interact_v_pd = 0.25 * alpha_pd_0 .* v_pd;

@time loss_spectrum!(psi_sub_pd,
    psi_incoherent_pd, interact_v_pd, dis_sp, elec)

for (ind, angle) in enumerate(angle_array)
    
    interact_v_or = (450.0/650.0) * alpha_or_0 .* circshift(v_struct_or["angle_"*string(angle)], (shift_or,0)) 
    interact_v = interact_v_or .+ interact_v_pd;

    @time loss_spectrum!(psi_sub_array,
    psi_incoherent_array,ind,interact_v, dis_sp, elec)

    @time loss_spectrum!(psi_sub_array_rect,
    psi_incoherent_array_rect,ind,interact_v_or, dis_sp, elec)
    
end


# plot_agg = [];

# gr(size=(850,300));

# for (ind, angle) in enumerate(angle_array)
    
    
#     @views p = heatmap(e_w, t_w .- 0.2, psi_incoherent_array[:,:,ind], c =:jet, 
#     aspect_ratio = 9/2.5,
#     xlims=[-4.5,4.5],
#     ylims=[-1,1.5],
#     colorbar=false,
#     xaxis=nothing,
#     yaxis=nothing,
#     right_margin = -0Plots.mm)
#     yflip!(true)
    
#     push!(plot_agg,p)
    
# end

# plot(plot_agg...,layout=(2,6))

plot_agg = [];

gr(size=(850,300));

p_pd = heatmap(e_w, t_w .- 0.2, psi_incoherent_pd, c =:jet, 
aspect_ratio = 9/2.5,
xlims=[-4.5,4.5],
ylims=[-1,1.5],
colorbar=false,
xaxis=nothing,
yaxis=nothing,
right_margin = -0Plots.mm)
yflip!(true)
# push!(plot_agg,p_pd)

for (ind, angle) in enumerate(angle_array)
    
    
    @views p_comb = heatmap(e_w, t_w .- 0.2, psi_incoherent_array[:,:,ind], c =:jet, 
    aspect_ratio = 9/2.5,
    xlims=[-4.5,4.5],
    ylims=[-1,1.5],
    colorbar=false,
    xaxis=nothing,
    yaxis=nothing,
    right_margin = -0Plots.mm)
    yflip!(true)
    
    @views p_or = heatmap(e_w, t_w .- 0.2, psi_incoherent_array_rect[:,:,ind], c =:jet, 
    aspect_ratio = 9/2.5,
    xlims=[-4.5,4.5],
    ylims=[-1,1.5],
    colorbar=false,
    xaxis=nothing,
    yaxis=nothing,
    right_margin = -0Plots.mm)
    yflip!(true)
    push!(plot_agg, p_comb, p_or, p_pd)
    
end


p_pd = heatmap(e_w, t_w .- 0.2, psi_sub_pd, c =:jet, 
aspect_ratio = 9/2.5,
xlims=[-4.5,4.5],
ylims=[-1,1.5],
colorbar=false,
xaxis=nothing,
yaxis=nothing,
right_margin = -0Plots.mm)
yflip!(true)


for (ind, angle) in enumerate(angle_array)
    
    
    @views p_comb = heatmap(e_w, t_w .- 0.2, psi_sub_array[:,:,ind], c =:jet, 
    aspect_ratio = 9/2.5,
    xlims=[-4.5,4.5],
    ylims=[-1,1.5],
    colorbar=false,
    xaxis=nothing,
    yaxis=nothing,
    right_margin = -0Plots.mm)
    yflip!(true)
    
    @views p_or = heatmap(e_w, t_w .- 0.2, psi_sub_array_rect[:,:,ind], c =:jet, 
    aspect_ratio = 9/2.5,
    xlims=[-4.5,4.5],
    ylims=[-1,1.5],
    colorbar=false,
    xaxis=nothing,
    yaxis=nothing,
    right_margin = -0Plots.mm)
    yflip!(true)
    push!(plot_agg, p_comb, p_or, p_pd)
    
end

plot(plot_agg...,layout=(length(2 * angle_array),3))
