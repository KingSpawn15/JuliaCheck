using Plots
using JLD2, FileIO, CodecZlib
using .cdem_julia: mod_laser, mod_physicalconstants, mod_utils, mod_customtypes, mod_discrete,
mod_cdem, mod_eels, mod_electron, mod_material
using ..mod_customtypes:Discretization
using ..mod_electron:Electron



base="./analysis_angles/saved-matrices/CLEO"
sub_base = "CLEO"
v_struct = load(base*".jld2");

t_c = dis_sp.t[dis_sp.t .> dis_sp.t0];
t_c_subsampled = t_c[1:numericalp.tc_subsampling:end];
w, e_w, t_w = mod_electron.energy_time_grid(elec,numericalp.subsampling_factor,dis_sp.energy, dis_sp.deltat);

angle_array = sort(vcat([0:20:180;],[45, 135]))

function loss_spectrum!(psi_sub_array::Array{Float64,3},
    psi_incoherent_array::Array{Float64,3},
    ind::Int64,
    interact_v::Array{Float64,2}, 
    discretization::Discretization,
    elec::Electron,
    w::Array{Float64,2}, t_w::Array{Float64,1}, e_w::Array{Float64,1})
    
    f_t = mod_eels.calculate_ft(discretization, interact_v , elec);
    
    psi = mod_eels.calculate_psi_coherent(dis_sp, elec, f_t);
    psi_sub = mod_eels.psi_subsampled(numericalp.subsampling_factor,psi, e_w);
    @views @inbounds psi_sub_array[:,:,ind] .= psi_sub
    @views @inbounds psi_incoherent_array[:,:,ind] .= mod_eels.incoherent_convolution_fast(psi_sub, w, t_w, e_w);
end


function plot_all_angles(array_alpha_pd_0::Array{Float64,1},
    array_alpha_or_0::Array{Float64,1},
    array_delay_pd::Array{Int64,1},
    v_struct::Dict{String, Any}, angle_array::Array{Int64, 1}, 
    elec::Electron, dis_sp::Discretization,  
    w::Array{Float64,2}, t_w::Array{Float64,1}, e_w::Array{Float64,1},
    sub_base::String)

    delay_or = -15;

    for alpha_pd_0 in array_alpha_pd_0
        for alpha_or_0 in array_alpha_or_0
            for delay_pd in array_delay_pd
            
                interact_v_pd = circshift(v_struct["photodember"], (delay_pd,0));
                psi_sub_array = zeros(length(t_w),length(e_w),length(angle_array));
                psi_incoherent_array = zeros(length(t_w),length(e_w),length(angle_array));
            
                plot_agg = [];
                gr(size=(850,300));
            
                for (ind, angle) in enumerate(angle_array)
                    
                    
                    interact_v = alpha_or_0 .* circshift(v_struct["angle_"*string(angle)], (delay_or,0)) .+
                    alpha_pd_0 .* interact_v_pd;
                    
                    @time loss_spectrum!(psi_sub_array,
                    psi_incoherent_array,ind,interact_v, dis_sp, elec, w, t_w, e_w)
            
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
                savefig("./analysis_angles/saved-plots/"*sub_base*"delay_pd_"*string(delay_pd)*
                "alpha_or_0"*string(alpha_or_0)*
                "alpha_pd_0"*string(alpha_pd_0)*".svg")

            end
        end
    end

end

array_alpha_pd_0 = [0.07, 0.05, 0.04]
array_alpha_or_0 = [20.0, 25.0, 31.2] * (1/4)
array_delay_pd = [12, 18, 24]

plot_all_angles(array_alpha_pd_0[1:1],
array_alpha_or_0[1:1],
array_delay_pd[1:1],
v_struct,
angle_array,
elec,
dis_sp,
w,
t_w,
e_w,
sub_base)

# plot_all_angles(array_alpha_pd_0,
# array_alpha_or_0,
# array_delay_pd,
# v_struct,
# angle_array,
# elec,
# dis_sp,
# w,
# t_w,
# e_w,
# sub_base)

# for alpha_pd_0 in [0.07, 0.05, 0.04]
# for alpha_or_0 in [20.0, 25.0, 31.2] * (1/4)
# for delay_pd in [12, 18, 24]

#     interact_v_pd = circshift(v_struct["photodember"], (delay_pd,0));
#     psi_sub_array = zeros(length(t_w),length(e_w),length(angle_array));
#     psi_incoherent_array = zeros(length(t_w),length(e_w),length(angle_array));

#     plot_agg = [];
#     gr(size=(850,300));

#     for (ind, angle) in enumerate(angle_array)
        
        
#         interact_v = alpha_or_0 .* circshift(v_struct["angle_"*string(angle)], (delay_or,0)) .+
#         alpha_pd_0 .* interact_v_pd;
        
#         @time loss_spectrum!(psi_sub_array,
#         psi_incoherent_array,ind,interact_v, dis_sp, elec, w, t_w, e_w)

#         @views p = heatmap(e_w, t_w .- 0.2, psi_incoherent_array[:,:,ind], c =:jet, 
#         aspect_ratio = 9/2.5,
#         xlims=[-4.5,4.5],
#         ylims=[-1,1.5],
#         colorbar=false,
#         xaxis=nothing,
#         yaxis=nothing,
#         right_margin = -0Plots.mm)
#         yflip!(true)
        
#         push!(plot_agg,p)
        
#     end

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

#     plot(plot_agg...,layout=(2,6))
#     savefig("./analysis_angles/saved-plots/"*sub_base*"delay_pd_"*string(delay_pd)*
#     "alpha_or_0"*string(alpha_or_0)*
#     "alpha_pd_0"*string(alpha_pd_0)*".svg")

# end
# end
# end


