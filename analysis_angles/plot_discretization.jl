using Plots
using JLD2, FileIO, CodecZlib
using .cdem_julia: mod_laser, mod_physicalconstants, mod_utils, mod_customtypes, mod_discrete,
mod_cdem, mod_eels, mod_electron, mod_material
using ..mod_customtypes:Discretization
using ..mod_electron:Electron

include("../src/utils.jl");


base="./analysis_angles/saved-matrices/cleo"

v_struct = load(base*".jld2");

# shift_pd = 18
# shift_or = -15
# alpha_pd_factor = 0.07
# alpha_or_factor = 3

# array_shift_pd = [10, 18, 23]
# array_shift_or = [15, 0 , -5 ]

shift_pd = 18
shift_or = -15
alpha_pd_factor = 0.07
alpha_or_factor = 3

array_shift_pd = [14, 18, 22]
# array_shift_or = [-15 ]

ang = 90
# interact_v_pd = circshift(v_struct["photodember"], (shift_pd,0));
# interact_v = alpha_or_0 .* circshift(v_struct["angle_"*string(ang)], (0,0))
# numericalp = mod_customtypes.NumericalParameters(;tc_subsampling = 5,subsampling_factor = 5)
w, e_w, t_w = mod_electron.energy_time_grid(elec,numericalp.subsampling_factor,dis_sp.energy, dis_sp.deltat);


psi_sub_array = zeros(length(t_w),length(e_w),length(array_shift_pd));
psi_incoherent_array =similar(psi_sub_array)
psi_sub_array_rect =similar(psi_sub_array)
psi_incoherent_array_rect =similar(psi_sub_array)
psi_sub_array_pd =similar(psi_sub_array)
psi_incoherent_array_pd =similar(psi_sub_array)

for (ind,shift_or) in enumerate(array_shift_or)
    
        
    v_pd = alpha_pd_factor * circshift(v_struct["photodember"],(shift_pd,0))
    v_or = alpha_or_factor * circshift(v_struct["angle_"*string(ang)],(shift_or,0))
    v_comb = v_pd .+ v_or
    
    @time loss_spectrum!(psi_sub_array,
    psi_incoherent_array,ind,v_comb, dis_sp, numericalp, elec, w, t_w, e_w)
    
    @time loss_spectrum!(psi_sub_array_rect,
    psi_incoherent_array_rect,ind,v_or, dis_sp,numericalp, elec ,w, t_w, e_w)

    @time loss_spectrum!(psi_sub_array_pd,
    psi_incoherent_array_pd,ind,v_pd, dis_sp, numericalp, elec, w, t_w, e_w)
    

end


# psi_sub_array ,
# psi_incoherent_array ,
# psi_sub_array_rect ,
# psi_incoherent_array_rect ,
# psi_sub_array_pd ,
# psi_incoherent_array_pd = aggregator(array_shift_pd, shift_or, alpha_pd_factor, alpha_or_factor, 
# dis_sp,
# elec, w, t_w, e_w)


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

base_sub = "shift_or "
titlep = plot(title= base_sub*"Dependence", grid = false, showaxis = false, bottom_margin = -50Plots.px)
l = @layout [A{0.01h}; (length(array_shift_pd),6)]
pp= plot(titlep,plot_agg_incoherent...,layout=l)

filename = "./analysis_angles/saved-plots/cleo"*string(ang)*".svg";
savefig(filename)
convert_svg_to_png(filename, 300)
    
    







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