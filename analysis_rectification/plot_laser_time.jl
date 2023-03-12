using LaTeXStrings
include("./init.jl");
include("../src/utils.jl");
#original values

mat = IndiumArsenide()

laser_pulse_time_fwhm = 40e-6;

base="analysis_rectification/saved-matrices/analysis_"
sub_base = "laser_pulse_time_fwhm_single"
base = base*sub_base

# array_laser_pulse_time_fwhm = laser_pulse_time_fwhm_o.*10.0.^(range(-3,stop=-2.5,length=3))

array_laser_pulse_time_fwhm  = [250.0, 650.0, 50.0]
alpha_pd_factor_array = [0.07, 0.07, 0.07]
alpha_or_factor_array = [31.2, 31.2, 11.2]

shift_pd = 18
shift_or = 0

v_struct = load(base*".jld2");

w, e_w, t_w = mod_electron.energy_time_grid(elec,numericalp.subsampling_factor,dis_sp.energy, dis_sp.deltat);

psi_sub_array = zeros(length(t_w),length(e_w),length(array_laser_pulse_time_fwhm));
psi_incoherent_array =similar(psi_sub_array)
psi_sub_array_rect =similar(psi_sub_array)
psi_incoherent_array_rect =similar(psi_sub_array)
psi_sub_array_pd =similar(psi_sub_array)
psi_incoherent_array_pd =similar(psi_sub_array)

for ind in 1:3

    v_pd = alpha_pd_factor_array[ind] * circshift(v_struct["photodember_"*string(ind)],(shift_pd,0))
    v_or = alpha_or_factor_array[ind] * circshift(v_struct["rectification_"*string(ind)],(shift_or,0))
    v_comb = v_pd .+ v_or

    @time loss_spectrum!(psi_sub_array,
    psi_incoherent_array,ind,v_comb, dis_sp, elec, w, t_w, e_w)

    @time loss_spectrum!(psi_sub_array_rect,
    psi_incoherent_array_rect,ind,v_or, dis_sp, elec, w, t_w, e_w)

    @time loss_spectrum!(psi_sub_array_pd,
    psi_incoherent_array_pd,ind,v_pd, dis_sp, elec, w, t_w, e_w)


end


# function my_heatmap(e_w::Array{Float64,1}, t_w::Array{Float64,1}, psi::Array{Float64,2})

#     p = heatmap(e_w, t_w, psi, c =:jet, 
#     aspect_ratio = 9/2.5,
#     xlims=[-4.5,4.5],
#     ylims=[-1,1.5],
#     colorbar=false,
#     xaxis=nothing,
#     yaxis=nothing,
#     right_margin = -0Plots.mm)
#     yflip!(true)

#     return p
# end


plot_agg_incoherent = [];
plot_agg_sub = [];
gr(size=(1200,650));


for ind in 1:3
    
    xl = ind==3;

    p_comb = my_heatmap_l(e_w, t_w .- 0.2, psi_incoherent_array[:,:,ind]; yl = true, xl=xl, rmargin=-5, lmargin = 10)
    p_or = my_heatmap(e_w, t_w .- 0.2, psi_incoherent_array_rect[:,:,ind]; xl=xl, rmargin=-5)
    p_pd = my_heatmap(e_w, t_w .- 0.2, psi_incoherent_array_pd[:,:,ind]; xl=xl, rmargin=10)

    push!(plot_agg_incoherent, p_comb, p_or, p_pd)

    p_comb = my_heatmap(e_w, t_w .- 0.2, psi_sub_array[:,:,ind]; xl=xl, rmargin=-5)
    p_or = my_heatmap(e_w, t_w .- 0.2, psi_sub_array_rect[:,:,ind]; xl=xl, rmargin=-5)
    p_pd = my_heatmap_o(e_w, t_w .- 0.2, psi_sub_array_pd[:,:,ind]; xl=xl)

    push!(plot_agg_incoherent, p_comb, p_or, p_pd)
end

titlep = plot(title= " Laser Time Dependence", grid = false, showaxis = false, bottom_margin = -50Plots.px)
l = @layout [A{0.01h}; (length(array_laser_pulse_time_fwhm),6)]
plot(titlep,plot_agg_incoherent...,plot_agg_sub...,layout=l)
savefig("./analysis_rectification/saved-plots/"*sub_base*".svg")


