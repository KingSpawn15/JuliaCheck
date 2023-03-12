include("./init.jl");
include("../src/utils.jl");
#original values
mat = IndiumArsenide()


base="analysis_rectification/saved-matrices/zmax_analysis_laser_spot_40_laser_time_50"

mod_laser.set_laser!(;laser=las,
laser_pulse_time_fwhm  = 50 * 1e-15,
laser_spot_fwhm = 40e-6,
theta_pol = 90.0*pi/180);

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


set_electron!(;electron=elec,
    electron_total_energy = 0.94,
    # electron_total_energy = 1.1,
    electron_total_time_fs = 360.0,
    electron_time_coherent_fwhm_fs = 50.0,
    electron_theta = -7*pi/180,
    electron_velocity_c = 0.7)

alpha_pd_factor_array = [0.07, 0.07, 0.07]
alpha_pd_factor_array.*=1
alpha_or_factor_array = [1.0, 1.0, 1.0]*1.2

shift_pd = 18
shift_or = -0


# y0_array = [-1e-7, -1e-6, -1e-5]
# x0_array = [-1e-7, 0.0, 1e-7]
# z_max_array = [20e-6, 30e-6, 50e-6]

# array_discretization_x0 = []
# array_discretization_y0 = []
# array_discretization_z_max = []

# # for y0 in y0_array
# #     push!(array_discretization, discretization_mod(;y0 = y0))
# # end

# for ind in 1 : 3
#     push!(array_discretization_z_max, discretization_mod(;z_max = z_max_array[ind]))
#     push!(array_discretization_x0, discretization_mod(;x0 = x0_array[ind]))
#     push!(array_discretization_y0, discretization_mod(;y0 = y0_array[ind]))
# end

# # dict = Dict("y0"=>array_discretization)
# dict = Dict("z_max"=>array_discretization_z_max,"x0"=>array_discretization_x0,
# "y0"=>array_discretization_y0)

z_max_array = [20e-6, 30e-6, 50e-6]


# dict = Dict("y0" => y0_array,
# )

dict = Dict("z_max" => z_max_array
)

function get_discretization(array_discretization::Vector{Any}, ss::Int64)::Discretization
    return array_discretization[ss]
end



function my_heatmap_o(e_w::Array{Float64,1}, t_w::Array{Float64,1}, psi::Array{Float64,2}; yl::Bool=false, xl::Bool=false)

    p = heatmap(e_w, t_w, psi, c =:jet, 
    aspect_ratio = 9/2.5,
    xlims=[-4.5,4.5],
    ylims=[-1,1.5],
    colorbar=false,
    xaxis=xl,
    yaxis=yl,
    yguidefontsize=12,
    xguidefontsize=8)
    yflip!(true)
    if (yl)
        ylabel!(L"\Delta t (ps)")
    end
    if (xl)
        xlabel!("Energy Loss (ev)")
    end

    return p
end

function my_heatmap_l(e_w::Array{Float64,1}, t_w::Array{Float64,1}, psi::Array{Float64,2}; yl::Bool=false, xl::Bool=false, 
    rmargin::Int64=-25, lmargin::Int64=5)

    p = heatmap(e_w, t_w, psi, c =:jet, 
    aspect_ratio = 9/2.5,
    xlims=[-4.5,4.5],
    ylims=[-1,1.5],
    colorbar=false,
    xaxis=xl,
    yaxis=yl,
    right_margin = eval(rmargin)Plots.mm,
    left_margin = eval(lmargin)Plots.mm,
    yguidefontsize=12,
    xguidefontsize=8)
    yflip!(true)
    if (yl)
        ylabel!(L"\Delta t (ps)")
    end
    if (xl)
        xlabel!("Energy Loss (ev)")
    end

    return p
end

function my_heatmap(e_w::Array{Float64,1}, t_w::Array{Float64,1}, psi::Array{Float64,2}; yl::Bool=false, xl::Bool=false, rmargin::Int64=-25)

    p = heatmap(e_w, t_w, psi, c =:jet, 
    aspect_ratio = 9/2.5,
    xlims=[-4.5,4.5],
    ylims=[-1,1.5],
    colorbar=false,
    xaxis=xl,
    yaxis=yl,
    right_margin = eval(rmargin)Plots.mm,
    yguidefontsize=12,
    xguidefontsize=8)
    yflip!(true)
    if (yl)
        ylabel!(L"\Delta t (ps)")
    end
    if (xl)
        xlabel!("Energy Loss (ev)")
    end

    return p
end


for field in keys(dict)

    base_sub = field
    base_i = base*"_"*base_sub
    v_struct = load(base_i*".jld2");

    
    w, e_w, t_w = mod_electron.energy_time_grid(elec,numericalp.subsampling_factor,dis_sp.energy, dis_sp.deltat);
    psi_sub_array = zeros(length(t_w),length(e_w),length(dict[field]));
    psi_incoherent_array =similar(psi_sub_array)
    psi_sub_array_rect =similar(psi_sub_array)
    psi_incoherent_array_rect =similar(psi_sub_array)
    psi_sub_array_pd =similar(psi_sub_array)
    psi_incoherent_array_pd =similar(psi_sub_array)
    
    for (ind,val) in enumerate(dict[field])
        
        ff = Symbol(field)
        dis = eval(:(discretization_mod(;$ff = $val)))

        v_pd = alpha_pd_factor_array[ind] * circshift(v_struct["photodember_"*string(ind)],(shift_pd,0))
        v_or = alpha_or_factor_array[ind] * circshift(v_struct["rectification_"*string(ind)],(shift_or,0))
        v_comb = v_pd .+ v_or

        @time loss_spectrum!(psi_sub_array,
        psi_incoherent_array,ind,v_comb, dis, elec, w, t_w, e_w)

        @time loss_spectrum!(psi_sub_array_rect,
        psi_incoherent_array_rect,ind,v_or, dis, elec ,w, t_w, e_w)

        @time loss_spectrum!(psi_sub_array_pd,
        psi_incoherent_array_pd,ind,v_pd, dis, elec, w, t_w, e_w)


    end


    plot_agg_incoherent = [];
    plot_agg_sub = [];
    gr(size=(1200,650));


    for (ind, _) in enumerate(dict[field])
        
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

    titlep = plot(title= base_sub*"Dependence", grid = false, showaxis = false, bottom_margin = -50Plots.px)
    l = @layout [A{0.01h}; (length(dict[field]),6)]
    plot(titlep,plot_agg_incoherent...,layout=l)
    savefig("./analysis_rectification/saved-plots/zmax_analysis_laser_spot_40_laser_time_50"*base_sub*".svg")



end






