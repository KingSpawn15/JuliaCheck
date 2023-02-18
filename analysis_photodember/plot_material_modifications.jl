include("./init.jl");
include("./utils.jl");
using LaTeXStrings
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

base="analysis_photodember/saved-matrices/materialanalysis_withcutoff"

array_alpha =[ 0.8 * orig_alpha, 1 * orig_alpha, 1.2 * orig_alpha]
array_kappa =[ 0.5 * orig_kappa , 1 * orig_kappa, 1.02 * orig_kappa]
array_gamma =[ 0.5 * orig_gamma, 1 * orig_gamma, 1.01 * orig_gamma]
array_me = [0.8 * orig_me, 1 * orig_me, 1.2 * orig_me]
array_mh = [0.8 * orig_mh, 1 * orig_mh, 1.02 * orig_mh]
array_hnew = [0.8 * orig_hnew, 1 * orig_hnew, 1.2 * orig_hnew]
array_alpha_gamma = [0.8 * orig_alpha_gamma, 1 * orig_alpha_gamma, 1.2 * orig_alpha_gamma]
array_n_eq = [0.98 * orig_n_eq, 1 * orig_n_eq, 1.2 * orig_n_eq]


mod_laser.set_laser!(;laser=las,
theta_pol = 90.0*pi/180);



alpha_pd_factor_array = [0.07, 0.07, 0.07]
alpha_pd_factor_array.*=1
alpha_or_factor_array = [31.2, 31.2, 31.2]

shift_pd = 18
shift_or = -15


############

dict = Dict("alpha"=>array_alpha ,
"kappa"=>array_kappa ,
"gamma"=>array_gamma ,
"me"=>array_me ,
"mh"=>array_mh ,
"hnew"=>array_hnew ,
"alpha_gamma"=>array_alpha_gamma ,
"n_eq"=>array_n_eq)

# dict = Dict("kappa"=>array_kappa ,
# "gamma"=>array_gamma ,
# "mh"=>array_mh ,
# "n_eq"=>array_n_eq)


# dict = Dict("alpha"=>array_alpha)

mod_laser.set_laser!(;laser=las,
theta_pol = 90.0*pi/180);


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

    for (ind,_) in enumerate(dict[field])
        
        v_pd = alpha_pd_factor_array[ind] * circshift(v_struct["photodember_"*string(ind)],(shift_pd,0))
        v_or = alpha_or_factor_array[ind] * circshift(v_struct["rectification_"*string(ind)],(shift_or,0))
        v_comb = v_pd .+ v_or

        @time loss_spectrum!(psi_sub_array,
        psi_incoherent_array,ind,v_comb, dis_sp, elec, w, t_w, e_w)

        @time loss_spectrum!(psi_sub_array_rect,
        psi_incoherent_array_rect,ind,v_or, dis_sp, elec ,w, t_w, e_w)

        @time loss_spectrum!(psi_sub_array_pd,
        psi_incoherent_array_pd,ind,v_pd, dis_sp, elec, w, t_w, e_w)


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
    savefig("./analysis_photodember/saved-plots/withcutoff_"*base_sub*".svg")



end






