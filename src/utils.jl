using LaTeXStrings
using DataFrames
using CSV

function loss_spectrum!(psi_sub_array::Array{Float64,3},
    psi_incoherent_array::Array{Float64,3},
    ind::Int64,
    interact_v::Array{Float64,2}, 
    discretization::Discretization,
    elec::Electron, w::Array{Float64,2}, t_w::Array{Float64,1}, e_w::Array{Float64,1})
        
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
    elec::Electron, w::Array{Float64,2}, t_w::Array{Float64,1}, e_w::Array{Float64,1})
        
    f_t = mod_eels.calculate_ft(discretization, interact_v , elec);
    
    psi = mod_eels.calculate_psi_coherent(dis_sp, elec, f_t);
    psi_sub = mod_eels.psi_subsampled(numericalp.subsampling_factor,psi, e_w);
    psi_sub_array .= psi_sub
    psi_incoherent_array .= mod_eels.incoherent_convolution_fast(psi_sub, w, t_w, e_w);
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

function save_to_database(directory::String, database_name::String, simulation_name::String, 
    dis_sp::Discretization, las::Laser, numericalp::NumericalParameters)
    

    if !isfile(directory*database_name*".csv")
        data = DataFrame()
        CSV.write(directory*database_name*".csv", data)
    end

    
    df = CSV.read(directory*database_name*".csv",DataFrame, header=true)
    # print(df)
    variable_names =  ["dis_x0",
    "dis_y0",
    "dis_d_xprime",
    "dis_d_yprime",
    "dis_d_zprime",
    "dis_xprime_max",
    "dis_yprime_max",
    "dis_zprime_max",
    "dis_d_z",
    "dis_zmax",
    "dis_z_max",
    "dis_t0",
    "dis_ddt",
    "dis_delay_max",
    "dis_fs",
    "dis_l",
    "las_pulse_energy_experiment",
    "las_pulse_energy_gain_factor",
    "las_laser_spot_fwhm",
    "las_laser_pulse_time_fwhm",
    "las_pulse_type",
    "numericalp_tc_subsampling",
    "numericalp_subsampling_factor"
    ]

    values = [dis_sp.x0,
    dis_sp.y0,
    dis_sp.xprime[2]-dis_sp.xprime[1],
    dis_sp.yprime[2]-dis_sp.yprime[1],
    dis_sp.zprime[2]-dis_sp.zprime[1],
    maximum(dis_sp.xprime),
    maximum(dis_sp.yprime),
    maximum(dis_sp.zprime),
    dis_sp.z[2]-dis_sp.z[1],
    maximum(dis_sp.z),
    dis_sp.z_max,
    dis_sp.t0,
    dis_sp.deltat[2]-dis_sp.deltat[1],
    maximum(dis_sp.deltat),
    1.0/(dis_sp.t[2]-dis_sp.t[1]),
    length(dis_sp.t)-1,
    las.pulse_energy_experiment,
    las.pulse_energy_gain_factor ,
    las.laser_spot_fwhm ,
    las.laser_pulse_time_fwhm ,
    las.pulse_type,
    numericalp.tc_subsampling,
    numericalp.subsampling_factor

    ]

    df[!, "Variable"] = variable_names
    df[!, simulation_name] = values

    # df = DataFrame(Variable = variable_names', Value = values');
    # pushfirst!(df, ("simulation",simulation_name), promote = true);
    CSV.write(directory*database_name*".csv", df) 
    # print(df)
end

