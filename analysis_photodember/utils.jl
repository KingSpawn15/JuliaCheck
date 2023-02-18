
    
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
