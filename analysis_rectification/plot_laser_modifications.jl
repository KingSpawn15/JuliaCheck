include("./init.jl");
using LaTeXStrings
#original values
mat = IndiumArsenide()

base="analysis_photodember/saved-matrices/laseranalysis_withoutcutoff"
v_struct = load(base*".jld2");

# array_pulse_energy_experiment = [0.5 * pulse_energy_experiment_o, pulse_energy_experiment_o, 2 * pulse_energy_experiment_o]
# array_pulse_energy_experiment = [1e-2 * pulse_energy_experiment_o, 1e-1 * pulse_energy_experiment_o, pulse_energy_experiment_o]

# array_pulse_energy_experiment = pulse_energy_experiment_o.*10.0.^(range(-3,stop=-2.5,length=3))

# array_pulse_energy_experiment = [1e-5 * pulse_energy_experiment_o, 
# 1e-4 * pulse_energy_experiment_o, 
# 1e-3 *pulse_energy_experiment_o]

# alpha_pd_0 =  1;
alpha_pd_0 =  .07;

w, e_w, t_w = mod_electron.energy_time_grid(elec,numericalp.subsampling_factor,dis_sp.energy, dis_sp.deltat);

psi_incoherent_array = zeros(length(t_w),length(e_w),
length(array_pulse_energy_experiment));

psi_sub_array = similar(psi_incoherent_array);

function loss_spectrum_single(interact_v::Array{Float64,2}, 
    discretization::Discretization,
    elec::Electron)
        
    f_t = mod_eels.calculate_ft(discretization, interact_v , elec);
    
    psi = mod_eels.calculate_psi_coherent(dis_sp, elec, f_t);
    psi_sub = mod_eels.psi_subsampled(numericalp.subsampling_factor,psi, e_w);
    
    psi_incoherent = mod_eels.incoherent_convolution_fast(psi_sub, w, t_w, e_w);

    return psi_sub, psi_incoherent
end

global ig = 1;

for (ind, pulse_energy_experiment) in enumerate(array_pulse_energy_experiment) 
    

    igl = ig
    v_pd = v_struct["photodember_pulse_energy_experiment"*string(ind)]
    v_pd .*= alpha_pd_0
    @time psi_sub, psi_incoherent = loss_spectrum_single(v_pd, dis_sp, elec)
    psi_sub_array[:,:,igl] .= psi_sub
    psi_incoherent_array[:,:,igl] .= psi_incoherent
    global ig += 1
end


"""
"""

plot_agg = [];
gr(;size=(600,350))

for ind in eachindex(psi_incoherent_array[1,1,:])
    print(ind)

    p_sub = heatmap(e_w, t_w, psi_sub_array[:,:,ind], c =:jet, 
    aspect_ratio = 9/2.5,
    xlims=[-4.5,4.5],
    ylims=[-1,1.5],
    colorbar=false,
    xaxis=(ind==size(psi_incoherent_array,3) ? true : false),
    # yaxis=false,
    ylabel=L"\Delta t (ps)",
    xlabel="Energy Loss (eV)",
    right_margin = 0Plots.mm,
    bottom_margin = -10Plots.mm);
    yflip!(true);

    p_inc = heatmap(e_w, t_w, psi_incoherent_array[:,:,ind], c =:jet, 
    aspect_ratio = 9/2.5,
    xlims=[-4.5,4.5],
    ylims=[-1,1.5],
    colorbar=false,
    xaxis=(ind==size(psi_incoherent_array,3) ? true : false),
    xlabel="Energy Loss (eV)",
    # yaxis=false,
    right_margin = 0Plots.mm,
    bottom_margin = -10Plots.mm);
    yflip!(true);

    push!(plot_agg,p_sub);
    push!(plot_agg,p_inc);
end

titlep = plot(title="Best Value \n (CLEO)", grid = false, showaxis = false, bottom_margin = -50Plots.px)
l = @layout [A{0.01h}; (size(psi_incoherent_array,3),2)]
plot(titlep, plot_agg[:]..., layout = l)

savefig("./analysis_photodember/saved-plots/"*"CLEO_withoutcutoff.svg")