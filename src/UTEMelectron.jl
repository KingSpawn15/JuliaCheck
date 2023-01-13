mutable struct Electron

    electron_total_energy::Float64
    electron_total_time::Float64
    electron_time_coherent_fwhm::Float64
    electron_theta::Float64
    electron_velocity::Float64
    electron_time_incoherent_sigma::Float64
    electron_energy_incoherent_sigma::Float64
    electron_time_coherent_sigma::Float64

    function Electron()
        return new()
    end

end

function ElectronSetup(utem_parameters::utem_params)
        
    utem = Electron() 
    utem.electron_total_energy = utem_parameters.electron_total_energy;
    utem.electron_total_time = utem_parameters.electron_total_time_fs * 1e-15;
    utem.electron_time_coherent_fwhm = utem_parameters.electron_time_coherent_fwhm_fs * 1e-15;
    utem.electron_theta = utem_parameters.electron_theta;
    utem.electron_velocity = c2msec(utem_parameters.electron_velocity_c);
    utem.electron_time_coherent_sigma = parameters_coherent(utem);
    utem.electron_time_incoherent_sigma,utem.electron_energy_incoherent_sigma = parameters_incoherent(utem);
    return utem

end

function parameters_coherent(utem_parameters::Electron)
            
    electron_time_coherent_sigma = utem_parameters.electron_time_coherent_fwhm./(2*sqrt(2*log(2)));#[s]
    return electron_time_coherent_sigma
    
end


function parameters_incoherent(utem_parameters::Electron)
    
    #Incoherent parameters

    electron_time_incoherent_sigma = utem_parameters.electron_total_time/(2*sqrt(2*log(2)));#[s]
    electron_energy_incoherent_sigma = utem_parameters.electron_total_energy/(2*sqrt(2*log(2)));#[eV]
    
    return electron_time_incoherent_sigma, electron_energy_incoherent_sigma

end

function energy_time_grid(utem_parameters::Electron, sub_sample_factor::Int, energy::Array{Float64}, deltat::Array{Float64})
    
    e_w = energy[1:sub_sample_factor:end];
    t_w = deltat*1e12;#[ps]
    
    sigma_t = utem_parameters.electron_time_incoherent_sigma*1e12;
    sigma_e= utem_parameters.electron_energy_incoherent_sigma;
    
    a = cos(utem_parameters.electron_theta)^2/(2*sigma_t^2) + sin(utem_parameters.electron_theta)^2/(2*sigma_e^2);
    b = (sin(2*utem_parameters.electron_theta)/4)*((1/sigma_e^2)-(1/sigma_t^2));
    c = sin(utem_parameters.electron_theta)^2/(2*sigma_t^2) + cos(utem_parameters.electron_theta)^2/(2*sigma_e^2);
    
    TW,EW = ndgrid_array(vec(t_w),vec(e_w));
    w = exp.(-(a*TW.^2 + 2*b*TW.*EW + c*EW.^2));
    # w = w';

    return w, e_w, t_w

end

function c2msec(velocity_c::Float64)

    consts = constants_fundamental()
    C = consts.C;
    velocity_msec = velocity_c * C;
    
    return velocity_msec

end