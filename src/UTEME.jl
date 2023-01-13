module ElectronMod
    
    using ..Utils
    using ..PhysicalContants
    export Electron, set_electron!, energy_time_grid

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


    function set_electron!(;electron::Union{Electron,Nothing}=nothing,
        electron_total_energy::Union{Float64,Nothing}=nothing,
        electron_total_time_fs::Union{Float64,Nothing}=nothing,
        electron_time_coherent_fwhm_fs::Union{Float64,Nothing}=nothing,
        electron_theta::Union{Float64,Nothing}=nothing,
        electron_velocity_c::Union{Float64,Nothing}=nothing)
        
        if isnothing(electron)
            electron = Electron()
        end
        
        if !isnothing(electron_total_energy)
            electron.electron_total_energy = electron_total_energy;
            electron.electron_energy_incoherent_sigma = Utils.calculate_sigma(electron.electron_total_energy)
        end

        if !isnothing(electron_total_time_fs)
            electron.electron_total_time = electron_total_time_fs * 1e-15;
            electron.electron_time_incoherent_sigma = Utils.calculate_sigma(electron.electron_total_time)
        end

        if !isnothing(electron_time_coherent_fwhm_fs)
            electron.electron_time_coherent_fwhm = electron_time_coherent_fwhm_fs * 1e-15;
            electron.electron_time_coherent_sigma = Utils.calculate_sigma(electron.electron_time_coherent_fwhm);
        end

        if !isnothing(electron_theta)
            electron.electron_theta = electron_theta;
        end

        if !isnothing(electron_velocity_c)
            electron.electron_velocity = c2msec(electron_velocity_c);
        end

        return electron

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

        # consts = constants_fundamental()
        # C = consts.C;
        velocity_msec = velocity_c * PhysicalContants.C;
        
        return velocity_msec

    end

end


