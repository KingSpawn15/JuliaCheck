mutable struct laser_params

    pulse_energy_experiment::Float64
    pulse_energy_gain_factor::Float64
    laser_spot_fwhm::Float64
    theta_pol::Float64
    laser_pulse_time_fwhm::Float64
    pulse_energy::Float64
    laser_spot_sigma::Float64
    laser_pulse_time_sigma::Float64

    function laser_params()
        return new()
    end

end

mutable struct discretization_params

    x0::Float64
    y0::Float64
    ddt::Float64
    delay_max::Float64
    fs::Float64
    l::Float64
    t0::Float64
    xprime_max::Float64
    d_xprime::Float64
    yprime_max::Float64
    d_yprime ::Float64
    zprime_max ::Float64
    d_zprime::Float64
    ddz::Float64
    zmax ::Float64
    z_max ::Float64

    function discretization_params()
        return new()
    end

end

mutable struct utem_params
    electron_total_energy::Float64
    electron_total_time_fs::Float64
    electron_time_coherent_fwhm_fs::Float64
    electron_theta::Float64
    electron_velocity_c::Float64

    function utem_params()
        return new()
    end

end

mutable struct numerical_params
    tc_subsampling::Int
    subsampling_factor::Int

    function numerical_params()
        return new()
    end

end

mutable struct default_params
    laser_parameters::laser_params
    discretization_parameters::discretization_params
    utem_parameters::utem_params
    numerical_parameters::numerical_params

    function default_params()
        return new()
    end
end
