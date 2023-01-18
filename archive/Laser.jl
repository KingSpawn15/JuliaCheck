module LaserMod

    using ..Utils
    export Laser, set_laser!

    mutable struct Laser

        pulse_energy_experiment::Float64
        pulse_energy_gain_factor::Float64
        laser_spot_fwhm::Float64
        theta_pol::Float64
        laser_pulse_time_fwhm::Float64
        pulse_energy::Float64
        laser_spot_sigma::Float64
        laser_pulse_time_sigma::Float64

        function Laser()
            return new()
        end

    end

    function set_laser!(;laser::Union{Laser,Nothing}=nothing,
        pulse_energy_experiment::Union{Float64,Nothing}=nothing,
        pulse_energy_gain_factor::Union{Float64,Nothing}=nothing,
        laser_spot_fwhm::Union{Float64,Nothing}=nothing,
        theta_pol::Union{Float64,Nothing}=nothing,
        laser_pulse_time_fwhm::Union{Float64,Nothing}=nothing)
        
        if isnothing(laser)
            laser = Laser()
        end

        if !isnothing(pulse_energy_experiment)
            laser.pulse_energy_experiment = pulse_energy_experiment
        end

        if !isnothing(pulse_energy_gain_factor)
            laser.pulse_energy_gain_factor = pulse_energy_gain_factor
        end
        
        if !isnothing(theta_pol)
            laser.theta_pol = theta_pol
        end

        if !isnothing(laser.pulse_energy_gain_factor) && !isnothing(pulse_energy_gain_factor)
            laser.pulse_energy = laser.pulse_energy_gain_factor * laser.pulse_energy_experiment
        end
        
        if !isnothing(laser_spot_fwhm)

            laser.laser_spot_fwhm = laser_spot_fwhm
            laser.laser_spot_sigma = Utils.calculate_sigma(laser.laser_spot_fwhm)
        
        end

        if !isnothing(laser_pulse_time_fwhm)
            laser.laser_pulse_time_fwhm = laser_pulse_time_fwhm
            laser.laser_pulse_time_sigma = Utils.calculate_sigma(laser.laser_pulse_time_fwhm)#[s]
        end

        return laser
    end

end


