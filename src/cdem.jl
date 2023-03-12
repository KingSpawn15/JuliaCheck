module cdem_julia
    export mod_laser, mod_physicalconstants, mod_utils, mod_customtypes, mod_discrete,
    mod_cdem, mod_eels, mod_electron, mod_material


    module mod_utils

        export calculate_sigma, index_calculator_2d

        function calculate_sigma(fwhm::Float64)
            return fwhm./(2*sqrt(2*log(2)))
        end

        function index_calculator_2d(k::Int, nrows::Int, ncols::Int)
        
            col = Int(floor(k / nrows)) + 1
            row = mod(k,nrows)
        
            if row == 0
                row = nrows
                col = col - 1
            end
        
            return row, col
        
        end

        function index_calculator_3d(k::Int, d1::Int, d2::Int, d3::Int)

            ind3 = Int(floor(k / d1 / d2)) + 1
            ind12 = mod(k , d1 * d2)
        
            if ind12 == 0
                ind3 = ind3 - 1
                ind12 = d1 * d2
            end 
        
            ind1, ind2 = index_calculator_2d(ind12, d1, d2)
        
            return ind1, ind2, ind3
        
        end
        
        function index_calculator_4d(k::Int, d1::Int, d2::Int, d3::Int, d4::Int)
        
            ind4 = Int(floor(k / d1 / d2 / d3)) + 1
            ind123 = mod(k , d1 * d2 * d3)
        
            if ind123 == 0
                ind4 = ind4 - 1
                ind123 = d1 * d2 * d3
            end 
        
            ind1, ind2 , ind3 = index_calculator_3d(ind123, d1, d2, d3)
        
            return ind1, ind2, ind3, ind4
            
        end

        
    end

    module mod_physicalconstants

        export   HBAR,EPSILON_0, C, M_E, Q_E, MU_0, ETA_0

        global const HBAR       = 6.62607004e-34/(2*pi)::Float64  # s J
        global const EPSILON_0  = 8.8541878128e-12::Float64  # V/m
        global const C = 299792458.0::Float64 
        global const M_E = 9.109e-31::Float64 
        global const Q_E = 1.60217662e-19::Float64  # coulomb::Float64 
        global const MU_0 = 1.25663706212e-6::Float64 #[H/m]::Float64 
        global const ETA_0 = 377.0::Float64 #[Ohm]

    end

    module mod_customtypes

        export Discretization, NumericalParameters, FactorsRect, RectificationParameters

        mutable struct Discretization
            x0::Float64
            y0::Float64    
            xprime::Array{Float64,1}
            yprime::Array{Float64,1}
            zprime::Array{Float64,1}
            z::Array{Float64,1}
            z_max::Float64
            ddt::Float64
            deltat::Array{Float64,1}
            t0::Float64
            t::Array{Float64,1}
            dt::Float64
            omega::Array{Float64,1}
            energy::Array{Float64,1}
            XPRIME::Array{Float64,4}
            YPRIME::Array{Float64,4}
            ZPRIME::Array{Float64,4}
            Z::Array{Float64,4}
        end

        mutable struct NumericalParameters
            tc_subsampling::Int
            subsampling_factor::Int
        
            function NumericalParameters(;tc_subsampling::Int64, subsampling_factor::Int64)
                return new(tc_subsampling, subsampling_factor)
            end
        
        end

        struct FactorsRect

            factor_rho::Array{Float64,4}
            factor_a::Array{Float64,4}
            factor_rho_Y0::Array{Float64,4}
            factor_rho_Z0::Array{Float64,4}

        end

        struct RectificationParameters
            laser_pulse_time_sigma::Float64
            mz0_ind::Int64
            pz0_ind::Int64
            factors_rect::FactorsRect
            t0::Float64
            y0_ind::Int64
            t_c_subsampled::Array{Float64,1}
            l_tc::Int64
            t_r::Array{Float64,4}
            electron_velocity::Float64
            xprime::Array{Float64,1}
            yprime::Array{Float64,1}
            zprime::Array{Float64,1}
        end

    end

    module mod_laser

        using ..mod_utils
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
                laser.laser_spot_sigma = mod_utils.calculate_sigma(laser.laser_spot_fwhm)
            
            end

            if !isnothing(laser_pulse_time_fwhm)
                laser.laser_pulse_time_fwhm = laser_pulse_time_fwhm
                laser.laser_pulse_time_sigma = mod_utils.calculate_sigma(laser.laser_pulse_time_fwhm)#[s]
            end

            return laser
        end

    end

    module mod_discrete

        using ..mod_customtypes:Discretization
        using ..mod_physicalconstants
        
        export discretization_setup, space_grid

        function discretization_setup(;x0::Float64, y0::Float64, 
            d_xprime::Float64, d_yprime::Float64, d_zprime::Float64, 
            xprime_max::Float64, yprime_max::Float64, zprime_max::Float64, 
            d_z::Float64, zmax::Float64, z_max::Float64,
            t0::Float64, ddt::Float64, delay_max::Float64, fs::Float64, l::Float64)

            xprime = -xprime_max : d_xprime : xprime_max
            yprime = 0.0 : d_yprime : yprime_max
            zprime = -zprime_max : d_zprime : zprime_max


            z = -zmax : d_z : zmax

            XPRIME, YPRIME, ZPRIME, Z = space_grid(xprime,yprime,zprime,z)

            deltat = -delay_max : ddt : delay_max
            t, omega, energy = incidence(fs, l); 
            dt = t[2] - t[1];

            return Discretization(x0, y0, xprime, yprime, zprime, z, z_max,
            ddt, deltat, t0, t, dt, omega, energy,
            XPRIME, YPRIME, ZPRIME, Z)

        end
        
        function space_grid(xprime::Union{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},Array{Float64,1}},
            yprime::Union{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},Array{Float64,1}},
            zprime::Union{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},Array{Float64,1}},
            z::Union{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},Array{Float64,1}})
            
            XPRIME = [i for i in xprime, _ in yprime, _ in zprime, _ in z]
            YPRIME = [i for _ in xprime, i in yprime, _ in zprime, _ in z]
            ZPRIME = [i for _ in xprime, _ in yprime, i in zprime, _ in z]
            Z = [i for _ in xprime, _ in yprime, _ in zprime, i in z]

            return XPRIME, YPRIME, ZPRIME, Z

        end

        function  incidence(fs::Float64, l::Float64)
        
            dt = 1/fs;#[s]
            domega = fs/l*(2*pi);#[Hz]
            t =  (-l/2:l/2).*dt;#[s
            omega = (-l/2:l/2)*domega;#[Hz
            energy = omega.*mod_physicalconstants.HBAR/mod_physicalconstants.Q_E;#[eV]
        
            return t, omega, energy
            
        end

    end

    module mod_material

        using ..mod_physicalconstants, ..mod_laser
        export Material, IndiumArsenide, excited_carriers, calculate_omega_y

        mutable struct Material 
            alpha::Float64
            kappa::Float64
            gamma::Float64
            me::Float64
            mh::Float64
            lambda::Float64
            hnew::Float64
            eg::Float64
            alpha_gamma::Float64
            n_eq::Float64
            m_eq::Float64
            epsilon_e::Float64
            v_t::Float64
            me0tilda::Float64
            d14::Float64
            gamma_factor::Float64
            phase ::Float64

            function Material()
                return new()
            end

        end

        function calculate_omega_y(mat::Material, n_exc::Float64,
            yprime_grid::Array{Float64,4}, gaussian_laser_spot::Array{Float64,4})
            
            Q_E = mod_physicalconstants.Q_E;
            
            omega_y = sqrt.( (Q_E^2/mat.kappa).*(n_exc.*gaussian_laser_spot.*exp.(-mat.alpha.*yprime_grid).*
            (1/mat.me0tilda+1/mat.mh) .+ mat.n_eq/mat.m_eq) .- (mat.gamma/2)^2 );
            
            return omega_y
        end


        function excited_carriers(laser::Laser, alpha::Float64, hnew::Float64)
                
            
            Q_E = mod_physicalconstants.Q_E;
            excited_volume = (pi/4)*laser.laser_spot_fwhm^2*alpha^(-1);#[m^3]
            n_exc = laser.pulse_energy/Q_E/hnew/excited_volume;#[m^-3]
            return n_exc
        
        end

        function IndiumArsenide()
            inas = Material()
            inas.alpha = 7e6;#[m^-1]
            
            EPSILON_0 = mod_physicalconstants.EPSILON_0;
            M_E = mod_physicalconstants.M_E;
            inas.kappa = 12.3*EPSILON_0;
            inas.gamma = 1.3e13;#[s^-1]
            inas.gamma_factor = 0.27;
            inas.phase = -0.81;
            inas.me = 0.022*M_E;#[kg]
            inas.mh = 0.6*M_E;#[kg]
            inas.lambda = 0.8;#[um]
            inas.hnew = 1.24./inas.lambda;#[eV]
            inas.eg = 0.354;#[eV]
            inas.alpha_gamma = 2.2;#[eV^-1]
            inas.n_eq = 1e17*1e6;#[m^-3]#From resistivity measurement!
            inas.m_eq = mass_equilibrium_carrier(inas); #p-type InAs
            inas.epsilon_e = photoexcited_electron_energy(inas);
            inas.v_t = velocity_t(inas);
            inas.me0tilda = photoelectron_mass(inas);
            inas.d14 = 237.6e-12;#[V/m]
            return inas
        end

        function mass_equilibrium_carrier(mat::Material)
            # %METHOD1 Summary of this method goes here
            # %   Detailed explanation goes here
            return mat.mh
        end

        function photoexcited_electron_energy(mat::Material)
            # %METHOD1 Summary of this method goes here
            # %   Detailed explanation goes here
            epsilon_e = (2*(mat.hnew-mat.eg)* mat.mh)./(mat.me+mat.mh+sqrt((mat.me+mat.mh)^2+4*mat.alpha_gamma.*(mat.hnew-mat.eg).*mat.me.*mat.mh));#[eV]
            return epsilon_e
        end


        function photoelectron_mass(mat::Material)
            # %METHOD1 Summary of this method goes here
            # %   Detailed explanation goes here
            
            # me0tilda = mat.me*3*(1+4*mat.alpha_gamma*mat.epsilon_e*(1+mat.alpha_gamma*mat.epsilon_e))^(3/2)
            # /(3+8*mat.alpha_gamma*mat.epsilon_e*(1+mat.alpha_gamma*mat.epsilon_e));#[kg]

            me0tilda = mat.me*3*(1+4*mat.alpha_gamma*mat.epsilon_e*(1+mat.alpha_gamma*mat.epsilon_e))^(3/2)/
            (3+8*mat.alpha_gamma*mat.epsilon_e*(1+mat.alpha_gamma*mat.epsilon_e));#[kg]
            
            return me0tilda

        end

        

        function velocity_t(mat::Material)
            # %METHOD1 Summary of this method goes here
            # %   Detailed explanation goes here
            
            # consts = constants_fundamental()
            Q_E = mod_physicalconstants.Q_E;
            
            v_t = sqrt((2*mat.epsilon_e.*Q_E.*(1+mat.alpha_gamma.*mat.epsilon_e))
            ./(3*mat.me.*(1+4*mat.alpha_gamma.*mat.epsilon_e.*(1+mat.alpha_gamma.*mat.epsilon_e))));#[m/s]

            return v_t

        end


    end

    module mod_electron
        
        using ..mod_utils
        using ..mod_physicalconstants
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
                electron.electron_energy_incoherent_sigma = mod_utils.calculate_sigma(electron.electron_total_energy)
            end

            if !isnothing(electron_total_time_fs)
                electron.electron_total_time = electron_total_time_fs * 1e-15;
                electron.electron_time_incoherent_sigma = mod_utils.calculate_sigma(electron.electron_total_time)
            end

            if !isnothing(electron_time_coherent_fwhm_fs)
                electron.electron_time_coherent_fwhm = electron_time_coherent_fwhm_fs * 1e-15;
                electron.electron_time_coherent_sigma = mod_utils.calculate_sigma(electron.electron_time_coherent_fwhm);
            end

            if !isnothing(electron_theta)
                electron.electron_theta = electron_theta;
            end

            if !isnothing(electron_velocity_c)
                electron.electron_velocity = c2msec(electron_velocity_c);
            end

            return electron

        end


        function energy_time_grid(utem_parameters::Electron, sub_sample_factor::Int, energy::Array{Float64,1}, deltat::Array{Float64,1})
            
            e_w = energy[1:sub_sample_factor:end];
            t_w = deltat*1e12;#[ps]
            
            sigma_t = utem_parameters.electron_time_incoherent_sigma*1e12;
            sigma_e= utem_parameters.electron_energy_incoherent_sigma;
            
            a = cos(utem_parameters.electron_theta)^2/(2*sigma_t^2) + sin(utem_parameters.electron_theta)^2/(2*sigma_e^2);
            b = (sin(2*utem_parameters.electron_theta)/4)*((1/sigma_e^2)-(1/sigma_t^2));
            c = sin(utem_parameters.electron_theta)^2/(2*sigma_t^2) + cos(utem_parameters.electron_theta)^2/(2*sigma_e^2);
        
            w = [exp(-(a*ti^2 + 2*b*ti.*ei + c*ei.^2)) for ti in t_w, ei in e_w]


            return w, e_w, t_w

        end

        function c2msec(velocity_c::Float64)


            velocity_msec = velocity_c * mod_physicalconstants.C;
            
            return velocity_msec

        end

    end

    module mod_cdem

        using Trapz, Interpolations
        
        using ..mod_physicalconstants, ..mod_customtypes, ..mod_material, ..mod_laser, ..mod_electron
        export interaction_potential_photodember, interaction_potential_rectification

        function interaction_potential_rectification(discretization::Discretization, material::Material,
            laser::Laser , electron::Electron, numerical_parameters::NumericalParameters)

        
            # # InAs parameters
            d14 = material.d14;
            alpha = material.alpha;

            # # Laser parameters

            laser_spot_sigma = laser.laser_spot_sigma;
            laser_pulse_time_fwhm = laser.laser_pulse_time_fwhm;#[s]
            laser_pulse_time_sigma = laser.laser_pulse_time_sigma;#[s]
            pulse_energy = laser.pulse_energy;
            e0_squared = (pulse_energy/laser_pulse_time_fwhm)*2*ETA_0;

            xprime = discretization.xprime;
            yprime = discretization.yprime;
            zprime = discretization.xprime;
            z_max = discretization.z_max;#[m] Sample +/-z boundary

            # # Down sample t to improve run speed

            t_c = discretization.t[discretization.t .> discretization.t0];
            t_c_subsampled = t_c[1:numerical_parameters.tc_subsampling:end];

            y0_ind = findfirst(yprime .>= 0);
            mz0_ind = findfirst(zprime .>= -z_max);
            pz0_ind = findfirst(zprime .>= z_max);

            # laser_xz = calculate_laser_xz(discretization,laser_spot_sigma);

            t0 = 0.2e-12;
            # dt = discretization.dt;

            interaction_v = zeros(length(t_c_subsampled),length(discretization.z));

            theta_pol = laser.theta_pol;
            electron_velocity = electron.electron_velocity;

            # t_r = calculate_t_r(discretization.x0, discretization.y0, discretization.XPRIME,
            # discretization.YPRIME, discretization.ZPRIME, discretization.Z)

            l_tc = length(t_c_subsampled)

            factors_rect = calculate_factors_rectification(calculate_green_kernel_rectification(discretization),
            calculate_laser_xz(discretization,laser_spot_sigma),
            d14,
            e0_squared,
            alpha,
            discretization,
            laser_spot_sigma,
            theta_pol)

            rect_params = RectificationParameters(laser_pulse_time_sigma,
            mz0_ind,
            pz0_ind,
            factors_rect,
            t0,
            y0_ind,
            t_c_subsampled,
            l_tc,
            calculate_t_r(discretization.x0, discretization.y0, discretization.XPRIME,
            discretization.YPRIME, discretization.ZPRIME, discretization.Z),
            electron_velocity,
            xprime,
            yprime,
            zprime)

            threaded_calculation_rectification!(interaction_v,rect_params)

            intf = linear_interpolation((vec(t_c_subsampled),[1:size(interaction_v)[2];]), interaction_v,extrapolation_bc=Flat());
            interaction_v = transpose([intf(x,y) for x in t_c, y in [1:size(interaction_v)[2];] ]);
            interaction_v[isnan.(interaction_v)] .= 0.0;
            interaction_v =  hcat(zeros(length(discretization.z),length(discretization.t)-length(t_c)),interaction_v);
            
            return t_c_subsampled, t_c, interaction_v

        end

        

        function calculate_factors_rectification(green_kernel::Array{Float64,4},
            laser_xz::Array{Float64,4},
            d14::Float64,
            e0_squared::Float64,
            alpha::Float64,
            discretization::Discretization,
            laser_spot_sigma::Float64,
            theta_pol::Float64)

            factor_rho = @. (1/(4*pi*EPSILON_0)).*green_kernel.*(1/sqrt(3.0)).*d14*e0_squared.*laser_xz.*exp.(-alpha.*discretization.YPRIME).*
            ((2*sqrt(2)/laser_spot_sigma^2).*(discretization.XPRIME.*sin(2*theta_pol)+discretization.ZPRIME.*cos(2*theta_pol)) .- alpha)
            factor_a=d14*e0_squared.*laser_xz.*exp.(-alpha.*discretization.YPRIME).*(2.0/sqrt(6.0)).*cos(2*theta_pol).*(MU_0/(4*pi)).*green_kernel
            factor_rho_Y0=(1/(4*pi*EPSILON_0))*(1/sqrt(3)).*d14*e0_squared.*laser_xz.*green_kernel
            factor_rho_Z0=-(1/(4.0*pi*EPSILON_0))*d14*e0_squared.*laser_xz.*exp.(-alpha.*discretization.YPRIME).*((2.0/sqrt(6.0)).*cos(2.0*theta_pol)).*green_kernel

            return FactorsRect(factor_rho,factor_a,factor_rho_Y0,factor_rho_Z0)

        end

        function calculate_laser_xz(discretization::Discretization, laser_spot_sigma::Float64)
            return @. exp.(-(discretization.XPRIME.^2 .+ discretization.ZPRIME.^2)./(laser_spot_sigma.^2)).*
            (discretization.ZPRIME .<= discretization.z_max) .* (discretization.ZPRIME .>= -discretization.z_max);
        end


        
        function threaded_calculation_rectification!(interaction_v::Array{Float64,2},
            rect_params::RectificationParameters)

            n = rect_params.l_tc
            Threads.@threads for time_ind in 1:n

                # if(mod(time_ind,10)==0)
                #     print("$(time_ind) out of $(n) \n")
                # end

                integrate_v_rectification!(interaction_v::Array{Float64,2},
                time_ind,
                rect_params.t_c_subsampled[time_ind],
                rect_params::RectificationParameters)

            end

        end
        

        
        function integrate_v_rectification!(interaction_v::Array{Float64,2},
            time_ind::Int64,
            t_c_subsampled_i::Float64,
            final_params::RectificationParameters)
            
            t_prime = calculate_t_prime(t_c_subsampled_i, final_params.t_r)
            
            laser_t = calculate_laser_t_triple(t_prime,
            final_params.t0,
            final_params.laser_pulse_time_sigma)

            dPhidA = calculate_internal(t_prime,
            final_params.t0,final_params.laser_pulse_time_sigma,
            final_params.factors_rect.factor_rho,final_params.factors_rect.factor_a,
            laser_t,final_params.electron_velocity)

            dPhi_Y0, dPhi_Z0 = calculate_boundaries(laser_t,final_params.y0_ind,final_params.mz0_ind, final_params.pz0_ind,
                final_params.factors_rect.factor_rho_Y0,final_params.factors_rect.factor_rho_Z0)

            @inbounds interaction_v[time_ind,:] .= trapz((final_params.xprime,final_params.yprime,final_params.zprime,:),(dPhidA)) .+
                trapz((final_params.xprime,final_params.zprime,:), dPhi_Y0) .+ trapz((final_params.xprime,final_params.yprime,:),dPhi_Z0)


        end


        
        function calculate_laser_t(t_prime::Array{Float64,4},
            t0::Float64,
            laser_pulse_time_sigma::Float64)::Array{Float64,4}
            return @. exp.(-(t_prime.-t0).^2 ./ (laser_pulse_time_sigma.^2));
        end

        function calculate_laser_t_triple(t_prime::Array{Float64,4},
            t0::Float64,
            laser_pulse_time_sigma::Float64)::Array{Float64,4}

            return  @. exp.(-(t_prime.-t0).^2 ./ (laser_pulse_time_sigma^2)) .+
            0.2.*exp.(-(t_prime.-t0).^2 ./ ((laser_pulse_time_sigma*2.5)^2)).+
            0.2.*exp.(-(t_prime .- 0.5e-12).^2 ./ ((laser_pulse_time_sigma*0.5)^2));
        end


        function calculate_t_prime(t_c_subsampled_i::Float64, t_r::Array{Float64,4})::Array{Float64,4}

            return  @. t_c_subsampled_i .- t_r;

        end



        function calculate_internal(t_prime::Array{Float64,4},t0::Float64,laser_pulse_time_sigma::Float64,
            factor_rho::Array{Float64,4},factor_a::Array{Float64,4},laser_t::Array{Float64,4}, electron_velocity::Float64)::Array{Float64,4}

            return @. factor_rho .* laser_t .- electron_velocity.*factor_a.*laser_t.*( -2.0 .*(t_prime .- t0)./laser_pulse_time_sigma.^2)
        
        end

        function calculate_boundaries(laser_t::Array{Float64,4},y0_ind::Int64,mz0_ind::Int64, pz0_ind::Int64,
            factor_rho_Y0::Array{Float64,4},factor_rho_Z0::Array{Float64,4})

            # rho_Y0 =  laser_t.*factor_rho_Y0;
            # rho_mZ0 = laser_t.*factor_rho_Z0;
            # rho_pZ0 = -rho_mZ0;

            
            # @views return rho_Y0[:,y0_ind,:,:], (rho_mZ0[:,:,mz0_ind,:] .+ rho_pZ0[:,:,pz0_ind,:])

            
            @views return laser_t[:,y0_ind,:,:].*factor_rho_Y0[:,y0_ind,:,:],
            laser_t[:,:,mz0_ind,:].*factor_rho_Z0[:,:,mz0_ind,:] .- 
            laser_t[:,:,pz0_ind,:].*factor_rho_Z0[:,:,pz0_ind,:];
        

            
            # @views return rho_Y0[:,y0_ind,:,:], (rho_mZ0[:,:,mz0_ind,:] .+ rho_pZ0[:,:,pz0_ind,:])


        end

        function interaction_potential_photodember(discretization::Discretization, material::Material, laser::Laser ,
            numerical_parameters::NumericalParameters; cutOff::Float64 = 0.0)
            # %METHOD1 Summary of this method goes here
            # %   Detailed explanation goes here

        
            EPSILON_0 = mod_physicalconstants.EPSILON_0;
            C = mod_physicalconstants.C;
            Q_E = mod_physicalconstants.Q_E;

            t_c = discretization.t[discretization.t .> discretization.t0];
            t_c_subsampled = t_c[1:numerical_parameters.tc_subsampling:end];

            z = discretization.z;

            interaction_v =  zeros(length(t_c_subsampled),length(discretization.z));
            # interaction_v = SharedArray{Float64,2}(length(t_c_subsampled),length(discretization.z))

            green_kernel  = calculate_green_kernel(discretization);
            
            gaussian_laser_spot = calculate_gaussian_laser_spot(discretization.XPRIME, discretization.ZPRIME,
            laser.laser_spot_sigma,discretization.zprime,discretization.z_max)

            l_tc = length(t_c_subsampled);

            n_exc = excited_carriers(laser, material.alpha, material.hnew);
            omega_y = calculate_omega_y(material, n_exc, discretization.YPRIME, gaussian_laser_spot);
            
            t_r = calculate_t_r(discretization.x0, discretization.y0, discretization.XPRIME,
                discretization.YPRIME, discretization.ZPRIME, discretization.Z)

            alpha = material.alpha;
            gamma = material.gamma;
            v_t = material.v_t;
            phase = material.phase;
            gamma_factor = material.gamma_factor;


            YPP_GK_g1 = @. (Q_E*alpha*v_t^2*
            n_exc.*gaussian_laser_spot.*exp.(-alpha.*discretization.YPRIME)./(omega_y.*(gamma^2 .+ 4*omega_y.^2))).*
            (discretization.y0 .- discretization.YPRIME).*green_kernel;
            
            threaded_calculation!(l_tc, interaction_v, discretization.xprime, discretization.yprime, discretization.zprime, 
            omega_y, t_c_subsampled, t_r, phase, gamma, gamma_factor, YPP_GK_g1; cut0ff=cutOff)

            intf = linear_interpolation((vec(t_c_subsampled),[1:size(interaction_v)[2];]), interaction_v,extrapolation_bc=Flat());
            interaction_v = transpose([intf(x,y) for x in t_c, y in [1:size(interaction_v)[2];] ]);
            interaction_v[isnan.(interaction_v)] .= 0.0;
            interaction_v =  hcat(zeros(length(z),length(discretization.t)-length(t_c)),interaction_v);
            
            return t_c_subsampled, t_c, interaction_v

        end

        function calculate_gaussian_laser_spot(XPRIME::Array{Float64,4}, ZPRIME::Array{Float64,4},
            laser_spot_sigma::Float64,zprime::Array{Float64,1},z_max::Float64)
            gaussian_laser_spot = @. exp.(-(XPRIME.^2+ZPRIME.^2)./(2*laser_spot_sigma.^2));
            gaussian_laser_spot[:,:,abs.(vec(zprime)) .> z_max,:] .= 0.0;
            return gaussian_laser_spot
        end

        function calculate_t_r(x0::Float64, y0::Float64, 
            XPRIME::Array{Float64,4}, YPRIME::Array{Float64,4}, 
            ZPRIME::Array{Float64,4}, Z::Array{Float64,4})
            
            
            C = mod_physicalconstants.C;
            return @. (1/C).*sqrt.((x0 .- XPRIME).^2+(y0 .- YPRIME).^2+(Z .- ZPRIME).^2);

        end

        function integrate_v!(interaction_v::Array{Float64,2}, ind::Int64, 
            dx::Array{Float64,1}, dy::Array{Float64,1}, dz::Array{Float64,1}, integrand::Array{Float64,4})

            
            EPSILON_0 = mod_physicalconstants.EPSILON_0;

            @inbounds interaction_v[ind,:] .= (1/(4*pi*EPSILON_0)).*trapz((dx,dy,dz,:),
            integrand)

        end

  
        function calculate_g2(omega_y::Array{Float64,4}, t_prime::Array{Float64,4}, phase::Float64, gamma::Float64, gamma_factor::Float64; cut0ff::Float64=0.0)
            
            g = @. (4.0 .* omega_y.*cos.(omega_y.*t_prime .+ phase).+2.0*gamma.*sin.(omega_y.*t_prime .+ phase)) .*
            exp.(-gamma.*gamma_factor.*t_prime./2.0);

            g[t_prime.<cut0ff] .= 0.0

            return g
        end

        function calculate_green_kernel_rectification(discretization::Discretization)
            return @. ((discretization.x0 .- discretization.XPRIME).^2 .+ 
                (discretization.y0 .- discretization.YPRIME).^2 .+ (discretization.Z .- discretization.ZPRIME).^2).^(-1/2);
        end

        function calculate_green_kernel(discretization::Discretization)
            # %METHOD1 Summary of this method goes here
            # %   Detailed explanation goes here
            return @. ((discretization.x0 .- discretization.XPRIME).^2 .+ 
                (discretization.y0 .- discretization.YPRIME).^2 .+ (discretization.Z .- discretization.ZPRIME).^2).^(-3/2);
        end

        function threaded_calculation!(n::Int, interaction_v::Array{Float64,2}, xprime::Array{Float64,1}, yprime::Array{Float64,1}, zprime::Array{Float64,1}, 
            omega_y::Array{Float64,4}, t_c_subsampled::Array{Float64,1}, t_r::Array{Float64,4}, phase::Float64, gamma::Float64, gamma_factor::Float64, YPP_GK_g1::Array{Float64,4};cut0ff::Float64=0.0)
            
            @inbounds Threads.@threads for time_ind in 1:n

                # if(mod(time_ind,10)==0)
                #     print("$(time_ind) out of $(n) \n")
                # end

                integrate_v!(interaction_v, time_ind, xprime,yprime,zprime,
                calculate_g2(omega_y, t_c_subsampled[time_ind] .- t_r, phase, gamma, gamma_factor; cut0ff = cut0ff).*YPP_GK_g1)

            end
        end

       

    end

    module mod_eels

        using ..mod_customtypes, ..mod_electron, ..mod_physicalconstants, ..mod_utils
        using FFTW
        # FFTW.set_num_threads(8)
        using ThreadsX 
        using Trapz

        export calculate_ft

        function calculate_ft(discretization::Discretization , interact_v::Array{Float64,2}, electron::Electron)
                    
            t = discretization.t;
            omega = discretization.omega;
            z = discretization.z;
            deltaz = discretization.z[2] - discretization.z[1];
            
            # beta = ft_beta(interact_v, t, omega, z, deltaz, electron.electron_velocity)

            return ft_cal(electron.electron_velocity, discretization.omega, discretization.t,
            ft_beta(interact_v, t, omega, z, deltaz, electron.electron_velocity))

        end

        function ft_cal(electron_velocity::Float64, omega::Array{Float64,1}, 
            t::Array{Float64,1}, beta::Array{ComplexF64,1})
            t_map = zeros(ComplexF64,size(t))
            ThreadsX.map!(tt-> ft_parameter(omega, tt, beta), t_map, t);
            t_map .*= ((- Q_E / (1.0im*HBAR*electron_velocity)));
            # return transpose(exp.(-circshift((-1.0 ./ (1.0im*HBAR*electron_velocity)).*Q_E.*t_map,(ceil(length(t)/2)))));
            return transpose(exp.(-circshift(t_map,(ceil(length(t)/2)))));
        end

        function ft_beta(interact_v::Array{Float64,2}, t::Array{Float64,1}, 
            omega::Array{Float64,1}, z::Array{Float64,1}, 
            deltaz::Float64, electron_velocity::Float64)::Array{ComplexF64,1}
            return vec(sum(fftshift(fft2(interact_v,length(t)),(2,))./maximum(omega).*[exp(-1.0im*omg*zz/electron_velocity) for zz in z, omg in omega], dims = 1).*deltaz)
        end

        # function ft_parameter(omega::Array{Float64,1}, tt::Float64, 
        #     beta::Array{ComplexF64,1})::Float64
        #     trapz(omega,2.0.*real.(exp.(1.0im.*omega.*tt).*beta))
        # end 

        function ft_parameter(omega::Array{Float64,1}, tt::Float64, 
            beta::Array{ComplexF64,1})::Float64
            trapz(omega,ft_integrand(omega, tt, beta))
        end 

        function ft_integrand(omega::Array{Float64,1}, tt::Float64, 
            beta::Array{ComplexF64,1})::Array{Float64,1}
            return 2.0.*real.(exp.(1.0im.*omega.*tt).*beta)
        end

        function incoherent_convolution_fast(psi::Array{Float64,2}, w::AbstractArray{Float64,2}, t_w::Array{Float64,1}, e_w::Array{Float64,1} ,
            w_cut_off_factor::Float64 = 0.01)

            ThreadsX.map!(x->isnan(x) ? 0.0 : x, psi,psi)
            psi_sum = zeros(size(psi));

            w_cutOff = w_cut_off_factor*maximum(w[:]);

            l_t_w = length(t_w)
            l_e_w = length(e_w)

            Threads.@threads for k = 1:l_t_w * l_e_w
                
                incoherent_circ!(psi_sum, w, w_cutOff, psi, t_w, e_w, k, l_t_w, l_e_w)
            
            end

            psi_sum./=trapz((:,e_w),psi_sum)
            return psi_sum./=maximum(psi[:])

        end

        function incoherent_convolution_fast(psi::SubArray{Float64, 2, Array{Float64, 3}}, w::AbstractArray{Float64,2}, t_w::Array{Float64,1}, e_w::Array{Float64,1} ,
            w_cut_off_factor::Float64 = 0.01)

            ThreadsX.map!(x->isnan(x) ? 0.0 : x, psi,psi)
            psi_sum = zeros(size(psi));

            w_cutOff = w_cut_off_factor*maximum(w[:]);

            l_t_w = length(t_w)
            l_e_w = length(e_w)

            Threads.@threads for k = 1:l_t_w * l_e_w
                
                incoherent_circ!(psi_sum, w, w_cutOff, psi, t_w, e_w, k, l_t_w, l_e_w)
            
            end

            psi_sum./=trapz((:,e_w),psi_sum)
            return psi_sum./=maximum(psi[:])

        end

        function incoherent_circ!(psi_sum::Array{Float64,2}, w::Array{Float64,2}, w_cutOff::Float64,
            psi::Array{Float64,2}, t_w::Array{Float64,1}, 
            e_w::Array{Float64,1}, k::Int, l_t_w::Int, l_e_w::Int)

            t_ind,e_ind = index_calculator_2d(k, l_t_w, l_e_w)

            return w[k]<w_cutOff ? psi_sum : psi_sum.+=w[k].*circshift(psi,(-ceil(length(t_w)/2)+t_ind,-ceil(length(e_w)/2)+e_ind))

        end

        function incoherent_convolution(psi::Array{Float64,2}, w::AbstractArray{Float64,2}, t_w::Array{Float64,1}, e_w::Array{Float64,1} ,
            w_cut_off_factor::Float64 = 0.01)

            psi_sum = zeros(size(psi));

            w_cutOff = w_cut_off_factor*maximum(w[:]);

            # % Example dont use such code. It becomes hard to modify
            # % w_cutOff = 0.01*max(w(:));
            for t_ind = 1:length(t_w)
                for e_ind = 1:length(e_w)
                    
                    if w[t_ind,e_ind] < w_cutOff
                        continue
                    end
                    
                    psi_sum = psi_sum + w[t_ind,e_ind].*circshift(circshift(psi,(-ceil(length(t_w)/2)+t_ind,0)),(0,-ceil(length(e_w)/2)+e_ind));
                    
                end
            end

            psi = psi_sum;
            psi_incoherent = psi./trapz((:,e_w),psi);
            psi_incoherent ./= maximum(psi_incoherent[:]);
            return psi_incoherent

        end

        function psi_subsampled(sub_sample_factor::Int64, psi::Array{Float64,2}, e_w::Array{Float64,1})
                    
            psi_sub = psi[:,1:sub_sample_factor:end];
            psi_sub = psi_sub./trapz((:,e_w),psi_sub);
            return psi_sub

        end

        function calculate_psi_coherent_old(discretization::Discretization,elec::Electron, f_t::Union{Matrix{ComplexF64},
            AbstractArray{ComplexF64}})
                    
            deltat = discretization.deltat;
            t = discretization.t;
            
            # delta_t_rep = repeat(deltat,outer = (1,length(t)));
            # t_rep = repeat(t,outer=(1,length(deltat)))';

            TT = repeat(t,outer=(1,length(deltat)))' .- repeat(deltat,outer = (1,length(t)));
            FT = repeat(f_t,outer = (length(deltat),1));
            

            # psi_coherent = FT.*exp.(-(t_rep.-delta_t_rep).^2.0 ./ (2.0*elec.electron_time_coherent_sigma^2));
            psi_coherent = FT.*exp.(-(TT).^2 ./ (2.0*elec.electron_time_coherent_sigma^2));
            psi_coherent = abs2.(fftshift(fft2(psi_coherent,length(t)),(2,)));
            # copy!(psi_coherent,abs2.(fftshift(fft2(psi_coherent,length(t)),(2,))));
           
            # psi_coherent = (abs.(psi_coherent)).^2;
            
            return psi_coherent./trapz((:,discretization.energy),psi_coherent);
            
            # return psi_coherent
            
        end

        function calculate_psi_coherent(discretization::Discretization,elec::Electron, f_t::Union{Matrix{ComplexF64},
            AbstractArray{ComplexF64,2}})
                    
            t_len = length(discretization.t);
            deltat_len = length(discretization.deltat);
            TT = repeat(discretization.t,outer=(1,deltat_len))' .- repeat(discretization.deltat,outer = (1,t_len));
            FT = repeat(f_t,outer = (deltat_len,1)).*exp.(-(TT).^2 ./ (2.0*elec.electron_time_coherent_sigma^2));

            psi_coherent_r = zeros(size(TT));
            abs2fftshift!(psi_coherent_r, FT, t_len)
            return psi_coherent_r./trapz((:,discretization.energy),psi_coherent_r);
  
            
        end
        
        function abs2fftshift!(psi_coherent_r::Array{Float64,2}, psi_coherent::Array{ComplexF64,2}, n::Int64)
            copy!(psi_coherent_r , abs2.(fftshift(fft2(psi_coherent,n),(2,))));
        end


        # function fft2(a::Union{AbstractArray{Float64},Matrix{ComplexF64}}, n::Int64=size(a,2))
        #     if n <= size(a,2)
        #         return fft(a[:,1:n],(2,))
        #     else
        #         A = zeros(size(a,1),n)
        #         @views A[:,1:size(a,2)] .= a
        #         return fft(A,(2,))
        #     end
                
        # end

        function fft2(a::Array{ComplexF64,2}, n::Int64=size(a,2))
            if n <= size(a,2)
                return fft(a[:,1:n],(2,))
            else
                A = zeros(size(a,1),n)
                @views A[:,1:size(a,2)] .= a
                return fft(A,(2,))
            end
                
        end

        function fft2(a::Array{Float64,2}, n::Int64=size(a,2))
            if n <= size(a,2)
                return fft(a[:,1:n],(2,))
            else
                A = zeros(size(a,1),n)
                @views A[:,1:size(a,2)] .= a
                return fft(A,(2,))
            end
                
        end




    end
end
