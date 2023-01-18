module DiscretizationS

    export Space, space_config, space_grid

    struct Space
        x0::Float64
        y0::Float64    
        xprime::Array{Float64}
        yprime::Array{Float64}
        zprime::Array{Float64}
        z::Array{Float64}
        z_max::Float64
    end

    function space_config(;x0::Float64, y0::Float64, 
        d_xprime::Float64, d_yprime::Float64, d_zprime::Float64, 
        xprime_max::Float64, yprime_max::Float64, zprime_max::Float64, 
        d_z::Float64, zmax::Float64, z_max::Float64)

        xprime = -xprime_max : d_xprime : xprime_max
        yprime = 0.0 : d_yprime : yprime_max
        zprime = -zprime_max : d_zprime : zprime_max

        z = -zmax : d_z : zmax

        return Space(x0, y0, xprime, yprime, zprime, z, z_max)

    end
    
    function space_grid(space::Space)
        
        XPRIME = zeros(length(space.xprime), length(space.yprime), length(space.zprime), length(space.z),)
        XPRIME = [i for i in space.xprime, _ in space.yprime, _ in space.zprime, _ in space.z]
        YPRIME = [i for _ in space.xprime, i in space.yprime, _ in space.zprime, _ in space.z]
        ZPRIME = [i for _ in space.xprime, _ in space.yprime, i in space.zprime, _ in space.z]
        Z = [i for _ in space.xprime, _ in space.yprime, _ in space.zprime, i in space.z]

        return XPRIME, YPRIME, ZPRIME, Z

    end
end

module DiscretizationEnergyTime

    using ..PhysicalContants
    export EnergyTime, energytime_config

    struct EnergyTime

        ddt::Float64
        deltat::Array{Float64}
        t0::Float64
        t::Array{Float64}
        dt::Float64
        omega::Array{Float64}
        energy::Array{Float64}

    end

    function energytime_config(;t0::Float64, ddt::Float64, delay_max::Float64, fs::Float64, l::Float64)
        
        deltat = -delay_max : ddt : delay_max
        t, omega, energy = incidence(fs, l); 
        dt = t[2] - t[1];

        return EnergyTime(ddt, deltat, t0, t, dt, omega, energy)
    
    end

    function  incidence(fs::Float64, l::Float64)
    
        dt = 1/fs;#[s]
        domega = fs/l*(2*pi);#[Hz]
        t =  (-l/2:l/2).*dt;#[s
        omega = (-l/2:l/2)*domega;#[Hz
        energy = omega.*PhysicalContants.HBAR/PhysicalContants.Q_E;#[eV]
    
        return t, omega, energy
        
    end
end

