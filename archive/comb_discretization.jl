module DiscretizationS

    export Space, space_config, space_grid

    struct Discretization
        x0::Float64
        y0::Float64    
        xprime::Array{Float64}
        yprime::Array{Float64}
        zprime::Array{Float64}
        z::Array{Float64}
        z_max::Float64

        ddt::Float64
        deltat::Array{Float64}
        t0::Float64
        t::Array{Float64}
        dt::Float64
        omega::Array{Float64}
        energy::Array{Float64}

        XPRIME::Array{Float64}
        YPRIME::Array{Float64}
        ZPRIME::Array{Float64}
        Z::Array{Float64}
    end

    function Discretization(;x0::Float64, y0::Float64, 
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
        ddt, deltat, t0, t, dt, omega, energy,XPRIME, YPRIME, ZPRIME, Z)

    end
    
    function space_grid(xprime::Array{Float64},
        yprime::Array{Float64},
        zprime::Array{Float64},
        z::Array{Float64})
        
        XPRIME = zeros(length(xprime), length(yprime), length(zprime), length(z),)
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
        energy = omega.*PhysicalContants.HBAR/PhysicalContants.Q_E;#[eV]
    
        return t, omega, energy
        
    end

end
