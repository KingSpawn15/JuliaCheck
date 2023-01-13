mutable struct Discretization

    ddt::Float64
    deltat::Array{Float64}
    z::Array{Float64}
    deltaz::Float64
    t::Array{Float64}
    omega::Array{Float64}
    energy::Array{Float64}
    x0::Float64
    y0::Float64
    t0::Float64
    xprime::Array{Float64}
    yprime::Array{Float64}
    zprime::Array{Float64}
    d_xprime::Float64
    d_yprime::Float64
    d_zprime::Float64
    XPRIME::Array{Float64}
    YPRIME::Array{Float64}
    ZPRIME::Array{Float64}
    Z::Array{Float64}
    z_max::Float64
    dt::Float64

    function Discretization(discretization_parameters::discretization_params)
        discretization = new()
        discretization.x0 = discretization_parameters.x0;
        discretization.y0 = discretization_parameters.y0;
        discretization.d_xprime = discretization_parameters.d_xprime;
        discretization.d_yprime = discretization_parameters.d_yprime;
        discretization.d_zprime = discretization_parameters.d_zprime;
        discretization.ddt, discretization.deltat = delay(discretization_parameters.ddt, discretization_parameters.delay_max);
        discretization.z , discretization.deltaz = zstep(discretization_parameters.ddz, discretization_parameters.zmax);
        discretization.t, discretization.omega, discretization.energy = incidence(discretization_parameters.fs, discretization_parameters.l); 
        discretization.dt = discretization.t[2] - discretization.t[1];
        discretization.t0 = discretization_parameters.t0;

        
        discretization.xprime , discretization.yprime , discretization.zprime , 
            discretization.XPRIME, discretization.YPRIME, discretization.ZPRIME, discretization.Z = 
        space_points_grid(discretization_parameters , discretization.z);
        discretization.z_max = discretization_parameters.z_max;
        return discretization
    end
    
end

function delay(ddt::Float64, delay_max::Float64)
    # %DISCRETIZATION_DT Summary of this function goes here
    # %   Detailed explanation goes here
    # % Time delay
    
    deltat = collect(transpose([-delay_max:ddt:delay_max;]))
    
    return ddt, deltat

end

function zstep(ddz::Float64, zmax::Float64)
    # %DISCRETIZATION_Z Summary of this function goes here
    # %   Detailed explanation goes here
    
    z = collect(transpose([-zmax : ddz : zmax;]))#[m]
    deltaz = z[2]-z[1];
    
    return z , deltaz
    
end

function  incidence(fs::Float64, l::Float64)
    
    consts = constants_fundamental()
    HBAR = consts.HBAR;
    Q_E = consts.Q_E;
    dt = 1/fs;#[s]
    domega = fs/l*(2*pi);#[Hz]
    t =  collect(transpose([-l/2:l/2;])).*dt;#[s
    omega = collect(transpose([-l/2:l/2;]))*domega;#[Hz
    energy = omega.*HBAR/Q_E;#[eV]

    return t, omega, energy
    
end

# function space_points_grid!(dis_obj::Discretization, discretization_parameters::discretization_params , z::Array{Float64})

#     dis_obj.xprime = collect(transpose([- discretization_parameters.xprime_max : discretization_parameters.d_xprime :  discretization_parameters.xprime_max ;]))
#     dis_obj.yprime = collect(transpose([0 : discretization_parameters.d_yprime :  discretization_parameters.yprime_max ;]))
#     dis_obj.zprime = collect(transpose([- discretization_parameters.zprime_max : discretization_parameters.d_zprime :  discretization_parameters.zprime_max ;]))

#     dis_obj.XPRIME,dis_obj.YPRIME, dis_obj.ZPRIME,dis_obj.Z = ndgrid_array(vec(dis_obj.xprime),vec(dis_obj.yprime),vec(dis_obj.zprime),vec(z));

#     return nothing
# end

function space_points_grid(discretization_parameters::discretization_params , z::Array{Float64})

    xprime = collect(transpose([- discretization_parameters.xprime_max : discretization_parameters.d_xprime :  discretization_parameters.xprime_max ;]))
    yprime = collect(transpose([0 : discretization_parameters.d_yprime :  discretization_parameters.yprime_max ;]))
    zprime = collect(transpose([- discretization_parameters.zprime_max : discretization_parameters.d_zprime :  discretization_parameters.zprime_max ;]))

    XPRIME,YPRIME, ZPRIME,Z = ndgrid_array(vec(xprime),vec(yprime),vec(zprime),vec(z));

    return xprime , yprime , zprime , XPRIME, YPRIME, ZPRIME, Z

end

function space_points_grid_2(discretization_parameters::discretization_params , z::Array{Float64})

    xprime = collect(transpose([- discretization_parameters.xprime_max : discretization_parameters.d_xprime :  discretization_parameters.xprime_max ;]))
    yprime = collect(transpose([0 : discretization_parameters.d_yprime :  discretization_parameters.yprime_max ;]))
    zprime = collect(transpose([- discretization_parameters.zprime_max : discretization_parameters.d_zprime :  discretization_parameters.zprime_max ;]))

    # XPRIME,YPRIME, ZPRIME,Z = ndgrid_array(vec(xprime),vec(yprime),vec(zprime),vec(z));

    XPRIME = [i for i in vec(xprime), _ in vec(yprime), _ in vec(zprime), _ in vec(z)]
    YPRIME = [i for _ in vec(xprime), i in vec(yprime), _ in vec(zprime), _ in vec(z)]
    ZPRIME = [i for _ in vec(xprime), _ in vec(yprime), i in vec(zprime), _ in vec(z)]
    Z = [i for _ in vec(xprime), _ in vec(yprime), _ in vec(zprime), i in vec(z)]

    return xprime , yprime , zprime , XPRIME, YPRIME, ZPRIME, Z
end

function space_points_grid_3(discretization_parameters::discretization_params)

    xprime = - discretization_parameters.xprime_max : discretization_parameters.d_xprime :  discretization_parameters.xprime_max 
    yprime = 0 : discretization_parameters.d_yprime :  discretization_parameters.yprime_max 
    zprime = - discretization_parameters.zprime_max : discretization_parameters.d_zprime :  discretization_parameters.zprime_max 
    z = - discretization_parameters.z_max : discretization_parameters.ddz : discretization_parameters.z_max
    # XPRIME,YPRIME, ZPRIME,Z = ndgrid_array(vec(xprime),vec(yprime),vec(zprime),vec(z));

    return xprime , yprime , zprime , z
end
