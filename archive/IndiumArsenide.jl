
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

function IndiumArsenide()
    inas = Material()
    inas.alpha = 7e6;#[m^-1]
    consts = constants_fundamental();
    EPSILON_0 = consts.EPSILON_0;
    M_E = consts.M_E;
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
    inas.m_eq = inas.mh;#p-type InAs
    inas.epsilon_e = photoexcited_electron_energy(inas);
    inas.v_t = velocity_t(inas);
    inas.me0tilda = photoelectron_mass(inas);
    inas.d14 = 237.6e-12;#[V/m]
    return inas
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

function calculate_omega_y(mat::Material, n_exc::Float64,
    yprime_grid::Array{Float64}, gaussian_laser_spot::Array{Float64})
    consts = constants_fundamental();
    Q_E = consts.Q_E;
    
    omega_y = sqrt.( (Q_E^2/mat.kappa).*(n_exc.*gaussian_laser_spot.*exp.(-mat.alpha.*yprime_grid).*
    (1/mat.me0tilda+1/mat.mh) .+ mat.n_eq/mat.m_eq) .- (mat.gamma/2)^2 );

    return omega_y
end

function velocity_t(mat::Material)
    # %METHOD1 Summary of this method goes here
    # %   Detailed explanation goes here
    
    consts = constants_fundamental()
    Q_E = consts.Q_E;
    
    v_t = sqrt((2*mat.epsilon_e.*Q_E.*(1+mat.alpha_gamma.*mat.epsilon_e))
    ./(3*mat.me.*(1+4*mat.alpha_gamma.*mat.epsilon_e.*(1+mat.alpha_gamma.*mat.epsilon_e))));#[m/s]

    return v_t

end