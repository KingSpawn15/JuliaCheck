mutable struct constants_fundamental

    HBAR::Float64       
    EPSILON_0::Float64
    C::Float64
    M_E::Float64
    Q_E::Float64
    MU_0::Float64
    ETA_0::Float64

    function constants_fundamental()
        h = 6.62607004e-34;
        consts = new()
        consts.HBAR       = h/(2*pi); # s J
        consts.EPSILON_0  = 8.8541878128e-12; # V/m
        consts.C = 299792458;
        consts.M_E = 9.109e-31;
        consts.Q_E = 1.60217662e-19; # coulomb;
        consts.MU_0 = 1.25663706212e-6;#[H/m];
        consts.ETA_0 = 377;#[Ohm]
        return consts
    end

end