module PhysicalContants

    export   HBAR,EPSILON_0, C, M_E, Q_E, MU_0, ETA_0

    global const HBAR       = 6.62607004e-34/(2*pi)::Float64  # s J
    global const EPSILON_0  = 8.8541878128e-12::Float64  # V/m
    global const C = 299792458.0::Float64 
    global const M_E = 9.109e-31::Float64 
    global const Q_E = 1.60217662e-19::Float64  # coulomb::Float64 
    global const MU_0 = 1.25663706212e-6::Float64 #[H/m]::Float64 
    global const ETA_0 = 377.0::Float64 #[Ohm]

end
