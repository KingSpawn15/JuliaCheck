function interaction_potential_photodember()
    return 0.0
end


function interaction_potential_photodember_optimized(discretization::Discretization, material::Material, laser::laser_params ,
     numerical_parameters::numerical_params)
# %METHOD1 Summary of this method goes here
# %   Detailed explanation goes here

    consts = constants_fundamental();
    EPSILON_0 = consts.EPSILON_0;
    C = consts.C;
    Q_E = consts.Q_E;

    t_c = discretization.t[discretization.t .> discretization.t0];
    t_c_subsampled = t_c[1:numerical_parameters.tc_subsampling:end];

    x0 = discretization.x0;
    y0 = discretization.y0;
    z = discretization.z;

    xprime = discretization.xprime;
    yprime = discretization.yprime;
    zprime = discretization.zprime;

    z_max = discretization.z_max;

    # XPRIME = discretization.XPRIME;
    # YPRIME = discretization.YPRIME;
    # ZPRIME= discretization.ZPRIME;

    # Z = discretization.Z;

    interaction_v =  zeros(length(t_c_subsampled),length(discretization.z));

    # green_kernel  = calculate_green_kernel(discretization);
    # gaussian_laser_spot = exp.(-(XPRIME.^2+ZPRIME.^2)./(2*laser.laser_spot_sigma.^2));
    # gaussian_laser_spot[:,:,abs.(vec(zprime)) .> discretization.z_max,:] .= 0;

    laser_spot_sigma = laser.laser_spot_sigma

    l_tc = length(t_c_subsampled);

    n_exc = excited_carriers(laser, material.alpha, material.hnew);
    # omega_y = calculate_omega_y(material, n_exc, YPRIME, gaussian_laser_spot);

    # t_r = (1/C)*sqrt.((x0 .- XPRIME).^2+(y0 .- YPRIME).^2+(Z .- ZPRIME).^2);

    alpha = material.alpha;
    gamma = material.gamma;
    v_t = material.v_t;
    phase = material.phase;
    gamma_factor = material.gamma_factor;

    # begin = '['; last = '] ... Photodember'; arrow = '>';
    # body = '='; empty = ' ';
    # disp("Calculating potential photodember ....");
    # %             omega_y = 0.05 * omega_y;
    
    for time_ind = 1:l_tc

        # body_disp = repmat(body , [1,fix(time_ind / length(t_c_subsampled) * 70)]);
        # empty_disp = repmat(empty , [1,70 - fix(time_ind / length(t_c_subsampled) * 70)]);
        # disp([begin, body_disp, arrow, empty_disp , last]);
        if(mod(time_ind,10)==0)
            print("$(time_ind) out of $(l_tc) \n")
        end
        # t_prime = t_c_subsampled[time_ind] .- t_r;

        # g = (Q_E*alpha*v_t^2*
        #     n_exc.*gaussian_laser_spot.*exp.(-alpha.*YPRIME)./(omega_y.*(gamma^2 .+ 4*omega_y.^2))).*
        #     (4*omega_y.*cos.(omega_y.*t_prime .+ phase)+2*gamma.*sin.(omega_y.*t_prime .+ phase)).*
        #     exp.(-gamma.*gamma_factor.*t_prime./2);
        # g[t_prime.<0] .= 0;

        # dv = green_kernel.*g.*(y0 .- YPRIME);

        # dv = [calculate_dv(x0, y0, xi, yi, zi, zzi,
        # z_max, laser_spot_sigma,
        # t_c_subsampled[time_ind], 
        # n_exc, 
        # material) for xi in vec(xprime), yi in vec(yprime), zi in vec(zprime), zzi in vec(z)]

        
        
        l_x = length(vec(xprime))
        l_y = length(vec(yprime))
        l_z = length(vec(zprime))
        l_zz = length(vec(z))
        
        function get_coords(k::Int)
            idx, idy, idz, idzz = index_calculator_4d(k, l_x, l_y, l_z, l_zz)
            xi = xprime[idx]; yi = yprime[idy]; zi = zprime[idz]; zzi = z[idzz];
            return (xi, yi, zi, zzi)
        end

        dv = zeros(length(vec(xprime)),length(vec(yprime))*length(vec(zprime))*length(vec(z)))

        t_c_subsampled_i = t_c_subsampled[time_ind]
        
        for k in 1:l_x*l_y*l_z*l_zz

            @inbounds    dv[k] = calculate_dv(x0, y0, get_coords(k),
                        z_max, laser_spot_sigma,
                        t_c_subsampled_i, 
                        n_exc, 
                        material)

        end

        dv = reshape(dv, l_x, l_y, l_z, l_zz)

        interaction_v[time_ind,:] .= (1/(4*pi*EPSILON_0)).*trapz((xprime,yprime,zprime,:),dv)

    end

    int_v = copy(interaction_v)
    @exfiltrate
    
    intf = linear_interpolation((vec(t_c_subsampled),[1:size(interaction_v)[2];]), interaction_v,extrapolation_bc=Flat());
    interaction_v = transpose([intf(x,y) for x in t_c, y in [1:size(interaction_v)[2];] ]);
    interaction_v[isnan.(interaction_v)] .= 0.0;
    interaction_v =  hcat(zeros(length(z),length(discretization.t)-length(t_c)),interaction_v);
    
    

    return t_c_subsampled, t_c, interaction_v

end

function calculate_dv(x0::Float64, y0::Float64, coords::Tuple{Float64, Float64, Float64,Float64},
    z_max::Float64, laser_spot_sigma::Float64,
    t_c_subsampled::Float64, 
    n_exc::Float64, 
    mat::Material)

    xi = coords[1]
    yi = coords[2]
    zi = coords[3]
    zzi = coords[4]

    @time gaussian_laser_spot = calculate_gaussian_laser_spot(xi, zi, z_max, laser_spot_sigma)
    @time t_prime = calculate_t_prime(t_c_subsampled, x0, y0, xi, yi, zzi, zi)
    @time omega_y = calculate_omega_y_optimized(mat, n_exc, yi, gaussian_laser_spot)
    
    alpha = mat.alpha
    gamma = mat.gamma
    phase = mat.phase
    v_t = mat.v_t
    gamma_factor = mat.gamma_factor

    consts = constants_fundamental();
    Q_E = consts.Q_E;

    @time g = (t_prime >= 0.0) ? (Q_E*alpha*v_t^2*
    n_exc*gaussian_laser_spot*exp(-alpha*yi)/(omega_y*(gamma^2 + 4.0*omega_y^2)))*
    (4.0*omega_y*cos(omega_y*t_prime + phase)+2.0*gamma*sin(omega_y*t_prime + phase))*
    exp(-gamma*gamma_factor*t_prime/2.0) : 0.0

    # g = 0.0;

    # if t_prime >= 0.0
    #     g = (Q_E*alpha*v_t^2*
    #         n_exc*gaussian_laser_spot*exp(-alpha*yi)/(omega_y*(gamma^2 + 4*omega_y^2)))*
    #         (4*omega_y*cos(omega_y*t_prime + phase)+2*gamma*sin(omega_y*t_prime + phase))*
    #         exp(-gamma*gamma_factor*t_prime/2);
    # end

    # g = (Q_E*alpha*v_t^2*
    #         n_exc*gaussian_laser_spot*exp(-alpha*yi)/(omega_y*(gamma^2 + 4*omega_y^2)))*
    #         (4*omega_y*cos(omega_y*t_prime + phase)+2*gamma*sin(omega_y*t_prime + phase))*
    #         exp(-gamma*gamma_factor*t_prime/2);
    # if t_prime < 0.0
    #     g = 0.0
    # end

    @timeit "final_mul" calculate_green_kernel_optimized(x0, y0, xi, yi, zi, zzi) * g * (y0 - yi)

end

function calculate_green_kernel_optimized(x0::Float64, y0::Float64, 
    xi::Float64, yi::Float64, zi::Float64, zzi::Float64)
    
    # %METHOD1 Summary of this method goes here
    # %   Detailed explanation goes here
    green_kernel = ((x0 - xi)^2 + (y0 - yi)^2 + (zzi - zi)^2)^(-3/2);
    return green_kernel

end

function calculate_gaussian_laser_spot(xi::Float64, zi::Float64, z_max::Float64, laser_spot_sigma::Float64)
    
    gaussian_laser_spot = exp(-(xi^2+zi^2)./(2*laser_spot_sigma^2));
    if (abs(zi) > z_max)
        gaussian_laser_spot = 0.0
    end
    
    return gaussian_laser_spot

end

function calculate_omega_y_optimized(mat::Material, n_exc::Float64,
    yi::Float64, gaussian_laser_spot::Float64)
    consts = constants_fundamental();
    Q_E = consts.Q_E;
    
    omega_y = sqrt( (Q_E^2/mat.kappa)*(n_exc*gaussian_laser_spot*exp(-mat.alpha.*yi)*
    (1/mat.me0tilda+1/mat.mh) + mat.n_eq/mat.m_eq) - (mat.gamma/2)^2 );

    return omega_y
end

function calculate_t_prime(t_c_subsampled::Float64, x0::Float64, y0::Float64, 
    xi::Float64, yi::Float64, zzi::Float64, zi::Float64)

    consts = constants_fundamental();
    C = consts.C;

    t_prime = t_c_subsampled - (1/C)*sqrt((x0 - xi)^2+(y0 - yi)^2+(zzi - zi)^2);

    return t_prime
end
