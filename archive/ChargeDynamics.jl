function interaction_potential_photodember(discretization::Discretization, material::Material, laser::laser_params ,
    numerical_parameters::numerical_params)
# %METHOD1 Summary of this method goes here
# %   Detailed explanation goes here

   consts = constants_fundamental();
   EPSILON_0 = consts.EPSILON_0;
   C = consts.C;
   Q_E = consts.Q_E;

   t_c = discretization.t[discretization.t .> discretization.t0];
   t_c_subsampled = t_c[1:numerical_parameters.tc_subsampling:end];

   z = discretization.z;

   interaction_v =  zeros(length(t_c_subsampled),length(discretization.z));

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


   YPP_GK_g1 = (Q_E*alpha*v_t^2*
   n_exc.*gaussian_laser_spot.*exp.(-alpha.*discretization.YPRIME)./(omega_y.*(gamma^2 .+ 4*omega_y.^2))).*
   (discretization.y0 .- discretization.YPRIME).*green_kernel;
   
   threaded_calculation!(l_tc, interaction_v, discretization.xprime, discretization.yprime, discretization.zprime, 
   omega_y, t_c_subsampled, t_r, phase, gamma, gamma_factor, YPP_GK_g1)

   intf = linear_interpolation((vec(t_c_subsampled),[1:size(interaction_v)[2];]), interaction_v,extrapolation_bc=Flat());
   interaction_v = transpose([intf(x,y) for x in t_c, y in [1:size(interaction_v)[2];] ]);
   interaction_v[isnan.(interaction_v)] .= 0.0;
   interaction_v =  hcat(zeros(length(z),length(discretization.t)-length(t_c)),interaction_v);
   
   return t_c_subsampled, t_c, interaction_v

end

function calculate_gaussian_laser_spot(XPRIME::Array{Float64}, ZPRIME::Array{Float64},
    laser_spot_sigma::Float64,zprime::Array{Float64},z_max::Float64)
    gaussian_laser_spot = exp.(-(XPRIME.^2+ZPRIME.^2)./(2*laser_spot_sigma.^2));
    gaussian_laser_spot[:,:,abs.(vec(zprime)) .> z_max,:] .= 0.0;
    return gaussian_laser_spot
end

function calculate_t_r(x0::Float64, y0::Float64, 
    XPRIME::Array{Float64}, YPRIME::Array{Float64}, 
    ZPRIME::Array{Float64}, Z::Array{Float64})
    consts = constants_fundamental();
    
    C = consts.C;
    return (1/C)*sqrt.((x0 .- XPRIME).^2+(y0 .- YPRIME).^2+(Z .- ZPRIME).^2);

end

function integrate_v!(interaction_v::Array{Float64}, ind::Int, 
    dx::Array{Float64}, dy::Array{Float64}, dz::Array{Float64}, integrand::Array{Float64})

    consts = constants_fundamental();
    EPSILON_0 = consts.EPSILON_0;

    interaction_v[ind,:] .= (1/(4*pi*EPSILON_0)).*trapz((dx,dy,dz,:),
    integrand)

end

function calculate_g2(omega_y::Array{Float64}, t_prime::Array{Float64}, phase::Float64, gamma::Float64, gamma_factor::Float64)

    g = (4.0 .* omega_y.*cos.(omega_y.*t_prime .+ phase).+2.0*gamma.*sin.(omega_y.*t_prime .+ phase)) .*
    exp.(-gamma.*gamma_factor.*t_prime./2.0);

    g[t_prime.<0] .= 0.0

    return g
end

function calculate_green_kernel(discretization::Discretization)
    # %METHOD1 Summary of this method goes here
    # %   Detailed explanation goes here
    return((discretization.x0 .- discretization.XPRIME).^2 .+ 
        (discretization.y0 .- discretization.YPRIME).^2 .+ (discretization.Z .- discretization.ZPRIME).^2).^(-3/2);
end

function threaded_calculation!(n::Int, interaction_v::Array{Float64}, xprime::Array{Float64}, yprime::Array{Float64}, zprime::Array{Float64}, 
    omega_y::Array{Float64}, t_c_subsampled::Array{Float64}, t_r::Array{Float64}, phase::Float64, gamma::Float64, gamma_factor::Float64, YPP_GK_g1::Array{Float64})
    
    Threads.@threads for time_ind in 1:n

        # if(mod(time_ind,10)==0)
        #     print("$(time_ind) out of $(n) \n")
        # end

        integrate_v!(interaction_v, time_ind, xprime,yprime,zprime,
        calculate_g2(omega_y, t_c_subsampled[time_ind] .- t_r, phase, gamma, gamma_factor).*YPP_GK_g1)

    end
end


