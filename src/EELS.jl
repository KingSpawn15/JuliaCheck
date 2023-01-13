function calculate_ft(discretization::Discretization , interact_v::AbstractArray{Float64}, electron::Electron)
            
    t = discretization.t;
    omega = discretization.omega;
    z = discretization.z;
    deltaz = discretization.deltaz;
    
    
    # %Time domain Fourier transform
    interact_v_fft = fftshift(fft2(interact_v,length(t)),(2,));
    interact_v_fft = interact_v_fft./maximum(omega);

    # % beta(omega)
    Omega = [x for zz in vec(z), x in vec(omega)]
    Z = [zz for zz in vec(z), x in vec(omega)]
    interact_v_z = interact_v_fft.*exp.(-1.0im*Omega.*Z/elec.electron_velocity);
    beta = sum(interact_v_z, dims = 1).*deltaz
    
    consts = constants_fundamental();
    HBAR = consts.HBAR;
    Q_E = consts.Q_E;

    multiplication_factor = (-1/(1.0im*HBAR*electron.electron_velocity))*Q_E;

    f_t_exp = [multiplication_factor.*trapz(omega,2*real(exp.(1.0im.*omega.*t[time_ind]).*beta)) for time_ind in 1:length(t)];
    f_t_exp = reduce(hcat,f_t_exp);
    f_t_exp_shift = circshift(f_t_exp,(0,ceil(length(t)/2)));
    f_t = exp.(-f_t_exp_shift);
    
    return f_t
end

function calculate_ft_fast(discretization::Discretization , interact_v::AbstractArray{Float64}, electron::Electron)
            
    t = discretization.t;
    omega = discretization.omega;
    z = discretization.z;
    deltaz = discretization.deltaz;
    
    # beta = sum(fftshift(fft2(interact_v,length(t)),(2,))./maximum(omega).*[exp(-1.0im*omg*zz/elec.electron_velocity) for zz in vec(z), omg in vec(omega)], dims = 1).*deltaz
    beta = ft_beta(interact_v, t, omega, z, deltaz)

    consts = constants_fundamental();
    HBAR = consts.HBAR;
    Q_E = consts.Q_E;

    # return exp.(-circshift(reduce(hcat,
    # (-1/(1.0im*HBAR*electron.electron_velocity))*Q_E.*ThreadsX.map(tt-> trapz(omega,2.0*real.(exp.(1.0im.*omega.*tt).*beta)), t)),(0,ceil(length(t)/2))));

    return exp.(-circshift(reduce(hcat,
    (-1/(1.0im*HBAR*electron.electron_velocity))*Q_E.*ThreadsX.map(tt-> ft_parameter(omega, tt, beta), t)),(0,ceil(length(t)/2))));

end

function ft_beta(interact_v::Array{Float64}, t::Array{Float64}, 
    omega::Array{Float64}, z::Array{Float64}, deltaz::Float64)
    return sum(fftshift(fft2(interact_v,length(t)),(2,))./maximum(omega).*[exp(-1.0im*omg*zz/elec.electron_velocity) for zz in vec(z), omg in vec(omega)], dims = 1).*deltaz
end

function ft_parameter(omega::Array{Float64}, tt::Float64, beta::Array{ComplexF64})
    trapz(omega,2.0*real.(exp.(1.0im.*omega.*tt).*beta))
end 


function incoherent_convolution_fast(psi::Array{Float64}, w::AbstractArray{Float64}, t_w::Array{Float64}, e_w::Array{Float64} ,
    w_cut_off_factor::Float64 = 0.01)

    psi_sum = zeros(size(psi));

    w_cutOff = w_cut_off_factor*maximum(w[:]);


    # % Example dont use such code. It becomes hard to modify
    # % w_cutOff = 0.01*max(w(:));
    l_t_w = length(t_w)
    l_e_w = length(e_w)

    Threads.@threads for k = 1:l_t_w * l_e_w
        
        incoherent_circ!(psi_sum, w, w_cutOff, psi, t_w, e_w, k, l_t_w, l_e_w)
    
    end

    # psi = psi_sum;
    # psi_incoherent = psi./trapz((:,e_w),psi);
    # psi_incoherent = psi_incoherent./maximum(psi_incoherent[:]);

    psi_sum./=trapz((:,e_w),psi_sum)
    return psi_sum./=maximum(psi[:])

    return psi_incoherent

end

function incoherent_circ!(psi_sum::Array{Float64}, w::Array{Float64}, w_cutOff::Float64,
    psi::Array{Float64}, t_w::Array{Float64}, 
    e_w::Array{Float64}, k::Int, l_t_w::Int, l_e_w::Int)

    t_ind,e_ind = index_calculator_2d(k, l_t_w, l_e_w)

    # return w[k]<w_cutOff ? psi_sum : psi_sum.+=w[k].*circshift(circshift(psi,(-ceil(length(t_w)/2)+t_ind,0)),(0,-ceil(length(e_w)/2)+e_ind))
    return w[k]<w_cutOff ? psi_sum : psi_sum.+=w[k].*circshift(psi,(-ceil(length(t_w)/2)+t_ind,-ceil(length(e_w)/2)+e_ind))

end

function incoherent_convolution(psi::Array{Float64}, w::AbstractArray{Float64}, t_w::Array{Float64}, e_w::Array{Float64} ,
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
    psi_incoherent = psi_incoherent./maximum(psi_incoherent[:]);
    return psi_incoherent

end

function psi_subsampled(sub_sample_factor, psi, e_w)
            
    psi_sub = psi[:,1:sub_sample_factor:end];
    psi_sub = psi_sub./trapz((:,e_w),psi_sub);
    return psi_sub

end

function calculate_psi_coherent(discretization::Discretization,electron::Electron, f_t::Matrix{ComplexF64})
            
    deltat = discretization.deltat;
    t = discretization.t;
    
    # delta_t_rep = repmat(deltat',1,length(t));
    # t_rep = repmat(t,length(deltat),1);
    
    delta_t_rep = repeat(deltat',outer = (1,length(t)));
    t_rep = repeat(t,outer=(length(deltat),1));
    
    FT = repeat(f_t,outer = (length(deltat),1));
    psi_coherent = FT.*exp.(-(t_rep-delta_t_rep).^2/(2*elec.electron_time_coherent_sigma^2));
    
    psi_coherent = fftshift(fft2(psi_coherent,length(t)),(2,));
    
    psi_coherent = (abs.(psi_coherent)).^2;
    psi_coherent = psi_coherent./trapz((:,dis.energy),psi_coherent);
    
    return psi_coherent
    
end

function fft2(a::Union{AbstractArray{Float64},Matrix{ComplexF64}}, n::Int64=size(a,2))
    if n <= size(a,2)
        return A = fft(a[:,1:n],(2,))
    else
        A = zeros(size(a,1),n)
        A[:,1:size(a,2)] .= a
        return fft(A,(2,))
    end
        
end
