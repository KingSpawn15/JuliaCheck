include("./init.jl");
using LaTeXStrings
mat = IndiumArsenide()
@time _, _, interaction_v_pd=mod_cdem.interaction_potential_photodember(dis_sp, mat, las, numericalp);
@time _, _, interaction_v_pd=mod_cdem.interaction_potential_photodember(dis_sp, mat, las, numericalp);
@time _, _, interaction_v_or=mod_cdem.interaction_potential_rectification(dis_sp, mat, las , elec, numericalp)
@time _, _, interaction_v_or=mod_cdem.interaction_potential_rectification(dis_sp, mat, las , elec, numericalp)


# include("./utils.jl")

# w, e_w, t_w = mod_electron.energy_time_grid(elec,numericalp.subsampling_factor,dis_sp.energy, dis_sp.deltat);

# psi_sub_array = zeros(length(t_w),length(e_w));
# psi_incoherent_array =similar(psi_sub_array);
# @time loss_spectrum!(psi_sub_array,psi_incoherent_array,interaction_v_pd, dis_sp, elec, w, t_w, e_w);

# using SharedArrays, Distributed

# function Test_distributed()

#     Z = SharedArray{Float64,2}(1000,1000)

#     @sync @distributed for p = 1:100000
#         parameter = rand()
#         Z[p] = parameter
#     end

#     return Z
# end

# function Test_serial()

#     Z = SharedArray{Float64,2}(1000,1000)
#     # Z =zeros(1000,1000)
#     for p = 1:100000
#         parameter = rand()
#         Z[p] = parameter
#     end

#     return Z
# end

# function Test_threads()

#     Z = SharedArray{Float64,2}(1000,1000)
#     # Z =zeros(1000,1000)
#     Threads.@threads for p = 1:100000
#         parameter = rand()
#         Z[p] = parameter
#     end

#     return Z
# end


# @time Test_distributed();
# @time Test_serial();
# @time Test_threads();

# a = rand(101,5);

# function trap1(a::Array{Float,2}, dx::Float64)
#     return dx.*(sum(a, dims = 1) .+ 0.5.*(a[1,:,:,:] .+  0.5.*a[end,:,:,:]));
# end

# dx = [0:0.1:10;]
# dy = [0:0.2:10;]
# a = rand(length(dx),length(dy));
# trapz((dx,:),a)