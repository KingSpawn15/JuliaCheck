include("./init.jl");

#original values
mat = IndiumArsenide()

# laser_pulse_time_fwhm = 650e-15;

base="analysis_rectification/saved-matrices/analysis_"
base = base*"laser_pulse_time_fwhm_single"

array_laser_pulse_time_fwhm  = [250.0, 650.0, 50.0]

function build_matrices_laser(array_laser_pulse_time_fwhm ::Array{Float64, 1}, base::String, las::Laser, 
    dis_sp::Discretization, mat::Material, numericalp::NumericalParameters, elec::Electron; startind::Int64 = 0)

    for (ind, laser_pulse_time_fwhm ) in enumerate(array_laser_pulse_time_fwhm ) 
        
        mod_laser.set_laser!(;laser=las,
        laser_pulse_time_fwhm  = laser_pulse_time_fwhm * 1e-15 ,
        theta_pol = 90.0*pi/180);
        
        println(las.laser_spot_sigma)
        @time _, _, interaction_v_pd=mod_cdem.interaction_potential_photodember(dis_sp, mat, las, numericalp)
        @time _, _, interaction_v_or=mod_cdem.interaction_potential_rectification(dis_sp, mat, las , elec, numericalp)

        jldopen(base*".jld2", "a+"; compress = true) do f
            f["photodember_"*string(ind + startind)] = interaction_v_pd
            f["rectification_"*string(ind  + startind)] = interaction_v_or
        end



    end

end

build_matrices_laser(array_laser_pulse_time_fwhm[1:1],base,las,dis_sp,mat,numericalp,elec;startind=0)
build_matrices_laser(array_laser_pulse_time_fwhm[2:3],base,las,dis_sp,mat,numericalp,elec;startind=1)



