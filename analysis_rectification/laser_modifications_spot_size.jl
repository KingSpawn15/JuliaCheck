include("./init.jl");

#original values
mat = IndiumArsenide()

laser_spot_fwhm = 40e-6;

base="analysis_rectification/saved-matrices/analysis_"
base = base*"laser_spot_fwhm"

array_laser_spot_fwhm  = [0.5 * laser_spot_fwhm, 1 * laser_spot_fwhm, 1.5 * laser_spot_fwhm]

function build_matrices_laser(array_laser_spot_fwhm ::Array{Float64, 1}, base::String, las::Laser, 
    dis_sp::Discretization, mat::Material, numericalp::NumericalParameters, elec::Electron; startind::Int64 = 0)

    for (ind, laser_spot_fwhm ) in enumerate(array_laser_spot_fwhm ) 
        
        mod_laser.set_laser!(;laser=las,
        laser_pulse_time_fwhm = 50 * 1e-15,
        laser_spot_fwhm  = laser_spot_fwhm ,
        pulse_energy_gain_factor = 0.014,
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

build_matrices_laser(array_laser_spot_fwhm[1:1],base,las,dis_sp,mat,numericalp,elec;startind=0)
build_matrices_laser(array_laser_spot_fwhm[2:3],base,las,dis_sp,mat,numericalp,elec;startind=1)


