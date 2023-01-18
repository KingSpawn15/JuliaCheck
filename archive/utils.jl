module Utils

    export calculate_sigma

    function calculate_sigma(fwhm::Float64)
        return fwhm./(2*sqrt(2*log(2)))
    end

end