function index_calculator_2d(k::Int, nrows::Int, ncols::Int)
    
    col = Int(floor(k / nrows)) + 1
    row = mod(k,nrows)

    if row == 0
        row = nrows
        col = col - 1
    end

    return row, col

end

function index_calculator_3d(k::Int, d1::Int, d2::Int, d3::Int)

    ind3 = Int(floor(k / d1 / d2)) + 1
    ind12 = mod(k , d1 * d2)

    if ind12 == 0
        ind3 = ind3 - 1
        ind12 = d1 * d2
    end 

    ind1, ind2 = index_calculator_2d(ind12, d1, d2)

    return ind1, ind2, ind3

end

function index_calculator_4d(k::Int, d1::Int, d2::Int, d3::Int, d4::Int)

    ind4 = Int(floor(k / d1 / d2 / d3)) + 1
    ind123 = mod(k , d1 * d2 * d3)

    if ind123 == 0
        ind4 = ind4 - 1
        ind123 = d1 * d2 * d3
    end 

    ind1, ind2 , ind3 = index_calculator_3d(ind123, d1, d2, d3)

    return ind1, ind2, ind3, ind4
    
end
