a = dis.xprime;
b = SVector{length(a)}(vec(a));

fn_a(a::AbstractArray{Float64}) = [bi*bi for bi in a, _ in a, _ in a, _ in a] 


c = [bi for bi in a, _ in a, _ in a, _ in a] 
fn_c(a::AbstractArray{Float64}) = a.*a

