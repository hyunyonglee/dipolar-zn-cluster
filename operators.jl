using LinearAlgebra

# Defining Qudit Z operator
function Zn(N::Int)
    op = zeros(ComplexF64, N, N)
    for i in 1:N
        op[i, i] = exp(2 * pi * im * (i - 1) / N)
    end
    return op
end

# Defining Qudit X operator
function Xn(N::Int)
    op = zeros(N, N)
    for i in 1:N-1
        op[i, i+1] = 1
    end
    op[N, 1] = 1
    return op
end

# Defining Qudit Y operator
function Yn(N::Int)
    op = adjoint(Zn(N)) * Xn(N) * adjoint(Zn(N))
    return op
end

