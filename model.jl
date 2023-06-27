using ITensors

include("operators.jl")

# Defining Qudit Z operator
function ITensors.op(::OpName"Z", ::SiteType"Qudit", N::Int)
    return Zn(N)
end

# Defining Qudit Z^dagger operator
function ITensors.op(::OpName"Z_dag", ::SiteType"Qudit", N::Int)
    return adjoint(Zn(N))
end

# Defining Qudit X operator
function ITensors.op(::OpName"X", ::SiteType"Qudit", N::Int)
    return Xn(N)
end

# Defining Qudit X_dag operator
function ITensors.op(::OpName"X_dag", ::SiteType"Qudit", N::Int)
    return adjoint(Xn(N))
end

# Defining Qudit Y operator
function ITensors.op(::OpName"Y", ::SiteType"Qudit", N::Int)
    return Yn(N)
end

# Defining Qudit Y_dag operator
function ITensors.op(::OpName"Y_dag", ::SiteType"Qudit", N::Int)
    return adjoint(Yn(N))
end

# Defining Dipolar ZN SPT Hamiltonian
function DIPOLAR_ZN_SPT(N, L, g)

    # Define sites
    sites = siteinds("Qudit", L; dim=N)

    # Defining MPO
    os = OpSum()

    # On-site terms
    for s = 1:L
        os .+= ((-1)^s) * g, "X", s
        os .+= ((-1)^s) * g, "X_dag", s
    end

    # Hopping
    for s = 1:(L-2)
        os .+= -1, "Z", s, "Y", s + 1, "Z", s + 2
        os .+= -1, "Z_dag", s, "Y_dag", s + 1, "Z_dag", s + 2
    end

    return MPO(os, sites; splitblocks=true), sites
end
