using Base.Filesystem
using ITensors

include("model.jl")
include("run_model.jl")


# Creating directory function
function create_directory(dir_path)
    if !isdir(dir_path)
        mkpath(dir_path)
    end
end

# main function
let

    # Defining params
    N = 3
    L = 100
    chi = 100
    nsweeps = 50

    # Defining directory path
    dir_path = pwd() * "/data/Z$(N)_L$(L)_chi_$(chi)"
    create_directory(dir_path * "/mps")
    create_directory(dir_path * "/logs")
    create_directory(dir_path * "/observables")

    # Defining g scan params
    g_min = 0.0
    g_max = 2.0
    g_step = 0.1

    # Defining params dictionary
    model_params = Dict{String,Any}(
        "N" => N,
        "L" => L
    )

    # Defining dmrg params dictionary
    dmrg_params = Dict{String,Any}(
        "nsweeps" => nsweeps,
        "maxdim" => min.([100, 200, 400, 800, 2000, 3000, chi], chi),
        "noise" => [1E-6, 1E-7, 1E-8, 0.0],
        "cutoff" => 1E-8 # 1E-9
    )

    # Loop over g
    for g = g_min:g_step:g_max

        # Define Hamiltonian
        H, sites = DIPOLAR_ZN_SPT(N, L, g)

        # Run DMRG
        psi, e, dmrg_observer = run_dmrg(H, dmrg_params, sites)

        # Calculate observables
        Xs = expect(psi, "X")
        Zs = expect(psi, "Z")

        # Write observables to the file
        writing_observable(Xs, dir_path * "/observables/Xs_g_$(g).txt")
        writing_observable(Zs, dir_path * "/observables/Zs_g_$(g).txt")

        # Calculate SvN
        b = Int(L/2)-1
        orthogonalize!(psi, b)
        U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
        SvN = 0.0
        for n = 1:dim(S, 1)
            p = S[n, n]^2
            SvN -= p * log(p)
        end
        @show SvN

        # Write mps
        writing_mps(psi, dir_path * "/mps/psi_g_$(g).h5")

        # Write logs
        model_params["g"] = g
        writing_logs(model_params, dmrg_params, dmrg_observer, dir_path * "/logs/log_g_$(g).txt")

        # Write observables to the file
        f1 = open(dir_path * "/observables.txt", "a+")
        @printf(f1, "%.3f %.8f %.8f %.8f %.8f\n", g, e, sum(abs.(Xs)) / L, sum(abs.(Zs)) / L, SvN)
        close(f1)

    end

end