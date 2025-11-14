# ----------------------------------------------------------
# Worker-side code (loaded by @everywhere include("worker.jl"))
# ----------------------------------------------------------

using ITensors
using ITensorMPS
using Printf
using ProgressMeter

# ----------------------------------------------------------
# Utility functions
# ----------------------------------------------------------

function semi_random(opts)
    # returns opts[1] or opts[2] randomly
    rand() < 0.5 ? opts[1] : opts[2]
end

function chop_small(x, tol=1e-12)
    if x isa Complex
        r = abs(real(x)) < tol ? 0.0 : real(x)
        i = abs(imag(x)) < tol ? 0.0 : imag(x)
        return complex(r, i)
    else
        return abs(x) < tol ? 0.0 : x
    end
end

# ----------------------------------------------------------
# Hamiltonian
# ----------------------------------------------------------

function TFIMDickeMPO(sites, w0::Float64, h::Float64, e::Float64, g::Float64)
    L = length(sites)
    hmpo = OpSum()

    # Boson energy
    hmpo += w0, "N", 1

    gfac = g / sqrt(L - 1)

    for j in 2:L
        hmpo += -h, "Sz", j
        hmpo += e * (-1)^j, "Sx", j
        hmpo += gfac, "Adag", 1, "Sx", j
        hmpo += gfac, "A",    1, "Sx", j
    end

    return MPO(hmpo, sites)
end

# ----------------------------------------------------------
# Worker task: compute g-scan at fixed e
# ----------------------------------------------------------

function compute_for_e(ee; NSpin, NBoson, ww0, hh)
    worker_file = joinpath(DATA_DIR, "TC_Dicke_worker$(myid()).txt")

    # updates the status of the worker; comment out if its spammy
    @info "Worker $(myid()): e=$ee"

    # Sites
    sites = [siteind("Boson", dim=NBoson+1)]
    append!(sites, [siteind("S=1/2") for _ in 1:NSpin])

    # g values
    g_min  = parse(Float64, get(ENV, "G_MIN",  "0.5"))
    g_step = parse(Float64, get(ENV, "G_STEP", "0.05"))
    g_max  = parse(Float64, get(ENV, "G_MAX",  "0.5"))
    gvals  = g_min:g_step:g_max

    for (i, gg) in enumerate(gvals)

        # initial state vector
        st = Vector{String}(undef, NSpin+1)
        st[1] = "0"
        for j in 2:(NSpin+1)
            st[j] = semi_random(["Up","Dn"])
        end

        psi0 = productMPS(sites, st)

        # Hamiltonian
        H = TFIMDickeMPO(sites, ww0, hh, ee, gg)

        # DMRG
        energy, psi = dmrg(H, psi0;
            nsweeps=20,
            maxdim=[10,20,40,80,160,200],
            cutoff=1E-12,
            noise=[1E-4,1E-6,1E-8,0.0],
            outputlevel=0
        )

        # observables
        avgN = expect(psi, "N"; sites=1)

        Sx = Sy = Sz = 0.0
        for j in 2:(NSpin+1)
            Sx += 0.5 * expect(psi, "S+ + S-"; sites=j)
            Sy += -0.5im * expect(psi, "S+ - S-"; sites=j)
            Sz += expect(psi, "Sz"; sites=j)
        end

        Sx /= NSpin
        Sy /= NSpin
        Sz /= NSpin

        # Write to file 
        open(worker_file, "a") do f
            println(f, "$ee \t $gg \t $(energy/NSpin) \t $(chop_small(avgN)) \t $(chop_small(Sx)) \t $(chop_small(Sy)) \t $(chop_small(Sz))")
        end 

    end # g-loop
end