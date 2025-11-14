using Distributed
using Plots
using ProgressMeter

# Number of workers 
nworkers = parse(Int, get(ENV, "NWORKERS", "4"))
addprocs(nworkers)

# Load worker code on ALL workers
@everywhere include("worker_TC-Dicke.jl")

# -----------------------------------------------
# Parameters for the scan: 
# Read from the Makefile; 
# default values included in case reading fails.
# -----------------------------------------------
NSpin  = parse(Int, get(ENV, "N_SPIN", "10"))
NBoson = parse(Int, get(ENV, "N_BOSON", "10"))
ww0    = parse(Float64, get(ENV, "W0", "1.0"))
hh     = parse(Float64, get(ENV, "HFIELD", "1.0"))

e_min  = parse(Float64, get(ENV, "E_MIN", "0.05"))
e_step = parse(Float64, get(ENV, "E_STEP", "0.05"))
e_max  = parse(Float64, get(ENV, "E_MAX", "0.05"))

e_values = e_min:e_step:e_max

data_dir = joinpath(@__DIR__, "data")
mkpath(data_dir)
data_file = joinpath(data_dir, "TC-Dicke_w0$(ww0)_h$(hh)_L$(NSpin)_Nb$(NBoson).txt")
@everywhere const DATA_DIR = $data_dir

# Prepare output file
open(data_file, "w") do f
    println(f, "# g-scan of the ground state of TFIM-Dicke model")
    println(f, "# e, g, E/N, <a†a>, <Sx>, <Sy>, <Sz>")
end

# -----------------------------------------------
# Run the distributed calculation
# -----------------------------------------------
results = @showprogress pmap(e_values) do ee
    compute_for_e(ee; NSpin=NSpin, NBoson=NBoson,
                  ww0=ww0, hh=hh)
end

# -----------------------------------------------------------
# Merge worker files into the final data_file (sorted by e, g)
# -----------------------------------------------------------
worker_files = filter(f -> occursin("_worker", f),
                      readdir(joinpath(@__DIR__, "data"), join=true))

println("Merging worker files...")

# Read all lines, parse into tuples (e, g, line)
all_rows = Vector{Tuple{Float64,Float64,String}}()

for wf in worker_files
    open(wf, "r") do fin
        for ln in eachline(fin)
            fields = split(ln)
            if length(fields) >= 2
                e = parse(Float64, fields[1])
                g = parse(Float64, fields[2])
                push!(all_rows, (e, g, ln))
            end
        end
    end
end

# Sort by e, then g
sort!(all_rows, by = x -> (x[1], x[2]))

# Rewrite final file with header + sorted data
open(data_file, "w") do fout
    println(fout, "# g-scan of the ground state of TFIM-Dicke model")
    println(fout, "# e, g, E/N, <a†a>, <Sx>, <Sy>, <Sz>")
    for (_, _, ln) in all_rows
        println(fout, ln)
    end
end

# Delete worker files
for wf in worker_files
    rm(wf; force=true)
end

# -----------------------------------------------------------
# Create density plot of abs(<a†a>) vs e and g
# -----------------------------------------------------------
e_vals = Float64[]
g_vals = Float64[]
n_vals = Float64[]

open(data_file, "r") do fin
    for ln in eachline(fin)
        startswith(ln, "#") && continue
        fields = split(ln)
        push!(e_vals, parse(Float64, fields[1]))
        push!(g_vals, parse(Float64, fields[2]))
        push!(n_vals, abs(parse(Float64, fields[4])))
    end
end


# Create heatmap (convert scattered data to matrix form)
uniq_e = sort(unique(e_vals))
uniq_g = sort(unique(g_vals))

zmat = fill(NaN, length(uniq_g), length(uniq_e))

for (ee, gg, nn) in zip(e_vals, g_vals, n_vals)
    i = findfirst(==(gg), uniq_g)
    j = findfirst(==(ee), uniq_e)
    zmat[i, j] = nn
end

plt = heatmap(uniq_e, uniq_g, zmat,
              xlabel="e", ylabel="g", colorbar_title="|<a†a>|",
              title="Heatmap of |<a†a>| vs e & g")

png_path = replace(data_file, ".txt" => "_density.png")
savefig(plt, png_path)
println("Saved density plot to $png_path")

println("Merged, sorted, and cleaned up worker files.")