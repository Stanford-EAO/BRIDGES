# rm(joinpath(DEPOT_PATH[1], "registries", "General"); force=true, recursive=true)
# import Pkg
# Pkg.Registry.add("General")
# Pkg.status()

# Pkg.add("Tables")
# Pkg.add("DataFrames")
# Pkg.add("CSV")
# Pkg.add("Clustering")
# Pkg.add("JuMP")
# Pkg.add("Distances")
# Pkg.add("Gurobi")
# Pkg.add("Tables")

# using CSV, Dates, Tables
# # creates output folder automatically
# function mk_out_dir()
#     timestamp = Dates.format(now(), "YYYYmmdd-HHMMSS")
#     dir_name = joinpath(@__DIR__, "output", "$timestamp")
#     @assert !ispath(dir_name) "File name already taken!"
#     mkpath(dir_name)
#     return dir_name
# end

# top_dir = mk_out_dir() * "/"
# t = Table(name = ["Alice", "Bob", "Charlie"], age = [25, 42, 37])

# CSV.write(top_dir * "test_file.csv", t,  writeheader=true)

################################################################################
versioninfo()
# import Pkg
# Pkg.add("DataFrames")
# Pkg.add("CSV")
# Pkg.add("Clustering")
# Pkg.add("JuMP")
# Pkg.add("Distances")
# Pkg.add("Gurobi")
# Pkg.add(Pkg.PackageSpec(name = "Gurobi", version = v"0.8"))
# Pkg.add("Tables")

print("Run Success!")