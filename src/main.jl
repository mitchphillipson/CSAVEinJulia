using JuMP, Ipopt, MPSGE, JLD2, DataFrames, CSV, GTAPdata

#using Pkg
#Pkg.develop(path="E:/work/Projects/2025-01_Julia/CSVtoDIC")
#Pkg.develop(path="E:/work/Projects/2025-01_Julia/GTAPdata")

io(joinpath(@__DIR__, "data"))

data = load("./IO.jld2")    #k = keys(data)
for (k, v) in data
    @eval $(Symbol(k)) = $v
end

model_generation_time = @elapsed include("model.jl")
df1 = DataFrame(run = "MG", runtime = model_generation_time)

solve!(MGE, cumulative_iteration_limit = 0)
benchmark = generate_report(MGE)
#println(benchmark)

solvetime = zeros(5)
n = length(solvetime)
for t ∈ 1:n
    for i ∈ set_fe, g ∈ set_g
    set_value!(rtfd[i, g, :USA], rtfd0[i, g, :USA]*2*(t-1)/(n-1))
    set_value!(rtfi[i, g, :USA], rtfi0[i, g, :USA]*2*(t-1)/(n-1))
    end
    solvetime[t] = @elapsed solve!(MGE; cumulative_iteration_limit = 1000, convergence_tolerance = 1e-8)
end
df2 = DataFrame(run = 1:length(solvetime), runtime = solvetime)

df = vcat(df1, df2)
#println(df)

path = joinpath(@__DIR__, "56x2_5.csv")
CSV.write(path, df)

df_results = generate_report(MGE)
println(df_results)

df_Y = filter(:var => x -> startswith(string(x), "Y["), df_results)
df_M = filter(:var => x -> startswith(string(x), "M["), df_results)
df_YM = vcat(df_Y, df_M)

df_P = filter(:var => x -> startswith(string(x), "P["), df_results)

path = joinpath(@__DIR__, "56x2_YM.csv")
CSV.write(path, df_YM)

#df_pf = filter(row -> startswith(string(row.var), "PF"), df)
#df_filtered = df[df.margin .> 0.001, :]
#println(df_filtered)
