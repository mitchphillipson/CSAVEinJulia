using GTAP9data, JuMP, Ipopt, MPSGE, JLD2, DataFrames, CSV

GTAP9data.source(joinpath(@__DIR__, "data"))

data = load("./IO.jld2")    #k = keys(data)
for (k, v) in data
    @eval const $(Symbol(k)) = $v
end

# Declare Vector similar to set declaration in GAMS
set_fe      = [:coa, :gas, :p_c]
set_elec    = [:ely]
set_ne      = setdiff(setdiff(set_i, set_fe), set_elec)
set_tr      = [:wtp, :atp, :otp]

MGE  = MPSGEModel()

@parameters(MGE, begin
    rtfd[i=set_i, g=set_g, r=set_r],    rtfd0[i, g, r], (description = "Firms' domestic tax rates")
    rtfi[i=set_i, g=set_g, r=set_r],    rtfi0[i, g, r], (description = "Firms' import tax rates")
    rtms[i=set_i, r=set_r, s=set_r],    rtms0[i, r, s], (description = "Import tax rates")
    rtxs[i=set_i, r=set_r, s=set_r],    rtxs0[i, r, s], (description = "Export subsidy rates")
    rto[g=set_g, r=set_r],              rto0[g, r],     (description = "Output subsidy rates")              
    rtf[f=set_f, i=set_i, r=set_r],     rtf0[f, i, r],  (description = "Primary factor tax rates")
end)

@sectors(MGE, begin
    Y[set_g, set_r],            (description = "Supply")
    M[set_i, set_r],            (description = "Imports")
    YT[set_i],                  (description = "Transportation services")
    E[set_i, set_r, set_r],     (description = "Subsidy and transport service included exports")
    A[set_i, set_g, set_r],     (description = "Armington good")
end)

@commodities(MGE, begin
    P[set_g, set_r],            (description = "Domestic output price")
    PM[set_i, set_r],           (description = "Import price")
    PT[set_i],                  (description = "Transportation services")
    PF[set_mf, set_r],          (description = "Non-sector-specific primary factor rent")
    PS[set_sf, set_g, set_r],   (description = "Sector-specific primary factor rent")  
    PX[set_i, set_r, set_r],    (description = "Price index for exports (include subsidy and transport service)")
    PA[set_i, set_g, set_r],    (description = "Price index for Armington good")
    PE[set_i, set_r],           (description = "Price index for exports (exclude subsidy and transport service)")
end)

@consumers(MGE, begin
    RA[set_r],                  (description = "Representative agent")
end)

for i ∈ set_i, g ∈ set_g, r ∈ set_r
    @production(MGE, A[i, g, r], [t = 0, s = esubd[i]], begin
        @output(PA[i, g, r],    vafm[i, g, r],  t)
        @input(P[i, r],         vdfm[i, g, r],  s,   taxes = [Tax(RA[r], rtfd[i, g, r])],   reference_price = 1 + rtfd0[i, g, r])
        @input(PM[i, r],        vifm[i, g, r],  s,   taxes = [Tax(RA[r], rtfi[i, g, r])],   reference_price = 1 + rtfi0[i, g, r])  
    end)
end

for g ∈ set_i, r ∈ set_r
    @production(MGE, Y[g, r], [t = etadx[g], s = esub[g], sn => s = esubn[g], sve => sn = esubve[g], sva => sve = esubva[g], sef => sve = esubef[g], sf => sef = esubf[g]], begin
        @output(P[g, r],        vhm[g, r], t, taxes = [Tax(RA[r], rto[g, r])], reference_price = 1-rto0[g, r])
        @output(PE[g, r],        vxm[g, r], t, taxes = [Tax(RA[r], rto[g, r])], reference_price = 1-rto0[g, r])    
        [@input(PA[i, g, r],     vafm[i, g, r], sf) for i ∈ set_fe]...
        [@input(PA[i, g, r],     vafm[i, g, r], sef) for i ∈ set_elec]...
        [@input(PA[i, g, r],     vafm[i, g, r], sn) for i ∈ set_ne]...
        [@input(PS[sf, g, r],   vfm[sf, g, r],  s, taxes = [Tax(RA[r], rtf[sf, g, r])],   reference_price = 1 + rtf0[sf, g, r])   for sf ∈ set_sf]...    
        [@input(PF[mf, r],      vfm[mf, g, r],  sva, taxes = [Tax(RA[r], rtf[mf, g, r])],   reference_price = 1 + rtf0[mf, g, r])   for mf ∈ set_mf]...    
    end)
end

for g ∈ set_cgi, r ∈ set_r
    @production(MGE, Y[g, r], [t = 0, s = esub[g], sn => s = esubn[g], sef => sn = esubef[g], sf => sef = esubf[g]], begin
        @output(P[g, r],        vom[g, r], t, taxes = [Tax(RA[r], rto[g, r])])
        [@input(PA[i, g, r],     vafm[i, g, r], sf) for i ∈ set_fe]...
        [@input(PA[i, g, r],     vafm[i, g, r], sef) for i ∈ set_elec]...
        [@input(PA[i, g, r],     vafm[i, g, r], sn) for i ∈ set_ne]...
    end)
end

for j ∈ set_i
    @production(MGE, YT[j], [t = 0, s = 1], begin
        @output(PT[j],          vtw[j],         t)
        [@input(PE[j, r],       vst[j, r],      s)   for r ∈ set_r]...
    end)
end

for i ∈ set_i, r ∈ set_r
    @production(MGE, M[i, r], [t = 0, s = esubm[i]], begin
        @output(PM[i, r],       vim[i, r],      t)
        [@input(PX[i, s, r], vxmd[i, s, r]*(1 - rtxs0[i, s, r]) + sum(vtwr[j, i, s, r] for j ∈ set_tr), s, taxes = [Tax(RA[r], rtms[i, s, r])], reference_price = pvtwr[i, s, r]) for s ∈ set_r]...
    end)
end

# vxmr = Dict((i, s, r) => vxmd[i, s, r]*(1 - rtxs0[i, s, r]) + sum(vtwr[j, i, s, r] for j ∈ set_tr)

for i ∈ set_i, s ∈ set_r, r ∈ set_r
    @production(MGE, E[i, s, r], [t = 0, s = 0], begin
        [@output(PX[i, s, r], vxmd[i, s, r]*(1 - rtxs0[i, s, r]) + sum(vtwr[j, i, s, r] for j ∈ set_tr), t)]...
        @input(PE[i, s],        vxmd[i, s, r], s,   taxes = [Tax(RA[s], -rtxs[i, s, r])],   reference_price = 1 - rtxs0[i, s, r])
        [@input(PT[j],          vtwr[j, i, s, r], s) for j ∈ set_i]...
    end)
end

for r ∈ set_r 
    @demand(MGE, RA[r], begin
        @final_demand(P[:c, r], vom[:c, r])
        @endowment(P[:c, :USA], vb[r])
        @endowment(P[:g, r], -vom[:g, r])
        @endowment(P[:i, r], -vom[:i, r])
        [@endowment(PF[f, r], evom[f, r]) for f ∈ set_mf]...
        [@endowment(PS[f, j, r], vfm[f, j, r]) for f ∈ set_sf, j ∈ set_i]...
    end)
end

fix(P[:c, :USA], 1)

solve!(MGE, cumulative_iteration_limit = 0)
benchmark = generate_report(MGE)
println(benchmark)

solvetime = zeros(5)
n = length(solvetime)
for t ∈ 1:n
    for i ∈ set_fe, g ∈ set_g
    set_value!(rtfd[i, g, :USA], rtfd0[i, g, :USA]*2*(t-1)/(n-1))
    set_value!(rtfi[i, g, :USA], rtfi0[i, g, :USA]*2*(t-1)/(n-1))
    end
    solvetime[t] = @elapsed solve!(MGE; cumulative_iteration_limit = 1000, convergence_tolerance = 1e-8)
end
df = DataFrame(run = 1:length(solvetime), runtime = solvetime)
println(df)

path = joinpath(@__DIR__, "56x2_5.csv")
CSV.write(path, df)

df = generate_report(MGE)
println(df)

df_pf = filter(row -> startswith(string(row.var), "PF"), df)

df_filtered = df[df.margin .> 0.001, :]
println(df_filtered)
