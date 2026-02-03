using Revise
using Dates, CSV, DataFrames, DimensionalData, Tables

using LagrangianERA5

########## Helpers ##########

function build_itps(interp_dict)
    itps = Dict{DateTime, Tuple{Any,Any,Any}}()
    for t in keys(interp_dict[:u])
        itps[t] = (
            interp_dict[:u][t],
            interp_dict[:v][t],
            interp_dict[:w][t]
        )
    end
    times = sort(collect(keys(itps)))
    return itps, times
end

function build_ps_fun(ps_itps)
    times_vec = collect(keys(ps_itps))
    return (lon, lat, t) -> begin
        nearest_time = findmin(abs.(Dates.value.(times_vec .- t)))[2]
        t_sel = times_vec[nearest_time]
        return ps_itps[t_sel](lon, lat)
    end
end

function count_crossings!(traj, key, n, i)
    if Bool(traj.extras[key][i]) && !Bool(traj.extras[key][i+1])
        n += 1
    end
end

function count_time!(traj, key, n, i)
    if Bool(traj.extras[key][i])
        n += 1
    end
end     

function integrate_bidirectional(
    traj_id::Int,
    x0::Tuple{Float64,Float64,Float64},
    times::Vector{DateTime},
    itps::Dict{DateTime, Tuple{Any,Any,Any}},
    start_idx::Int;
    extra_itps = nothing,
    vars = String[],
    ps_fun = nothing)

    # slice times around pivot
    times_back = times[1:start_idx]
    times_forw = times[start_idx:end]

    traj_back = LagrangianERA5.Integrator.integrate_trajectory(traj_id, x0, times_back, itps; direction=:backward)
    traj_forw = LagrangianERA5.Integrator.integrate_trajectory(traj_id, x0, times_forw, itps; direction=:forward)

    # reverse backward segment to chronological order
    rev_times = reverse(traj_back.time)
    rev_lon = reverse(traj_back.lon)
    rev_lat = reverse(traj_back.lat)
    rev_p = reverse(traj_back.p)

    if length(rev_times) > 1
        t_back_cut = rev_times[1:end-1]
        lon_back_cut = rev_lon[1:end-1]
        lat_back_cut = rev_lat[1:end-1]
        p_back_cut = rev_p[1:end-1]
    else
        t_back_cut = DateTime[]
        lon_back_cut = Float64[]
        lat_back_cut = Float64[]
        p_back_cut = Float64[]
    end

    t_forw = traj_forw.time
    lon_forw = traj_forw.lon
    lat_forw = traj_forw.lat
    p_forw = traj_forw.p

    times_comb = vcat(t_back_cut, t_forw)
    lons_comb  = vcat(lon_back_cut, lon_forw)
    lats_comb  = vcat(lat_back_cut, lat_forw)
    ps_comb    = vcat(p_back_cut, p_forw)

    pivot_idx = length(t_back_cut) + 1
    if pivot_idx >= 1 && pivot_idx <= length(times_comb)
        st_time = times_comb[pivot_idx]
        st_lon  = lons_comb[pivot_idx]
        st_lat  = lats_comb[pivot_idx]
        st_p    = ps_comb[pivot_idx]
    else
        st_time = length(times_comb) >= 1 ? times_comb[1] : DateTime(0)
        st_lon  = length(lons_comb)  >= 1 ? lons_comb[1]  : NaN
        st_lat  = length(lats_comb)  >= 1 ? lats_comb[1]  : NaN
        st_p    = length(ps_comb)    >= 1 ? ps_comb[1]    : NaN
    end

    traj = LagrangianERA5.TrajectoryModule.Trajectory(traj_id, st_time, st_lon, st_lat, st_p,
                                                        times_comb, lons_comb, lats_comb, ps_comb,
                                                        Float64[], Dict{String,Any}())


    if ps_fun !== nothing
        LagrangianERA5.TrajectoryModule.compute_altitude!(traj, ps_fun)
    end

    if extra_itps !== nothing
        
        LagrangianERA5.TrajectoryModule.interpolate_along_trajectory!(traj, extra_itps, vars)
        LagrangianERA5.TrajectoryModule.compute_rhi!(traj)
        LagrangianERA5.TrajectoryModule.compute_potential_temperature!(traj)
        LagrangianERA5.TrajectoryModule.compute_time_derivative!(traj, "t", "dTdt")
        LagrangianERA5.TrajectoryModule.compute_time_derivative!(traj, "q", "dqdt")
        LagrangianERA5.TrajectoryModule.compute_curvature!(traj)
        LagrangianERA5.TrajectoryModule.compute_cold_point!(traj)

        if haskey(traj.extras, "mask")
            mask = traj.extras["mask"]

            pidx = pivot_idx
            start_rhi_flag = (pidx >= 1 && pidx <= length(mask)) ? Bool(mask[pidx]) : false
            traj.extras["start_rhi"] = start_rhi_flag

            # Calculate quantities for the sample point being ice-supersaturated
            if start_rhi_flag

                j = pidx

                while j > 1 && Bool(mask[j-1])
                    j -= 1
                end

                k = pidx

                while k < length(mask) && Bool(mask[k+1])
                    k += 1
                end

                # j is the start index of the issr (based on 90%)
                # k is the end index of the issr (based on 90%)

                LagrangianERA5.TrajectoryModule.compute_cold_point2!(traj, j)

                rhi_start_time = times_comb[j]
                rhi_end_time = times_comb[k]
                # duration in hours between block start and block end
                dt_block = rhi_end_time - rhi_start_time
                duration_hours = Dates.value(dt_block) / 3_600_000.0

                traj.extras["rhi_start_time"] = rhi_start_time
                traj.extras["rhi_end_time"] = rhi_end_time
                traj.extras["rhi_duration_hours"] = duration_hours

                n_drop_100 = 0
                n_drop_110 = 0

                n_above_100 = 0
                n_above_110 = 0

                sum_rhi = 0
                sum_w = 0
                sum_dqdt = 0
                sum_curv = 0
                sum_d = 0
                count = 0

                for i in j:(k-1)

                    # Count the number of crossings (a measure of intermittency)
                    count_crossings!(traj, "mask110", n_drop_110, i)
                    count_crossings!(traj, "mask100", n_drop_100, i)

                    # Time above thresholds
                    count_time!(traj, "mask110", n_above_110, i)
                    count_time!(traj, "mask100", n_above_100, i)

                    # Calculate the averages over the ice supersaturated region
                    sum_rhi += traj.extras["rhi"][i]
                    sum_w += traj.extras["w"][i]
                    sum_dqdt += traj.extras["dqdt"][i]
                    sum_curv += traj.extras["curv"][i]
                    sum_d += traj.extras["d"][i]
                    count += 1
                end
                
                traj.extras["n_rhi100_drops"] = n_drop_100
                traj.extras["n_rhi110_drops"] = n_drop_110

                traj.extras["t_rhi100_hours"] = 0.1 * n_above_100
                traj.extras["t_rhi110_hours"] = 0.1 * n_above_110

                if k > j # i.e. most issrs

                    LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "w", j, 30)
                    LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "d", j, 30)
                    LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "dTdt", j, 30)
                    LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "dqdt", j, 30)
                    LagrangianERA5.TrajectoryModule.compute_deltaT!(traj, j, 30)

                    traj.extras["max_rhi_issr"] = maximum(traj.extras["rhi"][j:(k-1)])
                    traj.extras["max_w_issr"]   = maximum(traj.extras["w"][j:(k-1)])

                    traj.extras["avg_rhi_issr"]         = sum_rhi / count
                    traj.extras["avg_w_issr"]           = sum_w / count
                    traj.extras["avg_dqdt_issr"]        = sum_dqdt / count
                    traj.extras["avg_curvature_issr"]   = sum_curv / count
                    traj.extras["avg_divergence_issr"]  = sum_d / count

                    traj.extras["delta_q_issr"] = traj.extras["q"][k] - traj.extras["q"][j]

                    traj.extras["dt_cold_to_event"]         = (j - traj.extras["idx_cold_point"]) * 0.1
                    traj.extras["dt_cold_to_event_pre"]     = (j - traj.extras["idx_cold_point_before"]) * 0.1

                    traj.extras["init_w_issr"]         = traj.extras["w"][j]

                elseif k==j # i.e. issrs or length 1

                    LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "w", j, 30)
                    LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "d", j, 30)
                    LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "dTdt", j, 30)
                    LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "dqdt", j, 30)
                    #event integrated ascent
                    LagrangianERA5.TrajectoryModule.compute_deltaT!(traj, j, 30)

                    traj.extras["max_rhi_issr"]       = traj.extras["rhi"][j]
                    traj.extras["max_w_issr"]         = traj.extras["w"][j]
                    traj.extras["avg_rhi_issr"]       = traj.extras["rhi"][j]
                    traj.extras["avg_w_issr"]         = traj.extras["w"][j]
                    traj.extras["avg_dqdt_issr"]      = traj.extras["dqdt"][j]
                    traj.extras["avg_curvature_issr"] = traj.extras["curv"][j]
                    traj.extras["avg_divergence_issr"]= traj.extras["d"][j]

                    traj.extras["delta_q_issr"] = 0.0

                    traj.extras["dt_cold_to_event"]     = (j - traj.extras["idx_cold_point"]) * 0.1
                    traj.extras["dt_cold_to_event_pre"] = (j - traj.extras["idx_cold_point_before"]) * 0.1

                    traj.extras["init_w_issr"]          = traj.extras["w"][j]
                    
                end

            # Store the data for the case that the sample point was not ice supersaturated
            else

                LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "w", pidx, 30)
                LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "d", pidx, 30)
                LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "dTdt", pidx, 30)
                LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "dqdt", pidx, 30)
                LagrangianERA5.TrajectoryModule.compute_deltaT!(traj, pidx, 30)

                traj.extras["rhi_start_time"]       = missing
                traj.extras["rhi_end_time"]         = missing
                traj.extras["rhi_duration_hours"]   = missing

                traj.extras["max_rhi_issr"]         = missing
                traj.extras["max_w_issr"]           = missing
                traj.extras["avg_rhi_issr"]         = missing
                traj.extras["avg_w_issr"]           = missing
                traj.extras["avg_dqdt_issr"]        = missing
                traj.extras["avg_curvature_issr"]   = missing
                traj.extras["avg_divergence_issr"]  = missing

                traj.extras["delta_q_issr"]         = missing

                traj.extras["dt_cold_to_event"]     = (pidx - traj.extras["idx_cold_point"]) * 0.1
                traj.extras["dt_cold_to_event_pre"] = missing

            end
        end
    end

    return traj
end

function trajectory_to_df(traj::LagrangianERA5.TrajectoryModule.Trajectory)
    N = length(traj.time)

    z = length(traj.z) == N ? traj.z : fill(NaN, N)

    # core trajectory information
    df = DataFrame(
        traj_id = fill(traj.id, N),
        time    = traj.time,
        lon     = traj.lon,
        lat     = traj.lat,
        p       = traj.p,
        z       = z
    )

    # extras (vectors only)
    for (var, vals) in traj.extras
        if isa(vals, AbstractVector)
            df[!, var] = length(vals) == N ? vals : fill(missing, N)
        end
    end

    return df
end

function trajectory_to_metadata(traj::LagrangianERA5.TrajectoryModule.Trajectory)

    row = Dict{Symbol,Any}()

    # core metadata
    row[:id] = traj.id
    row[:start_time] = traj.start_time
    row[:start_lon] = traj.start_lon
    row[:start_lat] = traj.start_lat
    row[:start_p] = traj.start_p

    # extras (scalars only)
    for (k, v) in traj.extras
        if !(v isa AbstractVector)
            row[Symbol(k)] = v
        end
    end

    return DataFrame(row)
end

println("Lagrangian IceSupersaturated Tool (LIST)")
println("----------------------------------------\n")

########## CONFIG ##########

start_time = DateTime(2019,1,3)
end_time   = DateTime(2019,1,8)
n_time     = 31

lon_min   = -50.0
lon_max   = -10.0
lon_steps = 20

lat_min   = 35.0
lat_max   = 65.0
lat_steps = 20

pressure_level = 25000.0

vars = ["u", "v", "w", "t", "q", "pv", "d"]

###### GENERATE GRIDS ######

ics = LagrangianERA5.Experiments.generate_grid_ic(
    lon_min, lon_max, lon_steps, lat_min, lat_max, lat_steps, pressure_level;
    time_start = start_time, time_end = end_time, n_time = n_time
)

dates = collect(start_time-Hour(48):Day(1):end_time+Hour(48))
datadir = "/Users/lewisclark/Documents/python_work/lagrangian/data"

initial_conditions = [(i, ics[i]) for i in 1:length(ics)]

########## LOAD DATA ##########

print("Reading data...              ")

t = @elapsed begin 
stack = LagrangianERA5.GeoData.loadera5range(dates; datadir=datadir, variables=vars)
ps_da = LagrangianERA5.GeoData.load_surface_range(dates; datadir=datadir)
end

println("(time taken = $(round(t, digits=2)) s)")

##### BUILD INTERPOLATORS ######

print("Building interpolators...    ")

t = @elapsed begin
interp_dict = LagrangianERA5.Interpolation.build_interpolators(stack)
itps, times = build_itps(interp_dict)

ps_itps = LagrangianERA5.Interpolation.build_ps_interpolators(ps_da)
ps_fun = build_ps_fun(ps_itps)

extra_itps = LagrangianERA5.Interpolation.build_interpolators(stack)
end

println("(time taken = $(round(t, digits=2)) s)")
# snap requested start_time to available sample
start_idx = findmin(abs.(Dates.value.(times .- start_time)))[2]
start_time = times[start_idx]

########## INTEGRATE ##########

print("Integrating Trajectories...  ")

t = @elapsed begin
trajectories = Vector{LagrangianERA5.TrajectoryModule.Trajectory}()

for (traj_id, x0) in initial_conditions
    # x0 may be (lon,lat,p) or (lon,lat,p,start_time)
    if length(x0) == 4
        lon0, lat0, p0, ic_time = x0
        x0_pos = (lon0, lat0, p0)
        local_start_idx = findmin(abs.(Dates.value.(times .- ic_time)))[2]
    else
        x0_pos = x0
        local_start_idx = start_idx
    end

    traj = integrate_bidirectional(traj_id, x0_pos, times, itps, local_start_idx;
                                   extra_itps=extra_itps, vars=vars, ps_fun=ps_fun)
    push!(trajectories, traj)
end

end
println("(time taken = $(round(t, digits=2)) s)")

########## EXPORT ##########
print("Exporting Trajectories...    ")
t = @elapsed begin

dfs = [trajectory_to_df(traj) for traj in trajectories]
metadata = trajectory_to_metadata.(trajectories)
    
ensemble_df = vcat(dfs...; cols=:union)
metadata_df = vcat(metadata...; cols=:union)

CSV.write("ensemble.csv", ensemble_df)
CSV.write("ensemble_metadata.csv", metadata_df)

end
println("(time taken = $(round(t, digits=2)) s)")

println("Complete!")
