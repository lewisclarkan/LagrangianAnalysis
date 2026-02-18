using Revise
using Dates, CSV, DataFrames, DimensionalData, Tables

using LagrangianERA5

function round_to_nearest_hour(t::DateTime)
    minute_offset = Minute(minute(t))
    second_offset = Second(second(t))

    t0 = t - minute_offset - second_offset

    if minute(t) ≥ 30
        return t0 + Hour(1)
    else
        return t0
    end
end

function storm_relative_polar(lat, lon, lat_s, lon_s)

    # Great circle distance
    r = LagrangianERA5.PressureSystems.haversine_km(lat, lon, lat_s, lon_s)

    # Local Cartesian approximation
    dx = (lon - lon_s) * cosd(lat_s)
    dy = (lat - lat_s)

    θ = atan(dy, dx)   # radians

    return r, θ
end

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
    times_vec = sort!(collect(keys(ps_itps)))

    return function(lon, lat, t::DateTime)
        idx = LagrangianERA5.TrajectoryModule.searchsortednearest(times_vec, t)
        t_sel = times_vec[idx]
        return ps_itps[t_sel](lon, lat)
    end
end

function count_crossings(traj, key, n, i)
    if Bool(traj.extras[key][i]) && !Bool(traj.extras[key][i+1])
        n += 1
    end

    return n
end

function count_time(traj, key, n, i)
    if Bool(traj.extras[key][i])
        n += 1
    end

    return n
end     

function reverse_without_last(v)
    length(v) <= 1 && return similar(v,0)
    return reverse(@view v[1:end-1])
end

function integrate_bidirectional(
    traj_id     ::Int,
    x0          ::NTuple{3,Float64},
    times       ::AbstractVector{DateTime},
    itps        ::AbstractDict{DateTime, <:Tuple},
    start_idx   ::Int;
    extra_itps  ::AbstractDict{Symbol, <:AbstractDict{DateTime, <:Any}},
    vars        ::Vector{Symbol})

    # slice times around pivot
    times_back = times[1:start_idx]
    times_forw = times[start_idx:end]

    # perform the backward and forward integration
    traj_back = LagrangianERA5.Integrator.integrate_trajectory(traj_id, x0, times_back, itps; direction=:backward)
    traj_forw = LagrangianERA5.Integrator.integrate_trajectory(traj_id, x0, times_forw, itps; direction=:forward)

    # reverse backward segment to chronological order
    t_back      = reverse_without_last(traj_back.time)
    lon_back    = reverse_without_last(traj_back.lon)
    lat_back    = reverse_without_last(traj_back.lat)
    p_back      = reverse_without_last(traj_back.p)

    # combine backward and forward segments
    time_comb = vcat(t_back, traj_forw.time)
    lon_comb  = vcat(lon_back, traj_forw.lon)
    lat_comb  = vcat(lat_back, traj_forw.lat)
    p_comb    = vcat(p_back, traj_forw.p)

    # build trajectory struct
    traj = LagrangianERA5.TrajectoryModule.Trajectory(traj_id, time_comb, lon_comb, lat_comb, p_comb)

    # compute extra quantities along trajectory
    LagrangianERA5.TrajectoryModule.interpolate_along_trajectory!(traj, extra_itps, vars)
    LagrangianERA5.TrajectoryModule.compute_geopotential_height!(traj)
    LagrangianERA5.TrajectoryModule.compute_rhi!(traj)
    LagrangianERA5.TrajectoryModule.compute_potential_temperature!(traj)
    LagrangianERA5.TrajectoryModule.compute_adiabatic_dTdt!(traj)
    LagrangianERA5.TrajectoryModule.compute_time_derivative!(traj, :t, :dTdt)
    LagrangianERA5.TrajectoryModule.compute_time_derivative!(traj, :q, :dqdt)
    LagrangianERA5.TrajectoryModule.compute_curvature!(traj)
    LagrangianERA5.TrajectoryModule.compute_cold_point!(traj)
    LagrangianERA5.TrajectoryModule.compute_wind_vector!(traj)

    return traj

end

function formation_stats(traj::LagrangianERA5.TrajectoryModule.Trajectory, i_start::Int)

    window = (i_start - 18):(i_start - 1)
    tsec = [(Dates.value(t - traj.time[1])) / 1000.0 for t in traj.time]

    q0 = traj.q[window[1]]
    dq = traj.q[window[end]] - q0
    M  = dq / q0

    dq_si = traj.extras[:q_si][window[end]] - traj.extras[:q_si][window[1]]
    saturation_control = abs(dq) / (abs(dq) + abs(dq_si))

 
    T_avg  = LagrangianERA5.TrajectoryModule.mean(traj.t[window])
    dT     = traj.t[window[end]] - traj.t[window[1]]
    dtheta = traj.extras[:potential_t][window[end]] - traj.extras[:potential_t][window[1]]
    dT_adiab = sum(traj.extras[:dTdt_adiab][window]) * LagrangianERA5.TrajectoryModule.mean(diff(tsec[window]))

    T = - 2.5e6 / (461.5 * T_avg ^ 2) * (dT)

    ascent = sum(max.(traj.w[window], 0.0)) * LagrangianERA5.TrajectoryModule.mean(diff(tsec[window]))



    return (
        M = M,
        T = T,
        T_avg = T_avg,
        dT = dT,
        dT_adiab = dT_adiab,
        dtheta = dtheta,
        saturation_control = saturation_control,
        integrated_ascent = ascent
    )

    return 

end

function extract_issrs(traj::LagrangianERA5.TrajectoryModule.Trajectory)

    issr_bounds = LagrangianERA5.ISSRModule.detect_issrs(traj.extras[:rhi])
    events = Vector{LagrangianERA5.ISSRModule.ISSREvent}()

    for (k, (i_start, i_end)) in enumerate(issr_bounds)

        duration = i_end - i_start + 1
        mean_rhi = LagrangianERA5.TrajectoryModule.mean(traj.extras[:rhi][i_start:i_end])
        max_rhi  = maximum(traj.extras[:rhi][i_start:i_end])

        has_window = i_start > 18

        recent_event = any(prev_end >= i_start - 18 for (_, prev_end) in issr_bounds[1:k-1])

        if has_window
            stats = formation_stats(traj, i_start)
            push!(events, LagrangianERA5.ISSRModule.ISSREvent(
                traj.id, k, i_start, i_end, traj.time[i_start], traj.time[i_end], mean_rhi, max_rhi, true, recent_event,
                stats.M, stats.T, stats.T_avg, stats.dT, stats.dT_adiab, stats.dtheta, stats.saturation_control, stats.integrated_ascent
            ))
        else
            push!(events, LagrangianERA5.ISSRModule.ISSREvent(
                traj.id, k, i_start, i_end, traj.time[i_start], traj.time[i_end], mean_rhi, max_rhi, false, recent_event,
                missing, missing, missing, missing, missing, missing, missing, missing
            ))
        end

    end

    return events

end

function trajectory_to_df(traj::LagrangianERA5.TrajectoryModule.Trajectory)
    N = length(traj.time)

    df = DataFrame(
        traj_id = fill(traj.id, N),

        time    = traj.time,
        lon     = traj.lon,
        lat     = traj.lat,
        p       = traj.p,

        z       = traj.z,
        gh      = traj.gh,

        t       = traj.t,
        u       = traj.u,
        v       = traj.v,
        w       = traj.w,
        q       = traj.q,
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
            row[k] = v
        end
    end

    return DataFrame(row)
end

function tracks_to_dataframe(tracks)

    rows = []

    for tr in tracks
        for i in eachindex(tr.times)
            push!(rows, (
                storm_id = tr.id,
                time     = tr.times[i],
                lat      = tr.lats[i],
                lon      = tr.lons[i]
            ))
        end
    end

    return DataFrame(rows)
end

function system_to_dataframe(lows)

    return DataFrame(
        time = [low.time for low in lows],
        low_id = [low.id for low in lows],
        track_id = [low.track_id for low in lows],
        lat  = [low.lat  for low in lows],
        lon  = [low.lon  for low in lows]
    )
end

function basin_to_dataframe(lows)

    rows = Vector{NamedTuple}()

    for low in lows

        polygon = low.edge
        polygon === nothing && continue

        for (k, (lon, lat)) in enumerate(polygon)

            push!(rows, (
                time        = low.time,
                low_id      = low.id,
                vertex_id   = k,
                lon         = lon,
                lat         = lat
            ))

        end
    end

    return DataFrame(rows)
end

function build_storm_center_lookup(systems::Vector{LagrangianERA5.PressureSystems.PressureSystem})

    lookup = Dict{DateTime, Dict{Int,Tuple{Float64,Float64}}}()

    for sys in systems

        t = sys.time
        tid = sys.track_id

        inner = get!(lookup, t, Dict{Int,Tuple{Float64,Float64}}())
        inner[tid] = (sys.lat, sys.lon)
    end

    return lookup
end


function storm_relative_polar(lat_p, lon_p, lat_s, lon_s)

    # distance
    r = LagrangianERA5.PressureSystems.haversine_km(lat_p, lon_p, lat_s, lon_s)

    # bearing
    φ1 = deg2rad(lat_s)
    φ2 = deg2rad(lat_p)
    Δλ = deg2rad(lon_p - lon_s)

    y = sin(Δλ) * cos(φ2)
    x = cos(φ1)*sin(φ2) - sin(φ1)*cos(φ2)*cos(Δλ)

    θ = atan(y, x)   # radians

    return r, θ
end


function attach_roi_influence!(
    traj::LagrangianERA5.TrajectoryModule.Trajectory,
    rois_by_time::Dict{DateTime,Vector{LagrangianERA5.PressureSystems.ROI}},
    storm_center_by_time::Dict{DateTime,Dict{Int,Tuple{Float64,Float64}}},
    storm_type::Dict{Int,Symbol})

    n = length(traj.time)

    storm_id   = fill(-1, n)
    storm_r    = fill(NaN, n)
    storm_theta = fill(NaN, n)
    in_storm   = falses(n)
    in_low     = falses(n)
    in_high    = falses(n)

    current_hour = DateTime(0)
    current_rois = LagrangianERA5.PressureSystems.ROI[]
    current_centers = Dict{Int,Tuple{Float64,Float64}}()

    for i in 1:n

        # ---- causal hourly snap ----
        t_hour = floor(traj.time[i], Dates.Hour)

        if t_hour != current_hour
            current_hour = t_hour
            current_rois = get(rois_by_time, t_hour, LagrangianERA5.PressureSystems.ROI[])
            current_centers = get(storm_center_by_time, t_hour,
                                  Dict{Int,Tuple{Float64,Float64}}())
        end

        isempty(current_rois) && continue

        # ---- ROI membership ----
        roi = LagrangianERA5.PressureSystems.find_roi_membership(
            traj.lon[i],
            traj.lat[i],
            current_rois
        )

        if roi !== nothing

            storm_id[i] = roi
            in_storm[i] = true

            # ---- centre lookup ----
            lat_s, lon_s = current_centers[roi]

            # ---- polar coordinates ----
            r, θ = storm_relative_polar(
                traj.lat[i],
                traj.lon[i],
                lat_s,
                lon_s
            )

            storm_r[i] = r
            storm_theta[i] = θ

            # ---- storm type ----
            type = storm_type[roi]

            if type == :low
                in_low[i] = true
            elseif type == :high
                in_high[i] = true
            end
        end
    end

    traj.extras["storm_id"] = storm_id
    traj.extras["storm_r_km"] = storm_r
    traj.extras["storm_theta_rad"] = storm_theta
    traj.extras["in_storm"] = in_storm
    traj.extras["in_low"]   = in_low
    traj.extras["in_high"]  = in_high

    return nothing
end




########## CONFIG ##########

# Full set

start_time = DateTime(2019,1,3)
end_time   = DateTime(2019,1,8)
n_time     = 31

lon_min   = -50.0
lon_max   = -10.0
lon_steps = 30

lat_min   = 35.0
lat_max   = 65.0
lat_steps = 30

pressure_level = 25000.0

vars = [:u, :v, :w, :t, :q, :d, :z]

# Reduced set for testing

"""start_time = DateTime(2019,1,3)
end_time   = DateTime(2019,1,3)
n_time     = 1

lon_min   = -50.0
lon_max   = -10.0
lon_steps = 30

lat_min   = 35.0
lat_max   = 65.0
lat_steps = 30

pressure_level = 25000.0

vars = [:u, :v, :w, :t, :q, :z]"""

println("Lagrangian IceSupersaturated Tool (LIST)")
println("----------------------------------------")
println("Program start time: $(Dates.format(now(), "HH:MM:SS"))\n")

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

##### PRESSURE SYSTEMS #######

print("Resolving pressure systems...")

t = @elapsed begin 

lows = LagrangianERA5.PressureSystems.detect_systems(ps_da; mode=:low, threshold=1500.0)
highs = LagrangianERA5.PressureSystems.detect_systems(ps_da; mode=:high, threshold=1200.0)

low_tracks = LagrangianERA5.PressureSystems.track_systems(lows; kind=:low)
high_tracks = LagrangianERA5.PressureSystems.track_systems(highs, kind=:high)

systems = vcat(lows, highs)

rois = LagrangianERA5.PressureSystems.build_rois(systems)
rois_by_time = LagrangianERA5.PressureSystems.group_lows_by_time(rois)

storm_type = Dict(sys.track_id => sys.type for sys in systems)
storm_center_by_time = build_storm_center_lookup(systems)

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
all_events   = Vector{LagrangianERA5.ISSRModule.ISSREvent}()

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

    traj   = integrate_bidirectional(traj_id, x0_pos, times, itps, local_start_idx;
                                     extra_itps=extra_itps, vars=vars)
    events = extract_issrs(traj)

    attach_roi_influence!(
        traj,
        rois_by_time,
        storm_center_by_time,
        storm_type
    )

    push!(trajectories, traj)
    append!(all_events, events)
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
all_events_df = DataFrame(Tables.columntable(all_events))

df_lows  = system_to_dataframe(lows)
df_highs = system_to_dataframe(highs)

df_basin_lows  = basin_to_dataframe(lows)
df_basin_highs = basin_to_dataframe(highs)

df_low_tracks = tracks_to_dataframe(low_tracks)
df_high_tracks = tracks_to_dataframe(high_tracks)

CSV.write("pressure_lows.csv", df_lows)
CSV.write("pressure_highs.csv", df_highs)

CSV.write("lows_basins.csv", df_basin_lows)
CSV.write("highs_basins.csv", df_basin_highs)

CSV.write("low_tracks.csv", df_low_tracks)
CSV.write("high_tracks.csv", df_high_tracks)

CSV.write("ensemble.csv", ensemble_df)
CSV.write("ensemble_metadata.csv", metadata_df)
CSV.write("all_events.csv", all_events_df)

end
println("(time taken = $(round(t, digits=2)) s)")

println("Complete!")
