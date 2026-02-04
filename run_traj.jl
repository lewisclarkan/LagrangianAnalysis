using Revise
using Dates, CSV, DataFrames, DimensionalData, Tables

using LagrangianERA5

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
    LagrangianERA5.TrajectoryModule.compute_geopotential_height!(traj)
    LagrangianERA5.TrajectoryModule.interpolate_along_trajectory!(traj, extra_itps, vars)
    LagrangianERA5.TrajectoryModule.compute_rhi!(traj)
    LagrangianERA5.TrajectoryModule.compute_potential_temperature!(traj)
    LagrangianERA5.TrajectoryModule.compute_adiabatic_dTdt!(traj)
    LagrangianERA5.TrajectoryModule.compute_time_derivative!(traj, :t, :dTdt)
    LagrangianERA5.TrajectoryModule.compute_time_derivative!(traj, :q, :dqdt)
    LagrangianERA5.TrajectoryModule.compute_curvature!(traj)
    LagrangianERA5.TrajectoryModule.compute_cold_point!(traj)

    return traj

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

println("Lagrangian IceSupersaturated Tool (LIST)")
println("----------------------------------------\n")

########## CONFIG ##########

# Full set

"""start_time = DateTime(2019,1,3)
end_time   = DateTime(2019,1,8)
n_time     = 31

lon_min   = -50.0
lon_max   = -10.0
lon_steps = 30

lat_min   = 35.0
lat_max   = 65.0
lat_steps = 30

pressure_level = 25000.0

vars = [:u, :v, :w, :t, :q, :pv, :d]"""

# Reduced set for testing

start_time = DateTime(2019,1,3)
end_time   = DateTime(2019,1,3)
n_time     = 1

lon_min   = -50.0
lon_max   = -10.0
lon_steps = 30

lat_min   = 35.0
lat_max   = 65.0
lat_steps = 30

pressure_level = 25000.0

vars = [:u, :v, :w, :t, :q, :z]

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
                                   extra_itps=extra_itps, vars=vars)
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
