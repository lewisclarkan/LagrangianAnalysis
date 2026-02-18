using Revise
using Dates, CSV, DataFrames, DimensionalData, Tables, Plots

using LagrangianERA5

#### Functions (some from run_traj.jl) ####

function storm_relative_polar(lat, lon, lat_s, lon_s)

    # Great circle distance
    r = LagrangianERA5.PressureSystems.haversine_km(lat, lon, lat_s, lon_s)

    # Local Cartesian approximation
    dx = (lon - lon_s) * cosd(lat_s)
    dy = (lat - lat_s)

    θ = atan(dy, dx)   # radians

    return r, θ
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



############################# Test Code ########################

println("pressure_systems.jl                   test program")
println("--------------------------------------------------")
println("Program start time: $(Dates.format(now(), "HH:MM:SS"))\n")

start_time = DateTime(2019,1,3)
end_time   = DateTime(2019,1,5)

dates = collect(start_time:Day(1):end_time)
datadir = "/Users/lewisclark/Documents/python_work/lagrangian/data"

print("Reading data...              ")

t = @elapsed begin 
ps_da = LagrangianERA5.GeoData.load_surface_range(dates; datadir=datadir)
end

###############################################################

println("(time taken = $(round(t, digits=2)) s)")

print("Resolving weather systems... ")

t = @elapsed begin 
lows = LagrangianERA5.PressureSystems.detect_systems(ps_da; mode=:low, threshold=1500.0)
highs = LagrangianERA5.PressureSystems.detect_systems(ps_da; mode=:high, threshold=1200.0)

low_tracks = LagrangianERA5.PressureSystems.track_systems(lows; kind=:low)
high_tracks = LagrangianERA5.PressureSystems.track_systems(highs, kind=:high)
end

println("(time taken = $(round(t, digits=1)) s)")

###############################################################

print("Exporting data...            ")

t = @elapsed begin
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
end

println("(time taken = $(round(t, digits=2)) s)")


println("Complete!")