using Revise
using Dates, CSV, DataFrames, DimensionalData, Tables, Plots

using LagrangianERA5




#### Functions from run_traj ####

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

function lows_to_dataframe(lows)

    return DataFrame(
        time = [low.time for low in lows],
        low_id = [low.id for low in lows],
        lat  = [low.lat  for low in lows],
        lon  = [low.lon  for low in lows]
    )
end

using DataFrames
using CSV

function lows_to_basin_dataframe(lows)

    rows = Vector{NamedTuple}()

    for low in lows

        mask = low.mask
        sublon = low.lon_window
        sublat = low.lat_window

        for i in 1:length(sublat)
            for j in 1:length(sublon)

                if mask[i,j]

                    push!(rows, (
                        time = low.time,
                        low_id = low.id,
                        lon = sublon[j],
                        lat = sublat[i]
                    ))

                end
            end
        end
    end

    return DataFrame(rows)
end


#### Test code ####

println("        pressure_systems_test.jl          ")
println("----------------------------------------\n")

start_time = DateTime(2019,1,3)
end_time   = DateTime(2019,1,4)

dates = collect(start_time:Day(1):end_time)
datadir = "/Users/lewisclark/Documents/python_work/lagrangian/data"

print("Reading data...              ")

t = @elapsed begin 
ps_da = LagrangianERA5.GeoData.load_surface_range(dates; datadir=datadir)
end



println("(time taken = $(round(t, digits=2)) s)")

print("Resolving weather systems... ")

t = @elapsed begin 
lows = LagrangianERA5.PressureSystems.detect_lows(ps_da)
#lows = LagrangianERA5.PressureSystems.cluster_lows(lows, 150.0)
storms = LagrangianERA5.PressureSystems.track_storms(lows)
end

println("(time taken = $(round(t, digits=2)) s)")

df_storms = tracks_to_dataframe(storms)
df_basin = lows_to_basin_dataframe(lows)
df_lows = lows_to_dataframe(lows)

CSV.write("pressure_lows.csv", df_lows)
CSV.write("storm_tracks.csv", df_storms)
CSV.write("pressure_basins.csv", df_basin)


println("Complete!")