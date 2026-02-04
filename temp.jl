    """if haskey(traj.extras, "mask")
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
            #sum_d = 0
            count = 0

            for i in j:(k-1)

                # Count the number of crossings (a measure of intermittency)
                n_drop_110 = count_crossings(traj, "mask110", n_drop_110, i)
                n_drop_100 = count_crossings(traj, "mask100", n_drop_100, i)

                # Time above thresholds
                n_above_110 = count_time(traj, "mask110", n_above_110, i)
                n_above_100 = count_time(traj, "mask100", n_above_100, i)

                # Calculate the averages over the ice supersaturated region
                sum_rhi += traj.extras["rhi"][i]
                sum_w += traj.extras["w"][i]
                sum_dqdt += traj.extras["dqdt"][i]
                sum_curv += traj.extras["curv"][i]
                #sum_d += traj.extras["d"][i]
                count += 1
            end
            
            traj.extras["n_rhi100_drops"] = n_drop_100
            traj.extras["n_rhi110_drops"] = n_drop_110

            traj.extras["t_rhi100_hours"] = 0.1 * n_above_100
            traj.extras["t_rhi110_hours"] = 0.1 * n_above_110

            if k > j # i.e. most issrs

                LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "w", j, 30)
                #LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "d", j, 30)
                LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "dTdt", j, 30)
                LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "dqdt", j, 30)
                LagrangianERA5.TrajectoryModule.compute_deltaT!(traj, j, 30)

                traj.extras["max_rhi_issr"] = maximum(traj.extras["rhi"][j:(k-1)])
                traj.extras["max_w_issr"]   = maximum(traj.extras["w"][j:(k-1)])

                traj.extras["avg_rhi_issr"]         = sum_rhi / count
                traj.extras["avg_w_issr"]           = sum_w / count
                traj.extras["avg_dqdt_issr"]        = sum_dqdt / count
                traj.extras["avg_curvature_issr"]   = sum_curv / count
                #traj.extras["avg_divergence_issr"]  = sum_d / count

                traj.extras["delta_q_issr"] = traj.extras["q"][k] - traj.extras["q"][j]

                traj.extras["dt_cold_to_event"]         = (j - traj.extras["idx_cold_point"]) * 0.1
                traj.extras["dt_cold_to_event_pre"]     = (j - traj.extras["idx_cold_point_before"]) * 0.1

                traj.extras["init_w_issr"]         = traj.extras["w"][j]

            elseif k==j # i.e. issrs or length 1

                LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "w", j, 30)
                #LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "d", j, 30)
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
                #traj.extras["avg_divergence_issr"]= traj.extras["d"][j]

                traj.extras["delta_q_issr"] = 0.0

                traj.extras["dt_cold_to_event"]     = (j - traj.extras["idx_cold_point"]) * 0.1
                traj.extras["dt_cold_to_event_pre"] = (j - traj.extras["idx_cold_point_before"]) * 0.1

                traj.extras["init_w_issr"]          = traj.extras["w"][j]
                
            end

        # Store the data for the case that the sample point was not ice supersaturated
        else

            LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "w", pidx, 30)
            #LagrangianERA5.TrajectoryModule.compute_average_quantity_before_event!(traj, "d", pidx, 30)
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
            #traj.extras["avg_divergence_issr"]  = missing

            traj.extras["delta_q_issr"]         = missing

            traj.extras["dt_cold_to_event"]     = (pidx - traj.extras["idx_cold_point"]) * 0.1
            traj.extras["dt_cold_to_event_pre"] = missing

        end
    end"""