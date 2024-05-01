using Distributions

"""Helper function for trapezoidal rule"""
function trapz(x, y)
    return sum((x[2:end] - x[1:(end - 1)]) .* (y[2:end] + y[1:(end - 1)])) * 0.5
end

# Function adjusts depths for levee protection
function leveer(depth_ft_g, levee_height)
    post_levee = [] 
    for depth in depth_ft_g  #check each depth against the levee height
        if depth >= levee_height
            depth = depth # if above levee, depth remains
        else
            depth = 0 # if below levee, set depth to 0
        end
        push!(post_levee, depth) 
        #println(post_levee)
    end
    return post_levee 
end



# Run simulation given action, SOW, model paramters, and levee height

function run_sim(a::Action, sow::SOW, p::ModelParams, levee_height::Float64)
    # Calculate the cost of elevating the house bassed on action
    construction_cost = elevation_cost(p.house, a.Δh_ft)

    # Generate a range of storm surge heights
    storm_surges_ft = range(
        quantile(sow.surge_dist, 0.0005); stop=quantile(sow.surge_dist, 0.9995), length=130)

    # Map over the years to compute Expected Annual Damages (EADs)
    eads = map(p.years) do year
        # Get the sea level rise for the year
        slr_ft = sow.slr(year)

        # adjust storm surge depths for sea level rise
        # apply levee protection
        depth_ft_gauge = storm_surges_ft .+ slr_ft
        effective_depth_ft_gauge = leveer(depth_ft_gauge, levee_height)

        # Adjust for the house's elevation above the gauge and any additional elevation due to actions
        effective_depth_ft_house = effective_depth_ft_gauge .- (p.house.height_above_gauge_ft + a.Δh_ft)

        # Calculate the fraction of damages and their weighted value
        damages_frac = p.house.ddf.(effective_depth_ft_house) ./ 100
        weighted_damages = damages_frac .* pdf.(sow.surge_dist, storm_surges_ft)

        # Trapezoidal integration of weighted damages
        ead = trapz(storm_surges_ft, weighted_damages) * p.house.value_usd
        return ead
    end

    # Compute NPV of the damages 
    years_idx = p.years .- minimum(p.years)
    discount_fracs = (1 .- sow.discount_rate) .^ years_idx
    ead_npv = sum(eads .* discount_fracs)
    return -(ead_npv + construction_cost)
end


