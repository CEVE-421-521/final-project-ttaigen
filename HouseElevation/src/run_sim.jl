using Distributions

"""Helper function for trapezoidal rule"""
function trapz(x, y)
    return sum((x[2:end] - x[1:(end - 1)]) .* (y[2:end] + y[1:(end - 1)])) * 0.5
end


function leveer(depth_ft_g, levee_height)
    #for depth in depth_ft
    #println("range ", depth_ft_g)
    post_levee = []
    for depth in depth_ft_g  #for each depth in the range
        if depth >= levee_height
            depth = depth
        else
            depth = 0
        end
        push!(post_levee, depth)
        #println(post_levee)
    end
    # start = first(depth_ft_g)
    # #println("before function", start)s
    # stop = last(depth_ft_g)
    # if start >= levee_height
    #    start = start
    # else
    #     start = levee_height
    # end
    # #println("after function", start)
    # post_levee = range(start; stop=stop, length=130)
    #println("range after", post_levee)
    return post_levee 
end

"""
Run the model for a given action, state of the world (SOW), model parameters, and levee height

Expected Annual Damages are computed using the trapezoidal rule
"""
function run_sim(a::Action, sow::SOW, p::ModelParams, levee_height::Float64)
    # Calculate the cost of elevating the house
    construction_cost = elevation_cost(p.house, a.Δh_ft)

    # Generate storm surge heights
    storm_surges_ft = range(
        quantile(sow.surge_dist, 0.0005); stop=quantile(sow.surge_dist, 0.9995), length=130)

    # Map over the years to compute Expected Annual Damages (EADs)
    eads = map(p.years) do year
        # Get the sea level rise for the year
        slr_ft = sow.slr(year)

        # Compute depth at the gauge and adjust for the levee
        depth_ft_gauge = storm_surges_ft .+ slr_ft
        #println(depth_ft_gauge)
        #println.("Guage depth ", depth_ft_gauge)
        #effective_depth_ft_gauge = max.(0.0, depth_ft_gauge .- levee_height)
        #effective_depth_ft_gauge = depth_ft_gauge
        #println(depth_ft_gauge)

        # test_list = []
        # for i in range(1, 130)
        #     push!(test_list, i)
        # end
        # effective_depth_ft_gauge = test_list
        effective_depth_ft_gauge = leveer(depth_ft_gauge, levee_height)

        #println.("post-levee depth ", effective_depth_ft_gauge)

        # Adjust for the house's elevation above the gauge and any additional elevation due to actions
        effective_depth_ft_house = effective_depth_ft_gauge .- (p.house.height_above_gauge_ft + a.Δh_ft)

        # Calculate the fraction of damages and their weighted value
        damages_frac = p.house.ddf.(effective_depth_ft_house) ./ 100
        weighted_damages = damages_frac .* pdf.(sow.surge_dist, storm_surges_ft)

        # Trapezoidal integration of weighted damages
        ead = trapz(storm_surges_ft, weighted_damages) * p.house.value_usd
        return ead
    end

    # Compute NPV
    years_idx = p.years .- minimum(p.years)
    discount_fracs = (1 .- sow.discount_rate) .^ years_idx
    ead_npv = sum(eads .* discount_fracs)
    return -(ead_npv + construction_cost)
end


