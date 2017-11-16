#!/usr/bin/env python3
import wdt

# Compute cost field from image
cost_field = wdt.map_image_to_costs('images/ex2.png')
# Plot the cost field
wdt.plot(cost_field)
# Compute the distance transform from the cost field
distance_transform = wdt.get_weighted_distance_transform(cost_field)
# Plot the result
wdt.plot(distance_transform)
