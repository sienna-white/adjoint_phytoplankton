


"True values"
"pmax" => 0.05 * hr2s,           # maximum specific growth rate [1/hour]
"Li" => 0.005 * hr2s,             # specific loss rate [1/hour]

                
"first guess" values
"pmax" => 0.04 * hr2s,           # maximum specific growth rate [1/hour]
"Li" => 0.005 * hr2s,             # specific loss rate [1/hour]


June 23 / 2025
Trying to split up growth + loss in the adjoint model. Right now it's not working but hopefully soon!
Split them up, and added a penalty function following guidance from this paper: https://www.mdpi.com/2072-4292/15/1/148
