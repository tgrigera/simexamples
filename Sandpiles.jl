# Sandpiles.jl -- Two-dimensional BTW sandpile
#
# This file copyright (C) 2022 by Tomas S. Grigera.
# 
# This is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License (GPL) as
# published by the Free Software Foundation. You can use either
# version 3, or (at your option) any later version.
#
# This software is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.
#
# For details see the file LICENSE in the home directory.
#

"""
Simple implementation of the BTW sandpile model, roughly after

H. J. Jensen, _Self-Organized Criticality_, Cambridge University Press,
Cambridge (1998)

To use, create a Sandpile object and optionally a Sandpile_history,
and call `run` for the desired number of steps.
"""
module Sandpiles

export Sandpile,Sandpile_history, run!

using Random

"The `Sandpile` struct holds the instantaneous state of the sandpile"
mutable struct Sandpile
    L::Int
    z::Matrix{Int}
    zc::Int
    ztot::BigInt
    time::Int
    rng
end

"""
    Sandpile(L,zc::Int,rng=Random.GLOBAL_RNG)

Return a new sandpile with random heights, setting size `L`
and the critical height `zc`.  Optionally, an appropriately
seeded random number generator can be passed, which will be
used in the simulation.
"""
function Sandpile(L,zc::Int,rng=Random.GLOBAL_RNG)
    s=Sandpile(L,rand(rng,1:zc,L,L),zc,0,0,rng)
    s.ztot=sum(s.z)
    return s
end

"""
struct to hold historical information of the sandpile: average height
and avalanche characteristics (start, duration, size and energy).
Will be properly updated if passed to the `run` method.
"""
mutable struct Sandpile_history
    time::Vector{Int}
    zav::Vector{Float64}
    avalanche_start::Vector{Int}
    avalanche_duration::Vector{Int}
    avalanche_size::Vector{Float64}
    avalanche_energy::Vector{Int}
end

"Return an empty history object"
Sandpile_history()=Sandpile_history(Int[],Float64[],Int[],Int[],Float64[],Int[])

"""
Return a history object initialized with the average
height of the given sandpile, and empty avalanche
information
"""
Sandpile_history(pile::Sandpile)=Sandpile_history([pile.time],[pile.ztot/pile.L^2],Int[],Int[],Float64[],Int[])

"""
    run(pile::Sandpile,steps=1,history::Union{Sandpile_history,Nothing}=nothing)

Do at least `steps` steps of dynamical evolution for the sandpile
`pile`.  If at the requested number of steps an avalanche
is ongoing, continue until avalanche is finished and all
sites are subcritical.
"""
function run!(pile::Sandpile,steps=1,history::Union{Sandpile_history,Nothing}=nothing)
    step0=pile.time
    while pile.time<step0+steps
        ir=rand(pile.rng,1:pile.L,2)
        pile.z[ir...]+=1
        pile.ztot+=1
        # Choose which implementation of avalanches to use
        # if pile.z[ir...]>pile.zc avalanche(pile,history)
        if pile.z[ir...]>pile.zc avalanche(pile,(ir[1],ir[2]),history)
        else
            pile.time+=1
            if !isnothing(history)
                push!(history.time,pile.time)
                push!(history.zav,pile.ztot/pile.L^2)
            end
        end

    end
end

# Two implementations of avalanche dynamis

include("avalanche1.jl")
include("avalanche2.jl")


end
