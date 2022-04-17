"""
Simple implementation of the BTW sandpile model, roughly after

H. J. Jensen, _Self-Organized Criticality_, Cambridge University Press,
Cambridge (1998)

To use, create a Sandpile object and optionally a Sandpile_history,
and call `run` for the desired number of steps.
"""
module Sandpiles

export Sandpile,Sandpile_history,run

using Random

"struct to hold the instantaneous state of the sandpile"
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
"Return a history object initialized with the average
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
function run(pile::Sandpile,steps=1,history::Union{Sandpile_history,Nothing}=nothing)
    step0=pile.time
    while pile.time<step0+steps
        ir=rand(pile.rng,1:pile.L,2)
        pile.z[ir...]+=1
        pile.ztot+=1
        if pile.z[ir...]>pile.zc avalanche(pile,history)
        else
            pile.time+=1
            if !isnothing(history)
                push!(history.time,pile.time)
                push!(history.zav,pile.ztot/pile.L^2)
            end
        end

    end
end


"""
Develop an avalanche until no supercritial sites remain.  This
is called from `run` when some site becomes supercritical.
"""
function avalanche(pile::Sandpile,history::Union{Sandpile_history,Nothing}=nothing)
    more::Bool = true
    if !isnothing(history) push!(history.avalanche_start,pile.time) end
    avener=0
    avsites=zero(pile.z)
    while more
        z=copy(pile.z)

        # Bulk sites
        for jx = 2:pile.L-1, jy = 2:pile.L-1
            if z[jx, jy] > pile.zc
                avener += 1             # Increase energy count
                avsites[jx, jy] = 1     # Mark site as belonging to avalnche
                pile.z[jx, jy] -= 4     # Move grains
                pile.z[jx+1, jy] += 1
                pile.z[jx-1, jy] += 1
                pile.z[jx, jy+1] += 1
                pile.z[jx, jy-1] += 1
            end
        end
        # Border sites
        for jy=2:pile.L-1
            if z[1,jy]>pile.zc
                avener+=1
                avsites[1, jy] = 1
                pile.z[1,jy]-=3
                pile.z[1,jy+1]+=1
                pile.z[1,jy-1]+=1
                pile.z[2,jy]+=1
            end
            if z[pile.L,jy]>pile.zc
                avener+=1
                avsites[pile.L, jy] = 1
                pile.z[pile.L,jy]-=4
                pile.z[pile.L-1,jy]+=1
                pile.z[pile.L,jy+1]+=1
                pile.z[pile.L,jy-1]+=1
                pile.ztot-=1
            end
        end
        for jx=2:pile.L-1
            if z[jx,1]>pile.zc
                avener+=1
                avsites[jx,1] = 1
                pile.z[jx,1]-=3
                pile.z[jx-1,1]+=1
                pile.z[jx-1,1]+=1
                pile.z[jx,2]+=1
            end
            if z[jx,pile.L]>pile.zc
                avener+=1
                avsites[jx,pile.L] = 1
                pile.z[jx,pile.L]-=4
                pile.z[jx+1,pile.L]+=1
                pile.z[jx-1,pile.L]+=1
                pile.z[jx,pile.L-1]+=1
                pile.ztot-=1
            end
        end
        # Corners
        if z[1,1]>pile.zc
            avener+=1
            avsites[1,1] = 1
            pile.z[1,1]-=2
            pile.z[1,2]+=1
            pile.z[2,1]+=1
        end
        if z[1,pile.L]>pile.zc
            avener+=1
            avsites[1,pile.L] = 1
            pile.z[1,pile.L]-=3
            pile.z[2,pile.L]+=1
            pile.z[1,pile.L-1]+=1
            pile.ztot-=1
        end
        if z[pile.L,1]>pile.zc
            avener+=1
            avsites[pile.L,1] = 1
            pile.z[pile.L,1]-=3
            pile.z[pile.L-1,1]+=1
            pile.z[pile.L,2]+=1
            pile.ztot-=1
        end
        if z[pile.L,pile.L]>pile.zc
            avener+=1
            avsites[pile.L,pile.L] = 1
            pile.z[pile.L,pile.L]-=4
            pile.z[pile.L-1,pile.L]+=1
            pile.z[pile.L,pile.L-1]+=1
            pile.ztot-=2
        end

        pile.time+=1
        if !isnothing(history)
            push!(history.time,pile.time)
            push!(history.zav,pile.ztot/pile.L^2)
        end
        more=any(z->z>pile.zc,pile.z)
    end
    if !isnothing(history)
        push!(history.avalanche_duration,pile.time-history.avalanche_start[end])
        push!(history.avalanche_energy,avener)
        push!(history.avalanche_size,avsize(pile,avsites))
    end
end

function avsize(pile::Sandpile,avsites)
    CM=Float64[0, 0]
    N=0
    for ix=1:pile.L, iy=1:pile.L
        if avsites[ix,iy]>0
            N+=1
            CM+=[ix, iy]
        end
    end
    CM/=N
    l=0.
    for ix=1:pile.L, iy=1:pile.L
        if avsites[ix,iy]>0
            l+= sqrt( (ix-CM[1])^2 + (iy-CM[2])^2 )
        end
    end
    l/=N
    return l
end

end
