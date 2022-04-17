# avalanche1.jl -- Simplest implmentation avalanche dynamcs
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
#
# This is the simplest implementation of avalanche dynamcs.  It can be
# actually slightly faster that the more sophisticated implentation
# in avalanche2.jl for very small piles (less than 40x40 or so)
#

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
