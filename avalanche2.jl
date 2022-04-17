# avalanche1.jl -- Faster implementation of avalanche dynamcs
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
#
# A more sophisticated implementation of avalanche dynamics.  Noticeably
# faster than avalanche1.jl for piles larger than about 100x100
#

Siteset = Set{Tuple{Int,Int}}

"""
Develop an avalanche until no supercritial sites remain.  This
is called from `run!` when some site becomes supercritical.
"""
function avalanche(pile::Sandpile,scsite::Tuple{Int,Int},history::Union{Sandpile_history,Nothing}=nothing)
    if !isnothing(history) push!(history.avalanche_start,pile.time) end
    supercritical_sites::Siteset=Set([scsite])
    avener=0
    cluster=Set(supercritical_sites)
    while length(supercritical_sites)>0
        avener+=length(supercritical_sites)
        relax(pile,supercritical_sites)
        pile.time+=1
        if !isnothing(history)
            push!(history.time,pile.time)
            push!(history.zav,pile.ztot/pile.L^2)
        end
        update_sc!(pile,supercritical_sites)
        union!(cluster,supercritical_sites)
    end
    if !isnothing(history)
        push!(history.avalanche_duration,pile.time-history.avalanche_start[end])
        push!(history.avalanche_energy,avener)
        push!(history.avalanche_size,csize(cluster))
    end
end

"Compute size of cluster"
function csize(cluster::Set)
    CM=Float64[0, 0]
    N=length(cluster)
    for (ix,iy) in cluster CM+=[ix, iy] end
    CM/=N
    l=0.
    for (ix,iy) in cluster l+= sqrt( (ix-CM[1])^2 + (iy-CM[2])^2 ) end
    l/=N
    return l
end

"""
Relax all the sites in the `supercritical_sites` collection.
Sites are not checked again, it is assumend that `supercritical_sites`
holds the right list.  May generate new supercritical sites.
"""
function relax(pile::Sandpile,supercritical_sites)
    for (jx,jy) in supercritical_sites
        if jx==1
            if jy==1
                pile.z[1,1]-=2
                pile.z[1,2]+=1
                pile.z[2,1]+=1
            elseif jy==pile.L
                pile.z[1,pile.L]-=3
                pile.z[2,pile.L]+=1
                pile.z[1,pile.L-1]+=1
                pile.ztot-=1
            else
                pile.z[1,jy]-=3
                pile.z[1,jy+1]+=1
                pile.z[1,jy-1]+=1
                pile.z[2,jy]+=1
            end
        elseif jx==pile.L
            if jy==1
                pile.z[pile.L,1]-=3
                pile.z[pile.L-1,1]+=1
                pile.z[pile.L,2]+=1
                pile.ztot-=1
            elseif jy==pile.L
                pile.z[pile.L,pile.L]-=4
                pile.z[pile.L-1,pile.L]+=1
                pile.z[pile.L,pile.L-1]+=1
                pile.ztot-=2
            else
                pile.z[pile.L,jy]-=4
                pile.z[pile.L-1,jy]+=1
                pile.z[pile.L,jy+1]+=1
                pile.z[pile.L,jy-1]+=1
                pile.ztot-=1
            end
        elseif jy==1
            pile.z[jx,1]-=3
            pile.z[jx-1,1]+=1
            pile.z[jx-1,1]+=1
            pile.z[jx,2]+=1
        elseif jy==pile.L
            pile.z[jx,pile.L]-=4
            pile.z[jx+1,pile.L]+=1
            pile.z[jx-1,pile.L]+=1
            pile.z[jx,pile.L-1]+=1
            pile.ztot-=1
        else
            pile.z[jx, jy] -= 4
            pile.z[jx+1, jy] += 1
            pile.z[jx-1, jy] += 1
            pile.z[jx, jy+1] += 1
            pile.z[jx, jy-1] += 1
        end
    end
end

@inline check(pile,sc,jx,jy) =  if pile.z[jx,jy]>pile.zc union!(sc,[(jx,jy)]) end

"""
Update the list (a set actually) of supercritical sites.
"""
function update_sc!(pile::Sandpile,supercritical_sites::Siteset)
    sc=copy(supercritical_sites)
    empty!(supercritical_sites)
    for (ix,iy) in sc
        check(pile,supercritical_sites,ix,iy)
        if ix>1 check(pile,supercritical_sites,ix-1,iy) end
        if ix<pile.L check(pile,supercritical_sites,ix+1,iy) end
        if iy>1 check(pile,supercritical_sites,ix,iy-1) end
        if iy<pile.L check(pile,supercritical_sites,ix,iy+1) end
    end
end
