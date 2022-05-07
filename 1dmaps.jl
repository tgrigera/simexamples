#
# 1dmaps.jl --- Exploring some 1-d dynamical maps
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

# Define the maps

"One step of the logistic map"
lmap(x;r=1.) = r*x*(1-x)
"One step of the tent map"
tmap(x;r=1.) = x<0.5 ? r*x : r*(1-x)
"One step of the sine map"
smap(x;r=1.) = r*sin(π*x)

"""
    domap(map,x0,r;nsteps=100)

Return the time evolution of `map` during
`nsteps` steps starting from `x0`
"""
function domap(map,x0,r;nsteps=100)
    x=zeros(Float64,nsteps)
    x[1]=x0
    for i=2:nsteps x[i]=map(x[i-1],r=r) end
    return x
end

"""
    do_cobweb(map,x0,r;nsteps=100)

Return points drawing the Verhust or cobweb diagram
of `map`
"""
function do_cobweb(map,x0,r)
    nsteps=100
    xt=domap(map,x0,r,nsteps=nsteps)
    x=zeros(2*nsteps-1)
    y=zeros(2*nsteps-1)
    x[1]=xt[1]
    y[1]=0
    for i=1:length(xt)-1
        x[2i]=xt[i]
        y[2i]=xt[i+1]
        x[2i+1]=xt[i+1]
        y[2i+1]=xt[i+1]
    end
    return x,y
end

"Plot the cobweb for `map`in axis `ax`"
function plot_cobweb(map,x0,r,ax)
    x=collect(0:0.1:1)
    lines!(ax,x,x)
    lines!(x,lmap.(x,r=r))
    x,y=do_cobweb(map,x0,r)
    scatterlines!(x,y)
    return fig
end

using Random

function bifurcation_map(bmap;rmin,rmax,xmin,xmax,npoints,eqsteps)
    bddata=[]
    rs=range(rmin,rmax,length=npoints)
    for r in rs
        x=xmin+rand()*(xmax-xmin)
        for _=1:eqsteps x=bmap(x,r=r)  end
        xt=domap(bmap,x,r;nsteps=npoints)
        push!(bddata,map.(x->(r,x),xt)...)
    end
    return bddata
end

using GLMakie
GLMakie.activate!()

#
# Logistic map
#

# Plot evoluion

map_params=[ (0.8,1.) , (0.1,3), (0.75,3.8) ]

fig=Figure()
ax=fig[1,1]=Axis(fig,title="Logistic map trajectories")
for (x0,r) in map_params
    xt=domap(lmap,x0,r)
    scatterlines!(ax,xt)
end
save("./lm_dyn.jpg",fig)

# Plot cobweb

fig=Figure()
ax=fig[1,1]=Axis(fig,title="Logistic map cobweb")
plot_cobweb(lmap,map_params[1]...,ax)
plot_cobweb(lmap,map_params[2]...,ax)
plot_cobweb(lmap,map_params[3]...,ax)
save("./lm_cobweb.jpg",fig)

# Bifurcation diagram

bddata=bifurcation_map(lmap,rmin=3,rmax=4,xmin=0,xmax=1,npoints=1000,eqsteps=1000)
fig=Figure()
ax=fig[1,1]=Axis(fig,xlabel="r",ylabel="x",title="Logistic map bifurcation diagram")
scatter!(ax,first.(bddata),last.(bddata),markersize=0.1)

# Inset with zoom
bddatain=bifurcation_map(lmap,rmin=3.84,rmax=3.86,xmin=0,xmax=1,npoints=500,eqsteps=100000)
inax=fig[1,1]=Axis(fig,width=Relative(0.4),height=Relative(0.4),halign=0.1,valign=0.1,limits=(nothing,nothing,0.4,0.6))
scatter!(inax,first.(bddatain),last.(bddatain),markersize=0.1)

save("./lm_bifdiag.jpg",fig)


#
# Bifurcation diagram for other maps
#

bddata=bifurcation_map(tmap,rmin=0,rmax=2,xmin=0,xmax=1,npoints=1000,eqsteps=1000)
fig=Figure()
ax=fig[1,1]=Axis(fig,xlabel="r",ylabel="x",title="Tent map bifurcation diagram")
scatter!(ax,first.(bddata),last.(bddata),markersize=0.1)

# # Inset with zoom
# bddatain=bifurcation_map(lmap,rmin=3.84,rmax=3.87,xmin=0.4,xmax=0.6,npoints=500,eqsteps=10000)
# inax=fig[1,1]=Axis(fig,width=Relative(0.4),height=Relative(0.4),halign=0.1,valign=0.1,limits=(nothing,nothing,0.4,0.6))
# scatter!(inax,first.(bddatain),last.(bddatain),markersize=0.1)

save("./tm_bifdiag.jpg",fig)

bddata=bifurcation_map(smap,rmin=0,rmax=1,xmin=0,xmax=1,npoints=1000,eqsteps=1000)
fig=Figure()
ax=fig[1,1]=Axis(fig,xlabel="r",ylabel="x",title="Sine map bifurcation diagram")
scatter!(ax,first.(bddata),last.(bddata),markersize=0.1)

# # Inset with zoom
# bddatain=bifurcation_map(lmap,rmin=3.84,rmax=3.87,xmin=0.4,xmax=0.6,npoints=500,eqsteps=10000)
# inax=fig[1,1]=Axis(fig,width=Relative(0.4),height=Relative(0.4),halign=0.1,valign=0.1,limits=(nothing,nothing,0.4,0.6))
# scatter!(inax,first.(bddatain),last.(bddatain),markersize=0.1)

save("./sm_bifdiag.jpg",fig)
