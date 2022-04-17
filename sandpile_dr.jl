#
# sandpile_dr.jl --- An example driver for the Sandpiles model
#

# Add the location of Sandpiles.jl to the path so that
# the using statement below will find the module
#
push!(LOAD_PATH,homedir()*"/software/simexamples")

using Sandpiles
using Random

# Choose a RNG and seed
rng = Random.MersenneTwister(211)

# Set the critical height.  Note that sites overflow when
# height is GREATER than zc (some authors use greater than or
# equal to)
zc=3
# Set size (this simulation runs the 2-d sandpile in an LxL square lattice)
L=40

# Create the sandpile and history
pile=Sandpile(L,zc,rng)
his=Sandpile_history(pile)
# Evolve for some time
run(pile,100000,his)
# Start recording with a clean history for 
his_eq=Sandpile_history(pile)
run(pile,400000,his_eq)

using CairoMakie
CairoMakie.activate!()

# Plot evolution

fig1 = Figure()
ax1 = fig1[1,1] = Axis(fig1,xlabel="t",ylabel="<z>")
lines!(ax1,vcat(his.time,his_eq.time),vcat(his.zav,his_eq.zav))
display(fig1)

include("Histograms.jl")
histog=Histogram(100,0,500)
map(x->push!(histog,x),his_eq.avalanche_duration)
d,p=prob(histog)

# Plot histograms of duration, energy, and size

fig2 = Figure()
ax2 = fig2[1,1] = Axis(fig2,xlabel="log duration",ylabel="log P(duration)")
lines!(ax2,log10.(d),log10.(p))
display(fig2)

histog=Histogram(100,0,10000)
map(x->push!(histog,x),his_eq.avalanche_energy)
d,p=prob(histog)
fig3 = Figure()
ax3 = fig3[1,1] = Axis(fig3,xlabel="log energy",ylabel="log P(energy)")
lines!(ax3,log10.(d),log10.(p))
display(fig3)

histog=Histogram(100,0,15)
map(x->push!(histog,x),his_eq.avalanche_size)
d,p=prob(histog)
fig4 = Figure()
ax4 = fig4[1,1] = Axis(fig4,xlabel="log size",ylabel="log P(size)")
lines!(ax4,log10.(d),log10.(p))
display(fig4)
