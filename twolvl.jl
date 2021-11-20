#
# Metropolis MC of a two-level system
#

using Random;
rng = MersenneTwister(211);

#
# Define the energy levels
#
E=zeros(2)
E[1]=0.;
E[2]=1.;

#
# MC simulation of a single system
#
function TwoLevelSingle!(Erun,T)
    state=2
    Erun[1]=E[state]
    p = exp(-(E[2]-E[1])/T)
    for i in 2:length(Erun)
        state = state==2 ? 1 :
            ( rand(rng) < p ? 2 : 1 )
        Erun[i] = E[state]
    end
end

#
# Simulate Nsamp independent two-level systems at temperature T
# Returns average energy en Eav.  Initialize Eav to an array
# of the desired length of the runs (will be overwritten)
function TwoLevelRun!(Eav,Nsamp,T)
    E=zeros(Float64,2);
    len=length(Eav);
    Erun=zeros(Float64,len);
    for n in 1:Nsamp
        TwoLevelSingle!(Erun,T)
        Eav.+=Erun
    end
    Eav./=Nsamp
end

#
# Equilibrium energy at temperature T
#
function EqE(T)
    return (E[1]*exp(-E[1]/T)+E[2]*exp(-E[2]/T)) / (exp(-E[1]/T)+exp(-E[2]/T))
end

#
# Define set of temperatures to simulate, legnth of each run
# and number of samples
#
Temp=[10.,5.,1.,.5];
len=100;
Nsamp=10000;

#
# Run simulations and plot
using CairoMakie
GLMakie.activate!(type="svg")

fig=Figure()
ax=Axis(fig[1,1])
for T in Temp
    Eav=zeros(Float64,len);
    TwoLevelRun!(Eav,Nsamp,T);
    t=range(0,length(Eav)-1,step=1);
    lines!(ax,t,Eav,label="T = $T");
    lines!(ax,[t[1],t[end]],[EqE(T),EqE(T)]);
end
ax.title="Energy vs time, $Nsamp samples";
ax.xlabel="t";
ax.ylabel="E";

# Plot the exponential envelope
for T in Temp
    t=range(0,len-1,step=1);
    E∞=EqE(T);
    tau=T/( (E[2]-E[1]) );
    lines!(ax,t,E∞ .+ (E[2]-E∞).*exp.(-t/tau),linestyle=:dot);
end

axislegend();
fig
