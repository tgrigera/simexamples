"""
   Metropolis MC of a two-level system

Call `do_runs()` to simulate and plot for several temperatures,
or see `twolevel_single!` and `twolevel_run!`
"""
module TwoLevel

import Random

rng = Random.MersenneTwister(211);

#
# Define the energy levels
#
E=zeros(2)
E[1]=0.;
E[2]=1.;

"""
    twolevel_single!(Erun,T)

MC simulation of a two-level system at temperature `T`.  Initial state
is the excited state, and energy for each step is written to the array
`Erun`.  The initial size of `Erun` givis the munber of steps.
"""
function twolevel_single!(Erun::Array,T::Float64)
    state=2
    Erun[1]=E[state]
    p = exp(-(E[2]-E[1])/T)
    for i in 2:length(Erun)
        state = state==2 ? 1 :
            ( rand(rng) < p ? 2 : 1 )
        Erun[i] = E[state]
    end
end

"""
     twolevel_run!(Eav,Nsamp,T)

 Do `Nsamp` independent simulations of two-level systems at
temperature `T` (calling `twolevel_single!`).  Returns average energy
in `Eav`, which must be initialized to a Float64 array of the desired
length of the runs.
 """
function twolevel_run!(Eav,Nsamp,T)
    E=zeros(Float64,2);
    len=length(Eav);
    Erun=zeros(Float64,len);
    for n in 1:Nsamp
        twolevel_single!(Erun,T);
        Eav.+=Erun;
    end
    Eav./=Nsamp;
end


" Equilibrium energy at temperature T"
function EqE(T)
    return (E[1]*exp(-E[1]/T)+E[2]*exp(-E[2]/T)) / (exp(-E[1]/T)+exp(-E[2]/T))
end


using CairoMakie

function do_runs()

    #
    # Define set of temperatures to simulate, legnth of each run
    # and number of samples
    #
    Temp=[10.,5.,1.,.5];
    len=100;
    Nsamp=10000;

    #
    # Run simulations and plot
    CairoMakie.activate!(type="svg")

    fig=Figure()
    ax=Axis(fig[1,1])
    for T in Temp
        Eav=zeros(Float64,len);
        twolevel_run!(Eav,Nsamp,T);
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
end

end # module
