import Random
import Distributions
import Roots
using Roots
using Statistics

NFaces = 6;
seed=3232;

"""
    MEprobability(λ,face,Fval)

Compute the probability of finding face `face` of a dice according to
the maxent distributions (exponential) with parameter λ and obtaining
the value of each face using function  'Fval'.

Returns vector with probability of each face
"""
function MEprobability(λ,face,Fval)
    Fs=Fval.(1:NFaces);
    P=exp.(-λ* Fs);
    P/=sum(P);
    return P[face]
end

"""
    MEface_ave(λ,Fval)

Compute the average face value according to the maxent
distribuiton and `Fval` for face value function.
"""
function MEface_ave(λ,Fval)
    Fs=Fval.(1:NFaces);
    P=exp.(-λ* Fs);
    P/=sum(P);
    Fav=sum(Fs.*P);
    return Fav
    end
  

"""
   experiment_face_ave(Nsamp,seed,Fval,ddis)

Compute the average face value in a simulated experiment of `Nsamp`
dice rolls using `Fval` as the face value function and the discrete
random number distribution `ddis`
"""
function experiment_face_ave(Nsamp,Fval,ddis)
    f=rand(ddis,Nsamp)
    fav=sum(Fval.(f))/Nsamp
    return fav
    end

"""
    lambdadist(Nsamp,Fval,ddis)

Returns an array of 10000 values of λ each obtained by fixing the
average of the maxent distribution to the simulated experimental
average from `Nsamp` samples (obtained by calling
`experiment_face_ave()`.  `Fval` is the face-value function.

"""
function lambdadist(Nsamp,Fval,ddis)
    λ=zeros(10000)
    for i=1:10000
        av=experiment_face_ave(Nsamp,Fval,ddis);
        λ[i]=fzero(λ -> MEface_ave(λ,Fval)-av,0);
    end
    return λ
end

# Face value functions
Fval_linear = n -> n;
Fval_quadratic = n -> n^2;

"""
Obtain λ from maxent and plot
"""
function do_maxent(FacesProb)
    rng = Random.MersenneTwister(seed);
    ddis = Distributions.DiscreteNonParametric(1:NFaces,FacesProb)
    fig = Figure(res=1200);
    ax1 = Axis(fig[1,1]);
    barplot!(ax1,1:NFaces,FacesProb,label="Actual probability")

    println("Linear faces, exact face average ",sum(FacesProb.*Fval_linear.(1:NFaces)))
    l1=lambdadist(10,Fval_linear,ddis);
    λ=mean(l1)
    println("ME average ",MEface_ave(λ,Fval_linear));
    s=sqrt(var(l1))
    PP=MEprobability.(λ,1:NFaces,Fval_linear)
    Pmin=MEprobability.(λ-s,1:NFaces,Fval_linear)
    Pmax=MEprobability.(λ+s,1:NFaces,Fval_linear)
    err=abs.(Pmin.-Pmax)
    scatter!(ax1,1:NFaces,PP,label="ME prob with linear faces")
    errorbars!(ax1,1:NFaces,PP,err,color=:cyan)

    l1=lambdadist(10,Fval_quadratic,ddis);
    λ=mean(l1)
    println("Quadratic faces, exact face average ",sum(FacesProb.*Fval_quadratic.(1:NFaces)))
    println("Quadratic faces, ME average ",MEface_ave(λ,Fval_quadratic));
    s=sqrt(var(l1))
    PP=MEprobability.(λ,1:NFaces,Fval_quadratic)
    Pmin=MEprobability.(λ-s,1:NFaces,Fval_quadratic)
    Pmax=MEprobability.(λ+s,1:NFaces,Fval_quadratic)
    err=abs.(Pmin.-Pmax)
    scatter!(ax1,1:NFaces,PP,label="ME prob with quadratic faces")
    errorbars!(ax1,1:NFaces,PP,err,color=:yellow)

    axislegend();
    return fig
end


#FacesProb=[1,1,1,1,1,1];
FacesProb=[1,1,2,2,1,1];
#FacesProb=[2,2,1,1,2,2];
#FacesProb=[1,2,3,4,5,6];
FacesProb/=sum(FacesProb);

using CairoMakie
CairoMakie.activate!(type="svg");
fig=do_maxent(FacesProb);
save("ME.svg",fig)

