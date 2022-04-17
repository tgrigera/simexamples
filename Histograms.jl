# twolvl.jl
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

mutable struct Histogram
    min::Float64
    max::Float64
    nbins::Int
    delta::Float64
    ndata::Int
    out_below::Int
    out_above::Int
    counts::Vector{BigInt}
end

Histogram(nbins,min,max) = Histogram(min,max,nbins,(max-min)/nbins,0,0,0,zeros(BigInt,nbins))

import Base.push!

function push!(his::Histogram,datum)
    his.ndata+=1
    bin=floor(Int,(datum-his.min)/his.delta)+1;
    if bin<1 his.out_below+=1; return -1
    elseif bin>his.nbins his.out_above+=1; return -1
    else his.counts[bin]+=1 ; return bin end
end
    
outliers(his::Histogram)=his.out_below+his.out_above

binc(his::Histogram,bin)=his.min + (bin+0.5)*his.delta

prob(his::Histogram,bin)= his.counts[bin]/(his.ndata*his.delta)

area(his::Histogram)=(his.ndata-outliers(his))/his.ndata;

function median(his::Histogram)
    cc=his.out_below;
    m=1
    while cc<his.ndata/2
        cc+=his.counts[m]
        m+=1
    end
    return binc(his,m)
end

function prob(his::Histogram)
    x=collect(his.min+his.delta/2:his.delta:his.max)
    f = 1/(his.ndata*his.delta)
    return x,f*his.counts
end

function counts(his::Histogram)
    x=collect(his.min:his.delta:his.max)
    return x,his.counts
end
