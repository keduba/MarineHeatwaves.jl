# Expected output show/print vs display
# EventHW{T} / EventCS{T} ::MEvents{T}
Base.show(io::IO, m::MEvents) = print(io, m.maxes)
Base.show(io::IO, ::MIME"text/plain", m::MEvents) = print(io, "$(typeof(m)):\n ", m)

# MExtreme => MCS || MHW
function Base.show(io::IO, ::MIME"text/plain", m::MExtreme)
    pst = split(summary(m.temp), " ")[1]
    summary(io, m)
    print(io, " {", pst, "} :\n  ")
    # kinda works for a vector. 
    # needs to be modified for a matrix
    joy = String.(propertynames(m))
    padd = maximum([length(s) for s in joy])

    for fd in propertynames(m)
        print(io, rpad(String(fd), padd), ": ", getfield(m, fd), "\n  ") # right or left padding
    end
end

# temp = x X y, summary(m.temp)
# For. the print/show version.
#=
function Base.show(io::IO, m::MExtreme)
    pre = split(summary(m), "{")[1]
    pst = split(summary(m), ",")[2]
    print(io, pre, "{")
    summary(io, m.temp)
    print(io, ",", pst, "}")
    for f in propertynames(m)
        if f == :edtype
            print(io, getfield(m, f), ")")
        else
            print(io, getfield(m, f), ", " )
        end
    end
end
=#

# MWrapper (M[CH]Wrapper)
function Base.show(io::IO, ::MIME"text/plain", m::MWrapper)
    summary(io, m)
    print(io, ":\n  ")
    # kinda works for a vector. 
    # needs to be modified for a matrix
    joy = String.(propertynames(m))
    padd = maximum([length(s) for s in joy])

    for fd in propertynames(m)
        print(io, rpad(String(fd), padd), ": ", getfield(m, fd), "\n  ") # right or left padding
    end
end

# Events
function Base.show(io::IO, ::MIME"text/plain", m::Events)
    summary(io, m)
    print(io, ":\n  ")
    # kinda works for a vector.
    # needs to be modified for a matrix
    joy = String.(propertynames(m))
    padd = maximum([length(s) for s in joy])
    lp = length(propertynames(m))
    nm = propertynames(m)
    for i in 1:lp-1
        print(io, rpad(String(nm[i]), padd+1))  # right or left padding
        summary(io, getfield(m, nm[i]))
        print(io, ":\n", getfield(m, nm[i]))
        print(io, i < lp-1 ? "\n\n" : "")
    end
end

# MHCMetrics
function Base.show(io::IO, ::MIME"text/plain", m::MHCMetrics)
    summary(io, m)
    print(io, ":\n  ")
    # kinda works for a vector. 
    # needs to be modified for a matrix
    joy = String.(propertynames(m))
    padd = maximum([length(s) for s in joy])
    for fd in propertynames(m)
        print(io, rpad(String(fd), padd), ": ", getfield(m, fd), "\n  ") # right or left padding
    end
end

# MHCMetricsm
function Base.show(io::IO, ::MIME"text/plain", m::MHCMetricsm)
    summary(io, m)
    print(io, ":\n")
    # kinda works for a vector. 
    # needs to be modified for a matrix
    joy = String.(propertynames(m))
    padd = maximum([length(s) for s in joy])
    lp = length(propertynames(m))
    for (i, fd) in enumerate(propertynames(m))
        print(io, rpad(String(fd), padd+1))  # right or left padding
        summary(io, getfield(m, fd))
        print(io, ":\n", getfield(m, fd))
        print(io, i < lp ? "\n\n" : "")
    end
end


# MHWCSOutarray
function Base.show(io::IO, ::MIME"text/plain", m::MHWCSO)
    summary(io, m)
    print(io, ":\n  ")
    # kinda works for a vector. 
    # needs to be modified for a matrix
    joy = String.(propertynames(m))
    padd = maximum([length(s) for s in joy])
    lp = length(propertynames(m))
    for (i, fd) in enumerate(propertynames(m))
        print(io, rpad(String(fd), padd+1))  # right or left padding
        summary(io, getfield(m, fd))
        print(io, ":\n", getfield(m, fd))
        print(io, i < lp ? "\n\n" : "")
    end
end
