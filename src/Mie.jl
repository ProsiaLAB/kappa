module Mie

export deRooij1984

function deRooij1984(
    rad::Float64,
    lam::Float64,
    e1::Float64,
    e2::Float64,
    Csca::Float64,
    Cext::Float64,
    F11::Vector{Float64},
    F12::Vector{Float64},
    F33::Vector{Float64},
    F34::Vector{Float64},
    nangle::Int,
)
    nparts = 1
    develop = 0
    delta = 1e-8
    cutoff = 1e-8
    thmin = 180.0 * (1.0 - 0.5) / nangle
    thmax = 180.0 * (nangle - 0.5) / nangle
    step = (thmax - thmin) / (nangle - 1)
    wavel = lam
    Rem = e1
    fImm = e2
    nsubr = 1
    ngaur = 1
    idis = 0
    par1 = rad
    par2 = 0.0
    par3 = 0.0
    ndis = 1
    rdis = zeros(Float64, (1, 300))
    rdis[1, 1] = rad
    nwrdis = zeros(Float64, (1, 300))
    nwrdis[1, 1] = 1.0

    mie()

    for i = 1:nangle
        F11[i] = F[1, nangle - i + 1]
        F12[i] = F[2, nangle - i + 1]
        F33[i] = F[3, nangle - i + 1]
        F34[i] = F[4, nangle - i + 1]
    end
end

function mie()
    
end


end # module