module Fractal

using Printf


function mg_mixing(refrel::Complex{Float64}, f1::Float64)
    eps_1 = refrel * refrel
    eps_2 = complex(1.0)
    mg = eps_2 * (
        2.0 * f1 * (eps_1 - eps_2) + eps_1 + 2.0 * eps_2
    ) / (eps_1 + 2.0 * eps_2 - f1 * (eps_1 - eps_2))
    mgav = sqrt(mg)
    return mgav
end


"""
    lorenz_mie(x::Float64, refrel::Complex{Float64})

Calculate Lorenz-Mie scattering coefficients (an,bn) for a monomer particle.

Since monomer's size parameter is supposed not to be very large, 
We use simple Bohren & Huffman Mie algorithm is used. 
The original BHMIE code is taken from Bruce Draine's HP:
      https://www.astro.princeton.edu/~draine/scattering.html
although we have slightly modified it.
"""
function lorenz_mie(x::Float64, refrel::Complex{Float64})
    nmxx::Int64 = 150000

    y = refrel * x
    ymod = abs(y)
    xstop = x + 4.0 * x^(1.0 / 3.0) + 2.0
    nmx = round(max(xstop, ymod) + 15.0)
    if nmx > nmxx
        error("nmx > nmxx")
    end

    # Calculate logarithmic derivative D_n(mx)
    # by downward recurrence. Initial value is set as D(mx) = 0+0i at n=nmx
    d = zeros(Complex{Float64}, nmx)
    for n = 1:nmx-1
        en = nmx - n + 1
        enr = real(nmx - n + 1)
        d[nmx-n] = (enr / y) - (1.0 / (d[en] + enr / y))
    end

    psi0 = cos(x)
    psi1 = sin(x)
    chi0 = -sin(x)
    chi1 = cos(x)
    xi1 = complex(psi1, -chi1)

    a = zeros(Complex{Float64}, nmx)
    b = zeros(Complex{Float64}, nmx)

    for n = 1:nstop
        nr = real(n)
        # Calcualte psi and chi via upward recurrence
        psi = (2.0 * nr - 1.0) * psi1 / x - psi0
        chi = (2.0 * nr - 1.0) * chi1 / x - chi0
        xi = complex(psi, -chi)
        # Calculate the Lorenz-Mie coefficients an and bn
        a[n] = (d[n] / refrel + nr / x) * psi - psi1
        a[n] = a[n] / ((d[n] / refrel + nr / x) * xi - xi1)
        b[n] = (refrel * d[n] + nr / x) * psi - psi1
        b[n] = b[n] / ((refrel * d[n] + nr / x) * xi - xi1)
        # Prepare for the next iteration
        psi0 = psi1
        psi1 = psi
        chi0 = chi1
        chi1 = chi
        xi1 = complex(psi1, -chi1)
    end

    return a, b
end


"""
    mean_scat_t(lmd, R0, PN, df, k0, refrel, iqsca, iqcor, iqgeo, nang, iquiet)

Calculate the mean scattering properties of a fractal aggregate.

"""
function mean_scat_t(
    lmd::Float64,
    R0::Float64,
    PN::Float64,
    df::Float64,
    k0::Float64,
    refrel::Complex{Float64},
    iqsca::Int64,
    iqcor::Int64,
    iqgeo::Int64,
    nang::Int64,
    iquiet::Bool
)
    k = 2π / lmd
    Rg = R0 * (PN / k0)^(1.0 / df)
    Rc = sqrt(5.0 / 3.0) * Rg
    xg = k * Rg
    x0 = k * R0

    xstop = x0 + 4.0 * x0^(1.0 / 3.0) + 2.0
    nstop = round(xstop)
    numax = nstop
    nmax = nstop

    if !iquiet
        @printf("%f :: Wavelength (μm) \n", lmd)
        @printf("%f :: Monomer radius (μm) \n", R0)
        @printf("%f :: Radius of gyration (μm) \n", Rg)
        @printf("%f :: Characteristic radius (μm) \n", Rc)
        @printf("%f :: Size parameter of monomer (μm) \n", x0)
        @printf("%f :: Size parameter of an aggregate (μm) \n", xg)
        @printf("%8d :: Expansion order of the scattering field \n", nstop)
    end

    if iqsca != 1 && iqsca != 2 && iqsca != 3
        error("Method must be 1, 2 or 3")
    end

    if iqcor != 1 && iqcor != 2 && iqcor != 3
        error("Correlation function must be 1, 2 or 3")
    end

    if iqgeo != 1 && iqgeo != 2 && iqgeo != 3
        error("Geometric cross section must be 1, 2 or 3")
    end

    if nang <= 1
        error("Number of angles must be greater than 1")
    end

    if PN < 1.0
        error("Number of monomers must be greater than 1")
    end

    if df > 3.0
        error("Fractal dimension must be greater than 1")
    end

    if (numax + nmax) >= 500 && !iquiet
        println("WARNING: The truncation order of monomer's scattered light")
        println("         field exceeds the maximum value (=500).          ")
        println("         This may cause a code crush at computations of   ")
        println("         the spherical Bessel function of the first kind. ")
    end

    ff = PN * (R0 / Rc)^3.0
    mgmref = mg_mixing(refrel, ff)
    dphic = 2.0 * k * Rc * abs(mgmref - 1.0)
    dphi0 = 2.0 * x0 * abs(refrel - 1.0)
    dphi = max(dphic, dphi0)

    if dphi >= 1.0 && !iquiet
        println("WARNING: The phase shift by an aggregate exceeds unity.")
        println("         Output of scattering matrix elements are not  ")
        println("         physically guaranteed.")
    end


    ad = zeros(Complex{Float64}, 2, nstop)
    dd = zeros(Complex{Float64}, 2, nstop)

    an, bn = lorenz_mie(x0, refrel)




    c_ext = 0
    c_sca = 0
    c_abs = 0
    g_asym = 0
    dphi = 0
    angs = range(0, stop=2π, length=nang)
    smat = zeros(Float64, nang, 4)
    phase_function = zeros(Float64, nang)




    return c_ext, c_sca, c_abs, g_asym, dphi, angs, smat, phase_function
end




mean_scat_t(0.5, 0.1, 1.0, 2.0, 1.0, 1.0 + 0.0im, 1, 1, 1, 10, false)


end # module