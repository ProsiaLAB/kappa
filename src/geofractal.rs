use std::error::Error;
use std::f64::consts::PI;

use statrs::function::gamma::{gamma, gamma_ur};

/// Computing geometrical cross-section of randomly oriented
/// fractal dust aggregates based on a statistical distribution
/// model of monomer particles developed by Tazaki 2021+.
///
/// # Arguments
///
/// * `iqapp` - Select a method for solving radial and angular integration.
/// when computing the mean overlapping efficiency.
///
///         --------------------------------------------------------
///                   |    radial  (x)        |  angular (u)        |
///         ---------------------------------------------------------
///         iqapp = 1 |   numerical           |  numerical          |
///         iqapp = 2 |   numerical  ( D<=2 ) |  analytical         |
///                   |   analytical ( D> 2 ) |                     |
///         --------------------------------------------------------
///
pub fn get_geometric_cross_section(
    iqapp: i32,
    iqcon: i32,
    iqcor: i32,
    pn: f64,
    k0: f64,
    df: f64,
) -> Result<f64, Box<dyn Error>> {
    let g: f64;

    if pn <= 0.9999 {
        return Err("The number of monomers should be greater than 1".into());
    }

    if df < 0.9999 || df > 3.0001 {
        return Err("The fractal dimension should be between 1 and 3".into());
    }

    if iqcon != 1 && iqcon != 2 {
        return Err("iqcon should be 1 or 2".into());
    }

    if iqcor != 1 && iqcor != 2 && iqcor != 3 {
        return Err("iqcor should be 1, 2 or 3".into());
    }

    if iqapp != 1 && iqapp != 2 {
        return Err("iqapp should be 1 or 2".into());
    }

    let mut a = 1.0;

    if iqcon == 2 {
        let mut pnth = 11.0 * df - 8.5;
        pnth = pnth.min(8.0);
        if pn < pnth {
            g = gfit(pn);
            return Ok(g);
        }
        let sigmath = mean_overlap_efficiency(iqapp, iqcor, k0, df, pnth);
        a = (1.0 + (pnth - 1.0) * sigmath) * gfit(pnth);
    }

    let sigma = mean_overlap_efficiency(iqapp, iqcor, k0, df, pn);

    g = a / (1.0 + (pn - 1.0) * sigma);

    Ok(g)
}

/// Minato et al. (2006)'s fitting formula for `N < 16`
fn gfit(pn: f64) -> f64 {
    12.5 * pn.powf(-0.315) * (-2.53 / pn.powf(0.0920)).exp()
}

/// Compute the mean overlapping efficiency
///
/// This subroutine calculations the mean overlapping efficiency: sigma
///
/// For `iqcor = 1` (Gaussian cut-off model)
///     
///                    xmin       /ln(xmax)
///      sigma  = --------------  |   dlnx frad(x)*S(rho),
///               16*Gamma(df/2)  /ln(xmin)
///     
///            xmin = df*(k0/PN)^(2/df); a=(df-2)/2; rho = sqrt(x/xmin)
///     
/// For iqcor=2 [Exponential cut-off model]
///     
///                   xmin^2      /ln(xmax)
///      sigma  = --------------  |   dlnx frad(x)*S(rho),
///                16*Gamma(df)   /ln(xmin)
///     
///            xmin = 2*sqrt(df(df+1)/2)*(k0/PN)^(1.0/df); a=df-2; rho = x/xmin
///     
/// For iqcor=3 [Fractal dimension cut-off model]
///     
///               xmin^{2/df} /ln(xmax)
///      sigma  = ----------  |  dlnx frad(x)*S(rho),
///                  16      /ln(xmin)
///     
///            xmin = 2^(df-1)*k0/PN; a=(df-2)/df; rho = (x/xmin)^(1/df)
///     
/// where df is fractal dimension, R0 is the monomer radius, Rg is
/// the radius of gyration, frad(x) is the radial integrand function,
/// and S(rho) is the angular integral function.
///
/// frad(x)dlnx = x^{a}*exp(-x) dlnx = x^{a-1}exp(-x) dx
///     
///                16*rho^2  / (1/rho)
///      S(rho)  = --------  |  du fang(rho,u),
///                   pi     /0
///     
///     fang(rho,u) = [Arcsin{sqrt(1-rho^2u^2)} - rho*u*sqrt{1-rho^2u^2}]*u/sqrt(1-u^2)
///
/// Note that the angular integral function S(rho) ≈1 when rho >> 1.
fn mean_overlap_efficiency(iqapp: i32, iqcor: i32, k0: f64, df: f64, pn: f64) -> f64 {
    let xmax = 25.0;
    let nx = 1000;

    let aicgm: f64;
    let xmin: f64;
    let factor: f64;

    let mut sigma: f64;

    match iqcor {
        1 => {
            aicgm = 0.5 * (df - 2.0);
            xmin = df * (k0 / pn).powf(2.0 / df);
            factor = xmin / 16.0 / gamma(0.5 * df);
        }
        2 => {
            aicgm = df - 2.0;
            xmin = 2.0 * (0.5 * df * (df + 1.0)).sqrt() * (k0 / pn).powf(1.0 / df);
            factor = xmin * xmin / 16.0 / gamma(df);
        }
        3 => {
            aicgm = (df - 2.0) / df;
            xmin = 2.0f64.powf(df - 1.0) * k0 / pn;
            factor = xmin.powf(2.0 / df) / 16.0;
        }
        _ => unreachable!(),
    }

    if df > 2.0 && iqapp == 3 {
        sigma = factor * gamma(aicgm) * gamma_ur(aicgm, xmin);
    } else {
        let dlnx = (xmax - xmin).ln() / (nx - 1) as f64;
        let mut x: Vec<f64> = vec![0.0; nx];
        let mut sang: Vec<f64> = vec![1.0; nx];
        sigma = 0.0;
        for i in 0..nx {
            x[i] = (xmin.ln() + i as f64 * dlnx).exp();
            if iqapp == 1 {
                sang[i] = angular_integration(iqcor, x[i], xmin, df);
            }
        }
        for i in 0..nx - 1 {
            sigma += 0.5 * (frad(aicgm, x[i]) * sang[i] + frad(aicgm, x[i + 1]) * sang[i + 1]);
        }
        sigma *= factor * dlnx * factor;
    }

    sigma
}

/// Radial intgrand function in `ln(x)` space
fn frad(a: f64, x: f64) -> f64 {
    x.powf(a) * (-x).exp()
}

fn angular_integration(iqcor: i32, x: f64, xmin: f64, df: f64) -> f64 {
    let nmax_u = 1000;

    let mut sang: f64;

    let rho: f64;
    let mut u: Vec<f64> = vec![0.0; nmax_u];

    match iqcor {
        1 => {
            rho = (x / xmin).sqrt();
        }
        2 => {
            rho = x / xmin;
        }
        3 => {
            rho = (x / xmin).powf(1.0 / df);
        }
        _ => unreachable!(),
    }

    sang = 0.0;
    let umin = 0.0;
    let umax = 1.0 / rho;

    let du = (umax - umin) / (nmax_u - 1) as f64;
    for j in 0..nmax_u {
        u[j] = umin + du * j as f64;
    }
    for j in 0..nmax_u - 1 {
        sang += 0.5 * (fang(rho, u[j]) + fang(rho, u[j + 1])) * du;
    }
    sang = 16.0 * rho * rho * sang / PI;

    sang
}

/// This function calculates the integrand of the angular integral.
fn fang(rho: f64, u: f64) -> f64 {
    if (1.0 - rho * u).abs() <= 1e-10 {
        return 0.0;
    } else {
        ((1.0 - rho * rho * u * u).sqrt().asin() - rho * u * (1.0 - rho * rho * u * u).sqrt()) * u
            / (1.0 - u * u).sqrt()
    }
}
