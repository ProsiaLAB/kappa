use std::error::Error;
use std::f64::consts::PI;

use ndarray::s;
use ndarray::Array1;
use ndarray::Array2;
use statrs::function::gamma::ln_gamma;

pub fn de_rooij_1984(nangle: usize, lam: f64, f_11: Vec<f64>) {
    let nparts = 1;
    let develop = 0;
    let delta = 1e-8;
    let cutoff = 1e-8;

    let thmin = 180.0 * (1.0 - 0.5) / nangle as f64;
    let thmax = 180.0 * (nangle as f64 - 0.5) / nangle as f64;
    let step = (thmax - thmin) / (nangle as f64 - 1.0);

    mie();
}

fn mie() {
    todo!("Mie scattering");
}

fn get_integration_bounds(
    idis: usize,
    p1: f64,
    p2: f64,
    p3: f64,
) -> Result<(f64, f64), Box<dyn Error>> {
    let mut r: Vec<f64> = vec![0.0; 1];

    let sef: f64;
    let ref_0: f64;
    let rref: f64;

    match idis {
        0 => return Ok((p1, p2)),
        1 => {
            sef = 1.0 / (p2 + 3.0).sqrt();
            ref_0 = 1.0 / (sef * sef * p2);
            rref = ref_0;
        }
        2 => {
            ref_0 = p1;
            sef = p2.sqrt();
            rref = ref_0;
        }
        3 => {
            sef = p3.sqrt();
            ref_0 = p1.max(p2) + sef;
            rref = p1.min(p2);
        }
        4 => {
            sef = (p2.ln().powi(2).exp() - 1.0).sqrt();
            ref_0 = p1 * (1.0 + sef * sef).powf(0.4);
            rref = ref_0;
        }
        5 => {
            ref_0 = p1;
            sef = ref_0.sqrt();
            rref = ref_0;
        }
        6 => {
            return Ok((p2, p3));
        }
        7 => {
            ref_0 = p2;
            sef = 2.0 * ref_0;
            rref = 0.5 * ref_0;
        }
        8 => {
            ref_0 = (p1 / (p2 * p3)).powf(p3);
            sef = 2.0 * ref_0;
            rref = 0.5 * ref_0;
        }
        9 => {
            return Ok((p1, p2));
        }
        _ => return Err("Illegal size distribution index".into()),
    }

    r[0] = ref_0 + sef;
    let r0 = ref_0;

    get_size_distribution(idis, p1, p2);
}

fn get_size_distribution(
    idis: usize,
    p1: f64,
    p2: f64,
    p3: f64,
    r: Array1<f64>,
    iparts: usize,
    rdis: Array1<f64>,
    mut nwithr: Array1<f64>,
) {
    let alpha_0: f64;
    let alpha_1: f64;
    let b0: f64;
    let b1: f64;
    let b2: f64;
    let c0: f64;
    let log_c0: f64;
    let log_c1: f64;
    let log_c2: f64;
    let flogrg: f64;
    let flogsi: f64;
    let fac: f64;
    let rg: f64;
    let aperg: f64;

    let ndis_max = 300;
    let ndpart = 4;
    let nwrdis: Array2<f64> = Array2::zeros((ndpart, ndis_max));
    let nwrdisp: Array1<f64> = Array1::zeros(ndis_max);

    let rdis: Array2<f64> = Array2::zeros((ndpart, ndis_max));
    let rdisp: Array1<f64> = Array1::zeros(ndis_max);

    let y2: Array1<f64> = Array1::zeros(ndis_max);

    match idis {
        0 => {
            return;
        }
        1 => {
            // Two parameter gamms with alpha and b given
            alpha_0 = p1;
            b0 = p2;
            alpha_1 = alpha_0 + 1.0;
            log_c0 = alpha_1 * b0.ln() - ln_gamma(alpha_1);
            nwithr = (log_c0 + alpha_0 * r.ln() - b0 * r).exp();
        }
        2 => {
            // Two parameter gamma with `p1 = reff` and `p2 = veff` given
            alpha_0 = 1.0 / p2 - 3.0;
            b0 = 1.0 / (p1 * p2);
            alpha_1 = alpha_0 + 1.0;
            log_c0 = alpha_1 * b0.ln() - ln_gamma(alpha_1);
            nwithr = (log_c0 + alpha_0 * &r.ln() - b0 * &r).exp();
        }
        3 => {
            alpha_0 = 1.0 / p3 - 3.0;
            b1 = 1.0 / (p1 * p3);
            b2 = 1.0 / (p2 * p3);
            let gamlna = ln_gamma(alpha_0 + 1.0);
            log_c1 = (alpha_0 + 1.0) * b1.ln() - gamlna;
            log_c2 = (alpha_0 + 1.0) * b2.ln() - gamlna;
            nwithr = 0.5
                * ((log_c1 + alpha_0 * &r.ln() - b1 * &r).exp()
                    + (log_c2 + alpha_0 * &r.ln() - b2 * &r).exp());
        }
        4 => {
            flogrg = p1.ln();
            flogsi = p2.ln().abs();
            c0 = 1.0 / ((2.0 * PI).sqrt() * flogsi);
            fac = -0.5 / (flogsi * flogsi);
            nwithr = c0 * (fac * (r.ln() - flogrg).powi(2)).exp() / r;
        }
        5 => {
            rg = p1 / (1.0 + p2).powf(2.5);
            flogrg = rg.ln();
            flogsi = (1.0 + p2).ln().sqrt();
            c0 = 1.0 / ((2.0 * PI).sqrt() * flogsi);
            fac = -0.5 / (flogsi * flogsi);
            nwithr = c0 * (fac * (r.ln() - flogrg).powi(2)).exp() / r;
        }
        6 => {
            alpha_0 = p1;
            let rmin = p2;
            let rmax = p3;
            if (alpha_0 + 1.0).abs() < 1e-10 {
                c0 = 1.0 / (rmax / rmin).ln();
            } else {
                alpha_1 = alpha_0 - 1.0;
                c0 = alpha_1 * rmax.powf(alpha_1) / ((rmax / rmin).powf(alpha_1) - 1.0);
            }
            nwithr = r.mapv(|ri| {
                if (ri < rmax) && (ri > rmin) {
                    c0 * ri.powf(alpha_0)
                } else {
                    0.0
                }
            })
        }
        7 => {
            alpha_0 = p1;
            let rc = p2;
            let gamma = p3;
            b0 = alpha_0 / (gamma * rc.powf(gamma));
            aperg = (alpha_0 + 1.0) / gamma;
            log_c0 = gamma.ln() + aperg * b0.ln() - ln_gamma(aperg);
            nwithr = (log_c0 + alpha_0 * r.ln() - b0 * r.powf(gamma)).exp();
        }
        8 => {
            alpha_0 = p1;
            b0 = p2;
            let gamma = p3;
            aperg = (alpha_0 + 1.0) / gamma;
            log_c0 = gamma.ln() + aperg * b0.ln() - ln_gamma(aperg);
            nwithr = (log_c0 + alpha_0 * r.ln() - b0 * r.powf(gamma)).exp();
        }
        9 => {
            let rdisp: Vec<f64> = vec![0.0; ndis];
            let nwrdisp: Vec<f64> = vec![0.0; ndis];
            for k in 0..ndis {}
        }
        _ => {
            return;
        }
    }
}

/// Given the arrays `xv` and `yv` of lengh `n` containing a tabulated
/// function, i.e., `yv[i] = f(xv[i])`, with `xv[0] < xv[1] < ... < xv[n-1]`,
/// and given the first derivative values `yp1` and `ypn` for the first
/// derivative of the interpolating function at points `xv[0]` and `xv[n-1]`,
/// this function returns the second derivative values `y2` of length `n`
/// which contains the second derivatives of the interpolating function at
/// the tabulated points `xv[i]`.
///
/// If `yp1` and `ypn` are equal or larger than `0.99e99`, the function will
/// set the first derivatives at the boundaries to be zero.
///
/// Reference:
/// Numerical Recipes: The Art of Scientific Computing, 3rd ed.
/// Press et al., 2007
fn spline(xv: Vec<f64>, yv: Vec<f64>, n: usize, yp1: f64, ypn: f64) -> Vec<f64> {
    let mut y2: Vec<f64> = vec![0.0; n];
    let mut u: Vec<f64> = vec![0.0; n - 1];

    // The lower boundary condition is set either to be “natural”
    // or else to have a specified first derivative.
    if yp1 > 0.99e99 {
        y2[0] = 0.0;
        u[0] = 0.0;
    } else {
        y2[0] = -0.5;
        u[0] = (3.0 / (xv[1] - xv[0])) * ((yv[1] - yv[0]) / (xv[1] - xv[0]) - yp1);
    }

    // This is the decomposition loop of the tridiagonal algorithm.
    // `y2` and `u` are used for temporary storage of the decomposed factors.
    for i in 1..n - 1 {
        let sig = (xv[i] - xv[i - 1]) / (xv[i + 1] - xv[i - 1]);
        let p = sig * y2[i - 1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        u[i] =
            (yv[i + 1] - yv[i]) / (xv[i + 1] - xv[i]) - (yv[i] - yv[i - 1]) / (xv[i] - xv[i - 1]);
        u[i] = (6.0 * u[i] / (xv[i + 1] - xv[i - 1]) - sig * u[i - 1]) / p;
    }

    let qn: f64;
    let un: f64;

    // The upper boundary condition is set either to be “natural”
    // or else to have a specified first derivative.
    if ypn > 0.99e99 {
        qn = 0.5;
        un = 0.0;
    } else {
        qn = 0.5;
        un = (3.0 / (xv[n - 1] - xv[n - 2]))
            * (ypn - (yv[n - 1] - yv[n - 2]) / (xv[n - 1] - xv[n - 2]));
    }

    y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);

    // This is the backsubstitution loop of the tridiagonal algorithm.
    for k in (0..n - 1).rev() {
        y2[k] = y2[k] * y2[k + 1] + u[k];
    }

    y2
}

/// Given the arrays `xa` and `ya` which tabulate a function (with
/// the `xa[i]` in order), and given the array `y2a`, which is the
/// output from `spline`, and given a value of `x`, this function
/// returns a cubic-spline interpolated value `y`.
fn splint(
    xa: Vec<f64>,
    ya: Vec<f64>,
    y2a: Vec<f64>,
    n: usize,
    x: f64,
    nn: usize,
    ndmui: usize,
) -> Result<f64, Box<dyn Error>> {
    if nn != ndmui {
        return Err("Dimension mismatch in splint".into());
    }

    let mut klo = 0;
    let mut khi = n - 1;

    let mut k: usize;

    while khi - klo > 1 {
        k = (khi + klo) / 2;
        if xa[k] > x {
            khi = k;
        } else {
            klo = k;
        }
    }

    let h = xa[khi] - xa[klo];
    if h.abs() < 1e-10 {
        return Err("Bad xa input to splint".into());
    }

    let a = (xa[khi] - x) / h;
    let b = (x - xa[klo]) / h;
    let y = a * ya[klo]
        + b * ya[khi]
        + ((a.powi(3) - a) * y2a[klo] + (b.powi(3) - b) * y2a[khi]) * (h * h) / 6.0;

    Ok(y)
}
