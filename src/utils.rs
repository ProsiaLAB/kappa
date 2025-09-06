use extensions::types::RVector;
use ndarray::Array1;

use crate::opac::KappaConfig;

pub mod constants {
    //! Defines mathematical expressions commonly used when computing distribution
    //! values as constants
    //!
    //! This module is directly copied from the `statrs` crate.
    //! The original code is licensed under the MIT license.
    //! It was purely done to avoid a dependency on the `statrs` crate.
    //!
    //! The original code can be found here:
    //! <https://github.com/statrs-dev/statrs/blob/master/src/constants.rs>
    #![allow(clippy::excessive_precision)]

    /// Constant value for `ln(pi)`
    pub const LN_PI: f64 = 1.1447298858494001741434273513530587116472948129153;

    /// Constant value for `ln(2 * sqrt(e / pi))`
    pub const LN_2_SQRT_E_OVER_PI: f64 = 0.6207822376352452223455184457816472122518527279025978;

    /// Constant value for `2 * sqrt(e / pi)`
    pub const TWO_SQRT_E_OVER_PI: f64 = 1.8603827342052657173362492472666631120594218414085755;

    /// Default accuracy for `f64`, equivalent to `0.0 * F64_PREC`
    pub const DEFAULT_F64_ACC: f64 = 0.0000000000000011102230246251565;
}
pub mod spline {
    use anyhow::bail;
    use anyhow::Result;
    use extensions::types::RVector;
    use ndarray::Array1;

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
    /// # Reference
    /// - Numerical Recipes: The Art of Scientific Computing, 3rd ed.
    ///   Press et al., 2007
    pub fn spline(xv: &[f64], yv: &[f64], n: usize, yp1: f64, ypn: f64) -> RVector {
        let mut y2: RVector = Array1::zeros(n);
        let mut u: RVector = Array1::zeros(n - 1);

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
            u[i] = (yv[i + 1] - yv[i]) / (xv[i + 1] - xv[i])
                - (yv[i] - yv[i - 1]) / (xv[i] - xv[i - 1]);
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
    ///
    /// # Reference
    /// - Numerical Recipes: The Art of Scientific Computing, 3rd ed.
    ///   Press et al., 2007
    pub fn splint(xa: &[f64], ya: &[f64], y2a: &[f64], n: usize, x: f64) -> Result<f64> {
        let mut klo = 0;
        let mut khi = n - 1;

        while khi - klo > 1 {
            let k = (khi + klo) / 2;
            if xa[k] > x {
                khi = k;
            } else {
                klo = k;
            }
        }

        let h = xa[khi] - xa[klo];
        if h == 0.0 {
            bail!("Bad xv input to routine splint");
        }

        let a = (xa[khi] - x) / h;
        let b = (x - xa[klo]) / h;

        let y = a * ya[klo]
            + b * ya[khi]
            + ((a.powi(3) - a) * y2a[klo] + (b.powi(3) - b) * y2a[khi]) * (h.powi(2) / 6.0);

        Ok(y)
    }
}

pub mod legendre {
    use std::f64::consts::PI;

    use anyhow::bail;
    use anyhow::Result;
    use extensions::types::RVector;
    use ndarray::Array1;

    pub fn lpmns(m: usize, n: usize, x: f64) -> Result<(RVector, RVector)> {
        let mut pm: RVector = Array1::zeros(n);
        let mut pd: RVector = Array1::zeros(n);

        let mf = m as f64;

        if x.abs() == 1.0 {
            for k in 0..n {
                let kf = k as f64;
                match m {
                    0 => {
                        pm[k] = 1.0;
                        pd[k] = 0.5 * kf * (kf + 1.0);
                        if x < 0.0 {
                            pm[k] *= -1.0f64.powf(kf);
                            pd[k] *= -1.0f64.powf(kf + 1.0);
                        }
                    }
                    1 => {
                        pm[k] = x;
                        pd[k] = 0.5 * kf * (kf + 1.0) * x;
                        if x < 0.0 {
                            pm[k] *= -1.0f64.powf(kf);
                            pd[k] *= -1.0f64.powf(kf + 1.0);
                        }
                    }
                    2 => {
                        pm[k] = 0.5 * (3.0 * x.powi(2) - 1.0);
                        pd[k] = 0.5 * kf * (kf + 1.0) * (3.0 * x.powi(2) - 1.0);
                        if x < 0.0 {
                            pm[k] *= -1.0f64.powf(kf);
                            pd[k] *= -1.0f64.powf(kf + 1.0);
                        }
                    }
                    _ => {
                        bail!("m must be 0, 1, or 2");
                    }
                }
            }
        }

        let x0 = (1.0 - x * x).abs();
        let mut pm0 = 1.0;
        let mut pmk = pm0;
        for k in 0..m {
            let kf = k as f64;
            pmk = (2.0 * kf + 1.0) * x0.sqrt() * pm0;
            pm0 = pmk;
        }
        let mut pm1 = (2.0 * m as f64 + 1.0) * x * pm0;
        pm[m] = pmk;
        pm[m + 1] = pm1;

        for (k, val) in pm.iter_mut().enumerate().skip(m + 2).take(n - m - 2) {
            let kf = k as f64;
            *val = ((2.0 * kf - 1.0) * x * pm1 - (kf + mf - 1.0) * pmk) / (kf - mf);
            pmk = pm1;
            pm1 = *val;
        }

        pd[0] = ((1.0 - mf) * pm[1] - x * pm[0]) / (x * x - 1.0);

        for (k, val) in pd.iter_mut().enumerate().skip(1).take(n - 1) {
            let kf = k as f64;
            *val = (kf * x * pm[k] - (kf + mf) * pm[k - 1]) / (x * x - 1.0);
        }

        Ok((pm, pd))
    }

    pub fn lpn(n: usize, x: f64) -> (RVector, RVector) {
        let mut pn: RVector = Array1::zeros(n);
        let mut pd: RVector = Array1::zeros(n);

        pn[0] = 1.0;
        pn[1] = x;
        pd[0] = 0.0;
        pd[1] = 1.0;

        let mut p0 = 1.0;
        let mut p1 = x;

        for k in 2..n {
            let kf = k as f64;
            let pf = (2.0 * kf - 1.0) * kf * x * p1 - (kf - 1.0) / kf * p0;
            pn[k] = pf;
            if x.abs() == 1.0 {
                pd[k] = 0.5 * x.powf(kf + 1.0) * kf * (kf + 1.0);
            } else {
                pd[k] = kf * (p1 - x * pf) / (1.0 - x * x);
            }
            p0 = p1;
            p1 = pf;
        }
        (pn, pd)
    }

    pub fn gauss_legendre(x1: f64, x2: f64, n: usize) -> (RVector, RVector) {
        let eps = 1.0e-14;

        let mut pp = 0.0;
        let m = n.div_ceil(2);
        let xm = 0.5 * (x2 + x1);
        let xl = 0.5 * (x2 - x1);

        let mut x: RVector = Array1::zeros(n);
        let mut w: RVector = Array1::zeros(n);

        let nf = n as f64;

        for i in 0..m {
            let ir = i as f64;
            let mut z = (PI * (ir + 0.75) / (nf + 0.5)).cos();
            let mut z1 = 2.0 * z;
            while (z - z1).abs() > eps {
                let mut p1 = 1.0;
                let mut p2 = 0.0;
                for j in 0..n {
                    let jf = j as f64;
                    let p3 = p2;
                    p2 = p1;
                    p1 = ((2.0 * jf + 1.0) * z * p2 - jf * p3) / (jf + 1.0);
                }
                pp = nf * (z * p1 - p2) / (z * z - 1.0);
                z1 = z;
                z = z1 - p1 / pp;
            }
            x[i] = xm - xl * z;
            x[n - 1 - i] = xm + xl * z;
            w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
            w[n - 1 - i] = w[i];
        }
        (x, w)
    }
}

pub mod bessel {
    use anyhow::bail;
    use anyhow::Result;
    use extensions::types::{CVector, RVector};
    use ndarray::Array1;
    use num_complex::{Complex, Complex64};

    use crate::fractal::FractalGeometry;
    use crate::fractal::{FractalConfig, FractalCutoff};

    use super::special::two_point_correlation;

    pub enum BoundaryCondition {
        Jablonski,
        TazakiTanaka,
    }

    impl BoundaryCondition {
        fn coefficients(self) -> ([f64; 4], [f64; 4]) {
            match self {
                BoundaryCondition::Jablonski => {
                    let a = [-3.6693021e-8, -3.1745158e-5, 2.1567720e-2, 9.9123380e-1];
                    let b = [-9.9921351e-8, 8.7822303e-6, 1.0238752e-2, 3.7588265];
                    (a, b)
                }
                BoundaryCondition::TazakiTanaka => {
                    let a = [
                        1.69496268177237e-08,
                        -2.43299782942114e-05,
                        0.0158750501131321,
                        1.00672154148706,
                    ];
                    let b = [
                        2.49355951047228e-08,
                        -2.9387731648675e-05,
                        0.0135005554796179,
                        3.72312019844119,
                    ];
                    (a, b)
                }
            }
        }
    }

    ///  This subroutine performs integration of S_p(kRg) (Equation 31):
    ///  $$
    ///         S_p(k R_g) = \frac{\pi^2}{k^3} \int_0^{\infty} du \, u \, J_{p+1/2}(u) \, H_{p+1/2}^{(1)}(u) \, g(u/k)
    ///  $$   
    ///  where g(u) is the two-point correlation function (Equations 18 and 19):
    ///  $$
    ///         g(u) = \frac{1}{4 \pi R_g^3} \left( \frac{u}{R_g} \right)^{d_f - 3} \cdot fc \left( \frac{u}{R_g} \right)
    ///  $$
    ///  where fc is the cut-off function. By substituting g(u) into S_p, we have
    ///  $$
    ///         S_p(k R_g) = \frac{\pi}{4 X_g^{d_f}} \int_{u_{\min}}^{u_{\max}} du \, u^{d_f - 2} J_{p+1/2}(u) H_{p+1/2}^{(1)}(u) fc(u / X_g)
    ///  $$
    ///
    ///  where the integration range is approximated by the range `[u_min,u_max]`.
    ///  By using the spherical Bessel j_p(u) and Hankel functions of 1st kind h_p^(1)(u),
    ///  the Bessel function and the Hankel function are rewritten by
    ///  $$
    ///  J_{p+1/2} (u) = \sqrt{\frac{2u}{\pi}} \, j_p(u)
    ///  $$
    ///  
    ///  $$
    ///  H_{p+1/2}^{(1)}(u) = \sqrt{\frac{2u}{\pi}} \, h_p^{(1)}(u)
    ///  $$
    ///  
    ///  We have:
    ///  
    ///  $$
    ///  S_p(k R_g) = \frac{1}{2 X_g^{d_f}} \int_{u_{\min}}^{u_{\max}} du \, u^{d_f - 1} j_p(u) h_p^{(1)}(u) fc(u / X_g)
    ///  $$
    ///  
    ///  For the unitary condition of the two-point correlation function:
    ///  
    ///  $$
    ///  1 = \int_0^{\infty} dw \, 4\pi w^2 g(w)
    ///  $$
    ///  
    ///  If we take the integration variable as \( w = u/k \), then we obtain:
    ///  
    ///  $$
    ///  1 = \frac{1}{X_g^{d_f}} \int_{u_{\min}}^{u_{\max}} du \, u^{d_f - 1} fc(u / X_g) \quad \text{.... Eq. (*)}
    ///  $$
    ///  
    ///  The integration range `[umin,umax]` is determined as follows.
    ///  The integrand of Equation (*) is
    ///
    ///         (u/xg)^{df-1}fc(u/xg) du = (u/xg) ^{d_f}fc(u/xg) dlnu
    ///
    ///  u_max is chosen so that fc(u/xg) ~ exp[-eta1].
    ///       For iqcor=1 (Gauss)
    ///               u_max ~ 2 * xg * sqrt(eta1 / d_f)
    ///       For iqcor=1 (Exponential)
    ///               u_max ~ xg * eta1 * sqrt(2.0 /(d_f*(d_f+1))
    ///       For iqcor=1 (FLDIM)
    ///               u_max ~ xg * sqrt( 2.0 * eta1 ) ** (1.0/d_f)
    ///  I adopt eta1 = 25.0.
    ///
    ///  u_min is chosen so that (u_min/xg)^{d_f} ~ exp[-eta2], thus,
    ///
    ///               umin ~ xg * exp(-eta2/d_f)
    ///
    ///  where eta2 = 40.0.
    ///
    /// # Remarks
    /// [`BoundaryCondition::TazakiTanaka`] is used for the boundary condition and is
    /// recommended; although [`BoundaryCondition::Jablonski`] is also available.
    pub fn int_sph_bessel(fracc: &FractalConfig, x_g: f64, p: usize) -> Result<Complex64> {
        let pf = p as f64;

        let floor_val = 1e-30;
        let eta1 = 25.0;
        let eta2 = 40.0;

        let nn = 10000; // Number of grid points in the numerical integration of S_p(kRg)
        let nnf = nn as f64;

        let umax = match fracc.cutoff {
            FractalCutoff::Gaussian => 2.0 * x_g * (eta1 / fracc.df).sqrt(),
            FractalCutoff::Exponential => x_g * eta1 * (2.0 / (fracc.df * (fracc.df + 1.0))).sqrt(),
            FractalCutoff::FractalDimension => x_g * (2.0 * eta1).powf(1.0 / fracc.df),
        };

        let umin = x_g * (-eta2 / fracc.df).exp();
        let du = (umax - umin).powf(1.0 / (nnf - 1.0));

        let mut u: RVector = Array1::zeros(nn);
        let mut intg: CVector = Array1::zeros(nn);
        let mut intg_unit: RVector = Array1::zeros(nn);

        for (n, val) in u.iter_mut().enumerate().take(nn) {
            *val = umin * du.powi(n as i32);
        }

        let (a, b) = BoundaryCondition::TazakiTanaka.coefficients();

        let lnxa = a[0] * pf.powf(3.0) + a[1] * pf.powf(2.0) + a[2] * pf + a[3];
        let lnxb = b[0] * pf.powf(3.0) + b[1] * pf.powf(2.0) + b[2] * pf + b[3];

        for n in 0..nn {
            let isol = if u[n] < lnxa {
                1
            } else if u[n] > lnxb {
                3
            } else {
                2
            };
            let (sj, sy) = sph_bessel(p, u[n], isol)?;
            let jp = sj[p];
            let yp = sy[p];
            let hp = Complex::new(jp, yp);

            if jp * 0.0 != 0.0 || jp.is_nan() {
                bail!("Error in sph_bessel (jp)");
            } else if yp * 0.0 != 0.0 || yp.is_nan() {
                bail!("Error in sph_bessel (yp)");
            }

            intg[n] = u[n].powf(fracc.df - 1.0) * jp * hp * two_point_correlation(fracc, u[n], x_g);
            intg_unit[n] = u[n].powf(fracc.df - 1.0) * two_point_correlation(fracc, u[n], x_g);
        }

        // Use iterators to apply the trapezoidal rule
        let wa: Complex64 = intg
            .windows(2)
            .into_iter()
            .zip(u.windows(2))
            .map(|(intg_pair, u_pair)| {
                0.5 * (intg_pair[0] + intg_pair[1]) * (u_pair[1] - u_pair[0])
            })
            .sum::<Complex64>();

        let mut unitary: f64 = intg_unit
            .windows(2)
            .into_iter()
            .zip(u.windows(2))
            .map(|(intg_unit_pair, u_pair)| {
                0.5 * (intg_unit_pair[0] + intg_unit_pair[1]) * (u_pair[1] - u_pair[0])
            })
            .sum();

        unitary /= x_g.powf(fracc.df);
        let mut sp = 0.5 * wa / x_g.powf(fracc.df);
        let error = (1.0 - unitary).abs();

        if sp.re.abs() < floor_val {
            sp = Complex::new(floor_val, sp.im);
        }
        if sp.im.abs() < floor_val {
            sp = Complex::new(sp.re, -floor_val);
        }

        if fracc.geometry != FractalGeometry::Tazaki && error > 1.0e-3 {
            bail!("Error in int_sph_bessel: error = {}", error);
        }

        Ok(sp)
    }

    fn sph_bessel(m: usize, x: f64, isol: usize) -> Result<(RVector, RVector)> {
        let imax = 100; // truncation order of the series expansion
        let nwarmup = 100; // number of warm-up iterations

        let floor_val = 1.0e-70;
        let ceiling_val = 1.0e+70;

        let mut sj: RVector = Array1::zeros(m + 1);
        sj[0] = x.sin() / x;

        let mut sy: RVector = Array1::zeros(m + 1);
        sy[0] = -x.cos() / x;

        if m == 0 {
            return Ok((sj, sy));
        }

        match isol {
            1 => {
                // Series expansion
                let mut f_n = 1.0;
                for (n, val) in sj.iter_mut().enumerate().skip(1).take(m) {
                    let nf = n as f64;
                    f_n *= x / (2.0 * nf + 1.0);
                    let mut xi = 1.0;
                    let mut wa = 0.0;
                    for i in 1..imax {
                        let ir = i as f64;
                        xi = -x * x * xi / (2.0 * ir * (2.0 * ir + 2.0 * nf + 1.0));
                        wa += xi;
                        if (xi / wa).abs() <= floor_val {
                            break;
                        }
                    }
                    *val = f_n * (1.0 + wa);
                    if val.abs() <= floor_val {
                        break;
                    }
                }
            }
            2 => {
                // Downward recursion
                let mut k1 = 0.0;
                let mut k0 = 1.0;
                let mut k: RVector = Array1::zeros(m + nwarmup + 1);
                for n in (m + nwarmup..0).rev() {
                    let nf = n as f64;
                    k[n] = -k1 + (2.0 * nf + 3.0) * k0 / x;
                    k1 = k0;
                    k0 = k[n];
                }
                let s = sj[0] / k[0];
                for (n, val) in sj.iter_mut().enumerate().skip(1).take(m) {
                    *val = s * k[n];
                    if val.abs() <= floor_val {
                        break;
                    }
                }
            }
            3 => {
                // Upward recursion
                sj[1] = x.sin() / (x * x) - x.cos() / x;
                let mut k0 = sj[0];
                let mut k1 = sj[1];
                for n in 1..m {
                    let nf = n as f64;
                    sj[n + 1] = (2.0 * nf + 1.0) * k1 / x - k0;
                    k0 = k1;
                    k1 = sj[n + 1];
                }
            }
            _ => {
                bail!("Invalid isol value");
            }
        }

        // Spherical Bessel function of the second kind
        sy[1] = -x.cos() / (x * x) - x.sin() / x;
        let mut y0 = sy[0];
        let mut y1 = sy[1];
        for (n, val) in sy.iter_mut().enumerate().skip(2).take(m) {
            let nf = n as f64;
            *val = (2.0 * nf - 1.0) * y1 / x - y0;
            y0 = y1;
            y1 = *val;
            if val.abs() >= ceiling_val {
                break;
            }
        }

        Ok((sj, sy))
    }
}

pub mod linalg {
    use anyhow::{bail, Result};
    use extensions::types::{CMatrix, CVector, RMatrix, RVector, UVector};
    use ndarray::{Array1, Array2};
    use num_complex::Complex64;

    fn lu_decomposition(a: &mut RMatrix, n: usize) -> Result<UVector> {
        let mut d = 1.0;
        let mut imax = 0;
        let eps = 1e-20;

        let mut vv: RVector = Array1::zeros(1000);
        let mut indx: UVector = Array1::zeros(n);

        let mut aamax;

        for i in 0..n {
            aamax = 0.0;
            for j in 0..n {
                if a[[i, j]].abs() > aamax {
                    aamax = a[[i, j]].abs();
                }
            }
            if aamax == 0.0 {
                bail!("ERROR: Singular matrix in routine lu_decomposition");
            }
            vv[i] = 1.0 / aamax;
        }

        for j in 0..n {
            for i in 0..j {
                let mut sum = a[[i, j]];
                for k in 0..i {
                    sum -= a[[i, k]] * a[[k, j]];
                }
                a[[i, j]] = sum;
            }
            aamax = 0.0;
            for i in j..n {
                let mut sum = a[[i, j]];
                for k in 0..j {
                    sum -= a[[i, k]] * a[[k, j]];
                }
                a[[i, j]] = sum;
                let dum = vv[i] * sum.abs();
                if dum >= aamax {
                    imax = i;
                    aamax = dum;
                }
            }
            if j != imax {
                for k in 0..n {
                    let dum = a[[imax, k]];
                    a[[imax, k]] = a[[j, k]];
                    a[[j, k]] = dum;
                }
                d = -d;
                vv[imax] = vv[j];
            }
            indx[j] = imax;
            if a[[j, j]] == 0.0 {
                a[[j, j]] = eps;
            }
            if j != n - 1 {
                let dum = 1.0 / a[[j, j]];
                for i in j + 1..n {
                    a[[i, j]] *= dum;
                }
            }
        }

        Ok(indx)
    }

    fn lu_backsubstitution(a: &RMatrix, b: &mut RVector, indx: &UVector, n: usize) {
        let mut ii = 0;
        for i in 0..n {
            let ll = indx[i];
            let mut sum = b[ll];
            b[ll] = b[i];
            if ii != 0 {
                for j in ii..i {
                    sum -= a[[i, j]] * b[j];
                }
            } else if sum != 0.0 {
                ii = i + 1;
            }
            b[i] = sum;
        }
        for i in (0..n).rev() {
            let mut sum = b[i];
            for j in i + 1..n {
                sum -= a[[i, j]] * b[j];
            }
            b[i] = sum / a[[i, i]];
        }
    }

    pub fn complex_matrix_inverse(
        n: usize,
        np: usize,
        p: &CVector,
        t: &CMatrix,
    ) -> Result<CVector> {
        let mut a: RMatrix = Array2::zeros((np, np));
        let mut b: RVector = Array1::zeros(np);

        let mut r: CVector = Array1::zeros(np);

        for i in 0..n {
            for j in 0..n {
                a[[i, j]] = t[[i, j]].re;
                a[[i, j + n]] = -t[[i, j]].im;
                a[[i + n, j]] = t[[i, j]].im;
                a[[i + n, j + n]] = t[[i, j]].re;
            }
            b[i] = p[i].re;
            b[i + n] = p[i].im;
        }

        // LU decomposition
        let indx = lu_decomposition(&mut a, np)?;
        // LU backsubstitution
        lu_backsubstitution(&a, &mut b, &indx, np);

        for i in 0..n {
            r[i] = Complex64::new(b[i], b[i + n]);
        }

        Ok(r)
    }
}

pub mod gamma {
    //! Gamma function related utilities
    //!
    //! This module is directly copied from the `statrs` crate.
    //! The original code is licensed under the MIT license.
    //! It was purely done to avoid a dependency on the `statrs` crate.
    //!
    //! The original code can be found here:
    //! <https://github.com/statrs-dev/statrs/blob/master/src/function/gamma.rs>

    #![allow(clippy::excessive_precision)]
    use std::f64::consts::{E, PI};

    use anyhow::anyhow;
    use anyhow::Result;
    use approx::ulps_eq;

    use super::constants::{DEFAULT_F64_ACC, LN_2_SQRT_E_OVER_PI, LN_PI, TWO_SQRT_E_OVER_PI};
    use super::special::almost_eq;

    /// Auxiliary variable when evaluating the `gamma_ln` function
    const GAMMA_R: f64 = 10.900511;

    /// Polynomial coefficients for approximating the `gamma_ln` function
    const GAMMA_DK: &[f64] = &[
        2.48574089138753565546e-5,
        1.05142378581721974210,
        -3.45687097222016235469,
        4.51227709466894823700,
        -2.98285225323576655721,
        1.05639711577126713077,
        -1.95428773191645869583e-1,
        1.70970543404441224307e-2,
        -5.71926117404305781283e-4,
        4.63399473359905636708e-6,
        -2.71994908488607703910e-9,
    ];

    /// Computes the logarithm of the gamma function
    /// with an accuracy of 16 floating point digits.
    /// The implementation is derived from
    /// "An Analysis of the Lanczos Gamma Approximation",
    /// Glendon Ralph Pugh, 2004 p. 116
    pub fn ln_gamma(x: f64) -> f64 {
        if x < 0.5 {
            let s = GAMMA_DK
                .iter()
                .enumerate()
                .skip(1)
                .fold(GAMMA_DK[0], |s, t| s + t.1 / (t.0 as f64 - x));

            LN_PI
                - (PI * x).sin().ln()
                - s.ln()
                - LN_2_SQRT_E_OVER_PI
                - (0.5 - x) * ((0.5 - x + GAMMA_R) / E).ln()
        } else {
            let s = GAMMA_DK
                .iter()
                .enumerate()
                .skip(1)
                .fold(GAMMA_DK[0], |s, t| s + t.1 / (x + t.0 as f64 - 1.0));

            s.ln() + LN_2_SQRT_E_OVER_PI + (x - 0.5) * ((x - 0.5 + GAMMA_R) / E).ln()
        }
    }

    /// Computes the gamma function with an accuracy
    /// of 16 floating point digits. The implementation
    /// is derived from "An Analysis of the Lanczos Gamma Approximation",
    /// Glendon Ralph Pugh, 2004 p. 116
    pub fn gamma(x: f64) -> f64 {
        if x < 0.5 {
            let s = GAMMA_DK
                .iter()
                .enumerate()
                .skip(1)
                .fold(GAMMA_DK[0], |s, t| s + t.1 / (t.0 as f64 - x));

            PI / ((PI * x).sin() * s * TWO_SQRT_E_OVER_PI * ((0.5 - x + GAMMA_R) / E).powf(0.5 - x))
        } else {
            let s = GAMMA_DK
                .iter()
                .enumerate()
                .skip(1)
                .fold(GAMMA_DK[0], |s, t| s + t.1 / (x + t.0 as f64 - 1.0));

            s * TWO_SQRT_E_OVER_PI * ((x - 0.5 + GAMMA_R) / E).powf(x - 0.5)
        }
    }

    /// Computes the upper incomplete gamma function
    /// `Gamma(a,x) = int(exp(-t)t^(a-1), t=0..x) for a > 0, x > 0`
    /// where `a` is the argument for the gamma function and
    /// `x` is the lower intergral limit.
    ///
    /// # Panics
    ///
    /// if `a` or `x` are not in `(0, +inf)`
    pub fn gamma_ui(a: f64, x: f64) -> f64 {
        checked_gamma_ui(a, x).unwrap()
    }

    /// Computes the upper incomplete gamma function
    /// `Gamma(a,x) = int(exp(-t)t^(a-1), t=0..x) for a > 0, x > 0`
    /// where `a` is the argument for the gamma function and
    /// `x` is the lower intergral limit.
    ///
    /// # Errors
    ///
    /// if `a` or `x` are not in `(0, +inf)`
    pub fn checked_gamma_ui(a: f64, x: f64) -> Result<f64> {
        checked_gamma_ur(a, x).map(|x| x * gamma(a))
    }

    /// Computes the lower incomplete gamma function
    /// `gamma(a,x) = int(exp(-t)t^(a-1), t=0..x) for a > 0, x > 0`
    /// where `a` is the argument for the gamma function and `x`
    /// is the upper integral limit.
    ///
    ///
    /// # Panics
    ///
    /// if `a` or `x` are not in `(0, +inf)`
    pub fn gamma_li(a: f64, x: f64) -> f64 {
        checked_gamma_li(a, x).unwrap()
    }

    /// Computes the lower incomplete gamma function
    /// `gamma(a,x) = int(exp(-t)t^(a-1), t=0..x) for a > 0, x > 0`
    /// where `a` is the argument for the gamma function and `x`
    /// is the upper integral limit.
    ///
    ///
    /// # Errors
    ///
    /// if `a` or `x` are not in `(0, +inf)`
    pub fn checked_gamma_li(a: f64, x: f64) -> Result<f64> {
        checked_gamma_lr(a, x).map(|x| x * gamma(a))
    }

    /// Computes the upper incomplete regularized gamma function
    /// `Q(a,x) = 1 / Gamma(a) * int(exp(-t)t^(a-1), t=0..x) for a > 0, x > 0`
    /// where `a` is the argument for the gamma function and
    /// `x` is the lower integral limit.
    ///
    /// # Remarks
    ///
    /// Returns `f64::NAN` if either argument is `f64::NAN`
    ///
    /// # Panics
    ///
    /// if `a` or `x` are not in `(0, +inf)`
    pub fn gamma_ur(a: f64, x: f64) -> f64 {
        checked_gamma_ur(a, x).unwrap()
    }

    /// Computes the upper incomplete regularized gamma function
    /// `Q(a,x) = 1 / Gamma(a) * int(exp(-t)t^(a-1), t=0..x) for a > 0, x > 0`
    /// where `a` is the argument for the gamma function and
    /// `x` is the lower integral limit.
    ///
    /// # Remarks
    ///
    /// Returns `f64::NAN` if either argument is `f64::NAN`
    ///
    /// # Errors
    ///
    /// if `a` or `x` are not in `(0, +inf)`
    pub fn checked_gamma_ur(a: f64, x: f64) -> Result<f64> {
        if a.is_nan() || x.is_nan() {
            return Ok(f64::NAN);
        }
        if a <= 0.0 || a == f64::INFINITY {
            return Err(anyhow!("AInvalid"));
        }
        if x <= 0.0 || x == f64::INFINITY {
            return Err(anyhow!("XInvalid"));
        }

        let eps = 0.000000000000001;
        let big = 4503599627370496.0;
        let big_inv = 2.22044604925031308085e-16;

        if x < 1.0 || x <= a {
            return Ok(1.0 - gamma_lr(a, x));
        }

        let mut ax = a * x.ln() - x - ln_gamma(a);
        if ax < -709.78271289338399 {
            return if a < x { Ok(0.0) } else { Ok(1.0) };
        }

        ax = ax.exp();
        let mut y = 1.0 - a;
        let mut z = x + y + 1.0;
        let mut c = 0.0;
        let mut pkm2 = 1.0;
        let mut qkm2 = x;
        let mut pkm1 = x + 1.0;
        let mut qkm1 = z * x;
        let mut ans = pkm1 / qkm1;
        loop {
            y += 1.0;
            z += 2.0;
            c += 1.0;
            let yc = y * c;
            let pk = pkm1 * z - pkm2 * yc;
            let qk = qkm1 * z - qkm2 * yc;

            pkm2 = pkm1;
            pkm1 = pk;
            qkm2 = qkm1;
            qkm1 = qk;

            if pk.abs() > big {
                pkm2 *= big_inv;
                pkm1 *= big_inv;
                qkm2 *= big_inv;
                qkm1 *= big_inv;
            }

            if qk != 0.0 {
                let r = pk / qk;
                let t = ((ans - r) / r).abs();
                ans = r;

                if t <= eps {
                    break;
                }
            }
        }
        Ok(ans * ax)
    }

    /// Computes the lower incomplete regularized gamma function
    /// `P(a,x) = 1 / Gamma(a) * int(exp(-t)t^(a-1), t=0..x) for real a > 0, x > 0`
    /// where `a` is the argument for the gamma function and `x` is the upper
    /// integral limit.
    ///
    /// # Remarks
    ///
    /// Returns `f64::NAN` if either argument is `f64::NAN`
    ///
    /// # Panics
    ///
    /// if `a` or `x` are not in `(0, +inf)`
    pub fn gamma_lr(a: f64, x: f64) -> f64 {
        checked_gamma_lr(a, x).unwrap()
    }

    /// Computes the lower incomplete regularized gamma function
    /// `P(a,x) = 1 / Gamma(a) * int(exp(-t)t^(a-1), t=0..x) for real a > 0, x > 0`
    /// where `a` is the argument for the gamma function and `x` is the upper
    /// integral limit.
    ///
    /// # Remarks
    ///
    /// Returns `f64::NAN` if either argument is `f64::NAN`
    ///
    /// # Errors
    ///
    /// if `a` or `x` are not in `(0, +inf)`
    pub fn checked_gamma_lr(a: f64, x: f64) -> Result<f64> {
        if a.is_nan() || x.is_nan() {
            return Ok(f64::NAN);
        }
        if a <= 0.0 || a == f64::INFINITY {
            return Err(anyhow!("AInvalid"));
        }
        if x <= 0.0 || x == f64::INFINITY {
            return Err(anyhow!("XInvalid"));
        }

        let eps = 0.000000000000001;
        let big = 4503599627370496.0;
        let big_inv = 2.22044604925031308085e-16;

        if almost_eq(a, 0.0, DEFAULT_F64_ACC) {
            return Ok(1.0);
        }
        if almost_eq(x, 0.0, DEFAULT_F64_ACC) {
            return Ok(0.0);
        }

        let ax = a * x.ln() - x - ln_gamma(a);
        if ax < -709.78271289338399 {
            if a < x {
                return Ok(1.0);
            }
            return Ok(0.0);
        }
        if x <= 1.0 || x <= a {
            let mut r2 = a;
            let mut c2 = 1.0;
            let mut ans2 = 1.0;
            loop {
                r2 += 1.0;
                c2 *= x / r2;
                ans2 += c2;

                if c2 / ans2 <= eps {
                    break;
                }
            }
            return Ok(ax.exp() * ans2 / a);
        }

        let mut y = 1.0 - a;
        let mut z = x + y + 1.0;
        let mut c = 0;

        let mut p3 = 1.0;
        let mut q3 = x;
        let mut p2 = x + 1.0;
        let mut q2 = z * x;
        let mut ans = p2 / q2;

        loop {
            y += 1.0;
            z += 2.0;
            c += 1;
            let yc = y * f64::from(c);

            let p = p2 * z - p3 * yc;
            let q = q2 * z - q3 * yc;

            p3 = p2;
            p2 = p;
            q3 = q2;
            q2 = q;

            if p.abs() > big {
                p3 *= big_inv;
                p2 *= big_inv;
                q3 *= big_inv;
                q2 *= big_inv;
            }

            if q != 0.0 {
                let nextans = p / q;
                let error = ((ans - nextans) / nextans).abs();
                ans = nextans;

                if error <= eps {
                    break;
                }
            }
        }
        Ok(1.0 - ax.exp() * ans)
    }

    /// Computes the Digamma function which is defined as the derivative of
    /// the log of the gamma function. The implementation is based on
    /// "Algorithm AS 103", Jose Bernardo, Applied Statistics, Volume 25, Number 3
    /// 1976, pages 315 - 317
    pub fn digamma(x: f64) -> f64 {
        let c = 12.0;
        let d1 = -0.57721566490153286;
        let d2 = 1.6449340668482264365;
        let s = 1e-6;
        let s3 = 1.0 / 12.0;
        let s4 = 1.0 / 120.0;
        let s5 = 1.0 / 252.0;
        let s6 = 1.0 / 240.0;
        let s7 = 1.0 / 132.0;

        if x == f64::NEG_INFINITY || x.is_nan() {
            return f64::NAN;
        }
        if x <= 0.0 && ulps_eq!(x.floor(), x) {
            return f64::NEG_INFINITY;
        }
        if x < 0.0 {
            return digamma(1.0 - x) + PI / (-PI * x).tan();
        }
        if x <= s {
            return d1 - 1.0 / x + d2 * x;
        }

        let mut result = 0.0;
        let mut z = x;
        while z < c {
            result -= 1.0 / z;
            z += 1.0;
        }

        if z >= c {
            let mut r = 1.0 / z;
            result += z.ln() - 0.5 * r;
            r *= r;

            result -= r * (s3 - r * (s4 - r * (s5 - r * (s6 - r * s7))));
        }
        result
    }

    pub fn inv_digamma(x: f64) -> f64 {
        if x.is_nan() {
            return f64::NAN;
        }
        if x == f64::NEG_INFINITY {
            return 0.0;
        }
        if x == f64::INFINITY {
            return f64::INFINITY;
        }
        let mut y = x.exp();
        let mut i = 1.0;
        while i > 1e-15 {
            y += i * signum(x - digamma(y));
            i /= 2.0;
        }
        y
    }

    // modified signum that returns 0.0 if x == 0.0. Used
    // by inv_digamma, may consider extracting into a public
    // method
    fn signum(x: f64) -> f64 {
        if x == 0.0 {
            0.0
        } else {
            x.signum()
        }
    }
}

pub mod special {
    use std::f64::consts::PI;

    use approx::AbsDiffEq;

    use crate::fractal::{FractalConfig, FractalCutoff};

    use super::gamma::gamma;

    /// Compares if two floats are close via `approx::abs_diff_eq`
    /// using a maximum absolute difference (epsilon) of `acc`.
    ///
    /// This is directly copied from the `statrs` crate.
    /// The original code is licensed under the MIT license.
    /// It was purely done to avoid a dependency on the `statrs` crate.
    ///
    /// The original code can be found here:
    /// <https://github.com/statrs-dev/statrs/blob/master/src/prec.rs>
    pub fn almost_eq(a: f64, b: f64, acc: f64) -> bool {
        if a.is_infinite() && b.is_infinite() {
            return a == b;
        }
        a.abs_diff_eq(&b, acc)
    }

    pub fn two_point_correlation(fracc: &FractalConfig, u: f64, x: f64) -> f64 {
        match fracc.cutoff {
            FractalCutoff::Gaussian => {
                let c = 0.25 * fracc.df;
                (2.0 * c.powf(0.5 * fracc.df) / gamma(0.5 * fracc.df)) * (-c * u * u / x / x).exp()
            }
            FractalCutoff::Exponential => {
                let c = (0.5 * fracc.df * (fracc.df + 1.0)).sqrt();
                (c.powf(fracc.df) / gamma(fracc.df)) * (-c * u / x).exp()
            }
            FractalCutoff::FractalDimension => {
                let c = 0.5;
                (c * fracc.df) * (-c * u / x).exp()
            }
        }
    }

    pub fn output_structure() {
        todo!()
    }

    pub fn confluent_hypergeometric(mut a: f64, b: f64, mut x: f64) -> f64 {
        let mut a0 = a;
        let x0 = x;
        let mut y0 = 0.0;
        let mut y1 = 0.0;
        let mut hg = 0.0;
        if b == 0.0 || b == -(b.floor().abs()) {
            hg = 1e300;
        } else if a == 0.0 || x == 0.0 {
            hg = 1.0;
        } else if a == -1.0 {
            hg = 1.0 - x / b;
        } else if a == b {
            hg = x.exp();
        } else if (a - b) == 1.0 {
            hg = (1.0 + x / b) * x.exp();
        } else if a == 1.0 && b == 2.0 {
            hg = (x.exp() - 1.0) * x.exp();
        } else if a == a.floor() && a < 0.0 {
            let m = -a.floor() as i32;
            let mut r = 1.0;
            hg = 1.0;
            for k in 0..m {
                let kf = k as f64;
                r *= (a + kf) / kf / (b + kf) * x;
                hg += r;
            }
        }
        if hg != 0.0 {
            return hg;
        }

        if x < 0.0 {
            a = b - a;
            a0 = a;
            x = x.abs();
        }

        let nl = if a < 2.0 { 0 } else { 1 };
        let la = if a >= 2.0 { a as usize } else { 0 };
        let laf = la as f64;
        if a >= 2.0 {
            a = a - laf - 1.0;
        }

        for n in 0..=nl {
            if a0 >= 2.0 {
                a += 1.0;
            }
            if x <= 30.0 + b.abs() || a < 0.0 {
                hg = 1.0;
                let mut rg = 1.0;
                for j in 0..500 {
                    let jf = j as f64;
                    rg *= (a + jf) / (jf * (b + jf)) * x;
                    hg += rg;
                    if (hg / rg).abs() < 1e-15 {
                        break;
                    }
                }
            } else {
                let ta = gamma(a);
                let tb = gamma(b);
                let xg = b - a;
                let tba = gamma(xg);
                let mut sum1 = 1.0;
                let mut sum2 = 1.0;
                let mut r1 = 1.0;
                let mut r2 = 1.0;
                for i in 0..8 {
                    let ir = i as f64;
                    r1 = -r1 * (a + ir) * (a - b + ir + 1.0) / (x * (ir + 1.0));
                    r2 = -r2 * (b - a) * (a - ir - 1.0) / (x * (ir + 1.0));
                    sum1 += r1;
                    sum2 += r2;
                }
                let hg1 = tb / tba * x.powf(-a) * (PI * a).cos() * sum1;
                let hg2 = tb / ta * x.exp() * x.powf(a - b) * sum2;
                hg = hg1 + hg2;
            }
            if n == 0 {
                y0 = hg;
            } else if n == 1 {
                y1 = hg;
            }
        }
        if a0 >= 2.0 {
            for _ in 0..(la - 2) {
                hg = ((2.0 * a - b + x) * y1 + (b - a) * y0) / a;
                y0 = y1;
                y1 = hg;
                a += 1.0;
            }
        }
        if x0 < 0.0 {
            hg *= x0.exp()
        }

        hg
    }

    pub fn optic_limit() {
        todo!()
    }
}

pub fn regrid_lnk_data(
    l0: &[f64],
    n0: &[f64],
    k0: &[f64],
    lam: &RVector,
    loglog: bool,
) -> (RVector, RVector) {
    let n = lam.len();
    let n_actual = l0.len();
    let mut x0 = l0[0];
    let mut y01 = n0[0];
    let mut y02 = k0[0];
    // let wp = (1.0 - y01) / x0.powi(2);
    // let gamma = y02 / x0.powi(3);

    // The first block is space in grid that is before the first specified
    // wavelength. Simple extrapolation: keep same value
    let mut e1: RVector = Array1::zeros(n);
    let mut e2: RVector = Array1::zeros(n);
    let mut i = 0;

    while x0 >= lam[i] {
        e1[i] = y01;
        e2[i] = y02;
        i += 1;
        if i >= n {
            // All requested lambda are before first data point
            break;
        }
    }
    // Main interpolation loop
    for i0 in 0..n_actual {
        let x1 = l0[i0];
        let y11 = n0[i0];
        let y12 = k0[i0];

        while lam[i] <= x1 && lam[i] > x0 {
            // log-log interpolation between points
            let lgy11 = y11.log10();
            let lgy12 = y12.log10();
            let lgy01 = y01.log10();
            let lgy02 = y02.log10();
            let lgx0 = x0.log10();
            let lgx1 = x1.log10();
            let lggr = lam[i].log10();

            e1[i] = 10f64.powf(lgy11 + (lggr - lgx1) * (lgy01 - lgy11) / (lgx0 - lgx1));
            e2[i] = 10f64.powf(lgy12 + (lggr - lgx1) * (lgy02 - lgy12) / (lgx0 - lgx1));

            i += 1;
            if i >= n {
                break;
            }
        }

        x0 = x1;
        y01 = y11;
        y02 = y12;
    }

    // Extrapolation to long wavelengths
    if loglog {
        // Use log-log extrapolation
        let slope1 = (n0[n_actual - 1] / n0[n_actual - 2]).log10()
            / (l0[n_actual - 1] / l0[n_actual - 2]).log10();
        let sect1 = n0[n_actual - 2].log10();
        let slope2 = (k0[n_actual - 1] / k0[n_actual - 2]).log10()
            / (l0[n_actual - 1] / l0[n_actual - 2]).log10();
        let sect2 = k0[n_actual - 2].log10();

        for j in i..n {
            let dist = (lam[j] / l0[n_actual - 2]).log10();
            e1[j] = 10f64.powf(slope1 * dist + sect1);
            e2[j] = 10f64.powf(slope2 * dist + sect2);
        }
    } else {
        // Use the dielectric extrapolation, this is the default
        for j in i..n {
            e1[j] = e1[i - 1];
            e2[j] = e2[i - 1] * lam[i - 1] / lam[j];
        }
    }

    (e1, e2)
}

pub fn prepare_sparse(kpc: &mut KappaConfig) {
    let f = kpc.lam[kpc.nlam - 1] / kpc.lam[kpc.nlam - 2];

    for (valmin, valmax) in kpc.scatlammin.iter_mut().zip(kpc.scatlammax.iter_mut()) {
        if *valmax / *valmin < f {
            *valmin /= f.sqrt() / 1.0001;
            *valmax *= f.sqrt() / 1.0001;
        }
    }

    for (il, &lam_val) in kpc.lam.iter().enumerate() {
        if kpc
            .scatlammin
            .iter()
            .zip(kpc.scatlammax.iter())
            .any(|(&min_val, &max_val)| lam_val >= min_val && lam_val <= max_val)
        {
            kpc.sparse_indices.insert(il);
            if il > 0 {
                kpc.sparse_indices.insert(il - 1);
            }
            if il + 1 < kpc.nlam {
                kpc.sparse_indices.insert(il + 1);
            }
        }
    }
}

/// Compute the moments of the size distribution
///
/// The results are returned in `ameans`, an array of length 3:
/// $$
///                [\langle a \rangle,\quad \langle a^2 \rangle^{1/2},\quad \langle a^3 \rangle^{1/3}]
/// $$
/// If both mn and sig are nonzero and the product is positive, we use
/// the log-normal size distribution.  If not, we use the powerlaw.
pub fn get_sizedis_moments(kpc: &KappaConfig) -> [f64; 3] {
    let ns = 1000;

    let aminlog = kpc.amin.log10();
    let amaxlog = kpc.amax.log10();
    let pow = -kpc.apow;

    if ((kpc.amax - kpc.amin) / kpc.amin).abs() < 1e-6 {
        [kpc.amin, kpc.amin, kpc.amin]
    } else {
        let mut tot = [0.0; 3];
        let mut totn = 0.0;
        let mut ameans = [0.0; 3];
        for is in 0..ns {
            let isf = is as f64;
            let nsf = ns as f64;
            let r = 10.0f64.powf(aminlog + (amaxlog - aminlog) * isf / nsf);
            let nr = if (kpc.amean * kpc.asigma).abs() > 0.0 {
                // normal or log-normal size distribution
                let expo = if kpc.asigma > 0.0 {
                    0.5 * ((r - kpc.amean) / kpc.asigma).powi(2)
                } else {
                    0.5 * ((r / kpc.amean).ln() / kpc.asigma).powi(2)
                };
                if expo > 99.0 {
                    0.0
                } else {
                    -expo.exp()
                }
            } else {
                r.powf(pow + 1.0)
            };
            totn += nr;
            tot[0] += nr * r;
            tot[1] += nr * r.powi(2);
            tot[2] += nr * r.powi(3);
        }
        if totn == 0.0 {
            totn = 1.0;
        }
        ameans[0] = tot[0] / totn;
        ameans[1] = (tot[1] / totn).sqrt();
        ameans[2] = (tot[2] / totn).powf(1.0 / 3.0);
        ameans
    }
}
