pub mod spline {
    use anyhow::bail;
    use anyhow::Result;
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
    pub fn spline(xv: &[f64], yv: &[f64], n: usize, yp1: f64, ypn: f64) -> Vec<f64> {
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

    pub fn lpmns(m: usize, n: usize, x: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        let mut pm = vec![0.0; n];
        let mut pd = vec![0.0; n];

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

    pub fn lpn(n: usize, x: f64) -> (Vec<f64>, Vec<f64>) {
        let mut pn = vec![0.0; n];
        let mut pd = vec![0.0; n];

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

    pub fn gauss_legendre(x1: f64, x2: f64, n: usize) -> (Vec<f64>, Vec<f64>) {
        let eps = 1.0e-14;

        let mut pp = 0.0;
        let m = (n + 1) / 2;
        let xm = 0.5 * (x2 + x1);
        let xl = 0.5 * (x2 - x1);

        let mut x: Vec<f64> = vec![0.0; n];
        let mut w: Vec<f64> = vec![0.0; n];

        let nf = n as f64;

        for i in 0..m {
            let ir = i as f64;
            let mut z = PI * (ir + 0.75) / (nf + 0.5);
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
    ///     
    ///                    pi^2     /Infinity
    ///       S_p(k*Rg) =  ----  *  |  du  u * J_{p+1/2}(u) * H_{p+1/2}^(1)(u) * g(u/k),
    ///                     k^3     /0
    ///     
    ///  where g(u) is the two-point correlation function (Equations 18 and 19):
    ///     
    ///                    1       / u  \ ^{d_f-3}       / u  \
    ///       g(u) =   ---------  |----- |       *  fc  | ---- |,
    ///                4*pi*Rg^3   \ Rg /                \ Rg /
    ///     
    ///  where fc is the cut-off function. By substituting g(u) into S_p, we have
    ///     
    ///                      pi   /u_max
    ///       S_p(k*Rg) = ------ *|  du  u^{df-2} * J_{p+1/2}(u) * H_{p+1/2}^(1)(u) * fc(u/xg),
    ///                    4Xg^df /u_min
    ///     
    ///  where the integration range is approximated by the range [u_min,u_max].
    ///  By using the spherical Bessel j_p(u) and Hankel functions of 1st kind h_p^(1)(u),
    ///  the Bessel function and the Hankel function are rewritten by
    ///
    ///               J_{p+1/2}    (u) = sqrt(2u/pi) * j_p(u)
    ///               H_{p+1/2}^(1)(u) = sqrt(2u/pi) * h_p^(1)(u)
    ///
    ///  we have
    ///     
    ///                      1        /u_max
    ///       S_p(k*Rg) =  ------  *  |  du  u^{df-1} * j_{p}(u) * h_{p}^(1)(u) * fc(u/xg),
    ///                    2*Xg^df    /u_min
    ///     
    ///  For the unitary condition of the two-point correlation function is
    ///     
    ///                     /Infinity
    ///       1         =   |  dw  4*pi*w^2 g(w),
    ///                     /0
    ///     
    ///  If we take the integration variable as w = u/k, then we obtain
    ///     
    ///                     1     /u_max
    ///       1         = ------  |  du u^{df-1} fc(u/xg),    .... Eq. (*)
    ///                    Xg^df  /u_min
    ///     
    ///     
    ///  The integration range [umin,umax] is determined as follows.
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

        let mut u = vec![0.0; nn];
        let mut intg: Vec<Complex64> = vec![Complex::new(0.0, 0.0); nn];
        let mut intg_unit = vec![0.0; nn];

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
            .zip(u.windows(2))
            .map(|(intg_pair, u_pair)| {
                0.5 * (intg_pair[0] + intg_pair[1]) * (u_pair[1] - u_pair[0])
            })
            .sum::<Complex64>();

        let mut unitary: f64 = intg_unit
            .windows(2)
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

    fn sph_bessel(m: usize, x: f64, isol: usize) -> Result<(Vec<f64>, Vec<f64>)> {
        let imax = 100; // truncation order of the series expansion
        let nwarmup = 100; // number of warm-up iterations

        let floor_val = 1.0e-70;
        let ceiling_val = 1.0e+70;

        let mut sj = vec![0.0; m + 1];
        sj[0] = x.sin() / x;

        let mut sy = vec![0.0; m + 1];
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
                let mut k = vec![0.0; m + nwarmup + 1];
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
    pub fn lu_decomposition() {
        todo!()
    }

    pub fn lu_backsubstitution() {
        todo!()
    }
    pub fn complex_matrix_inverse() {
        todo!()
    }
}

pub mod special {
    use statrs::function::gamma::gamma;

    use crate::fractal::{FractalConfig, FractalCutoff};

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

    pub fn confluent_hypergeometric() {
        todo!()
    }

    pub fn geometric_cross_secrtion() {
        todo!()
    }

    pub fn optic_limit() {
        todo!()
    }
}
