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
    pub fn int_sph_bessel() {
        todo!()
    }

    pub fn sph_bessel() {
        todo!()
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
    pub fn two_point_correlation() {
        todo!()
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
