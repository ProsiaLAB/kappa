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
    pub fn lpmns() {
        todo!()
    }
    pub fn lpn() {
        todo!()
    }

    pub fn gauss_legendre() {
        todo!()
    }

    pub fn gauss_legendre_alt() {
        todo!()
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
