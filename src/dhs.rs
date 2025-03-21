//! Distribution of Hollow Spheres (DHS) for irregular grains and
//! low-porosity aggregates.
//! In order to simulate irregularities in grains (irregular shapes, or the properties of
//! low-porosity aggregates), `kappa` averages the opacities of grains with an inner
//! empty region, over a range of volume fractions of this inner region between 0
//! and `fmax`. The subroutine used to compute the opacities and scattering matrix
//! elements for these structures is [`toon_ackerman_1981`]
//! ([Toon & Ackerman 1981](https://doi.org/10.1364/AO.20.003657)). For speed, you
//! can set xlim 1e3 or so to set a limit for the size parameter (`x = 2πa/λ`) where
//! `kappa` switches from DHS to [`crate::mie::de_rooij_1984`].

use std::f64::consts::PI;

use anyhow::anyhow;
use anyhow::Result;
use ndarray::Array2;
use num_complex::ComplexFloat;
use num_complex::{Complex, Complex64};

pub struct DHSConfig {
    r_core: f64,
    r_shell: f64,
    wave_number: f64,
    r_indsh: Complex64,
    r_indco: Complex64,
    mu: Vec<f64>,
    numang: usize,
    max_angle: usize,
}

pub struct DHSResult {
    q_ext: f64,
    q_sca: f64,
    q_bs: f64,
    g_qsc: f64,
    m0: Array2<f64>,
    m1: Array2<f64>,
    s10: Array2<f64>,
    d10: Array2<f64>,
}

// pub enum DHSError {
//     InvalidWaveNumber,
//     InvalidShellRadius,
//     InvalidCoreRadius,
//     InvalidRefractiveIndex,
//     NotEnoughAngles,
//     TooManyAngles,
//     InvalidAngle,
//     InsufficentDimensions,
// }

/// This function computes electromagnetic scattering by a
/// stratified sphere, i.e., a particle consisting of a
/// spherical core surrounded by a spherical shell.  The surrounding
/// medium is assumed to have refractive index unity.  The formulas,
/// manipulated to avoid the ill-conditioning that plagued earlier
/// formulations, were published in:
/// [Toon, O. and T. Ackerman, Applied Optics 20, 3657 (1981)](
/// https://doi.org/10.1364/AO.20.003657).
///
/// # Note
/// The refractive index entering into this routine use
/// a convention where the imaginary part has a different singn than
/// what is used in modern books.
pub fn toon_ackerman_1981(dhsc: &DHSConfig) -> Result<DHSResult> {
    let ll = 300000;
    let mxang = 1440;

    let tol = 1e-6;
    let c_0 = Complex::new(0.0, 0.0);
    let c_i = Complex::new(0.0, 1.0);

    let mut w: Array2<Complex64> = Array2::zeros((3, ll));
    let mut acap: Vec<Complex64> = vec![Complex::new(0.0, 0.0); ll];

    let x_shell = dhsc.r_shell * dhsc.wave_number;
    let x_core = dhsc.r_core * dhsc.wave_number;

    let mut t: Vec<f64> = vec![0.0; 5];
    t[0] = x_shell * dhsc.r_indsh.abs();
    let mut nmx_1 = (1.1 * t[0]) as usize;
    let mut nmx_2 = t[0] as usize;

    if nmx_1 <= 150 {
        nmx_1 = 150;
        nmx_2 = 135;
    }

    if dhsc.wave_number <= 0.0 {
        return Err(anyhow!("InvalidWaveNumber"));
    }

    if dhsc.r_shell <= 0.0 {
        return Err(anyhow!("InvalidShellRadius"));
    }

    if dhsc.r_core <= 0.0 || dhsc.r_core > dhsc.r_shell {
        return Err(anyhow!("InvalidCoreRadius"));
    }

    if dhsc.r_indsh.re <= 0.0 || dhsc.r_indsh.im > 0.0 {
        return Err(anyhow!("InvalidRefractiveIndex"));
    }
    if dhsc.r_indco.re <= 0.0 || dhsc.r_indco.im > 0.0 {
        return Err(anyhow!("InvalidRefractiveIndex"));
    }

    if dhsc.numang > mxang || dhsc.numang > dhsc.max_angle {
        return Err(anyhow!("`numang` is too large."));
    }
    if nmx_1 + 1 > ll {
        return Err(anyhow!("`nmx_1` is too large."));
    }

    for j in 0..dhsc.numang {
        if dhsc.mu[j] < -tol || dhsc.mu[j] > 1.0 + tol {
            return Err(anyhow!("`mu` is out of bounds."));
        }
    }

    let k1 = dhsc.r_indco * dhsc.wave_number;
    let k2 = dhsc.r_indsh * dhsc.wave_number;
    let k3 = Complex::new(1.0, 0.0) * dhsc.wave_number;

    let mut z: Vec<Complex64> = vec![c_0; 4];
    z[0] = dhsc.r_indsh * x_shell;
    z[1] = Complex::new(1.0, 0.0) * x_shell;
    z[2] = dhsc.r_indco * x_core;
    z[3] = dhsc.r_indsh * x_core;

    let x0 = z[0].re;
    let y0 = z[0].im;

    let x3 = z[3].re;
    let y3 = z[3].im;

    let rx = 1.0 / x_shell;

    let rrfx = 1.0 / (dhsc.r_indsh * x_shell);

    for nn in (0..nmx_1).rev() {
        let nnf = nn as f64;
        acap[nn] = (nnf + 1.0) * rrfx - (1.0 / ((nnf + 1.0) * rrfx + acap[nn + 1]));
        for m in 0..3 {
            w[(m, nn)] = (nnf + 1.0) * z[m + 1] - (1.0 / ((nnf + 1.0) * z[m + 1] + w[(m, nn + 1)]));
        }
    }

    let mut si2tht: Vec<f64> = vec![0.0; dhsc.numang];
    let mut pi: Array2<f64> = Array2::zeros((dhsc.numang, 3));
    let mut tau: Array2<f64> = Array2::zeros((dhsc.numang, 3));

    for j in 0..dhsc.numang {
        si2tht[j] = 1.0 - dhsc.mu[j] * dhsc.mu[j];
        pi[[j, 0]] = 0.0;
        pi[[j, 1]] = 1.0;
        tau[[j, 0]] = 0.0;
        tau[[j, 1]] = dhsc.mu[j];
    }

    t[0] = x_shell.cos();
    t[1] = x_shell.sin();

    let mut wm1 = Complex::new(t[0], -t[1]);
    let mut wfn: Vec<Complex64> = vec![Complex::new(0.0, 0.0); 2];
    let mut ta: Vec<f64> = vec![0.0; 4];

    wfn[0] = Complex::new(t[1], t[0]);
    ta[0] = t[1];
    ta[1] = t[0];
    wfn[1] = rx * wfn[0] - wm1;
    ta[2] = wfn[1].re;
    ta[3] = wfn[1].im;

    let mut dhsr = DHSResult {
        q_ext: 0.0,
        q_sca: 0.0,
        q_bs: 0.0,
        g_qsc: 0.0,
        m0: Array2::zeros((dhsc.numang, 2)),
        m1: Array2::zeros((dhsc.numang, 2)),
        s10: Array2::zeros((dhsc.numang, 2)),
        d10: Array2::zeros((dhsc.numang, 2)),
    };

    let mut n = 0;
    let sin_x0 = x0.sin();
    let sin_x3 = x3.sin();
    let cos_x0 = x0.cos();
    let cos_x3 = x3.cos();

    let e_y0 = y0.exp();
    let e2_y0 = e_y0 * e_y0;
    let e_y3 = y3.exp();
    let ey0_my3 = (y0 - y3).exp();
    let ey0_py3 = e_y0 * e_y3;

    let aa = sin_x3 * (ey0_py3 + ey0_my3);
    let bb = cos_x3 * (ey0_py3 - ey0_my3);
    let cc = sin_x0 * (e2_y0 + 1.0);
    let dd = cos_x0 * (e2_y0 - 1.0);

    let denom = 1.0 + e2_y0 * (4.0 * sin_x0 * sin_x0 - 2.0 + e2_y0);
    let mut temp = Complex::new((aa * cc + bb * dd) / denom, (bb * cc - aa * dd) / denom);
    temp *= (acap[n] + 1.0 / z[0]) / (w[[2, n]] + 1.0 / z[3]);
    let mut tempsq = temp * temp;

    let mut p23_h23 = 0.5 + Complex::new(sin_x3 * sin_x3 - 0.5, cos_x3 * sin_x3) * e_y3;
    let mut p23_h20 =
        0.5 * Complex::new(
            sin_x0 * sin_x3 - cos_x0 * cos_x3,
            sin_x0 * cos_x3 + cos_x0 * sin_x3,
        ) * ey0_py3
            + 0.5
                * Complex::new(
                    sin_x0 * sin_x3 + cos_x0 * cos_x3,
                    -sin_x0 * cos_x3 + cos_x0 * sin_x3,
                )
                * ey0_my3;

    let mut d_h0 = z[0] / (1.0 + c_i * z[0]) - 1.0 / z[0];
    let mut d_h1 = z[1] / (1.0 + c_i * z[1]) - 1.0 / z[1];
    let mut d_h3 = z[3] / (1.0 + c_i * z[3]) - 1.0 / z[3];

    p23_h23 /= (d_h3 + 1.0 / z[3]) * (w[[2, n]] + 1.0 / z[3]);
    p23_h20 /= (d_h0 + 1.0 / z[0]) * (w[[2, n]] + 1.0 / z[3]);

    let mut u: Vec<Complex64> = vec![Complex::new(0.0, 0.0); 8];

    u[0] = k3 * acap[n] - k2 * w[[0, n]];
    u[1] = k3 * acap[n] - k2 * d_h1;
    u[2] = k2 * acap[n] - k3 * w[[0, n]];
    u[3] = k2 * acap[n] - k3 * d_h1;
    u[4] = k1 * w[[2, n]] - k2 * w[[1, n]];
    u[5] = k2 * w[[2, n]] - k1 * w[[1, n]];
    u[6] = -c_i * (temp * p23_h20 - p23_h23);
    u[7] = ta[2] / wfn[1];

    let mut acoe = u[7] * (u[0] * u[4] * u[6] + k1 * u[0] - tempsq * k3 * u[4])
        / (u[1] * u[4] * u[6] + k1 * u[1] - tempsq * k3 * u[4]);
    let mut bcoe = u[7] * (u[2] * u[5] * u[6] + k2 * u[2] - tempsq * k2 * u[5])
        / (u[3] * u[5] * u[6] + k2 * u[3] - tempsq * k2 * u[5]);

    let mut acoem_0 = acoe;
    let mut bcoem_0 = bcoe;
    let mut are = acoe.re;
    let mut aim = acoe.im;
    let mut bre = bcoe.re;
    let mut bim = bcoe.im;

    let mut dqext = 3.0 * (are + bre);
    let mut dqsca = 3.0 * (are * are + aim * aim + bre * bre + bim * bim);
    let mut dgqsc = 0.0;
    let mut sback = 3.0 * (acoe - bcoe);
    let mut rmm = 1.0;

    let mut ac = 1.5 * acoe;
    let mut bc = 1.5 * bcoe;

    let mut s0: Array2<Complex64> = Array2::zeros((3, 3));
    let mut s1: Array2<Complex64> = Array2::zeros((3, 3));

    for j in 0..dhsc.numang {
        s0[[j, 0]] = ac * pi[[j, 1]] + bc * tau[[j, 1]];
        s0[[j, 1]] = ac * pi[[j, 1]] - bc * tau[[j, 1]];
        s1[[j, 0]] = bc * pi[[j, 1]] + ac * tau[[j, 1]];
        s1[[j, 1]] = bc * pi[[j, 1]] - ac * tau[[j, 1]];
    }

    n = 1;

    t[3] = 1.5;

    while t[3] >= 1e-14 {
        let nf = n as f64;
        t[0] = 2.0 * nf - 1.0;
        t[1] = nf - 1.0;

        for j in 0..dhsc.numang {
            pi[[j, 2]] = (t[0] * pi[[j, 1]] * dhsc.mu[j] - nf * pi[[j, 0]]) / t[1];
            tau[[j, 2]] = dhsc.mu[j] * (pi[[j, 2]] - pi[[j, 0]]) - t[0] * si2tht[j] * pi[[j, 1]]
                + tau[[j, 0]];
        }

        wm1 = wfn[0];
        wfn[0] = wfn[1];
        wfn[1] = t[0] * rx * wfn[0] - wm1;
        ta[0] = wfn[0].re;
        ta[1] = wfn[0].im;
        ta[2] = wfn[1].re;
        ta[3] = wfn[1].im;

        d_h0 = -(nf) / z[0] + 1.0 / (nf / z[0] - d_h0);
        d_h1 = -(nf) / z[1] + 1.0 / (nf / z[1] - d_h1);
        d_h3 = -(nf) / z[3] + 1.0 / (nf / z[3] - d_h3);
        p23_h23 /= (d_h3 + nf / z[3]) * (w[[2, n]] + nf / z[3]);
        p23_h20 /= (d_h0 + nf / z[0]) * (w[[2, n]] + nf / z[3]);
        temp *= (acap[n] + nf / z[0]) / (w[[2, n]] + nf / z[3]);
        tempsq = temp * temp;

        u[0] = k3 * acap[n] - k2 * w[[0, n]];
        u[1] = k3 * acap[n] - k2 * d_h1;
        u[2] = k2 * acap[n] - k3 * w[[0, n]];
        u[3] = k2 * acap[n] - k3 * d_h1;
        u[4] = k1 * w[[2, n]] - k2 * w[[1, n]];
        u[5] = k2 * w[[2, n]] - k1 * w[[1, n]];
        u[6] = -c_i * (temp * p23_h20 - p23_h23);
        u[7] = ta[2] / wfn[1];

        acoe = u[7] * (u[0] * u[4] * u[6] + k1 * u[0] - tempsq * k3 * u[4])
            / (u[1] * u[4] * u[6] + k1 * u[1] - tempsq * k3 * u[4]);
        bcoe = u[7] * (u[2] * u[5] * u[6] + k2 * u[2] - tempsq * k2 * u[5])
            / (u[3] * u[5] * u[6] + k2 * u[3] - tempsq * k2 * u[5]);

        are = acoe.re;
        aim = acoe.im;
        bre = bcoe.re;
        bim = bcoe.im;

        let am_0_re = acoem_0.re;
        let am_0_im = acoem_0.im;
        let bm_0_re = bcoem_0.re;
        let bm_0_im = bcoem_0.im;

        t[3] = (2.0 * (nf) - 1.0) / ((nf) * (nf - 1.0));
        t[1] = (nf - 1.0) / (nf + 1.0) / (nf);
        dgqsc += t[1] * (am_0_re * are + am_0_im * aim + bm_0_re * bre + bm_0_im * bim)
            + t[3] * (am_0_re * bm_0_re + am_0_im * bm_0_im);

        t[2] = 2.0 * (nf) + 1.0;
        dqext += t[2] * (are + bre);
        t[3] = are * are + aim * aim + bre * bre + bim * bim;
        dqsca += t[2] * t[3];
        rmm = -rmm;
        sback += t[2] * rmm * (acoe - bcoe);

        t[1] = (nf) / (nf + 1.0);
        t[0] = t[2] / t[1];

        ac = t[0] * acoe;
        bc = t[0] * bcoe;

        for j in 0..dhsc.numang {
            s0[[j, 0]] += ac * pi[[j, 2]] + bc * tau[[j, 2]];
            s1[[j, 0]] += bc * pi[[j, 2]] + ac * tau[[j, 2]];
        }

        if n % 2 == 0 {
            for j in 0..dhsc.numang {
                s0[[j, 1]] -= ac * pi[[j, 2]] + bc * tau[[j, 2]];
                s1[[j, 1]] -= bc * pi[[j, 2]] + ac * tau[[j, 2]];
            }
        } else {
            for j in 0..dhsc.numang {
                s0[[j, 1]] += ac * pi[[j, 2]] - bc * tau[[j, 2]];
                s1[[j, 1]] += bc * pi[[j, 2]] - ac * tau[[j, 2]];
            }
        }

        n += 1;
        if n > nmx_2 {
            return Err(anyhow!("InsufficentDimensions"));
        }

        for j in 0..dhsc.numang {
            pi[[j, 0]] = pi[[j, 1]];
            pi[[j, 1]] = pi[[j, 2]];
            tau[[j, 0]] = tau[[j, 1]];
            tau[[j, 1]] = tau[[j, 2]];
        }

        acoem_0 = acoe;
        bcoem_0 = bcoe;
    }

    for j in 0..dhsc.numang {
        for k in 0..1 {
            dhsr.m0[[j, k]] = (s0[[j, k]].re).powi(2) + (s0[[j, k]].im).powi(2);
            dhsr.m1[[j, k]] = (s1[[j, k]].re).powi(2) + (s1[[j, k]].im).powi(2);
            dhsr.s10[[j, k]] = s0[[j, k]].re * s1[[j, k]].re + s0[[j, k]].im * s1[[j, k]].im;
            dhsr.d10[[j, k]] = s1[[j, k]].re * s0[[j, k]].im - s0[[j, k]].re * s1[[j, k]].im;
        }
    }

    t[0] = 2.0 * rx * rx;
    dhsr.q_ext = t[0] * dqext;
    dhsr.q_sca = t[0] * dqsca;
    dhsr.g_qsc = 2.0 * t[0] * dgqsc;
    sback *= 0.5;
    dhsr.q_bs = ((sback.re).powi(2) + (sback.im).powi(2)) / (PI * x_shell * x_shell);

    Ok(dhsr)
}
