use ndarray::Array1;
use std::error::Error;

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

    get_size_distribution();
}

fn get_size_distribution() {
    todo!("Size distribution");
}
