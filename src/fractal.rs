//! Computes light scattering properties of randomly oriented
//! fractal dust aggregates by means of the modified mean field theory developed
//! in Tazaki & Tanaka (2018). This code is also capable of computing the light
//! scattering solution based on the Rayleigh-Gans-Debye theory and the Mean field
//! theory.

pub fn mean_scattering() {
    maxwell_garnett_mixing();
    structure_factor_integration();
    lorenz_mie();
    todo!()
}

fn maxwell_garnett_mixing() {
    todo!()
}

fn structure_factor_integration() {
    structure_factor_fn();
    todo!()
}

fn structure_factor_fn() {
    todo!()
}

fn lorenz_mie() {
    renormalize();
    todo!()
}

fn renormalize() {
    todo!()
}
