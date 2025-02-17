pub struct Porosity {
    pub core: f64,
    pub mantle: f64,
}

pub struct GrainSize {
    pub amin: f64,
    pub amax: f64,
    pub sd: String,
    pub na: i32,
}

pub struct WavelengthGrid {
    pub lmin: f64,
    pub lmax: f64,
    pub nlam: i32,
}

pub struct Config {
    pub porosity: Porosity,
    pub dhs_fmax: Option<f64>,
    pub mmf_size_a0: Option<f64>,
    pub mmf_df_or_fill: Option<f64>,
    pub wavelength_grid: WavelengthGrid,
}
