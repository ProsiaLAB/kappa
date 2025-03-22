use ndarray::{Array1, Array2, Array3};
use num_complex::Complex;

pub type RVector = Array1<f64>;
pub type CVector = Array1<Complex<f64>>;
pub type RMatrix = Array2<f64>;
pub type CMatrix = Array2<Complex<f64>>;
pub type RTensor = Array3<f64>;
pub type CTensor = Array3<Complex<f64>>;
pub type UVector = Array1<usize>;
