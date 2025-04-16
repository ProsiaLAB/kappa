use ndarray::ArrayView1;
use ndarray::{Array1, Array2, Array3};
use num_complex::Complex;

/// A fixed-length array of real ([f64]) numbers.
pub type RVector = Array1<f64>;

/// A view of [RVector].
pub type RVecView<'a> = ArrayView1<'a, f64>;

/// A fixed-length array of complex ([f64]) numbers.
pub type CVector = Array1<Complex<f64>>;

/// A 2-dimensional array (matrix) of real ([f64]) numbers.
pub type RMatrix = Array2<f64>;

/// A 2-dimensional array (matrix) of complex ([f64]) numbers.
pub type CMatrix = Array2<Complex<f64>>;

/// A 3-dimensional array (tensor) of real ([f64]) numbers.
pub type RTensor = Array3<f64>;

/// A 3-dimensional array (tensor) of complex ([f64]) numbers.
pub type CTensor = Array3<Complex<f64>>;

/// A fixed-length array of unsigned integers.
pub type UVector = Array1<usize>;

/// A fixed-length array of booleans.
pub type BVector = Array1<bool>;
