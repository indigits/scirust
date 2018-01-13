use std::cmp;
use std::vec::Vec;

use interpolation::Interpolation;

#[doc="Piece-wise Linear Interpolation.

# Remarks

Supports both differentiation and integration.
"]
pub struct LinearSpline {
    x: Vec<f64>,
    c0: Vec<f64>,
    c1: Vec<f64>,
    indefinite_integral: Vec<f64>,
}

impl LinearSpline {

    #[doc="x - Sample points (N+1), sorted ascending.
    c0 - Sample values (N or N+1) at the corresponding points; intercept, zero order coefficients.
    c1 - Slopes (N) at the sample points (first order coefficients): N.
    "]
    pub fn new(x: &Vec<f64>, c0: &Vec<f64>, c1: &Vec<f64>) -> LinearSpline {
        if (x.len() != c0.len() + 1 && x.len() != c0.len()) || x.len() != c1.len() + 1  {
            panic!("All vectors must have the same dimensionality.");
        }

        if x.len() < 2 {
            panic!("The given array is too small. It must be at least {} long.", 2);
        }

        let result = LinearSpline {
            x: x.clone(), c0: c0.clone(), c1: c1.clone(), indefinite_integral: vec![]
        };

        LinearSpline { indefinite_integral: result.compute_indefinite_integral(), .. result }
    }

    #[doc="Create a linear spline interpolation from a set of (x,y) value pairs,
    sorted ascendingly by x.
    "]
    pub fn interpolate_sorted(x: &Vec<f64>, y: &Vec<f64>) -> LinearSpline {
        if x.len() != y.len() {
            panic!("All vectors must have the same dimensionality.");
        }

        if x.len() < 2 {
            panic!("The given array is too small. It must be at least {} long.", 2);
        }

        let mut c1 = vec![0f64; x.len() - 1];
        for i in 0..c1.len() {
            c1[i] = (y[i + 1] - y[i])/(x[i + 1] - x[i]);
        }

        LinearSpline::new(x, y, &c1)
    }

    #[doc="Create a linear spline interpolation from an unsorted set of (x,y) value pairs.
    WARNING: Works in-place and can thus causes the data array to be reordered.
    "]
    pub fn interpolate_inplace(x: &Vec<f64>, y: &Vec<f64>) -> LinearSpline {
        if x.len() != y.len() {
            panic!("All vectors must have the same dimensionality.");
        }

        // TODO sort x and y zipped with x
        LinearSpline::interpolate_sorted(x, y)
    }

    #[doc="Create a linear spline interpolation from an unsorted set of (x,y) value pairs.
    "]
    pub fn interpolate<A : Iterator<Item=f64>>(mut x: A, mut y: A) -> LinearSpline {
        LinearSpline::interpolate_inplace(&x.collect(), &y.collect())
    }

    fn compute_indefinite_integral(&self) -> Vec<f64> {
        let mut integral = vec![0f64; self.c1.len()];
        for i in 0..integral.len() {
            let w = self.x[i + 1] - self.x[i];
            integral[i + 1] = integral[i] + w*(self.c0[i] + w*self.c1[i]/2f64);
        }

        integral
    }

    #[doc="Find the index of the greatest sample point smaller than t,
    or the left index of the closest segment for extrapolation.
    "]
    fn left_segment_index(&self, t: f64) -> usize {
        let location = self.x.binary_search_by(|v| {
            v.partial_cmp(&t).expect("Couldn't compare values.")
        });
        let index = match location {
            Ok(i) if i < 0 => !i - 1,
            Ok(i) => i,
            Err(i) => i,
        };

        cmp::min(cmp::max(index, 0), self.x.len() - 2)
    }
    
}

impl Interpolation for LinearSpline {

    fn supports_differentiation(&self) -> bool { true }

    fn supports_integration(&self) -> bool { true }

    fn interpolate(&self, t: f64) -> f64 {
        let k = self.left_segment_index(t);
        self.c0[k] + (t - self.x[k])*self.c1[k]
    }

    fn differentiate(&self, t: f64) -> f64 {
        let k = self.left_segment_index(t);
        self.c1[k]
    }

    fn differentiate2(&self, t: f64) -> f64 {
        0f64
    }

    fn integrate(&self, t: f64) -> f64 {
        let k = self.left_segment_index(t);
        let x = t - self.x[k];
        return self.indefinite_integral[k] + x*(self.c0[k] + x*self.c1[k]/2f64);
    }

    fn integrate_interval(&self, a: f64, b: f64) -> f64 {
        self.integrate(b) - self.integrate(a)
    }

}