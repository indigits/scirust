use interpolation::Interpolation;
use linear_spline::LinearSpline;

#[doc="Interpolation Factory.
"]
struct Intepolate;

impl Intepolate {
    
    #[doc="Create a piecewise linear interpolation based on arbitrary points.

    # Remarks

    If your data is already sorted in vectors, consider to use
    LinearSpline.interpolate_sorted instead, which is more efficient.
    "]
    pub fn linear<A : Iterator<Item=f64>>(mut points: A, mut values: A) -> Box<Interpolation> {
        Box::new(LinearSpline::interpolate(points, values))
    }

}