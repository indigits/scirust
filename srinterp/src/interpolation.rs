#[doc="Defines the basic interface for interpolation methods
"]
pub trait Interpolation {

    #[doc="Gets a value indicating whether the algorithm supports
    differentiation (interpolated derivative).
    "]
    fn supports_differentiation(&self) -> bool;

    #[doc="Gets a value indicating whether the algorithm supports
    integration (interpolated quadrature).
    "]
    fn supports_integration(&self) -> bool;

    #[doc="Interpolate at point t.
    "]
    fn interpolate(&self, t: f64) -> f64;

    #[doc="Differentiate at point t.
    "]
    fn differentiate(&self, t: f64) -> f64;

    #[doc="Differentiate twice at point t.
    "]
    fn differentiate2(&self, t: f64) -> f64;

    #[doc="Indefinite integral at point t.
    "]
    fn integrate(&self, t: f64) -> f64;

    #[doc="Definite integral between points a and b.
    "]
    fn integrate_interval(&self, a: f64, b: f64) -> f64;
    
}
