

// external imports
use num::complex::Complex;
use num::traits::{Float};
use std::fmt::Debug;


// local imports
use algebra::ops::Division;
use algebra::structure::semigroup::{SemiGroupAddPartial, 
    SemiGroupMulPartial};
use algebra::structure::quasigroup::{QuasiGroupAddPartial};
use algebra::structure::monoid::{CommutativeMonoidAddPartial
, CommutativeMonoidMulPartial};
use algebra::structure::group::{CommutativeGroupAddPartial};
use algebra::structure::integral_domain::{IntegralDomainPartial};
use algebra::structure::field::{FieldPartial};


use algebra::ops::Recip;





impl <T : Float> Division for Complex<T> {

}

impl<T: Float> Recip for Complex<T> {

    type Output = Self;
    #[inline]
    fn recip(self) -> Complex<T>{
        self.inv()
    }    
}

impl <T: Float + Debug> SemiGroupAddPartial for Complex<T>  {

}

impl <T: Float + Debug> SemiGroupMulPartial for Complex<T>  {

}


impl <T: Float + Debug> CommutativeMonoidAddPartial for  Complex<T>  {

}

impl <T: Float + Debug> CommutativeMonoidMulPartial for  Complex<T>  {

}

impl <T: Float + Debug> QuasiGroupAddPartial for Complex<T>  {

}

impl <T: Float + Debug> CommutativeGroupAddPartial for Complex<T>  {

}


impl <T: Float + Debug> IntegralDomainPartial for Complex<T>  {

}

impl <T: Float + Debug> FieldPartial for Complex<T>  {

}

#[cfg(test)]
mod tests {
    use num::complex::{Complex, Complex32};
    use algebra::structure::*;

    #[test]
    fn test_complex_traits() {
        let c : Complex32 = Complex { re: 0.0, im: 0.0 };
        is_magma_base(&c);
        
        is_magma_add_partial(&c);
        //is_magma_add(&c);
        is_magma_mul_partial(&c);
        //is_magma_mul(&c);
        
        is_quasigroup_add_partial(&c);
        //is_quasigroup_add(&c);
        //is_quasigroup_mul_partial(&c);
        //is_quasigroup_mul(&c);
        
        is_semigroup_add_partial(&c);
        //is_semigroup_add(&c);
        is_semigroup_mul_partial(&c);
        //is_semigroup_mul(&c);
        
        is_loop_add_partial(&c);
        //is_loop_add(&c);
        //is_loop_mul_partial(&c);
        //is_loop_mul(&c);


        is_monoid_add_partial(&c);
        //is_monoid_add(&c);
        is_monoid_mul_partial(&c);
        //is_monoid_mul(&c);

        is_group_add_partial(&c);
        //is_group_add(&c);
        //is_group_mul_partial(&c);
        //is_group_mul(&c);

        is_commutative_group_add_partial(&c);
        //is_commutative_group_add(&c);
        //is_commutative_group_mul_partial(&c);
        //is_commutative_group_mul(&c);


        is_ring_partial(&c);
        //is_ring(&c);

        is_commutative_ring_partial(&c);
        //is_commutative_ring(&c);

        is_integral_domain_partial(&c);
        //is_integral_domain_partial(&c);


        is_field_partial(&c);
        //is_field(&c);
    }

}

