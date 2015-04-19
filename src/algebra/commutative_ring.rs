#![doc="Defines the commutative ring algebraic structure.

A commutative ring is a ring where the multiplication
operation is commutative.

A commutative ring is a set R equipped with binary operations + and * 
satisfying the following nine axioms, called the ring axioms:

R is an abelian group under addition, meaning:

* (a + b) + c = a + (b + c) for all a, b, c in R (+ is associative).
* There is an element 0 in R such that a + 0 = a and 0 + a = a for all a in R (0 is the additive identity).
* For each a in R there exists −a in R such that a + (−a) = (−a) + a = 0 (−a is the additive inverse of a).
* a + b = b + a for all a, b in R (+ is commutative).


R is a monoid under multiplication, meaning:

* (a * b) * c = a * (b * c) for all a, b, c in R (* is associative).
* There is an element 1 in R such that a * 1 = a and 1 * a = a for all a in R (1 is the multiplicative identity).[2]
* a * b = b * a for all a, b in R (* is commutative).

Multiplication distributes over addition:

* a * (b + c) = (a * b) + (a * c) for all a, b, c in R (left distributivity).
* (b + c) * a = (b * a) + (c * a) for all a, b, c in R (right distributivity).


References:

* http://en.wikipedia.org/wiki/Algebraic_structure
* http://en.wikipedia.org/wiki/Group_(mathematics)
* http://en.wikipedia.org/wiki/Monoid
* http://en.wikipedia.org/wiki/Ring_(mathematics)

"]

use algebra::monoid::{CommutativeMonoidMulPartial, CommutativeMonoidMul};
use algebra::ring::{RingPartial, Ring};


/// Commutative ring with partial equivalence
pub trait  CommutativeRingPartial : RingPartial 
    + CommutativeMonoidMulPartial
{

    fn prop_multiplication_is_commutative(a : Self, b : Self) -> bool {
        CommutativeMonoidMulPartial::prop_is_commutative(a, b)
    }


    fn check_all_properties(a: Self, b: Self, c : Self)-> bool {
        if !RingPartial::check_all_properties(
            a.clone(), b.clone(), c.clone()) {
            return false;
        }
        if !CommutativeRingPartial::prop_multiplication_is_commutative(a.clone(), b.clone()){
            return false;
        }
        if !CommutativeRingPartial::prop_multiplication_is_commutative(b.clone(), c.clone()){
            return false;
        }
        if !CommutativeRingPartial::prop_multiplication_is_commutative(a.clone(), c.clone()){
            return false;
        }
        true
    }
}

impl<T> CommutativeRingPartial for T where
    T: RingPartial + CommutativeMonoidMulPartial
{}

///////////////////////////////////////////////////////////

/// Commutative ring with full equivalence
pub trait CommutativeRing : CommutativeRingPartial 
    + Ring + CommutativeMonoidMul{
    fn check_all_properties(a: Self, b: Self, c : Self)-> bool {
        CommutativeRingPartial::check_all_properties(a, b, c)
    }
}


impl<T> CommutativeRing for T where
    T: Ring + CommutativeMonoidMul
{

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_commutative_ring_partial() {
        assert!(CommutativeRingPartial::check_all_properties(2i8, 3i8, 1i8));
        assert!(CommutativeRingPartial::check_all_properties(2i16, 3i16, 1i16));
        assert!(CommutativeRingPartial::check_all_properties(2i32, 3i32, 1i32));
        assert!(CommutativeRingPartial::check_all_properties(2i64, 3i64, 1i64));
        assert!(CommutativeRingPartial::check_all_properties(2f32, 3f32, 1f32));
        assert!(CommutativeRingPartial::check_all_properties(2f64, 3f64, 1f64));
    }

    #[test]
    fn test_commutative_ring() {
        assert!(CommutativeRing::check_all_properties(2i8, 3i8, 1i8));
        assert!(CommutativeRing::check_all_properties(2i16, 3i16, 1i16));
        assert!(CommutativeRing::check_all_properties(2i32, 3i32, 1i32));
        assert!(CommutativeRing::check_all_properties(2i64, 3i64, 1i64));
        // following cannot work
        //assert!(Ring::check_all_properties(2f64, 3f64, 1f64));
        //assert!(Ring::check_all_properties(2u64, 3u64, 1u64));
    }

}
