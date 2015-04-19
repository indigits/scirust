#![doc="Defines the ring algebraic structure.


Briefly, a ring is an abelian group with a second binary operation 
that is distributive over the abelian group operation and is associative. 
The abelian group operation is called ``addition`` and the 
second binary operation is called ``multiplication`` in analogy with the integers. 
One familiar example of a ring is the set of integers.


A ring is a set R equipped with binary operations + and * 
satisfying the following eight axioms, called the ring axioms:

R is an abelian group under addition, meaning:

* (a + b) + c = a + (b + c) for all a, b, c in R (+ is associative).
* There is an element 0 in R such that a + 0 = a and 0 + a = a for all a in R (0 is the additive identity).
* For each a in R there exists −a in R such that a + (−a) = (−a) + a = 0 (−a is the additive inverse of a).
* a + b = b + a for all a, b in R (+ is commutative).


R is a monoid under multiplication, meaning:

* (a * b) * c = a * (b * c) for all a, b, c in R (* is associative).
* There is an element 1 in R such that a * 1 = a and 1 * a = a for all a in R (1 is the multiplicative identity).[2]

Multiplication distributes over addition:

* a * (b + c) = (a * b) + (a * c) for all a, b, c in R (left distributivity).
* (b + c) * a = (b * a) + (c * a) for all a, b, c in R (right distributivity).


References:

* http://en.wikipedia.org/wiki/Algebraic_structure
* http://en.wikipedia.org/wiki/Group_(mathematics)
* http://en.wikipedia.org/wiki/Monoid
* http://en.wikipedia.org/wiki/Ring_(mathematics)

"]

// std imports


// local imports
use algebra::zero::Zero;
use algebra::one::One;
use algebra::group::{CommutativeGroupAddPartial, CommutativeGroupAdd};
use algebra::monoid::{MonoidMulPartial, MonoidMul};
use algebra::semigroup::{SemiGroupAddPartial, SemiGroupMulPartial};


///////////////////////////////////////////////////////////


/// Ring with partial equivalence
pub trait RingPartial : CommutativeGroupAddPartial 
    + MonoidMulPartial
{
    fn prop_addition_is_associative(a : Self, b : Self, c : Self) -> bool {
        SemiGroupAddPartial::prop_is_associative(a, b, c)
    }

    fn prop_additive_identity(a : Self) -> bool{
        let z : Self = Zero::zero();
        if a != a.clone() + z.clone(){
            return false;
        }
        if a != z.clone() + a.clone(){
            return false;
        }
        true
    }

    fn prop_additive_inverse(a : Self) -> bool{
        let z : Self = Zero::zero();
        let b = -a.clone();
        if a.clone() + b.clone() != z{
            return false;
        }
        if b.clone() + a.clone() != z{
            return false;
        }
        true
    }

    fn prop_addition_is_commutative(a : Self, b : Self) -> bool {
        CommutativeGroupAddPartial::prop_is_commutative(a, b)
    }


    fn prop_multiplication_is_associative(a : Self, b : Self, c : Self) -> bool {
        SemiGroupMulPartial::prop_is_associative(a, b, c)
    }

    fn prop_multiplicative_identity(a : Self) -> bool{
        let o : Self = One::one();
        if a != a.clone() * o.clone(){
            return false;
        }
        if a != o.clone() * a.clone(){
            return false;
        }
        true
    }

    fn prop_is_distributive(a : Self, b : Self, c : Self) -> bool {
        let bc = b.clone() + c.clone();
        let lhs = a.clone()  * bc.clone();
        let ab = a.clone() * b.clone();
        let ac = a.clone() * c.clone();
        let rhs = ab.clone() + ac.clone();
        if lhs != rhs{
            return false;
        }
        let lhs = bc.clone() * a.clone();
        let ba = b.clone() * a.clone();
        let ca = c.clone() * a.clone();
        let rhs = ba + ca;
        if lhs != rhs{
            return false;
        }
        true
    }

    fn check_all_properties(a: Self, b: Self, c : Self)-> bool {
        // check addition associativity
        if !RingPartial::prop_addition_is_associative(
            a.clone(), b.clone(), c.clone()) {
            return false;
        }
        

        if !RingPartial::prop_additive_identity(a.clone()){
            return false;
        }
        if !RingPartial::prop_additive_identity(b.clone()){
            return false;
        }
        if !RingPartial::prop_additive_identity(c.clone()){
            return false;
        }

        
        if !RingPartial::prop_additive_inverse(a.clone()){
            return false;
        }
        if !RingPartial::prop_additive_inverse(b.clone()){
            return false;
        }
        if !RingPartial::prop_additive_inverse(c.clone()){
            return false;
        }
        
        if !RingPartial::prop_addition_is_commutative(a.clone(), b.clone()){
            return false;
        }
        if !RingPartial::prop_addition_is_commutative(b.clone(), c.clone()){
            return false;
        }
        if !RingPartial::prop_addition_is_commutative(a.clone(), c.clone()){
            return false;
        }
        
        if !RingPartial::prop_multiplication_is_associative(
            a.clone(), b.clone(), c.clone()) {
            return false;
        }

        if !RingPartial::prop_multiplicative_identity(a.clone()){
            return false;
        }
        if !RingPartial::prop_multiplicative_identity(b.clone()){
            return false;
        }
        if !RingPartial::prop_multiplicative_identity(c.clone()){
            return false;
        }

        if !RingPartial::prop_is_distributive(
            a.clone(), b.clone(), c.clone()) {
            return false;
        }

        true
    }
}

impl<T> RingPartial for T where
    T: CommutativeGroupAddPartial + MonoidMulPartial
{}


///////////////////////////////////////////////////////////

/// Ring with full equivalence
pub trait Ring : RingPartial + CommutativeGroupAdd 
    + MonoidMul
{
    fn check_all_properties(a: Self, b: Self, c : Self)-> bool {
        // check addition associativity
        if !RingPartial::prop_addition_is_associative(
            a.clone(), b.clone(), c.clone()) {
            return false;
        }
        

        if !RingPartial::prop_additive_identity(a.clone()){
            return false;
        }
        if !RingPartial::prop_additive_identity(b.clone()){
            return false;
        }
        if !RingPartial::prop_additive_identity(c.clone()){
            return false;
        }

        
        if !RingPartial::prop_additive_inverse(a.clone()){
            return false;
        }
        if !RingPartial::prop_additive_inverse(b.clone()){
            return false;
        }
        if !RingPartial::prop_additive_inverse(c.clone()){
            return false;
        }
        
        if !RingPartial::prop_addition_is_commutative(a.clone(), b.clone()){
            return false;
        }
        if !RingPartial::prop_addition_is_commutative(b.clone(), c.clone()){
            return false;
        }
        if !RingPartial::prop_addition_is_commutative(a.clone(), c.clone()){
            return false;
        }
        
        if !RingPartial::prop_multiplication_is_associative(
            a.clone(), b.clone(), c.clone()) {
            return false;
        }

        if !RingPartial::prop_multiplicative_identity(a.clone()){
            return false;
        }
        if !RingPartial::prop_multiplicative_identity(b.clone()){
            return false;
        }
        if !RingPartial::prop_multiplicative_identity(c.clone()){
            return false;
        }

        if !RingPartial::prop_is_distributive(
            a.clone(), b.clone(), c.clone()) {
            return false;
        }
        true
    }

}

impl<T> Ring for T where
    T: RingPartial + CommutativeGroupAdd + MonoidMul
{}


/*

pub trait SemiRing{

}

pub trait NonAssociativeRing{

}

pub trait LieRing{

}

pub trait BooleanRing{

}
*/


#[cfg(test)]
mod tests {
    use super::*;



    #[test]
    fn test_ring_partial() {
        assert!(RingPartial::check_all_properties(2i8, 3i8, 1i8));
        assert!(RingPartial::check_all_properties(2i16, 3i16, 1i16));
        assert!(RingPartial::check_all_properties(2i32, 3i32, 1i32));
        assert!(RingPartial::check_all_properties(2i64, 3i64, 1i64));
        assert!(RingPartial::check_all_properties(2f32, 3f32, 1f32));
        assert!(RingPartial::check_all_properties(2f64, 3f64, 1f64));
    }


    #[test]
    fn test_ring() {
        assert!(Ring::check_all_properties(2i8, 3i8, 1i8));
        assert!(Ring::check_all_properties(2i16, 3i16, 1i16));
        assert!(Ring::check_all_properties(2i32, 3i32, 1i32));
        assert!(Ring::check_all_properties(2i64, 3i64, 1i64));
        // following cannot work
        //assert!(Ring::check_all_properties(2f64, 3f64, 1f64));
        //assert!(Ring::check_all_properties(2u64, 3u64, 1u64));
    }




}