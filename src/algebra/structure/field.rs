#![doc="Defines the field algebraic structure.



A field is a set together with two operations, 
usually called addition and multiplication, 
and denoted by + and *, respectively, 
such that the following axioms hold; 
subtraction and division are defined in terms of the 
inverse operations of addition and multiplication:


Closure of F under addition and multiplication

* For all a, b in F, both a + b and a * b are in F 
  (or more formally, + and * are binary operations on F).


Associativity of addition and multiplication

* For all a, b, and c in F, the following equalities hold: 
  a + (b + c) = (a + b) + c and a * (b * c) = (a * b) * c.


Commutativity of addition and multiplication

* For all a and b in F, the following equalities hold: 
  a + b = b + a and a * b = b * a.

Existence of additive and multiplicative identity elements

* There exists an element of F, called the additive identity element 
  and denoted by 0, such that for all a in F, a + 0 = a. 
* Likewise, there is an element, called the multiplicative identity element 
  and denoted by 1, such that for all a in F, a * 1 = a. 
* To exclude the trivial ring, the additive identity and the multiplicative identity 
  are required to be distinct.


Existence of additive inverses and multiplicative inverses

* For every a in F, there exists an element -a in F, 
  such that a + (-a) = 0. 
* Similarly, for any a in F other than 0, there exists an element 
  b in F, such that a * b = 1. 

Distributivity of multiplication over addition

* For all a, b and c in F, the following equality holds: 
  a * (b + c) = (a * b) + (a  * c).


"]
// std imports
use std::ops::Div;

// external imports
use num::traits::{Zero, One};

// local imports
use algebra::ops::Recip;
use algebra::structure::integral_domain::{IntegralDomainPartial, IntegralDomain};
use algebra::structure::semigroup::{SemiGroupAdd, SemiGroupAddPartial, 
    SemiGroupMulPartial};
use algebra::structure::group::{CommutativeGroupAdd, CommutativeGroupAddPartial};
use algebra::structure::monoid::{CommutativeMonoidMulPartial, CommutativeMonoidMul};


/// Marker trait for fields with partial equivalence
pub trait FieldPartial : IntegralDomainPartial
        + Div<Output=Self>
        + Recip<Output=Self>
{

    fn prop_addition_is_associative(a : Self, b : Self, c : Self) -> bool {
        SemiGroupAddPartial::prop_is_associative(a, b, c)
    }

    fn prop_multiplication_is_associative(a : Self, b : Self, c : Self) -> bool {
        SemiGroupMulPartial::prop_is_associative(a, b, c)
    }

    fn prop_addition_is_commutative(a : Self, b : Self) -> bool {
        CommutativeGroupAddPartial::prop_is_commutative(a, b)
    }

    fn prop_multiplication_is_commutative(a : Self, b : Self) -> bool {
        CommutativeMonoidMulPartial::prop_is_commutative(a, b)
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

    #[allow(unused_variables)]
    fn prop_distinct_zero_one(a : Self) -> bool{
        let o : Self = One::one();
        let z : Self = Zero::zero();
        o != z
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

    fn prop_multiplicative_inverse(a : Self) -> bool{
        let z : Self = Zero::zero();
        if z == a {
            return true;
        }
        let b = a.clone().recip();
        let o : Self = One::one();

        if a.clone() *  b.clone() != o{
            return false;
        }
        if b.clone() * a.clone() != o{
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
        if !FieldPartial::prop_addition_is_associative(
            a.clone(), b.clone(), c.clone()) {
            return false;
        }
        

        if !FieldPartial::prop_additive_identity(a.clone()){
            return false;
        }
        if !FieldPartial::prop_additive_identity(b.clone()){
            return false;
        }
        if !FieldPartial::prop_additive_identity(c.clone()){
            return false;
        }

        
        if !FieldPartial::prop_additive_inverse(a.clone()){
            return false;
        }
        if !FieldPartial::prop_additive_inverse(b.clone()){
            return false;
        }
        if !FieldPartial::prop_additive_inverse(c.clone()){
            return false;
        }
        
        if !FieldPartial::prop_multiplicative_inverse(a.clone()){
            return false;
        }
        if !FieldPartial::prop_multiplicative_inverse(b.clone()){
            return false;
        }
        if !FieldPartial::prop_multiplicative_inverse(c.clone()){
            return false;
        }
        
        if !FieldPartial::prop_addition_is_commutative(a.clone(), b.clone()){
            return false;
        }
        if !FieldPartial::prop_addition_is_commutative(b.clone(), c.clone()){
            return false;
        }
        if !FieldPartial::prop_addition_is_commutative(a.clone(), c.clone()){
            return false;
        }
        
        if !FieldPartial::prop_multiplication_is_commutative(a.clone(), b.clone()){
            return false;
        }
        if !FieldPartial::prop_multiplication_is_commutative(b.clone(), c.clone()){
            return false;
        }
        if !FieldPartial::prop_multiplication_is_commutative(a.clone(), c.clone()){
            return false;
        }
        
        if !FieldPartial::prop_multiplication_is_associative(
            a.clone(), b.clone(), c.clone()) {
            return false;
        }

        if !FieldPartial::prop_multiplicative_identity(a.clone()){
            return false;
        }
        if !FieldPartial::prop_multiplicative_identity(b.clone()){
            return false;
        }
        if !FieldPartial::prop_multiplicative_identity(c.clone()){
            return false;
        }

        if !FieldPartial::prop_is_distributive(
            a.clone(), b.clone(), c.clone()) {
            return false;
        }

        if !FieldPartial::prop_distinct_zero_one(a.clone()){
            return false;
        }

        true
    }
}



impl FieldPartial for f32  {}
impl FieldPartial for f64  {}


/// Marker trait for fields with full equivalence
pub trait Field : FieldPartial + IntegralDomain {
    fn prop_addition_is_associative(a : Self, b : Self, c : Self) -> bool {
        SemiGroupAdd::prop_is_associative(a, b, c)
    }

    fn prop_multiplication_is_associative(a : Self, b : Self, c : Self) -> bool {
        // TODO : using SemiGroupMul doesn't work here. Why? 
        SemiGroupMulPartial::prop_is_associative(a, b, c)
    }

    fn prop_addition_is_commutative(a : Self, b : Self) -> bool {
        CommutativeGroupAdd::prop_is_commutative(a, b)
    }

    fn prop_multiplication_is_commutative(a : Self, b : Self) -> bool {
        CommutativeMonoidMul::prop_is_commutative(a, b)
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

     #[allow(unused_variables)]
    fn prop_distinct_zero_one(a : Self) -> bool{
        let o : Self = One::one();
        let z : Self = Zero::zero();
        o != z
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

    fn prop_multiplicative_inverse(a : Self) -> bool{
        let z : Self = Zero::zero();
        if z == a {
            return true;
        }
        let b = a.clone().recip();
        let o : Self = One::one();

        if a.clone() *  b.clone() != o{
            return false;
        }
        if b.clone() * a.clone() != o{
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
        if !Field::prop_addition_is_associative(
            a.clone(), b.clone(), c.clone()) {
            return false;
        }
        

        if !Field::prop_additive_identity(a.clone()){
            return false;
        }
        if !Field::prop_additive_identity(b.clone()){
            return false;
        }
        if !Field::prop_additive_identity(c.clone()){
            return false;
        }

        
        if !Field::prop_additive_inverse(a.clone()){
            return false;
        }
        if !Field::prop_additive_inverse(b.clone()){
            return false;
        }
        if !Field::prop_additive_inverse(c.clone()){
            return false;
        }
        
        if !Field::prop_multiplicative_inverse(a.clone()){
            return false;
        }
        if !Field::prop_multiplicative_inverse(b.clone()){
            return false;
        }
        if !Field::prop_multiplicative_inverse(c.clone()){
            return false;
        }
        
        if !Field::prop_addition_is_commutative(a.clone(), b.clone()){
            return false;
        }
        if !Field::prop_addition_is_commutative(b.clone(), c.clone()){
            return false;
        }
        if !Field::prop_addition_is_commutative(a.clone(), c.clone()){
            return false;
        }
        
        if !Field::prop_multiplication_is_commutative(a.clone(), b.clone()){
            return false;
        }
        if !Field::prop_multiplication_is_commutative(b.clone(), c.clone()){
            return false;
        }
        if !Field::prop_multiplication_is_commutative(a.clone(), c.clone()){
            return false;
        }
        
        if !Field::prop_multiplication_is_associative(
            a.clone(), b.clone(), c.clone()) {
            return false;
        }

        if !Field::prop_multiplicative_identity(a.clone()){
            return false;
        }
        if !Field::prop_multiplicative_identity(b.clone()){
            return false;
        }
        if !Field::prop_multiplicative_identity(c.clone()){
            return false;
        }

        if !Field::prop_is_distributive(
            a.clone(), b.clone(), c.clone()) {
            return false;
        }

        if !Field::prop_distinct_zero_one(a.clone()){
            return false;
        }

        true
    }

}

