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

use algebra::integral_domain::{IntegralDomainPartial, IntegralDomain};

/// Marker trait for fields with partial equivalence
pub trait FieldPartial : IntegralDomainPartial{

}


/// Marker trait for fields with full equivalence
pub trait Field : IntegralDomain {

}

