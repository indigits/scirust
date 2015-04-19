#![doc="Defines the integral domain algebraic structure.

There are a number of equivalent definitions of integral domain:

* An integral domain is a nonzero commutative ring in which the product of any two nonzero elements is nonzero.
* An integral domain is a nonzero commutative ring with no nonzero zero divisors.
* An integral domain is a commutative ring in which the zero ideal {0} is a prime ideal.
* An integral domain is a commutative ring for which every non-zero element is cancellable under multiplication.
* An integral domain is a ring for which the set of nonzero elements is a commutative monoid under multiplication (because the monoid is closed under multiplication).
* An integral domain is a ring that is (isomorphic to) a subring of a field. (This implies it is a nonzero commutative ring.)
* An integral domain is a nonzero commutative ring in which for every nonzero element r, the function that maps each element x of the ring to the product xr is injective. Elements r with this property are called regular, so it is equivalent to require that every nonzero element of the ring be regular.


References:

* http://en.wikipedia.org/wiki/Algebraic_structure
* http://en.wikipedia.org/wiki/Ring_(mathematics)
* http://en.wikipedia.org/wiki/Integral_domain

"]

// local imports
use algebra::commutative_ring::{CommutativeRingPartial, CommutativeRing};

/// Marker trait for integral domains with partial equivalence
pub trait IntegralDomainPartial : CommutativeRingPartial {

}

/// Marker trait for integral domains with full equivalence
pub trait IntegralDomain : IntegralDomainPartial + CommutativeRing{

}
