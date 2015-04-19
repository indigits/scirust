#![doc="Defines generic traits for algebraic structures on which
SciRust library works.


Supporting modules

* ``ops``:  Some operations not defined in std


Algebraic structures

* Magma : closure (set with an operation)
* Quasi-group: closure, division (a set where division is always possible)
* Loop: closure, identity, inverse (Quasigroup with identity)
* Semigroup: closure, associativity
* Monoid: closure, associativity, identity
* Commutative monoid: closure, associativity, identity, commutativity (Monoid with commutativity)
* Group: closure, associativity, identity, inverse  
* Commutative group: closure, associativity, identity, inverse, commutativity
* Ring : commutative group under addition, monoid under multiplication, distributive
* Commutative Ring: commutative group under addition, commutative monoid under multiplication, distributive (ring + commutative monoid under multiplication)

Modules:

* ``magma``: Magma
* ``quasigroup`` : Quasi-group
* ``semigroup``: Semi-group
* ``loop_``: Loop
* ``monoid``: Monoid and commutative monoid
* ``group``: Group  and commutative group
* ``ring``: Ring

References:

* http://en.wikipedia.org/wiki/Algebraic_structure
* http://en.wikipedia.org/wiki/Magma_(algebra)
* http://en.wikipedia.org/wiki/Quasigroup
* http://en.wikipedia.org/wiki/Semigroup
* http://en.wikipedia.org/wiki/Monoid
* http://en.wikipedia.org/wiki/Group_(mathematics)
* http://en.wikipedia.org/wiki/Abelian_group
* http://en.wikipedia.org/wiki/Ring_(mathematics)
* http://en.wikipedia.org/wiki/Integral_domain
* http://en.wikipedia.org/wiki/Integrally_closed_domain
* http://en.wikipedia.org/wiki/Unique_factorization_domain
* http://en.wikipedia.org/wiki/Principal_ideal_domain
* http://en.wikipedia.org/wiki/Euclidean_domain
* http://en.wikipedia.org/wiki/Field_(mathematics)
* http://en.wikipedia.org/wiki/Lattice_(order)


Similar libraries

* http://hackage.haskell.org/package/numeric-prelude



Items on the agenda


"]
// Re-exporting symbols
pub use self::number::Number;
pub use self::number::num_range;
pub use self::signed::Signed;
pub use self::zero::Zero;
pub use self::one::One;
pub use self::entry::Entry;


//pub use self::complex::{Complex, Complex32, Complex64};

pub mod number;
pub mod entry;
pub mod zero;
pub mod one;
pub mod signed;
pub mod ops;

pub mod structure {
pub mod magma;
pub mod quasigroup;
pub mod semigroup;
pub mod loop_;
pub mod monoid;
pub mod group;
pub mod ring;
pub mod commutative_ring;
pub mod integral_domain;
pub mod field;

}
//pub mod complex;
