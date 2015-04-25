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

pub use self::magma::{MagmaBase, MagmaAddPartial, MagmaAdd, 
    MagmaMulPartial, MagmaMul};

pub use self::quasigroup::{QuasiGroupAddPartial, QuasiGroupAdd, 
    QuasiGroupMulPartial, QuasiGroupMul};

pub use self::semigroup::{SemiGroupAddPartial, SemiGroupAdd, 
    SemiGroupMulPartial, SemiGroupMul};


pub use self::loop_::{LoopAddPartial, LoopAdd, 
    LoopMulPartial, LoopMul};


pub use self::monoid::{MonoidAddPartial, MonoidAdd, 
    MonoidMulPartial, MonoidMul};
pub use self::monoid::{CommutativeMonoidAddPartial, CommutativeMonoidAdd, 
    CommutativeMonoidMulPartial, CommutativeMonoidMul};

pub use self::group::{GroupAddPartial, GroupAdd, 
    GroupMulPartial, GroupMul};
pub use self::group::{CommutativeGroupAddPartial, CommutativeGroupAdd, 
    CommutativeGroupMulPartial, CommutativeGroupMul};

pub use self::ring::{RingPartial, Ring};
pub use self::commutative_ring::{CommutativeRingPartial, CommutativeRing};
    

pub use self::integral_domain::{IntegralDomainPartial, IntegralDomain};

pub use self::field::{FieldPartial, Field};


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


#[allow(unused_variables)]
fn noop<T>(v : T){}

pub fn is_magma_base<T : MagmaBase>(a : &T){noop(a)}

pub fn is_magma_add_partial<T>(a : &T) where T : MagmaAddPartial {noop(a)}
pub fn is_magma_add<T>(a : &T) where T : MagmaAdd {noop(a)}
pub fn is_magma_mul_partial<T>(a : &T) where T : MagmaMulPartial {noop(a)}
pub fn is_magma_mul<T>(a : &T) where T : MagmaMul {noop(a)}

pub fn is_quasigroup_add_partial<T>(a : &T) where T : QuasiGroupAddPartial {noop(a)}
pub fn is_quasigroup_add<T>(a : &T) where T : QuasiGroupAdd {noop(a)}
pub fn is_quasigroup_mul_partial<T>(a : &T) where T : QuasiGroupMulPartial {noop(a)}
pub fn is_quasigroup_mul<T>(a : &T) where T : QuasiGroupMul {noop(a)}


pub fn is_semigroup_add_partial<T>(a : &T) where T : SemiGroupAddPartial {noop(a)}
pub fn is_semigroup_add<T>(a : &T) where T : SemiGroupAdd {noop(a)}
pub fn is_semigroup_mul_partial<T>(a : &T) where T : SemiGroupMulPartial {noop(a)}
pub fn is_semigroup_mul<T>(a : &T) where T : SemiGroupMul {noop(a)}

pub fn is_loop_add_partial<T>(a : &T) where T : LoopAddPartial {noop(a)}
pub fn is_loop_add<T>(a : &T) where T : LoopAdd {noop(a)}
pub fn is_loop_mul_partial<T>(a : &T) where T : LoopMulPartial {noop(a)}
pub fn is_loop_mul<T>(a : &T) where T : LoopMul {noop(a)}

pub fn is_monoid_add_partial<T>(a : &T) where T : MonoidAddPartial {noop(a)}
pub fn is_monoid_add<T>(a : &T) where T : MonoidAdd {noop(a)}
pub fn is_monoid_mul_partial<T>(a : &T) where T : MonoidMulPartial {noop(a)}
pub fn is_monoid_mul<T>(a : &T) where T : MonoidMul {noop(a)}

pub fn is_group_add_partial<T>(a : &T) where T : GroupAddPartial {noop(a)}
pub fn is_group_add<T>(a : &T) where T : GroupAdd {noop(a)}
pub fn is_group_mul_partial<T>(a : &T) where T : GroupMulPartial {noop(a)}
pub fn is_group_mul<T>(a : &T) where T : GroupMul {noop(a)}

pub fn is_commutative_group_add_partial<T>(a : &T) where T : CommutativeGroupAddPartial {noop(a)}
pub fn is_commutative_group_add<T>(a : &T) where T : CommutativeGroupAdd {noop(a)}
pub fn is_commutative_group_mul_partial<T>(a : &T) where T : CommutativeGroupMulPartial {noop(a)}
pub fn is_commutative_group_mul<T>(a : &T) where T : CommutativeGroupMul {noop(a)}

pub fn is_ring_partial<T>(a : &T) where T : RingPartial {noop(a)}
pub fn is_ring<T>(a : &T) where T : Ring {noop(a)}

pub fn is_commutative_ring_partial<T>(a : &T) where T : CommutativeRingPartial {noop(a)}
pub fn is_commutative_ring<T>(a : &T) where T : CommutativeRing {noop(a)}

pub fn is_integral_domain_partial<T>(a : &T) where T : IntegralDomainPartial {noop(a)}
pub fn is_integral_domain<T>(a : &T) where T : IntegralDomain {noop(a)}


pub fn is_field_partial<T>(a : &T) where T : FieldPartial {noop(a)}
pub fn is_field<T>(a : &T) where T : Field {noop(a)}
