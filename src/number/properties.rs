#![doc="Properties of algebraic structures
"]

#[allow(unused_variables)]
fn noop<T>(v : T){

}

/// A unary operation defined on a type T
pub trait UnaryOp<T>{
    fn un_op(value : &T) -> T;
}

/// A binary operation defined on a type T
pub trait BinaryOp<T> {
    fn bin_op(lhs : &T, rhs: &T) -> T;

}

pub trait Commutativity<T, Operation: BinaryOp<T>>{
    fn is_commutative(lhs : &T, rhs: &T)-> bool{
        //let v1 : T = Operation::bin_op(lhs, rhs);
        //let v2 : T = Operation::bin_op(rhs, lhs);
        //v1 == v2
        noop(lhs);
        noop(rhs);
        true
    }
}


/// An algebraic structure which is closed
pub trait Closure{

}

/// An algebraic structure which is associative
pub trait Associativity{

}

/// An algebraic structure which has an additive identity element
pub trait AdditiveIdentity{

}


/// An algebraic structure which has a multiplicative identity element
pub trait MultiplicativeIdentity{

}




pub trait Divisibility{

}


pub trait SemiCategory{

}

pub trait Category{

}

pub trait Groupoid{

}

pub trait Magma{

}

pub trait QuasiGroup{

}

pub trait Loop{

}

pub trait SemiGroup{

}

pub trait Monoid{

}

pub trait Group{

}

pub trait CommutativeGroup{

}

pub trait SemiRing{

}

pub trait NonAssociativeRing{

}


pub trait Ring{

}

pub trait CommutativeRing{

}

pub trait LieRing{

}

pub trait BooleanRing{

}


pub trait IntegralDomain{

}


pub trait Field{

}

pub trait SkewLattice{

}

pub trait Lattice{

}

pub trait BoundedLattice{

}


