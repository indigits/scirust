#![doc="Provides functions to create special matrices.
"]

// std imports
use std::num;


// local imports
use matrix::matrix::{Matrix, MatrixI64, MatrixF64};
use matrix::error::*;
use matrix::element::MatrixElt;

#[doc="Returns a hadamard matrix of size n x n


n must be a power of 2.
"] 
pub fn hadamard(n : uint) -> Result<MatrixF64, MatrixError>{
    if !num::is_power_of_two(n){
        return Err(IsNotPowerOfTwo);
    }
    let mut m : MatrixF64 = Matrix::new(n, n);
    // Take the log of n with respect to 2.
    let order = n.trailing_zeros();
    // We are going to use Sylvester's construction
    // http://en.wikipedia.org/wiki/Hadamard_matrix


    // Let's fill the first level Hadamard matrix
    m.set(0, 0, 1.0);
    for  o in range(0, order){
        // We will construct four views.
        let size = num::pow(2u, o);
        // top left block
        let tl = m.view(0, 0, size, size);
        // top right block
        let mut tr = m.view(0, size, size, size);
        // bottom left block
        let mut bl = m.view(size, 0, size, size);
        // bottom right block
        let mut br = m.view(size, size, size, size);
        tr.copy_from(&tl);
        bl.copy_from(&tl);
        br.copy_scaled_from(&tl, -1.0);
    }
    Ok(m)
}


#[doc="Returns a matrix whose entries are picked up from
a range in column wise order.

# Remarks

Say you wish to create a matrix of 100 elements. If you
provide a range of 80 numbers, the first 80 entries in the
matrix will be filled by the numbers from the range. The
remaining 20 entries will be filled with zeros. On the
other hand, if you provide a range with more than 100 numbers,
then only the first 100 numbers will be used to fill the matrix
(off course in column major order). The remaining numbers
in the range will not be used. They will also not be generated.


# Examples

Constructing a 4x4 matrix of floating point numbers:

        use scirust::matrix::MatrixF64;
        use scirust::matrix::special::from_range;
        let start  = 0.0;
        let stop = 16.0;
        let m : MatrixF64 = from_range(4, 4, start, stop);
        for i in range(0, 16){
            let c = i >> 2;
            let r = i & 3;
            assert_eq!(m.get(r, c), i as f64);
        }


"]
pub fn from_range<T:MatrixElt+PartialOrd+Clone+ToPrimitive>(rows : uint, cols : uint, 
    start : T, stop : T )-> Matrix<T> {
    let m : Matrix<T> = Matrix::from_iter(rows, cols, range(start, stop));
    m 
}

#[doc="Returns a 64-bit integer matrix whose entries are picked up from
a range in column wise order.

See from_range function  for further discussion.


# Examples

    use scirust::matrix::special::from_range_i64;

    let m = from_range_i64(4, 4, 0, 16);
    for i in range(0, 16){
        let c = i >> 2;
        let r = i & 3;
        assert_eq!(m.get(r, c), i as i64);
    }


"]pub fn from_range_i64(rows : uint, cols : uint, 
    start : i64, stop : i64)->MatrixI64 {
    let m : MatrixI64 = Matrix::from_iter(rows, cols, range(start, stop));
    m 
}


#[cfg(test)]
mod test{
    use super::*;
    use matrix::matrix::*;

    #[test]
    fn test_hadamard(){
        let m = hadamard(1).unwrap();
        assert_eq!(m.num_cells(), 1);
        assert_eq!(m.get(0,0), 1.0);
        let m = hadamard(2).unwrap();
        assert_eq!(m.num_cells(), 4);
        assert_eq!(m.get(0,0), 1.0);
        assert_eq!(m.get(0,1), 1.0);
        assert_eq!(m.get(1,0), 1.0);
        assert_eq!(m.get(1,1), -1.0);
        let m = hadamard(4).unwrap();
        assert!(m.is_square());
    }

    #[test]
    fn test_range(){
        let m = from_range_i64(4, 4, 0, 16);
        for i in range(0, 16){
            let c = i >> 2;
            let r = i & 3;
            assert_eq!(m.get(r, c), i as i64);
        }
        let start  = 0.0;
        let stop = 16.0;
        let m : MatrixF64 = from_range(4, 4, start, stop);
        for i in range(0, 16){
            let c = i >> 2;
            let r = i & 3;
            assert_eq!(m.get(r, c), i as f64);
        }
    }
}
