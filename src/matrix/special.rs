#![doc="Provides functions to create special matrices.
"]

// std imports
use std::num;


// local imports
use matrix::matrix::{Matrix, MatrixF64};
use matrix::error::*;

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


#[cfg(test)]
mod test{
    use super::*;

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
}
