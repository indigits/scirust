#![doc="Provides functions to create different types of matrices.
"]

// std imports

// external imports
use num::{Num, One, Zero};


// local imports
use algebra::structure::{MagmaBase, CommutativeMonoidAddPartial, FieldPartial};
use matrix::matrix::{Matrix, 
    MatrixI8, MatrixI16, MatrixI32, MatrixI64,
    MatrixU8, MatrixU16, MatrixU32, MatrixU64,
    MatrixF32, MatrixF64,
    MatrixC32, MatrixC64};
use matrix::traits::{Shape};
use error::SRError;

// complex numbers
use num::complex::{Complex32, Complex64};


#[doc="Returns a Hadamard matrix of size n x n


n must be a power of 2.
"] 
pub fn hadamard(n : usize) -> Result<MatrixF64, SRError>{
    if !n.is_power_of_two(){
        return Err(SRError::IsNotPowerOfTwo);
    }
    let mut m : MatrixF64 = Matrix::new(n, n);
    // Take the log of n with respect to 2.
    let order = n.trailing_zeros();
    // We are going to use Sylvester's construction
    // http://en.wikipedia.org/wiki/Hadamard_matrix


    // Let's fill the first level Hadamard matrix
    m.set(0, 0, 1.0);
    for  o in 0..order{
        // We will construct four views.
        let size : usize = 2i32.pow(o) as usize;
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


#[doc="Returns a Hilbert matrix.
"]
pub fn hilbert(n : usize) -> MatrixF64{
    let mut m : MatrixF64 = Matrix::new(n, n);
    for r in 0..n{
        for c in 0..n{
            let l = (r + c + 1) as f64;
            m.set(r, c, 1.0 / l);
        }
    }
    m
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


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

        use scirust::api::{MatrixF64, from_range_cw, Shape};
        let start  = 0.0;
        let stop = 16.0;
        let m : MatrixF64 = from_range_cw(4, 4, start, stop);
        for i in 0..16{
            let c = i >> 2;
            let r = i & 3;
            assert_eq!(m.get(r, c).unwrap(), i as f64);
        }


"]
pub fn from_range_cw<T:MagmaBase+Num>(rows : usize, cols : usize, 
    start : T, stop : T )-> Matrix<T> {
    // TODO this is not working.
    //let m : Matrix<T> = Matrix::from_iter_cw(rows, cols, range);
    let mut m : Matrix<T> = Matrix::new(rows, cols);
    let mut cur = start;
    'outer: for c in 0..cols{
        for r in 0..rows{
            m.set(r, c, cur);
            cur = cur + One::one();
            if cur == stop {
                break 'outer;
            }
        }
    }
    m 
}


#[doc="Returns a 64-bit floating point matrix whose entries are
picked up from a range in column wise order.
"]
#[inline]
pub fn from_range_cw_f64(rows : usize, cols : usize, 
    start : f64, stop : f64)->MatrixF64 {
    from_range_cw(rows, cols, start, stop)
}


#[doc="Returns an 8-bit signed integer matrix whose entries are
picked up from a range in column wise order.
"]
#[inline]
pub fn from_range_cw_i8(rows : usize, cols : usize, 
    start : i8, stop : i8)->MatrixI8 {
    from_range_cw(rows, cols, start, stop)
}


#[doc="Returns an 16-bit signed integer matrix whose entries are
picked up from a range in column wise order.
"]
#[inline]
pub fn from_range_cw_i16(rows : usize, cols : usize, 
    start : i16, stop : i16)->MatrixI16 {
    from_range_cw(rows, cols, start, stop)
}

#[doc="Returns an 32-bit signed integer matrix whose entries are
picked up from a range in column wise order.
"]
#[inline]
pub fn from_range_cw_i32(rows : usize, cols : usize, 
    start : i32, stop : i32)->MatrixI32 {
    from_range_cw(rows, cols, start, stop)
}


#[doc="Returns an 64-bit signed integer matrix whose entries are
picked up from a range in column wise order.

See from_range_cw function  for further discussion.


# Examples

    use scirust::api::{from_range_cw_i64, Shape};

    let m = from_range_cw_i64(4, 4, 0, 16);
    for i in 0..16{
        let c = i >> 2;
        let r = i & 3;
        assert_eq!(m.get(r, c).unwrap(), i as i64);
    }
"]
#[inline]
pub fn from_range_cw_i64(rows : usize, cols : usize, 
    start : i64, stop : i64)->MatrixI64 {
    from_range_cw(rows, cols, start, stop)
}

#[doc="Returns an 8-bit unsigned integer matrix whose entries are
picked up from a range in column wise order.
"]
#[inline]
pub fn from_range_cw_u8(rows : usize, cols : usize, 
    start : u8, stop : u8)->MatrixU8 {
    from_range_cw(rows, cols, start, stop)
}


#[doc="Returns an 16-bit unsigned integer matrix whose entries are
picked up from a range in column wise order.
"]
#[inline]
pub fn from_range_cw_u16(rows : usize, cols : usize, 
    start : u16, stop : u16)->MatrixU16 {
    from_range_cw(rows, cols, start, stop)
}

#[doc="Returns an 32-bit unsigned integer matrix whose entries are
picked up from a range in column wise order.
"]
#[inline]
pub fn from_range_cw_u32(rows : usize, cols : usize, 
    start : u32, stop : u32)->MatrixU32 {
    from_range_cw(rows, cols, start, stop)
}


#[doc="Returns an 64-bit unsigned integer matrix whose entries are
picked up from a range in column wise order.

See from_range_cw function  for further discussion.


# Examples

    use scirust::api::{from_range_cw_u64, Shape};

    let m = from_range_cw_u64(4, 4, 0, 16);
    for i in 0..16{
        let c = i >> 2;
        let r = i & 3;
        assert_eq!(m.get(r, c).unwrap(), i as u64);
    }
"]
#[inline]
pub fn from_range_cw_u64(rows : usize, cols : usize, 
    start : u64, stop : u64)->MatrixU64 {
    from_range_cw(rows, cols, start, stop)
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


#[doc="Returns a matrix whose entries are picked up from
a range in row wise order.
"]
pub fn from_range_rw<T:MagmaBase+Num>(rows : usize, cols : usize, 
    start : T, stop : T )-> Matrix<T> {
    // TODO make it work.
    //let m : Matrix<T> = Matrix::from_iter_rw(rows, cols, start..stop);
    let mut m : Matrix<T> = Matrix::new(rows, cols);
    let mut cur = start;
    'outer: for r in 0..rows{
        for c in 0..cols{
            m.set(r, c, cur);
            cur = cur + One::one();
            if cur == stop {
                break 'outer;
            }
        }
    }
    m 
}


#[doc="Returns an 8-bit signed integer matrix whose entries are
picked up from a range in row wise order.
"]
#[inline]
pub fn from_range_rw_i8(rows : usize, cols : usize, 
    start : i8, stop : i8)->MatrixI8 {
    from_range_rw(rows, cols, start, stop)
}


#[doc="Returns an 16-bit signed integer matrix whose entries are
picked up from a range in row wise order.
"]
#[inline]
pub fn from_range_rw_i16(rows : usize, cols : usize, 
    start : i16, stop : i16)->MatrixI16 {
    from_range_rw(rows, cols, start, stop)
}

#[doc="Returns an 32-bit signed integer matrix whose entries are
picked up from a range in row wise order.
"]
#[inline]
pub fn from_range_rw_i32(rows : usize, cols : usize, 
    start : i32, stop : i32)->MatrixI32 {
    from_range_rw(rows, cols, start, stop)
}


#[doc="Returns an 64-bit signed integer matrix whose entries are
picked up from a range in row wise order.

See from_range_rw function  for further discussion.


# Examples

    use scirust::api::{from_range_rw_i64, Shape};

    let m = from_range_rw_i64(4, 4, 0, 16);
    for i in 0..16{
        let r = i >> 2;
        let c = i & 3;
        assert_eq!(m.get(r, c).unwrap(), i as i64);
    }
"]
#[inline]
pub fn from_range_rw_i64(rows : usize, cols : usize, 
    start : i64, stop : i64)->MatrixI64 {
    from_range_rw(rows, cols, start, stop)
}

#[doc="Returns an 8-bit unsigned integer matrix whose entries are
picked up from a range in row wise order.
"]
#[inline]
pub fn from_range_rw_u8(rows : usize, cols : usize, 
    start : u8, stop : u8)->MatrixU8 {
    from_range_rw(rows, cols, start, stop)
}


#[doc="Returns an 16-bit unsigned integer matrix whose entries are
picked up from a range in row wise order.
"]
#[inline]
pub fn from_range_rw_u16(rows : usize, cols : usize, 
    start : u16, stop : u16)->MatrixU16 {
    from_range_rw(rows, cols, start, stop)
}

#[doc="Returns an 32-bit unsigned integer matrix whose entries are
picked up from a range in row wise order.
"]
#[inline]
pub fn from_range_rw_u32(rows : usize, cols : usize, 
    start : u32, stop : u32)->MatrixU32 {
    from_range_rw(rows, cols, start, stop)
}


#[doc="Returns an 64-bit unsigned integer matrix whose entries are
picked up from a range in row wise order.

See from_range_rw function  for further discussion.


# Examples

    use scirust::api::{from_range_rw_u64, Shape};

    let m = from_range_rw_u64(4, 4, 0, 16);
    for i in 0..16{
        let r = i >> 2;
        let c = i & 3;
        assert_eq!(m.get(r, c).unwrap(), i as u64);
    }
"]
#[inline]
pub fn from_range_rw_u64(rows : usize, cols : usize, 
    start : u64, stop : u64)->MatrixU64 {
    from_range_rw(rows, cols, start, stop)
}

#[doc="Returns a 64-bit floating point matrix whose entries are
picked up from a range in row wise order.
"]
#[inline]
pub fn from_range_rw_f64(rows : usize, cols : usize, 
    start : f64, stop : f64)->MatrixF64 {
    from_range_rw(rows, cols, start, stop)
}



#[doc="Returns a 32-bit floating point matrix whose entries are
picked up from a range in row wise order.
"]
#[inline]
pub fn from_range_rw_f32(rows : usize, cols : usize, 
    start : f32, stop : f32)->MatrixF32 {
    from_range_rw(rows, cols, start, stop)
}



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////



#[doc="Returns an 8-bit unsigned integer matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_cw_u8(rows : usize, cols : usize, values: &[u8])->MatrixU8 {
    Matrix::from_slice_cw(rows, cols, values)
}

#[doc="Returns a 16-bit unsigned integer matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_cw_u16(rows : usize, cols : usize, values: &[u16])->MatrixU16 {
    Matrix::from_slice_cw(rows, cols, values)
}

#[doc="Returns a 32-bit unsigned integer matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_cw_u32(rows : usize, cols : usize, values: &[u32])->MatrixU32 {
    Matrix::from_slice_cw(rows, cols, values)
}

#[doc="Returns a 64-bit unsigned integer matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_cw_u64(rows : usize, cols : usize, values: &[u64])->MatrixU64 {
    Matrix::from_slice_cw(rows, cols, values)
}

#[doc="Returns an 8-bit signed integer matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_cw_i8(rows : usize, cols : usize, values: &[i8])->MatrixI8 {
    Matrix::from_slice_cw(rows, cols, values)
}

#[doc="Returns a 16-bit signed integer matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_cw_i16(rows : usize, cols : usize, values: &[i16])->MatrixI16 {
    Matrix::from_slice_cw(rows, cols, values)
}

#[doc="Returns a 32-bit signed integer matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_cw_i32(rows : usize, cols : usize, values: &[i32])->MatrixI32 {
    Matrix::from_slice_cw(rows, cols, values)
}

#[doc="Returns a 64-bit signed integer matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_cw_i64(rows : usize, cols : usize, values: &[i64])->MatrixI64 {
    Matrix::from_slice_cw(rows, cols, values)
}

#[doc="Returns a 32-bit float matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_cw_f32(rows : usize, cols : usize, values: &[f32])->MatrixF32 {
    Matrix::from_slice_cw(rows, cols, values)
}


#[doc="Returns a 64-bit float matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_cw_f64(rows : usize, cols : usize, values: &[f64])->MatrixF64 {
    Matrix::from_slice_cw(rows, cols, values)
}


#[doc="Returns a 32-bit complex matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_cw_c32(rows : usize, cols : usize, values: &[Complex32])->MatrixC32 {
    Matrix::from_slice_cw(rows, cols, values)
}


#[doc="Returns a 64-bit complex matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_cw_c64(rows : usize, cols : usize, values: &[Complex64])->MatrixC64 {
    Matrix::from_slice_cw(rows, cols, values)
}

#[doc="Returns an 8-bit unsigned integer matrix whose entries are
picked up from a slice in row wise order.
"]
#[inline]
pub fn matrix_rw_u8(rows : usize, cols : usize, values: &[u8])->MatrixU8 {
    Matrix::from_slice_rw(rows, cols, values)
}

#[doc="Returns a 16-bit unsigned integer matrix whose entries are
picked up from a slice in row wise order.
"]
#[inline]
pub fn matrix_rw_u16(rows : usize, cols : usize, values: &[u16])->MatrixU16 {
    Matrix::from_slice_rw(rows, cols, values)
}

#[doc="Returns a 32-bit unsigned integer matrix whose entries are
picked up from a slice in row wise order.
"]
#[inline]
pub fn matrix_rw_u32(rows : usize, cols : usize, values: &[u32])->MatrixU32 {
    Matrix::from_slice_rw(rows, cols, values)
}

#[doc="Returns a 64-bit unsigned integer matrix whose entries are
picked up from a slice in row wise order.
"]
#[inline]
pub fn matrix_rw_u64(rows : usize, cols : usize, values: &[u64])->MatrixU64 {
    Matrix::from_slice_rw(rows, cols, values)
}

#[doc="Returns an 8-bit signed integer matrix whose entries are
picked up from a slice in row wise order.
"]
#[inline]
pub fn matrix_rw_i8(rows : usize, cols : usize, values: &[i8])->MatrixI8 {
    Matrix::from_slice_rw(rows, cols, values)
}

#[doc="Returns a 16-bit signed integer matrix whose entries are
picked up from a slice in row wise order.
"]
#[inline]
pub fn matrix_rw_i16(rows : usize, cols : usize, values: &[i16])->MatrixI16 {
    Matrix::from_slice_rw(rows, cols, values)
}

#[doc="Returns a 32-bit signed integer matrix whose entries are
picked up from a slice in row wise order.
"]
#[inline]
pub fn matrix_rw_i32(rows : usize, cols : usize, values: &[i32])->MatrixI32 {
    Matrix::from_slice_rw(rows, cols, values)
}

#[doc="Returns a 64-bit signed integer matrix whose entries are
picked up from a slice in row wise order.
"]
#[inline]
pub fn matrix_rw_i64(rows : usize, cols : usize, values: &[i64])->MatrixI64 {
    Matrix::from_slice_rw(rows, cols, values)
}

#[doc="Returns a 32-bit float matrix whose entries are
picked up from a slice in row wise order.
"]
#[inline]
pub fn matrix_rw_f32(rows : usize, cols : usize, values: &[f32])->MatrixF32 {
    Matrix::from_slice_rw(rows, cols, values)
}


#[doc="Returns a 64-bit float matrix whose entries are
picked up from a slice in row wise order.
"]
#[inline]
pub fn matrix_rw_f64(rows : usize, cols : usize, values: &[f64])->MatrixF64 {
    Matrix::from_slice_rw(rows, cols, values)
}

/*
#[doc="Returns a 32-bit complex matrix whose entries are
picked up from a slice in row wise order.
"]
#[inline]
pub fn matrix_rw_c32(rows : usize, cols : usize, values: &[Complex32])->MatrixC32 {
    Matrix::from_slice_rw(rows, cols, values)
}


#[doc="Returns a 64-bit complex matrix whose entries are
picked up from a slice in row wise order.
"]
#[inline]
pub fn matrix_rw_c64(rows : usize, cols : usize, values: &[Complex64])->MatrixC64 {
    Matrix::from_slice_rw(rows, cols, values)
}
*/

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


#[doc="Returns a column vector with entries from a slice.
"]
#[inline]
pub fn col_vector<T:CommutativeMonoidAddPartial>(values: &[T])-> Matrix<T> {
    let m : Matrix<T> = Matrix::from_slice_cw(values.len(), 1, values);
    m 
}

#[doc="Returns a column vector with entries from an iterator.
"]
#[inline]
pub fn col_vector_from_iter<T:CommutativeMonoidAddPartial, A : Iterator<Item=T>>(
    values: A,
    len : usize)-> Matrix<T> {
    let m : Matrix<T> = Matrix::from_iter_rw(len, 1, values);
    m 
}



#[doc="Returns a 8-bit unsigned int column vector with entries from a slice.
"]
#[inline]
pub fn vector_u8(values: &[u8])->MatrixU8 {
    col_vector(values)
}

#[doc="Returns a 16-bit unsigned int column vector with entries from a slice.
"]
#[inline]
pub fn vector_u16(values: &[u16])->MatrixU16 {
    col_vector(values)
}

#[doc="Returns a 32-bit unsigned int column vector with entries from a slice.
"]
#[inline]
pub fn vector_u32(values: &[u32])->MatrixU32 {
    col_vector(values)
}


#[doc="Returns a 64-bit unsigned int column vector with entries from a slice.
"]
#[inline]
pub fn vector_u64(values: &[u64])->MatrixU64 {
    col_vector(values)
}

#[doc="Returns an 8-bit signed int column vector with entries from a slice.
"]
#[inline]
pub fn vector_i8(values: &[i8])->MatrixI8 {
    col_vector(values)
}


#[doc="Returns a 16-bit signed int column vector with entries from a slice.
"]
#[inline]
pub fn vector_i16(values: &[i16])->MatrixI16 {
    col_vector(values)
}

#[doc="Returns a 32-bit signed int column vector with entries from a slice.
"]
#[inline]
pub fn vector_i32(values: &[i32])->MatrixI32 {
    col_vector(values)
}


#[doc="Returns a 64-bit signed int column vector with entries from a slice.
"]
#[inline]
pub fn vector_i64(values: &[i64])->MatrixI64 {
    col_vector(values)
}


#[doc="Returns a 32-bit float column vector with entries from a slice.
"]
#[inline]
pub fn vector_f32(values: &[f32])->MatrixF32 {
    col_vector(values)
}



#[doc="Returns a 64-bit float column vector with entries from a slice.
"]
#[inline]
pub fn vector_f64(values: &[f64])->MatrixF64 {
    col_vector(values)
}

//pub fn vec_f64<T:Number+ToPrimitive>(values: &[T])->MatrixF64 {
//    let n = values.len();
//    let iter = values.iter().map(|&x| x.to_f64().unwrap());
//    col_vector_from_iter(iter, n)
//}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
//              Elementary row operation matrices
//
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

#[doc="Returns elementary matrix which can 
exchange rows i and j
on left multiplication.
"]
pub fn ero_switch<T:FieldPartial>(n : usize, 
    i  : usize, 
    j : usize)-> Matrix<T> {
    debug_assert! (i  < n);
    debug_assert! (j  < n);
    let mut m : Matrix<T> = Matrix::identity(n, n);
    let z : T = Zero::zero();
    let o : T = One::one();
    m.set(i, i, z);
    m.set(j, j, z);
    m.set(i, j, o);
    m.set(j, i, o);
    m
}

#[doc="Returns elementary matrix which can scale
a particular row by a factor on left multiplication.
"]
pub fn ero_scale<T:FieldPartial>(n : usize, 
    r  : usize, 
    scale : T)-> Matrix<T> {

    let mut m : Matrix<T> = Matrix::identity(n, n);
    m.set(r,r, scale);
    m
}


#[doc="Returns elementary matrix which can scale
a particular row by a factor and add it to
another row
on left multiplication.
r_i = r_i + k * r_j

"]
pub fn ero_scale_add<T:FieldPartial>(n : usize, 
    i  : usize, 
    j : usize, 
    scale : T)-> Matrix<T> {

    let mut m : Matrix<T> = Matrix::identity(n, n);
    m.set(i, j, scale);
    m
}



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
//              Unit tests follow
//
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod test{
    use api::*;

    #[test]
    fn test_hadamard(){
        let m = hadamard(1).unwrap();
        assert_eq!(m.num_cells(), 1);
        assert_eq!(m.get(0,0).unwrap(), 1.0);
        let m = hadamard(2).unwrap();
        assert_eq!(m.num_cells(), 4);
        assert_eq!(m.get(0,0).unwrap(), 1.0);
        assert_eq!(m.get(0,1).unwrap(), 1.0);
        assert_eq!(m.get(1,0).unwrap(), 1.0);
        assert_eq!(m.get(1,1).unwrap(), -1.0);
        let m = hadamard(4).unwrap();
        assert!(m.is_square());
    }

    #[test]
    fn test_range_cw_functions(){
        let m = from_range_cw_i64(4, 4, 0, 16);
        for i in 0..16{
            let c = i >> 2;
            let r = i & 3;
            assert_eq!(m.get(r, c).unwrap(), i as i64);
        }
        let start  = 0.0;
        let stop = 16.0;
        let m : MatrixF64 = from_range_cw(4, 4, start, stop);
        for i in 0..16{
            let c = i >> 2;
            let r = i & 3;
            assert_eq!(m.get(r, c).unwrap(), i as f64);
        }
        let m = from_range_cw_f64(4, 4, start, stop);
        for i in 0..16{
            let c = i >> 2;
            let r = i & 3;
            assert_eq!(m.get(r, c).unwrap(), i as f64);
        }
    }

    #[test]
    fn test_range_rw_functions(){
        let start  = 0.0;
        let stop = 16.0;
        let m : MatrixF64 = from_range_rw(4, 4, start, stop);
        for i in 0..16{
            let r = i >> 2;
            let c = i & 3;
            assert_eq!(m.get(r, c).unwrap(), i as f64);
        }
        let m = from_range_rw_i8(4, 4, 0, 16);
        for i in 0..16{
            let r = i >> 2;
            let c = i & 3;
            assert_eq!(m.get(r, c).unwrap(), i as i8);
        }
        let m = from_range_rw_i16(4, 4, 0, 16);
        for i in 0..16{
            let r = i >> 2;
            let c = i & 3;
            assert_eq!(m.get(r, c).unwrap(), i as i16);
        }
        let m = from_range_rw_i32(4, 4, 0, 16);
        for i in 0..16{
            let r = i >> 2;
            let c = i & 3;
            assert_eq!(m.get(r, c).unwrap(), i as i32);
        }
        let m = from_range_rw_i64(4, 4, 0, 16);
        for i in 0..16{
            let r = i >> 2;
            let c = i & 3;
            assert_eq!(m.get(r, c).unwrap(), i as i64);
        }


        let m = from_range_rw_u8(4, 4, 0, 16);
        for i in 0..16{
            let r = i >> 2;
            let c = i & 3;
            assert_eq!(m.get(r, c).unwrap(), i as u8);
        }
        let m = from_range_rw_u16(4, 4, 0, 16);
        for i in 0..16{
            let r = i >> 2;
            let c = i & 3;
            assert_eq!(m.get(r, c).unwrap(), i as u16);
        }
        let m = from_range_rw_u32(4, 4, 0, 16);
        for i in 0..16{
            let r = i >> 2;
            let c = i & 3;
            assert_eq!(m.get(r, c).unwrap(), i as u32);
        }
        let m = from_range_rw_u64(4, 4, 0, 16);
        for i in 0..16{
            let r = i >> 2;
            let c = i & 3;
            assert_eq!(m.get(r, c).unwrap(), i as u64);
        }


        let m = from_range_rw_f64(4, 4, start, stop);
        for i in 0..16{
            let r = i >> 2;
            let c = i & 3;
            assert_eq!(m.get(r, c).unwrap(), i as f64);
        }

        let m = from_range_rw_f32(4, 4, 0.0, 16.0);
        for i in 0..16{
            let r = i >> 2;
            let c = i & 3;
            assert_eq!(m.get(r, c).unwrap(), i as f32);
        }

    }


    #[test]
    fn test_matrix_type_functions(){
        let m = matrix_cw_u8(2,2, &[1,2,3,4]);
        for i in 0..4{
            let c = i >> 1;
            let r = i & 1;            
            assert_eq!(m.get(r, c).unwrap(), (i + 1) as u8);
        }
        let m = matrix_cw_u16(2,2, &[1,2,3,4]);
        for i in 0..4{
            let c = i >> 1;
            let r = i & 1;            
            assert_eq!(m.get(r, c).unwrap(), (i + 1) as u16);
        }
        let m = matrix_cw_u32(2,2, &[1,2,3,4]);
        for i in 0..4{
            let c = i >> 1;
            let r = i & 1;            
            assert_eq!(m.get(r, c).unwrap(), (i + 1) as u32);
        }
        let m = matrix_cw_u64(2,2, &[1,2,3,4]);
        for i in 0..4{
            let c = i >> 1;
            let r = i & 1;            
            assert_eq!(m.get(r, c).unwrap(), (i + 1) as u64);
        }
        let m = matrix_cw_i8(2,2, &[1,2,3,4]);
        for i in 0..4{
            let c = i >> 1;
            let r = i & 1;            
            assert_eq!(m.get(r, c).unwrap(), (i + 1) as i8);
        }
        let m = matrix_cw_i16(2,2, &[1,2,3,4]);
        for i in 0..4{
            let c = i >> 1;
            let r = i & 1;            
            assert_eq!(m.get(r, c).unwrap(), (i + 1) as i16);
        }
        let m = matrix_cw_i32(2,2, &[1,2,3,4]);
        for i in 0..4{
            let c = i >> 1;
            let r = i & 1;            
            assert_eq!(m.get(r, c).unwrap(), (i + 1) as i32);
        }
        let m = matrix_cw_i64(2,2, &[1,2,3,4]);
        for i in 0..4{
            let c = i >> 1;
            let r = i & 1;            
            assert_eq!(m.get(r, c).unwrap(), (i + 1) as i64);
        }
        let m = matrix_cw_f64(2,2, &[1.0,2.0,3.0,4.0]);
        for i in 0..4{
            let c = i >> 1;
            let r = i & 1;            
            assert_eq!(m.get(r, c).unwrap(), (i + 1) as f64);
        }

        //  We will now test row wise construction functions.


        let m = matrix_rw_u8(2,2, &[1,2,3,4]);
        for i in 0..4{
            let r = i >> 1;
            let c = i & 1;            
            assert_eq!(m.get(r, c).unwrap(), (i + 1) as u8);
        }
        let m = matrix_rw_u16(2,2, &[1,2,3,4]);
        for i in 0..4{
            let r = i >> 1;
            let c = i & 1;            
            assert_eq!(m.get(r, c).unwrap(), (i + 1) as u16);
        }
        let m = matrix_rw_u32(2,2, &[1,2,3,4]);
        for i in 0..4{
            let r = i >> 1;
            let c = i & 1;            
            assert_eq!(m.get(r, c).unwrap(), (i + 1) as u32);
        }
        let m = matrix_rw_u64(2,2, &[1,2,3,4]);
        for i in 0..4{
            let r = i >> 1;
            let c = i & 1;            
            assert_eq!(m.get(r, c).unwrap(), (i + 1) as u64);
        }
        let m = matrix_rw_i8(2,2, &[1,2,3,4]);
        for i in 0..4{
            let r = i >> 1;
            let c = i & 1;            
            assert_eq!(m.get(r, c).unwrap(), (i + 1) as i8);
        }
        let m = matrix_rw_i16(2,2, &[1,2,3,4]);
        for i in 0..4{
            let r = i >> 1;
            let c = i & 1;            
            assert_eq!(m.get(r, c).unwrap(), (i + 1) as i16);
        }
        let m = matrix_rw_i32(2,2, &[1,2,3,4]);
        for i in 0..4{
            let r = i >> 1;
            let c = i & 1;            
            assert_eq!(m.get(r, c).unwrap(), (i + 1) as i32);
        }
        let m = matrix_rw_i64(2,2, &[1,2,3,4]);
        for i in 0..4{
            let r = i >> 1;
            let c = i & 1;            
            assert_eq!(m.get(r, c).unwrap(), (i + 1) as i64);
        }
        let m = matrix_rw_f64(2,2, &[1.0,2.0,3.0,4.0]);
        for i in 0..4{
            let r = i >> 1;
            let c = i & 1;            
            assert_eq!(m.get(r, c).unwrap(), (i + 1) as f64);
        }
    }

    #[test]
    fn test_vector_type_functions(){
        let v = vector_u8(&[1,2,3,4]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in 0..4{
            assert_eq!(v.get(i, 0).unwrap(), (i + 1) as u8);
        }

        let v = vector_u16(&[1,2,3,4]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in 0..4{
            assert_eq!(v.get(i, 0).unwrap(), (i + 1) as u16);
        }

        let v = vector_u32(&[1,2,3,4]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in 0..4{
            assert_eq!(v.get(i, 0).unwrap(), (i + 1) as u32);
        }

        let v = vector_u64(&[1,2,3,4]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in 0..4{
            assert_eq!(v.get(i, 0).unwrap(), (i + 1) as u64);
        }

        let v = vector_i8(&[1,2,3,4]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in 0..4{
            assert_eq!(v.get(i, 0).unwrap(), (i + 1) as i8);
        }

        let v = vector_i16(&[1,2,3,4]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in 0..4{
            assert_eq!(v.get(i, 0).unwrap(), (i + 1) as i16);
        }

        let v = vector_i32(&[1,2,3,4]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in 0..4{
            assert_eq!(v.get(i, 0).unwrap(), (i + 1) as i32);
        }

        let v = vector_i64(&[1,2,3,4]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in 0..4{
            assert_eq!(v.get(i, 0).unwrap(), (i + 1) as i64);
        }

        let v = vector_f32(&[1.,2.,3.,4.]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in 0..4{
            assert_eq!(v.get(i, 0).unwrap(), (i + 1) as f32);
        }

        let v = vector_f64(&[1.,2.,3.,4.]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in 0..4{
            assert_eq!(v.get(i, 0).unwrap(), (i + 1) as f64);
        }

    }

    #[test]
    fn test_ero_switch_scale(){
        let eswitch : MatrixF64 = ero_switch(4, 1, 3);
        let escale : MatrixF64 = ero_scale(4, 2, 2.0);
        let mut m = matrix_rw_f64(4,4, &[0., 1., 2., 3., 
            4., 5., 6., 7.,
            8., 9., 10., 11.,
            12., 13., 14., 15.]);
        // Carry out transformation through multiplying
        // elementary matrices 
        let m2 = &eswitch * &m;
        let m3 = &escale * &m2;
        // Do ERO operations directly.
        m.ero_switch(1, 3);
        m.ero_scale(2, 2.0);
        println!("eswitch: {:?}", eswitch);
        println!("escale: {:?}", escale);
        println!("m2: {:?}", m2);
        println!("m3: {:?}", m3);
        assert_eq!(m3, m);
    }


    #[test]
    fn test_ero_scale_add(){
        let mut m = matrix_rw_f64(4,4, &[0., 1., 2., 3., 
            4., 5., 6., 7.,
            8., 9., 10., 11.,
            12., 13., 14., 15.]);
        let esa : MatrixF64 = ero_scale_add(4, 1, 2, 3.0);
        println!("esa: {:?}", esa);
        let m2 = &esa * &m;
        println!("m2: {:?}", m2);
        m.ero_scale_add(1, 2, 3.);
        assert_eq!(m2, m);
    }
}
