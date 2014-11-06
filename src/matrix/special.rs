#![doc="Provides functions to create special matrices.
"]

// std imports
use std::num;


// local imports
use matrix::matrix::*;
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
        use scirust::matrix::from_range;
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


#[doc="Returns a 64-bit floating point matrix whose entries are
picked up from a range in column wise order.
"]
#[inline]
pub fn from_range_f64(rows : uint, cols : uint, 
    start : f64, stop : f64)->MatrixF64 {
    from_range(rows, cols, start, stop)
}


#[doc="Returns an 8-bit signed integer matrix whose entries are
picked up from a range in column wise order.
"]
#[inline]
pub fn from_range_i8(rows : uint, cols : uint, 
    start : i8, stop : i8)->MatrixI8 {
    from_range(rows, cols, start, stop)
}


#[doc="Returns an 16-bit signed integer matrix whose entries are
picked up from a range in column wise order.
"]
#[inline]
pub fn from_range_i16(rows : uint, cols : uint, 
    start : i16, stop : i16)->MatrixI16 {
    from_range(rows, cols, start, stop)
}

#[doc="Returns an 32-bit signed integer matrix whose entries are
picked up from a range in column wise order.
"]
#[inline]
pub fn from_range_i32(rows : uint, cols : uint, 
    start : i32, stop : i32)->MatrixI32 {
    from_range(rows, cols, start, stop)
}


#[doc="Returns an 64-bit signed integer matrix whose entries are
picked up from a range in column wise order.

See from_range function  for further discussion.


# Examples

    use scirust::matrix::from_range_i64;

    let m = from_range_i64(4, 4, 0, 16);
    for i in range(0, 16){
        let c = i >> 2;
        let r = i & 3;
        assert_eq!(m.get(r, c), i as i64);
    }
"]
#[inline]
pub fn from_range_i64(rows : uint, cols : uint, 
    start : i64, stop : i64)->MatrixI64 {
    from_range(rows, cols, start, stop)
}

#[doc="Returns an 8-bit unsigned integer matrix whose entries are
picked up from a range in column wise order.
"]
#[inline]
pub fn from_range_u8(rows : uint, cols : uint, 
    start : u8, stop : u8)->MatrixU8 {
    from_range(rows, cols, start, stop)
}


#[doc="Returns an 16-bit unsigned integer matrix whose entries are
picked up from a range in column wise order.
"]
#[inline]
pub fn from_range_u16(rows : uint, cols : uint, 
    start : u16, stop : u16)->MatrixU16 {
    from_range(rows, cols, start, stop)
}

#[doc="Returns an 32-bit unsigned integer matrix whose entries are
picked up from a range in column wise order.
"]
#[inline]
pub fn from_range_u32(rows : uint, cols : uint, 
    start : u32, stop : u32)->MatrixU32 {
    from_range(rows, cols, start, stop)
}


#[doc="Returns an 64-bit unsigned integer matrix whose entries are
picked up from a range in column wise order.

See from_range function  for further discussion.


# Examples

    use scirust::matrix::from_range_u64;

    let m = from_range_u64(4, 4, 0, 16);
    for i in range(0, 16){
        let c = i >> 2;
        let r = i & 3;
        assert_eq!(m.get(r, c), i as u64);
    }
"]
#[inline]
pub fn from_range_u64(rows : uint, cols : uint, 
    start : u64, stop : u64)->MatrixU64 {
    from_range(rows, cols, start, stop)
}


#[doc="Returns an 8-bit unsigned integer matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_u8(rows : uint, cols : uint, values: &[u8])->MatrixU8 {
    Matrix::from_slice(rows, cols, values)
}

#[doc="Returns a 16-bit unsigned integer matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_u16(rows : uint, cols : uint, values: &[u16])->MatrixU16 {
    Matrix::from_slice(rows, cols, values)
}

#[doc="Returns a 32-bit unsigned integer matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_u32(rows : uint, cols : uint, values: &[u32])->MatrixU32 {
    Matrix::from_slice(rows, cols, values)
}

#[doc="Returns a 64-bit unsigned integer matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_u64(rows : uint, cols : uint, values: &[u64])->MatrixU64 {
    Matrix::from_slice(rows, cols, values)
}

#[doc="Returns an 8-bit signed integer matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_i8(rows : uint, cols : uint, values: &[i8])->MatrixI8 {
    Matrix::from_slice(rows, cols, values)
}

#[doc="Returns a 16-bit signed integer matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_i16(rows : uint, cols : uint, values: &[i16])->MatrixI16 {
    Matrix::from_slice(rows, cols, values)
}

#[doc="Returns a 32-bit signed integer matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_i32(rows : uint, cols : uint, values: &[i32])->MatrixI32 {
    Matrix::from_slice(rows, cols, values)
}

#[doc="Returns a 64-bit signed integer matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_i64(rows : uint, cols : uint, values: &[i64])->MatrixI64 {
    Matrix::from_slice(rows, cols, values)
}

#[doc="Returns a 32-bit float matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_f32(rows : uint, cols : uint, values: &[f32])->MatrixF32 {
    Matrix::from_slice(rows, cols, values)
}


#[doc="Returns a 64-bit float matrix whose entries are
picked up from a slice in column wise order.
"]
#[inline]
pub fn matrix_f64(rows : uint, cols : uint, values: &[f64])->MatrixF64 {
    Matrix::from_slice(rows, cols, values)
}


#[doc="Returns a column vector with entries from a slice.
"]
#[inline]
pub fn col_vector<T:MatrixElt>(values: &[T])-> Matrix<T> {
    let m : Matrix<T> = Matrix::from_slice(values.len(), 1, values);
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
        let m = from_range_f64(4, 4, start, stop);
        for i in range(0, 16){
            let c = i >> 2;
            let r = i & 3;
            assert_eq!(m.get(r, c), i as f64);
        }
    }

    #[test]
    fn test_matrix_type_functions(){
        let m = matrix_u8(2,2, [1,2,3,4]);
        for i in range(0, 4){
            let c = i >> 1;
            let r = i & 1;            
            assert_eq!(m.get(r, c), (i + 1) as u8);
        }
        let m = matrix_u16(2,2, [1,2,3,4]);
        for i in range(0, 4){
            let c = i >> 1;
            let r = i & 1;            
            assert_eq!(m.get(r, c), (i + 1) as u16);
        }
        let m = matrix_u32(2,2, [1,2,3,4]);
        for i in range(0, 4){
            let c = i >> 1;
            let r = i & 1;            
            assert_eq!(m.get(r, c), (i + 1) as u32);
        }
        let m = matrix_u64(2,2, [1,2,3,4]);
        for i in range(0, 4){
            let c = i >> 1;
            let r = i & 1;            
            assert_eq!(m.get(r, c), (i + 1) as u64);
        }
        let m = matrix_i8(2,2, [1,2,3,4]);
        for i in range(0, 4){
            let c = i >> 1;
            let r = i & 1;            
            assert_eq!(m.get(r, c), (i + 1) as i8);
        }
        let m = matrix_i16(2,2, [1,2,3,4]);
        for i in range(0, 4){
            let c = i >> 1;
            let r = i & 1;            
            assert_eq!(m.get(r, c), (i + 1) as i16);
        }
        let m = matrix_i32(2,2, [1,2,3,4]);
        for i in range(0, 4){
            let c = i >> 1;
            let r = i & 1;            
            assert_eq!(m.get(r, c), (i + 1) as i32);
        }
        let m = matrix_i64(2,2, [1,2,3,4]);
        for i in range(0, 4){
            let c = i >> 1;
            let r = i & 1;            
            assert_eq!(m.get(r, c), (i + 1) as i64);
        }
        let m = matrix_f64(2,2, [1.0,2.0,3.0,4.0]);
        for i in range(0, 4){
            let c = i >> 1;
            let r = i & 1;            
            assert_eq!(m.get(r, c), (i + 1) as f64);
        }
    }

    #[test]
    fn test_vector_type_functions(){
        let v = vector_u8([1,2,3,4]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in range(0, 4){
            assert_eq!(v.get(i, 0), (i + 1) as u8);
        }

        let v = vector_u16([1,2,3,4]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in range(0, 4){
            assert_eq!(v.get(i, 0), (i + 1) as u16);
        }

        let v = vector_u32([1,2,3,4]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in range(0, 4){
            assert_eq!(v.get(i, 0), (i + 1) as u32);
        }

        let v = vector_u64([1,2,3,4]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in range(0, 4){
            assert_eq!(v.get(i, 0), (i + 1) as u64);
        }

        let v = vector_i8([1,2,3,4]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in range(0, 4){
            assert_eq!(v.get(i, 0), (i + 1) as i8);
        }

        let v = vector_i16([1,2,3,4]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in range(0, 4){
            assert_eq!(v.get(i, 0), (i + 1) as i16);
        }

        let v = vector_i32([1,2,3,4]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in range(0, 4){
            assert_eq!(v.get(i, 0), (i + 1) as i32);
        }

        let v = vector_i64([1,2,3,4]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in range(0, 4){
            assert_eq!(v.get(i, 0), (i + 1) as i64);
        }

        let v = vector_f32([1.,2.,3.,4.]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in range(0, 4){
            assert_eq!(v.get(i, 0), (i + 1) as f32);
        }

        let v = vector_f64([1.,2.,3.,4.]);
        assert!(v.is_vector());
        assert!(v.is_col());
        for i in range(0, 4){
            assert_eq!(v.get(i, 0), (i + 1) as f64);
        }

    }

}
