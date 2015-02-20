#![doc="Provides simple methods to convert slices 
of numbers from one type to another
"]

// std imports
pub use std::raw::Slice as RawSlice;
use std::slice::Iter;

pub trait SliceTypeConversions{

    fn as_f64_vector(&mut self) -> Vec<f64>;
}


impl SliceTypeConversions for RawSlice<int>{

    fn as_f64_vector(&mut self) -> Vec<f64>{
        let mut result = Vec::new();
        let p = self.data;
        for i in range(0, self.len){
            let v = unsafe {*p.offset(i as isize)};
            result.push(v as f64);
        }
        result
    }
}


impl<'a> SliceTypeConversions for Iter<'a, isize> {

    fn as_f64_vector(&mut self) -> Vec<f64>{
        let mut result = Vec::new();
        for v in *self{
            result.push(*v as f64);
        }
        result
    }
}



/******************************************************
 *
 *   Unit tests follow.
 *
 *******************************************************/

 #[cfg(test)]
mod test{
    use super::*;

    #[test]
    fn test_slice_conv_1(){
        let a = [1i,2,3,4];
        let v = a.iter().as_f64_vector();
        assert_eq!(v, vec![1., 2., 3., 4.]);
        let v = a.as_slice().iter().as_f64_vector();
        assert_eq!(v, vec![1., 2., 3., 4.]);
    }
 
}