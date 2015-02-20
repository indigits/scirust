#![doc="Sorting algorithms
"]

pub mod quick_sort;
pub mod insertion_sort;


/// Returns whether the slice is sorted in ascending order or not.
pub fn is_ascending_slice<T : PartialOrd>(data : & [T]) -> bool {
    let n = data.len();
    if n == 0 {
        return true;
    }
    for i in range(0, n-1){
        if data[i] > data[i + 1]{
            return false;
        }
    }
    true
}

/// Returns whether the slice is sorted in descending order or not.
pub fn is_descending_slice<T : PartialOrd>(data : & [T]) -> bool {
    let n = data.len();
    if n == 0 {
        return true;
    }
    for i in range(0, n-1){
        if data[i] < data[i + 1]{
            return false;
        }
    }
    true
}

/// Returns whether the data in buffer is sorted in ascending order or not.
pub unsafe fn is_ascending_buffer<T : PartialOrd>(data : *const T, n : usize) -> bool {    
    if n == 0 {
        return true;
    }
    let mut ptr = data;
    for _ in range(0, n-1){
        if *ptr > *ptr.offset(1i){
            return false;
        }
        ptr = ptr.offset(1i);
    }
    true
}

/// Returns whether the data in buffer is sorted in descending order or not.
pub unsafe fn is_descending_buffer<T : PartialOrd>(data : *const T, n : usize) -> bool {    
    if n == 0 {
        return true;
    }
    let mut ptr = data;
    for _ in range(0, n-1){
        if *ptr < *ptr.offset(1i){
            return false;
        }
        ptr = ptr.offset(1i);
    }
    true
}


/// Returns whether the data in strided buffer is sorted in ascending order or not
/// every stride-th entry is considered in the buffer
/// assumes that the buffer as at least n* stride entries
pub unsafe fn is_ascending_buffer_strided<T : PartialOrd>(data : *const T, 
    n : usize, 
    stride : usize) -> bool {    
    if n == 0 {
        return true;
    }
    let mut ptr = data;
    let stride = stride as isize;
    for _ in range(0, n-1){
        if *ptr > *ptr.offset(stride){
            return false;
        }
        ptr = ptr.offset(stride);
    }
    true
}

/// Returns whether the data in strided buffer is sorted in ascending order or not
/// every stride-th entry is considered in the buffer
/// assumes that the buffer as at least n* stride entries
pub unsafe fn is_descending_buffer_strided<T : PartialOrd>(data : *const T, 
    n : usize, 
    stride : usize) -> bool {    
    if n == 0 {
        return true;
    }
    let mut ptr = data;
    let stride = stride as isize;
    for _ in range(0, n-1){
        if *ptr < *ptr.offset(stride){
            return false;
        }
        ptr = ptr.offset(stride);
    }
    true
}


/******************************************************
 *
 *   Unit tests
 *
 *******************************************************/


#[cfg(test)]
mod test{
    use super::*;

    #[test]
    fn test_is_ascending_slice(){
        let v = vec!(1i, 2, 3, 4);
        assert!(is_ascending_slice(v.as_slice()));
        let v = vec!(1i);
        assert!(is_ascending_slice(v.as_slice()));
        let v : Vec<int> = vec!();
        assert!(is_ascending_slice(v.as_slice()));
        let v = vec!(2i, 1);
        assert!(!is_ascending_slice(v.as_slice()));
        let v = vec!(3i, 2, 1);
        assert!(!is_ascending_slice(v.as_slice()));
        let v = vec!(1i, 3, 2);
        assert!(!is_ascending_slice(v.as_slice()));
        let v = vec!(2i, 1, 3);
        assert!(!is_ascending_slice(v.as_slice()));
    }


    #[test]
    fn test_is_ascending_buffer(){
        unsafe {
            let v = vec!(1i, 2, 3, 4);
            assert!(is_ascending_buffer(v.as_ptr(), v.len()));
            let v = vec!(1i);
            assert!(is_ascending_buffer(v.as_ptr(), v.len()));
            let v : Vec<int> = vec!();
            assert!(is_ascending_buffer(v.as_ptr(), v.len()));
            let v = vec!(2i, 1);
            assert!(!is_ascending_buffer(v.as_ptr(), v.len()));
            let v = vec!(3i, 2, 1);
            assert!(!is_ascending_buffer(v.as_ptr(), v.len()));
            let v = vec!(1i, 3, 2);
            assert!(!is_ascending_buffer(v.as_ptr(), v.len()));
            let v = vec!(2i, 1, 3);
            assert!(!is_ascending_buffer(v.as_ptr(), v.len()));
        }
    }

    #[test]
    fn test_is_descending_slice(){
        let v = vec!(4i, 3, 2, 1);
        assert!(is_descending_slice(v.as_slice()));
        let v = vec!(1i);
        assert!(is_descending_slice(v.as_slice()));
        let v : Vec<int> = vec!();
        assert!(is_descending_slice(v.as_slice()));
        let v = vec!(1i, 2);
        assert!(!is_descending_slice(v.as_slice()));
        let v = vec!(1i, 2, 3);
        assert!(!is_descending_slice(v.as_slice()));
        let v = vec!(2i, 3, 1);
        assert!(!is_descending_slice(v.as_slice()));
        let v = vec!(3i, 1, 2);
        assert!(!is_descending_slice(v.as_slice()));
    }


    #[test]
    fn test_is_descending_buffer(){
        unsafe {
            let v = vec!(4i, 3, 2, 1);
            assert!(is_descending_buffer(v.as_ptr(), v.len()));
            let v = vec!(1i);
            assert!(is_descending_buffer(v.as_ptr(), v.len()));
            let v : Vec<int> = vec!();
            assert!(is_descending_buffer(v.as_ptr(), v.len()));
            let v = vec!(1i, 2);
            assert!(!is_descending_buffer(v.as_ptr(), v.len()));
            let v = vec!(1i, 2, 3);
            assert!(!is_descending_buffer(v.as_ptr(), v.len()));
            let v = vec!(2i, 3, 1);
            assert!(!is_descending_buffer(v.as_ptr(), v.len()));
            let v = vec!(3i, 1, 2);
            assert!(!is_descending_buffer(v.as_ptr(), v.len()));
        }
    }

    #[test]
    fn test_is_ascending_buffer_strided(){
        unsafe {
            let v = vec!(4i, 3, 2, 1);
            assert!(!is_ascending_buffer_strided(v.as_ptr(), v.len() / 2, 2));
            let v = vec!(4i, 3, 5, 1);
            assert!(is_ascending_buffer_strided(v.as_ptr(), v.len() / 2, 2));
            let v = vec!(4i, 1, 5, 2, 3, 7);
            assert!(!is_ascending_buffer_strided(v.as_ptr(), v.len() / 2, 2));
            let v = vec!(4i, 3, 5, 2, 7, 1);
            assert!(is_ascending_buffer_strided(v.as_ptr(), v.len() / 2, 2));
        }
    }

}


/******************************************************
 *
 *   Bench marks
 *
 *******************************************************/


#[cfg(test)]
mod bench{
    extern crate test;
    use self::test::Bencher;
    use super::*;
    #[bench]
    fn bench_is_ascending_slice(b: &mut Bencher){
        let v = Vec::from_fn(10000, |idx| idx * 2);        
        b.iter(|| {
            assert!(is_ascending_slice(v.as_slice()));
                });
    }
    #[bench]
    fn bench_is_ascending_buffer(b: &mut Bencher){
        let v = Vec::from_fn(10000, |idx| idx * 2);        
        b.iter(|| unsafe {
            assert!(is_ascending_buffer(v.as_ptr(), v.len()));
                });
    }
    #[bench]
    fn bench_is_descending_slice(b: &mut Bencher){
        let v = Vec::from_fn(10000, |idx| -(idx as isize) * 2i);        
        b.iter(|| {
            assert!(is_descending_slice(v.as_slice()));
                });
    }
    #[bench]
    fn bench_is_descending_buffer(b: &mut Bencher){
        let v = Vec::from_fn(10000, |idx| -(idx as isize) * 2i);        
        b.iter(|| unsafe {
            assert!(is_descending_buffer(v.as_ptr(), v.len()));
                });
    }
}



