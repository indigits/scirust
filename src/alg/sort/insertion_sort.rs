#![doc="Insertion sort algorithms
"]

// std imports
use std::ptr;

/// Performs insertion sort on a slice of type T
pub fn insertion_sort_slice<T : PartialOrd>(data : &mut [T]){
    let n  = data.len();
    for j in range(1, n){
        // we insert data[j] into the sorted sequence 
        //data[0...j-1]
        let mut i = j -1;
        while data[i] > data[i+1]{
            data.swap(i + 1, i);
            if i == 0{
                break;
            }
            i -= 1;
        }
    }
}

/// Performs insertion sort on a buffer of data
pub unsafe fn insertion_sort_buffer<T : PartialOrd>(data : *mut T, n : uint){
    for j in range(1, n){
        // we insert data[j] into the sorted sequence 
        //data[0...j-1]
        let mut i = (j -1) as int;
        let mut di = data.offset(i);
        let mut dj = data.offset(i+1);
        while *di > *dj {
            ptr::swap(di, dj);
            if i == 0{
                break;
            }
            i -= 1;
            di = di.offset(-1);
            dj = dj.offset(-1);
        }
    }
}

/******************************************************
 *
 *   Unit tests
 *
 *******************************************************/


#[cfg(test)]
mod test{
    extern crate quickcheck_macros;
    use super::*;
    use alg::sort::*; 

    #[test]
    fn test_insertion_sort_slice_1() {
        let mut x : [int, ..5] = [1, 2, 3, 4, 5];
        insertion_sort_slice(&mut x);
        assert!(is_ascending_slice(&x));
    }

    #[test]
    fn test_insertion_sort_slice_2() {
        let mut x : [int, ..5] = [5,4,3,2,1];
        insertion_sort_slice(&mut x);
        assert!(is_ascending_slice(&x));
    }

    #[test]
    fn test_insertion_sort_buffer_1() {
        let mut x : [int, ..5] = [5,5,3,3,1];
        unsafe {
            insertion_sort_buffer(&mut x[0], x.len());
        }
        assert!(is_ascending_slice(&x));
    }

    //#[quickcheck]
    //fn test_insertion_sort_slice(mut xs: Vec<int>) -> bool {
    //    insertion_sort_slice(xs.as_mut_slice());
    //    is_ascending_slice(xs.as_slice())
    //}

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
    use std::rand;
    use std::rand::Rng;

    #[bench]
    fn bench_insertion_sort_slice_reverse_data(b: &mut Bencher){
        let mut v = Vec::from_fn(10000, |idx| (20000 - idx));        
        b.iter(|| {
            insertion_sort_slice(v.as_mut_slice());
        });
    }

    #[bench]
    fn bench_insertion_sort_slice_random_data(b: &mut Bencher){
        // create a task-local Random Number Generator
        let mut rng = rand::task_rng();
        let mut v: Vec<uint> = rng.gen_iter::<uint>().take(10000).collect();
        b.iter(|| {
            insertion_sort_slice(v.as_mut_slice());
        });
    }

    #[bench]
    fn bench_insertion_sort_buffer_reverse_data(b: &mut Bencher){
        let mut v = Vec::from_fn(10000, |idx| (20000 - idx));        
        b.iter(|| {
            unsafe {
                insertion_sort_buffer(v.as_mut_ptr(), v.len());
            }
        });
    }
}


