#![doc="Quick sort algorithms
"]

// std imports
use std::ptr;
use std::fmt;


/// Quick sort algorithm on a slice of data
pub fn quick_sort_slice<T: PartialOrd>(xs: &mut [T]){
    if xs.len() <= 1 {
        return;
    }
    // Compute the pivot as median between left-right and center.
    let pivot = {
        let (l, r) = (0, xs.len() - 1);
        // identify the index of middle element
        let m = l + ((r - l) / 2);
        let (left, middle, right) = (&xs[l], &xs[m], &xs[r]);
        if middle >= left && middle <= right {
            m
        } else if left >= middle && left <= right {
            l
        } else {
            r
        }
    };
    // Partition the array in two parts around the pivot
    let lasti = xs.len() - 1;
    // Send the pivot to the end.
    xs.swap(lasti, pivot);
    // now start from the beginning
    let (mut i, mut pivot) = (0, 0);
    while i < lasti {
        if xs[i] <= xs[lasti] {
            xs.swap(i, pivot);
            pivot = pivot + 1;
        }
        i = i + 1;
    }
    xs.swap(pivot, lasti);
    // Now run quick sort recursively
    // we split original slice into two slices
    // the left slide [0, pivot]
    quick_sort_slice(xs.slice_to_mut(pivot));
    // the right slice [pivot + 1, end]
    quick_sort_slice(xs.slice_from_mut(pivot+1));
}


// Helper function for debugging
pub unsafe fn  print_arr<T: PartialOrd + fmt::Show>(data : *mut T, n : uint){
    for i in range(0, n){
        print!("{} ", *data.offset(i as int));
    }
    println!("");
}

/// Quick sort algorithm on a buffer of data
pub unsafe fn quick_sort_buffer<T: PartialOrd>(data : *mut T, n : uint){
    if n <= 1 {
        return;
    }
    let right = (n - 1) as int;
    // pointer to the last entry
    let pright = data.offset(right);
    let mut pleft = data;
    // Compute the pivot as median between left-right and center.
    let pivot = {
        // identify the index of middle element
        let mid = right / 2;
        let pmid = data.offset(mid);
        if *pmid >= *pleft && *pmid <= *pright {
            mid
        } else if *pleft >= *pmid && *pleft <= *pright {
            0i
        } else {
            right
        }
    };
    // Partition the array in two parts around the pivot
    // Send the pivot to the end.
    let mut ppivot = data.offset(pivot);
    //println!("swapping pivot: {}, right: {}", *ppivot, *pright);
    ptr::swap(pright, ppivot);
    //print_arr(data, n);
    // now start from the beginning
    let (mut i, mut pivot) = (0i, 0i);
    ppivot = data.offset(pivot);
    while i < right {
        if *pleft <= *pright {
            //println!("swapping left: {}, pivot: {}", *pleft, *ppivot);
            ptr::swap(pleft, ppivot);
            //print_arr(data, n);
            pivot = pivot + 1;
            ppivot = ppivot.offset(1);
        }
        i = i + 1;
        pleft = pleft.offset(1);
    }
    //println!("swapping right: {}, pivot: {}", *pright, *ppivot);
    ptr::swap(ppivot, pright);
    //print_arr(data, n);
    // Now run quick sort recursively
    // we split original slice into two slices
    // the left slide [0, pivot]
    quick_sort_buffer(data, pivot as uint);
    // the right slice [pivot + 1, end]
    quick_sort_buffer(data.offset(pivot + 1i), n - (pivot as uint) - 1);
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
    use std::rand;
    use std::rand::Rng;

    #[test]
    fn test_quick_sort_slice_1() {
        let mut x : [int, ..5] = [1, 2, 3, 4, 5];
        quick_sort_slice(&mut x);
        assert!(is_ascending_slice(&x));
    }

    #[test]
    fn test_quick_sort_slice_2() {
        let mut x : [int, ..5] = [5,4,3,2,1];
        quick_sort_slice(&mut x);
        assert!(is_ascending_slice(&x));
    }

    #[test]
    fn test_quick_sort_buffer_1() {
        let mut x : [int, ..5] = [1, 2, 3, 4, 5];
        println!("{}", x);
        unsafe {quick_sort_buffer(&mut x[0], x.len());}
        println!("{}", x);
        assert!(is_ascending_slice(&x));
    }

    #[test]
    fn test_quick_sort_buffer_2() {
        let mut x : [int, ..5] = [5,4,3,2,1];
        unsafe { quick_sort_buffer(&mut x[0], x.len());}
        assert!(is_ascending_slice(&x));
    }

    #[test]
    fn test_quick_sort_buffer_reverse_large() {
        let mut v = Vec::from_fn(10000, |idx| (20000 - idx));        
        unsafe { quick_sort_buffer(v.as_mut_ptr(), v.len());}
        assert!(is_ascending_slice(v.as_slice()));
    }

    #[test]
    fn test_quick_sort_buffer_random_data() {
        let mut rng = rand::task_rng();
        let mut v: Vec<uint> = rng.gen_iter::<uint>().take(10000).collect();
        unsafe { quick_sort_buffer(v.as_mut_ptr(), v.len());}
        assert!(is_ascending_slice(v.as_slice()));
    }

    // #[test]
    // fn test_insertion_sort_buffer_1() {
    //     let mut x : [int, ..5] = [5,5,3,3,1];
    //     unsafe {
    //         quick_sort_buffer(&mut x[0], x.len());
    //     }
    //     assert!(is_ascending_slice(&x));
    // }

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
    fn bench_quick_sort_slice_reverse_data(b: &mut Bencher){
        let mut v = Vec::from_fn(10000, |idx| (20000 - idx));        
        b.iter(|| {
            quick_sort_slice(v.as_mut_slice());
        });
    }

    #[bench]
    fn bench_quick_sort_slice_random_data(b: &mut Bencher){
        // create a task-local Random Number Generator
        let mut rng = rand::task_rng();
        let mut v: Vec<uint> = rng.gen_iter::<uint>().take(10000).collect();
        b.iter(|| {
            quick_sort_slice(v.as_mut_slice());
        });
    }


    #[bench]
    fn bench_quick_sort_buffer_reverse_data(b: &mut Bencher){
        let mut v = Vec::from_fn(10000, |idx| (20000 - idx));        
        b.iter(|| {
            unsafe {
                quick_sort_buffer(v.as_mut_ptr(), v.len());
            }
        });
    }

    #[bench]
    fn bench_quick_sort_buffer_random_data(b: &mut Bencher){
        // create a task-local Random Number Generator
        let mut rng = rand::task_rng();
        let mut v: Vec<uint> = rng.gen_iter::<uint>().take(10000).collect();
        b.iter(|| {
            unsafe {
                quick_sort_buffer(v.as_mut_ptr(), v.len());
            }
        });
    }


}


