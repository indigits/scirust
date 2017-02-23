
// std imports
use std::mem;
use alloc::heap::{deallocate};



/***
Private helper functions follow
*/

// /// Deallocates a buffer of memory
// #[inline]
// pub unsafe fn dealloc<T>(ptr: *mut T, len: usize) {
//     if mem::size_of::<T>() != 0 {
//         deallocate(ptr as *mut u8,
//                    len * mem::size_of::<T>(),
//                    mem::align_of::<T>())
//     }
// }



