
// std imports
use std::mem;
use std::rt::heap::{deallocate};



/***
Private helper functions follow
*/

/// Deallocates a buffer of memory
#[inline]
pub unsafe fn dealloc<T>(ptr: *mut T, len: usize) {
    if mem::size_of::<T>() != 0 {
        deallocate(ptr as *mut u8,
                   len * mem::size_of::<T>(),
                   mem::min_align_of::<T>())
    }
}



