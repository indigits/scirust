/// This module provides useful discrete mathematics functions

// std imports

// external imports

// local imports


#[inline]
pub fn mod_n (x : isize, n : isize) -> usize {
    let x = x % n;
    if x < 0{
        (x + n.abs()) as usize
    }else{
        x as usize
    }
}

/******************************************************
 *
 *   Unit tests follow.
 *
 *******************************************************/


#[cfg(test)]
mod test{
    #[test]
    fn test_modulo_n(){
        assert_eq!(super::mod_n(4, 5), 4);
        assert_eq!(super::mod_n(6, 5), 1);
        assert_eq!(super::mod_n(-3, 5), 2);
        assert_eq!(super::mod_n(-5, -5), 0);
        assert_eq!(super::mod_n(-6, -5), 4);
    }
}
