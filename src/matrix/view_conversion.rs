
// local imports
use matrix::view::MatrixView;
use matrix::traits::*;
use algebra::structure::FieldPartial;

/// Implements matrix conversion API
impl <'a, T:FieldPartial> Conversion<T> for MatrixView<'a, T> {
    /// Converts the view to vector from standard library
    fn to_std_vec(&self) -> Vec<T> {
        let mut vec: Vec<T> = Vec::with_capacity(self.num_cells());
        // We iterate over elements in matrix and push in the vector
        let ptr = self.matrix().as_ptr();
        for c in 0..self.num_cols(){
            for r in 0..self.num_rows(){
                let offset = self.cell_to_offset(r, c);
                vec.push(unsafe{*ptr.offset(offset)});
            }
        } 
        vec
    }
}


#[cfg(test)]
mod test{
    use matrix::traits::*;
    use matrix::constructors::*;

    #[test]
    fn test_view_to_scalar(){
        let m = matrix_rw_i32(3, 3, [1, 2, 3,
            4, 5, 6,
            7, 8, 9].as_slice());
        let v = m.view(2,1, 1, 1);
        assert_eq!(v.to_scalar(), 8);
    }

}


