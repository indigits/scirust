// local imports
use matrix::matrix::Matrix;
use matrix::traits::*;
use algebra::structure::FieldPartial;

/// Implements matrix conversion API
impl <T:FieldPartial> Conversion<T> for Matrix<T> {
    /// Converts the matrix to vector from standard library
    fn to_std_vec(&self) -> Vec<T> {
        let mut vec: Vec<T> = Vec::with_capacity(self.num_cells());
        // We iterate over elements in matrix and push in the vector
        let ptr = self.as_ptr();
        for c in 0..self.num_cols(){
            for r in 0..self.num_rows(){
                let offset = self.cell_to_offset(r, c);
                vec.push(unsafe{*ptr.offset(offset)});
            }
        } 
        vec
    }


}
