
// local imports
use matrix::view::MatrixView;
use matrix::traits::*;
use matrix::element::MatrixElt;

/// Implements matrix conversion API
impl <'a, T:MatrixElt> Conversion<T> for MatrixView<'a, T> {
    /// Converts the view to vector from standard library
    fn to_std_vec(&self) -> Vec<T> {
        let mut vec: Vec<T> = Vec::with_capacity(self.num_cells());
        // We iterate over elements in matrix and push in the vector
        let ptr = self.matrix().as_ptr();
        for c in range(0, self.num_cols()){
            for r in range (0, self.num_rows()){
                let offset = self.cell_to_offset(r, c);
                vec.push(unsafe{*ptr.offset(offset)});
            }
        } 
        vec
    }
}


