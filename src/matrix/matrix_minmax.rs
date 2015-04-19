

// local imports
use matrix::matrix::Matrix;
use matrix::traits::*;
use algebra::Number;

/// Implements matrix extraction API
impl <T:Number+PartialOrd> MinMax<T> for Matrix<T> {

    /// Returns a column vector consisting of maximum over each row
    fn max_row_wise(&self) -> Matrix<T>{
        // Pick the first column
        let mut result = self.col(0);
        let pd = result.as_mut_ptr();
        let ps = self.as_ptr();
        for r  in range(0, self.num_rows()){
            let dst_offset = result.cell_to_offset(r, 0);
            for c in range (1, self.num_cols()){
                let src_offset = self.cell_to_offset(r, c);
                unsafe{
                    let s = *ps.offset(src_offset);
                    let d = *pd.offset(dst_offset);
                    *pd.offset(dst_offset) = if s > d { s } else {d};

                }
            }
        }
        result
    }

    /// Returns a column vector consisting of minimum over each row
    fn min_row_wise(&self) -> Matrix<T>{
        // Pick the first column
        let mut result = self.col(0);
        let pd = result.as_mut_ptr();
        let ps = self.as_ptr();
        for r  in range(0, self.num_rows()){
            let dst_offset = result.cell_to_offset(r, 0);
            for c in range (1, self.num_cols()){
                let src_offset = self.cell_to_offset(r, c);
                unsafe{
                    let s = *ps.offset(src_offset);
                    let d = *pd.offset(dst_offset);
                    *pd.offset(dst_offset) = if s < d { s } else {d};

                }
            }
        }
        result
    }

   /// Returns a row vector consisting of maximum over each column
    fn max_col_wise(&self) -> Matrix<T>{
        // Pick the first row
        let mut result = self.row(0);
        let pd = result.as_mut_ptr();
        let ps = self.as_ptr();
        for c in range (0, self.num_cols()){
            let dst_offset = result.cell_to_offset(0, c);
            for r  in range(1, self.num_rows()){
                let src_offset = self.cell_to_offset(r, c);
                unsafe{
                    let s = *ps.offset(src_offset);
                    let d = *pd.offset(dst_offset);
                    *pd.offset(dst_offset) = if s > d { s } else {d};

                }
            }
        }
        result
    }



   /// Returns a row vector consisting of minimum over each column
    fn min_col_wise(&self) -> Matrix<T>{
        // Pick the first row
        let mut result = self.row(0);
        let pd = result.as_mut_ptr();
        let ps = self.as_ptr();
        for c in range (0, self.num_cols()){
            let dst_offset = result.cell_to_offset(0, c);
            for r  in range(1, self.num_rows()){
                let src_offset = self.cell_to_offset(r, c);
                unsafe{
                    let s = *ps.offset(src_offset);
                    let d = *pd.offset(dst_offset);
                    *pd.offset(dst_offset) = if s < d { s } else {d};
                }
            }
        }
        result
    }


}
