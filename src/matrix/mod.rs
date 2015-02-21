#![doc="Fundamental matrix structures
"]
pub mod constructors;
pub mod iter;
pub mod matrix;
//pub mod matrix_conversion;
//pub mod matrix_minmax;
pub mod random;
pub mod traits;
pub mod vector;
pub mod view;
pub mod view_conversion;
pub mod view_minmax;
pub mod triangular_matrix;
pub mod eo {
    //pub use self::eo_traits;
    //pub use self::eo_matrix;
    //pub use self::eo_view;
    pub mod eo_traits;
    pub mod eo_matrix;
    pub mod eo_view;
}


pub mod update{

    pub mod traits;
    pub mod matrix_updates;
}


pub mod transpose{
    pub mod traits;
    pub mod matrix_transpose;
}

pub mod extract{
    pub mod traits;
    pub mod matrix_extract;
    pub mod view_extract;
}


