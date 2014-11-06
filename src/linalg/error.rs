// std imports
use std::fmt;


/// Errors in algorithms related to linear systems
pub enum LinearSystemError{

    /// There is no solution to the system of equations
    NoSolution,
    /// There are infinite solutions to the system of equations.
    InfiniteSolutions
}


impl LinearSystemError{

    /// Converts enum values to string representation
    pub fn to_string(&self) -> String {
        match *self{
            NoSolution => format!("No solution"),
            InfiniteSolutions => format!("Infinite solutions"),
        }
    }
}

impl fmt::Show for LinearSystemError {
    /// Formatting of enum values
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(write!(f, "{}", self.to_string()));
        Ok(())
    }
}
