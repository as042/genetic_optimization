pub mod sim;
pub mod genome;
pub mod chromosome;
pub mod gene;

pub mod prelude {
    pub use crate::{sim::*, genome::*, chromosome::*, gene::*};
}