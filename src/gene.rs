#[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd)]
pub struct Gene {
    pub value: f32,
}

impl Gene {
    /// Creates a new `Gene`.
    #[inline]
    pub fn new(value: f32) -> Self {
        Self { 
            value
        }
    }
}