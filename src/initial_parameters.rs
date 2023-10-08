pub(crate) type Float = f32;
pub(crate) const G: Float = 1.184e-4;
pub(crate) const ROCK_DENSITY: Float = 2.803e12; //5 g/cm^3 = 5000 kg/m^3
pub(crate) const DIMENSIONALITY: usize = 2;

pub(crate) struct InitialParameters {
    pub(crate) body_count: u32,
    pub(crate) total_mass: Float,
    pub(crate) stellar_mass_fraction: Float,
    pub(crate) position_variance: Float,
    pub(crate) velocity_variance: Float,
}

impl Default for InitialParameters {
    fn default() -> Self {
        Self {
            body_count: 100,
            total_mass: 1.,
            stellar_mass_fraction: 0.5,
            position_variance: 1.,
            velocity_variance: 1.,
        }
    }
}
