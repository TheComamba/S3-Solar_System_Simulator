use super::units::Float;

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
            body_count: 2,
            total_mass: 333_000.,
            stellar_mass_fraction: 0.9,
            position_variance: 3.,
            velocity_variance: 0.,
        }
    }
}
