use crate::{
    body::Body,
    initial_parameters::{Float, InitialParameters, DIMENSIONALITY},
};

struct System {
    central_body: Body,
    bodies: Vec<Body>,
}

impl System {
    fn new(params: InitialParameters) -> System {
        let central_body = Body::new(
            vec![0.; DIMENSIONALITY],
            vec![0.; DIMENSIONALITY],
            params.total_mass * params.stellar_mass_fraction,
        );
        let bodies = (0..params.body_count)
            .map(|_| {
                let position = vec![0.; DIMENSIONALITY];
                let velocity = vec![0.; DIMENSIONALITY];
                let mass = params.total_mass * (1. - params.stellar_mass_fraction)
                    / params.body_count as Float;
                Body::new(position, velocity, mass)
            })
            .collect();
        System {
            central_body,
            bodies,
        }
    }
}
