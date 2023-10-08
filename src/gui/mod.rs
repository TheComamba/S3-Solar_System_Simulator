use iced::{widget::Row, Sandbox};

use crate::sim::{initial_parameters::InitialParameters, system::StellarSystem};

pub(crate) struct Gui {
    system: StellarSystem,
}

impl Sandbox for Gui {
    type Message = ();

    fn new() -> Self {
        let params = InitialParameters::default();
        let system = StellarSystem::new(params);
        Gui { system }
    }

    fn title(&self) -> String {
        String::from("S3 - Solar System Simulator")
    }

    fn update(&mut self, _message: Self::Message) {}

    fn view(&self) -> iced::Element<'_, Self::Message> {
        Row::new().into()
    }
}
