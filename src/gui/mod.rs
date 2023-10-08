use std::vec;

use iced::{
    widget::{canvas, Row},
    Sandbox,
};

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

struct CanvasState {
    cache: canvas::Cache,
    system: StellarSystem,
}

impl CanvasState {
    fn new(system: StellarSystem) -> CanvasState {
        CanvasState {
            cache: canvas::Cache::default(),
            system,
        }
    }
}

impl<GuiMessage> canvas::Program<GuiMessage> for CanvasState {
    type State = ();

    fn draw(
        &self,
        _state: &Self::State,
        renderer: &iced::Renderer,
        _theme: &iced::theme::Theme,
        bounds: iced::Rectangle,
        _cursor: iced::mouse::Cursor,
    ) -> Vec<canvas::Geometry> {
        let visualisation = self.cache.draw(renderer, bounds.size(), |frame| {});
        vec![visualisation]
    }
}

enum GuiMessage {}
