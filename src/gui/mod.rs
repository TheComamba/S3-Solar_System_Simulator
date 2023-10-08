use std::vec;

use iced::{
    widget::{
        canvas::{self, Path},
        Button, Row, Text,
    },
    Color, Length, Sandbox, Size,
};

use crate::sim::{initial_parameters::InitialParameters, system::StellarSystem};

pub(crate) struct Gui {
    canvas_state: CanvasState,
}

impl Sandbox for Gui {
    type Message = ();

    fn new() -> Self {
        Gui {
            canvas_state: CanvasState::new(),
        }
    }

    fn title(&self) -> String {
        String::from("S3 - Solar System Simulator")
    }

    fn update(&mut self, _message: Self::Message) {}

    fn view(&self) -> iced::Element<'_, Self::Message> {
        let evolve_button = Button::new(Text::new("Evolve"));
        let canvas = iced::widget::canvas(&self.canvas_state)
            .width(Length::Fill)
            .height(Length::Fill);
        Row::new().push(evolve_button).push(canvas).into()
    }
}

struct CanvasState {
    background_cache: canvas::Cache,
    bodies_cache: canvas::Cache,
    system: StellarSystem,
}

impl CanvasState {
    fn new() -> CanvasState {
        let params = InitialParameters::default();
        let system = StellarSystem::new(params);
        CanvasState {
            background_cache: canvas::Cache::default(),
            bodies_cache: canvas::Cache::default(),
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
        let background = self
            .background_cache
            .draw(renderer, bounds.size(), |frame| {
                frame.fill_rectangle(frame.center(), Size::UNIT, Color::BLACK)
            });
        let bodies = self.bodies_cache.draw(renderer, bounds.size(), |frame| {
            let bodies = Path::new(|path_builder| {
                for body in &self.system.bodies {
                    let pos =
                        frame.center() + iced::Vector::new(body.position[0], body.position[1]);
                    path_builder.circle(pos, body.radius());
                }
            });
            frame.fill(&bodies, Color::BLACK);
        });
        vec![background, bodies]
    }
}

enum GuiMessage {}
