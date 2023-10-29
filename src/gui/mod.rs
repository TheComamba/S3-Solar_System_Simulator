use std::vec;

use iced::{
    widget::{
        canvas::{self, Path},
        Button, Column, Text,
    },
    Color, Length, Sandbox, Size,
};

use crate::sim::{
    body::Body, initial_parameters::InitialParameters, system::StellarSystem, units::Float,
};

pub(crate) struct Gui {
    canvas_state: CanvasState,
}

impl Sandbox for Gui {
    type Message = GuiMessage;

    fn new() -> Self {
        Gui {
            canvas_state: CanvasState::new(),
        }
    }

    fn title(&self) -> String {
        String::from("S3 - Solar System Simulator")
    }

    fn update(&mut self, message: Self::Message) {
        match message {
            GuiMessage::Evolve => {
                self.canvas_state.system.evolve_for(1e1);
                self.canvas_state.bodies_cache.clear();
                // println!("\n{:?}\n", self.canvas_state.system);
            }
        };
    }

    fn view(&self) -> iced::Element<'_, Self::Message> {
        let evolve_button = Button::new(Text::new("Evolve")).on_press(GuiMessage::Evolve);
        let canvas = iced::widget::canvas(&self.canvas_state)
            .width(Length::Fill)
            .height(Length::Fill);
        Column::new()
            .push(evolve_button)
            .push(canvas)
            .align_items(iced::Alignment::Center)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }
}

struct CanvasState {
    background_cache: canvas::Cache,
    bodies_cache: canvas::Cache,
    system: StellarSystem,
    body_size_factor: Float,
    system_size_factor: Float,
}

impl CanvasState {
    fn new() -> CanvasState {
        let params = InitialParameters::default();
        let system = StellarSystem::new(params);
        CanvasState {
            background_cache: canvas::Cache::default(),
            bodies_cache: canvas::Cache::default(),
            system,
            body_size_factor: 1e5,
            system_size_factor: 3e2,
        }
    }

    fn body_display_radius(&self, body: &Body) -> f32 {
        (body.radius() * self.body_size_factor).log(2.)
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
                    let radius = self.body_display_radius(body);
                    let pos = frame.center()
                        + iced::Vector::new(
                            body.position[0] * self.system_size_factor,
                            body.position[1] * self.system_size_factor,
                        );
                    path_builder.circle(pos, radius);
                }
            });
            frame.fill(&bodies, Color::BLACK);
        });
        vec![bodies]
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) enum GuiMessage {
    Evolve,
}
