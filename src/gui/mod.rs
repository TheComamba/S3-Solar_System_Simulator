use std::vec;

use iced::{
    widget::{
        canvas::{self, Path},
        Button, Column, Row, Text,
    },
    Color, Length, Sandbox, Size,
};

use crate::sim::{
    body::Body,
    initial_parameters::InitialParameters,
    system::StellarSystem,
    units::{Float, TIME_TO_SECONDS},
};

pub(crate) struct Gui {
    canvas_state: CanvasState,
    time_step: Float,
}

impl Sandbox for Gui {
    type Message = GuiMessage;

    fn new() -> Self {
        Gui {
            canvas_state: CanvasState::new(),
            time_step: 1e0,
        }
    }

    fn title(&self) -> String {
        String::from("S3 - Solar System Simulator")
    }

    fn update(&mut self, message: Self::Message) {
        match message {
            GuiMessage::Evolve => {
                self.canvas_state.system.evolve_for(self.time_step);
                // println!("\n{:?}\n", self.canvas_state.system);
            }
            GuiMessage::ChangeTimeStep(step) => self.time_step = step,
            GuiMessage::ChangeBodySize(fac) => self.canvas_state.body_size_factor = fac,
            GuiMessage::ChangeSystemSize(fac) => self.canvas_state.system_size_factor = fac,
        };
        self.canvas_state.bodies_cache.clear();
    }

    fn view(&self) -> iced::Element<'_, Self::Message> {
        let header = Row::new()
            .push(self.status_block())
            .push(self.control_block());
        let canvas = iced::widget::canvas(&self.canvas_state)
            .width(Length::Fill)
            .height(Length::Fill);
        Column::new()
            .push(header)
            .push(canvas)
            .align_items(iced::Alignment::Center)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }
}

impl Gui {
    fn time_format(time_in_seconds: Float) -> String {
        if time_in_seconds < 1e2 {
            format!("{:.2} s", time_in_seconds)
        } else if time_in_seconds < 1e4 {
            format!("{:.2} min", time_in_seconds / 60.)
        } else if time_in_seconds < 1e5 {
            format!("{:.2} h", time_in_seconds / 3600.)
        } else if time_in_seconds < 1e8 {
            format!("{:.2} d", time_in_seconds / 86400.)
        } else if time_in_seconds < 1e11 {
            format!("{:.2} y", time_in_seconds / 31536000.)
        } else {
            format!("{:.2} ky", time_in_seconds / 31536000. / 1e3)
        }
    }

    fn status_block(&self) -> iced::Element<'_, GuiMessage> {
        let time_text = Self::time_format(self.canvas_state.system.current_time * TIME_TO_SECONDS);
        let time_text = Text::new(format!("Evolution time: {}", time_text));
        Column::new().push(time_text).into()
    }

    fn value_control_block(
        &self,
        value_name: &'static str,
        value: String,
        decrease: GuiMessage,
        increase: GuiMessage,
    ) -> iced::Element<'_, GuiMessage> {
        let decrease_time_step_button = Button::new(Text::new("<<")).on_press(decrease);
        let increase_time_step_button = Button::new(Text::new(">>")).on_press(increase);
        let time_step_block = Row::new()
            .push(Text::new(value_name))
            .push(decrease_time_step_button)
            .push(Text::new(value))
            .push(increase_time_step_button);
        time_step_block.into()
    }

    fn control_block(&self) -> iced::Element<'_, GuiMessage> {
        let time_step_block = self.value_control_block(
            "Time step: ",
            Self::time_format(self.time_step * TIME_TO_SECONDS),
            GuiMessage::ChangeTimeStep(0.5 * self.time_step),
            GuiMessage::ChangeTimeStep(2. * self.time_step),
        );
        let system_size_block = self.value_control_block(
            "System size factor: ",
            format!("{}", self.canvas_state.system_size_factor),
            GuiMessage::ChangeSystemSize(0.5 * self.canvas_state.system_size_factor),
            GuiMessage::ChangeSystemSize(2. * self.canvas_state.system_size_factor),
        );
        let body_size_block = self.value_control_block(
            "Body size factor: ",
            format!("{}", self.canvas_state.body_size_factor),
            GuiMessage::ChangeBodySize(0.5 * self.canvas_state.body_size_factor),
            GuiMessage::ChangeBodySize(2. * self.canvas_state.body_size_factor),
        );
        let evolve_button = Button::new(Text::new("Evolve")).on_press(GuiMessage::Evolve);
        Column::new()
            .push(time_step_block)
            .push(system_size_block)
            .push(body_size_block)
            .push(evolve_button)
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
    ChangeTimeStep(Float),
    ChangeBodySize(Float),
    ChangeSystemSize(Float),
}
