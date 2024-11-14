use egui::{Slider, Ui};
use macroquad::prelude::*;

pub struct Grid {
    spacing: f32,
    color: Color,
}

impl Grid {
    pub fn new() -> Grid {
        let mut color = ORANGE;
        color.a = 0.5;
        Grid {
            spacing: 100.0,
            color,
        }
    }

    pub fn draw(&self) {
        let spacing = self.spacing;
        let axes_color = self.color;
        let other_color = self.color;
        let sc_width = screen_width();
        let sc_height = screen_height();
        let slices = ((sc_width / spacing) as u32).max((sc_height / spacing) as u32);
        for i in 0..slices + 1 {
            let color = if i == 0 { axes_color } else { other_color };

            draw_line_3d(
                vec3(i as f32 * spacing, 0., 0.),
                vec3(i as f32 * spacing, slices as f32 * spacing, 0.),
                color,
            );
            draw_line_3d(
                vec3(0., i as f32 * spacing, 0.),
                vec3(slices as f32 * spacing, i as f32 * spacing, 0.),
                color,
            );
        }
    }

    pub fn mouse_position(&self) -> (f32, f32) {
        let mut res = mouse_position();
        if is_key_down(KeyCode::LeftShift) || is_key_down(KeyCode::RightShift) {
            res.0 += 0.5 * self.spacing;
            res.1 += 0.5 * self.spacing;
            res.0 -= res.0 % self.spacing;
            res.1 -= res.1 % self.spacing;
        }
        res
    }

    pub fn ui(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            ui.label("spacing");
            ui.add(Slider::new(&mut self.spacing, 5.0..=500.0));
        });
        let mut color = [self.color.r, self.color.g, self.color.b, self.color.a];
        ui.color_edit_button_rgba_unmultiplied(&mut color);
        self.color.r = color[0];
        self.color.g = color[1];
        self.color.b = color[2];
        self.color.a = color[3];
    }
}

impl Default for Grid {
    fn default() -> Self {
        Self::new()
    }
}
