"""Pillow presentation layer and FFmpeg frame stream.

Every frame is derived from absolute simulation time, so the renderer has no
integration drift and produces the same physical state at any frame rate.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import subprocess

from PIL import Image, ImageDraw, ImageFont

from moon_libration.orbit import LunarOrbit, OrbitalState, TAU

BG, PANEL, PANEL_EDGE = "#07111f", "#101e31", "#29415f"
WHITE, MUTED, DIM = "#f2f5f9", "#a7b2c2", "#617086"
EARTH_BLUE, EARTH_LAND = "#3a8dde", "#72b879"
MOON, MOON_DARK = "#b8bdc7", "#858c98"
CYAN, AMBER, GREEN, YELLOW = "#52d6e8", "#ffca58", "#63d28d", "#ffe06b"


def _font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    names = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf",
    ]
    for name in names:
        if Path(name).exists():
            return ImageFont.truetype(name, size)
    return ImageFont.load_default(size=size)


def _point(origin: tuple[float, float], length: float, angle: float) -> tuple[float, float]:
    return origin[0] + length * math.cos(angle), origin[1] - length * math.sin(angle)


def _line(draw: ImageDraw.ImageDraw, points, fill, width: int, dash: int | None = None) -> None:
    if dash is None:
        draw.line(points, fill=fill, width=width)
        return
    (x1, y1), (x2, y2) = points
    length = math.hypot(x2 - x1, y2 - y1)
    if length == 0:
        return
    ux, uy = (x2 - x1) / length, (y2 - y1) / length
    cursor = 0.0
    while cursor < length:
        end = min(cursor + dash, length)
        draw.line((x1 + ux*cursor, y1 + uy*cursor, x1 + ux*end, y1 + uy*end), fill=fill, width=width)
        cursor += dash * 1.75


def _arrow(draw: ImageDraw.ImageDraw, start, end, fill, width: int, head: float = 14) -> None:
    draw.line((start, end), fill=fill, width=width)
    angle = math.atan2(-(end[1]-start[1]), end[0]-start[0])
    draw.polygon((end, _point(end, head, angle+math.pi*.82), _point(end, head, angle-math.pi*.82)), fill=fill)


@dataclass
class Layout:
    width: int
    height: int

    @property
    def sx(self) -> float:
        return self.width / 1920

    @property
    def sy(self) -> float:
        return self.height / 1080

    def p(self, x: float, y: float) -> tuple[int, int]:
        return round(x*self.sx), round(y*self.sy)

    def n(self, value: float) -> int:
        return max(1, round(value*min(self.sx, self.sy)))


class VideoRenderer:
    def __init__(self, width: int, height: int, fps: int, start_apsis: str):
        self.layout, self.fps, self.start = Layout(width, height), fps, start_apsis
        self.orbit = LunarOrbit()
        self.font_title = _font(self.layout.n(43), True)
        self.font_subtitle = _font(self.layout.n(23))
        self.font_heading = _font(self.layout.n(24), True)
        self.font_label = _font(self.layout.n(19))
        self.font_value = _font(self.layout.n(20))
        self.font_small = _font(self.layout.n(16))
        self.font_libration = _font(self.layout.n(45), True)
        self.earth = self.layout.p(620, 575)
        self.a_px = self.layout.n(400)
        self.moon_radius_px = self.layout.n(32)
        # Preserve the real body-size ratio while exaggerating both bodies
        # relative to their displayed separation.
        self.earth_radius_px = self.moon_radius_px * (6371.0 / 1737.4)
        self.e_px = self.a_px*self.orbit.eccentricity
        # In focus-centered coordinates x=a(cos E-e), so the ellipse center is
        # ae to the left of Earth, not to the right.
        self.ellipse_center = (self.earth[0]-self.e_px, self.earth[1])
        self.b_px = self.a_px*math.sqrt(1-self.orbit.eccentricity**2)

    def position(self, s: OrbitalState) -> tuple[float, float]:
        return (self.earth[0]+self.a_px*s.x_km/self.orbit.semi_major_axis_km,
                self.earth[1]-self.a_px*s.y_km/self.orbit.semi_major_axis_km)

    def frame(self, elapsed_days: float, intro_alpha: float = 1.0) -> Image.Image:
        s = self.orbit.state(elapsed_days, self.start)
        image = Image.new("RGB", (self.layout.width, self.layout.height), BG)
        draw = ImageDraw.Draw(image)
        self._headings(draw); self._orbit(draw); self._earth(draw)
        pos = self.position(s)
        self._bearing_and_velocity(draw, pos, s); self._moon(draw, pos, s.spin_angle)
        self._body_marker(draw, pos, s); self._telemetry(draw, s); self._legend(draw)
        if intro_alpha > 0:
            self._caption(image, intro_alpha)
        return image

    def _headings(self, draw: ImageDraw.ImageDraw) -> None:
        draw.text(self.layout.p(64, 34), "WHY THE MOON APPEARS TO ROCK", font=self.font_title, fill=WHITE)
        draw.text(self.layout.p(66, 91), "Uniform rotation + nonuniform orbital motion = longitudinal libration", font=self.font_subtitle, fill=MUTED)

    def _orbit(self, draw: ImageDraw.ImageDraw) -> None:
        cx, cy = self.ellipse_center
        draw.ellipse((cx-self.a_px, cy-self.b_px, cx+self.a_px, cy+self.b_px), outline="#43526a", width=self.layout.n(2))
        _line(draw, ((cx-self.a_px, cy), (cx+self.a_px, cy)), "#26384e", self.layout.n(1), self.layout.n(11))
        for name, x in (("APOGEE", cx-self.a_px), ("PERIGEE", cx+self.a_px)):
            r = self.layout.n(4); draw.ellipse((x-r, cy-r, x+r, cy+r), fill=MUTED)
            box = draw.textbbox((0, 0), name, font=self.font_small)
            draw.text((x-(box[2]-box[0])/2, cy+self.layout.n(14)), name, font=self.font_small, fill=DIM)

    def _earth(self, draw: ImageDraw.ImageDraw) -> None:
        x, y = self.earth; r = self.earth_radius_px
        draw.ellipse((x-r, y-r, x+r, y+r), fill=EARTH_BLUE, outline="#a9d8ff", width=self.layout.n(2))
        draw.ellipse((x-r*.55, y-r*.35, x-r*.05, y-r*.03), fill=EARTH_LAND)
        draw.ellipse((x+r*.05, y-r*.05, x+r*.52, y+r*.50), fill=EARTH_LAND)
        text = "Earth (at focus)"; box = draw.textbbox((0, 0), text, font=self.font_label)
        draw.text((x-(box[2]-box[0])/2, y+r+self.layout.n(11)), text, font=self.font_label, fill=EARTH_BLUE)

    def _moon(self, draw: ImageDraw.ImageDraw, pos, angle: float) -> None:
        x, y = pos; r = self.moon_radius_px
        draw.ellipse((x-r, y-r, x+r, y+r), fill=MOON, outline=WHITE, width=self.layout.n(2))
        for radial, theta, radius in ((.38, 2.4, 6), (.48, -.45, 4), (.51, .85, 3)):
            cx, cy = _point(pos, r*radial, angle+theta); rr = self.layout.n(radius)
            draw.ellipse((cx-rr, cy-rr, cx+rr, cy+rr), fill=MOON_DARK)

    def _bearing_and_velocity(self, draw: ImageDraw.ImageDraw, pos, s: OrbitalState) -> None:
        # The instantaneous Earth-bearing is shown literally center-to-center.
        _line(draw, (pos, self.earth), CYAN, self.layout.n(3), self.layout.n(8))
        e = self.orbit.eccentricity
        vx, vy = -math.sin(s.eccentric_anomaly), math.sqrt(1-e*e)*math.cos(s.eccentric_anomaly)
        direction = math.atan2(vy, vx)
        length = self.layout.n(57+34*(s.speed_km_s-.96)/.12)
        _arrow(draw, _point(pos, self.layout.n(37), direction), _point(pos, length, direction), GREEN, self.layout.n(4), self.layout.n(12))

    def _body_marker(self, draw: ImageDraw.ImageDraw, pos, s: OrbitalState) -> None:
        marker_length = self.a_px * 0.5
        _arrow(draw, pos, _point(pos, marker_length, s.spin_angle), AMBER, self.layout.n(6), self.layout.n(15))
        if abs(s.libration_rad) > math.radians(.12):
            steps = max(4, int(abs(s.libration_rad)*30))
            pts = [_point(pos, marker_length*.82, s.spin_angle+s.libration_rad*i/steps) for i in range(steps+1)]
            draw.line(pts, fill=YELLOW, width=self.layout.n(4))

    def _telemetry(self, draw: ImageDraw.ImageDraw, s: OrbitalState) -> None:
        x0, y0 = self.layout.p(1374, 165); x1, y1 = self.layout.p(1870, 935)
        draw.rounded_rectangle((x0, y0, x1, y1), radius=self.layout.n(18), fill=PANEL, outline=PANEL_EDGE, width=self.layout.n(2))
        draw.text(self.layout.p(1420, 202), "ORBITAL STATE", font=self.font_heading, fill=CYAN)
        draw.line((*self.layout.p(1420, 245), *self.layout.p(1823, 245)), fill=PANEL_EDGE, width=self.layout.n(1))
        rows = [
            ("Elapsed time", f"{s.elapsed_days:,.2f} d"),
            ("Orbit complete", f"{100*s.orbit_fraction:,.1f}%"),
            ("Earth–Moon distance", f"{s.distance_km:,.0f} km"),
            ("Orbital speed", f"{s.speed_km_s:.3f} km/s"),
            ("Uniform rotation", f"{math.degrees(TAU*s.orbit_fraction)%360:.1f}°"),
            ("True anomaly", f"{math.degrees(s.true_anomaly)%360:.1f}°"),
        ]
        for i, (label, value) in enumerate(rows):
            y = 285+i*65; draw.text(self.layout.p(1420, y), label, font=self.font_label, fill=MUTED)
            box = draw.textbbox((0, 0), value, font=self.font_value)
            draw.text((self.layout.p(1823, y)[0]-(box[2]-box[0]), self.layout.p(0, y)[1]), value, font=self.font_value, fill=WHITE)
        draw.text(self.layout.p(1420, 700), "LONGITUDINAL LIBRATION", font=self.font_heading, fill=YELLOW)
        value = f"{s.libration_deg:+.2f}°"; box = draw.textbbox((0, 0), value, font=self.font_libration)
        draw.text((self.layout.p(1622, 742)[0]-(box[2]-box[0])/2, self.layout.p(0, 742)[1]), value, font=self.font_libration, fill=YELLOW)
        draw.text(self.layout.p(1535, 808), "+ east   ·   − west", font=self.font_small, fill=MUTED)

    def _legend(self, draw: ImageDraw.ImageDraw) -> None:
        for i, (color, label, dashed) in enumerate(((AMBER, "Moon-fixed rotation angle", False), (CYAN, "Center-to-center direction", True), (GREEN, "Orbital velocity", False))):
            y = 862+i*28
            _line(draw, (self.layout.p(1422, y), self.layout.p(1470, y)), color, self.layout.n(4), self.layout.n(7) if dashed else None)
            draw.text(self.layout.p(1485, y-11), label, font=self.font_small, fill=MUTED)

    def _caption(self, image: Image.Image, alpha: float) -> None:
        overlay = Image.new("RGBA", image.size, (0, 0, 0, 0)); draw = ImageDraw.Draw(overlay)
        text = f"Starting at {self.start} · the lunar marker initially points at Earth"
        box = draw.textbbox((0, 0), text, font=self.font_subtitle)
        x, y = self.layout.width*.36-(box[2]-box[0])/2, self.layout.height-self.layout.n(55)
        draw.rounded_rectangle((x-self.layout.n(22), y-self.layout.n(9), x+(box[2]-box[0])+self.layout.n(22), y+(box[3]-box[1])+self.layout.n(10)), radius=self.layout.n(12), fill=(8, 18, 31, round(220*alpha)))
        draw.text((x, y), text, font=self.font_subtitle, fill=(190, 201, 216, round(255*alpha)))
        image.paste(overlay, (0, 0), overlay)

    def render(self, output: Path, duration: float, hold: float, ffmpeg: str = "ffmpeg") -> None:
        total_frames = round((hold+duration+hold)*self.fps)
        command = [ffmpeg, "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24",
                   "-s", f"{self.layout.width}x{self.layout.height}", "-r", str(self.fps), "-i", "-",
                   "-an", "-c:v", "libx264", "-preset", "medium", "-crf", "17",
                   "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(output)]
        process = subprocess.Popen(command, stdin=subprocess.PIPE)
        assert process.stdin is not None
        try:
            for index in range(total_frames):
                wall = index/self.fps
                if wall < hold:
                    elapsed, alpha = 0.0, 1.0
                elif wall < hold+duration:
                    elapsed = self.orbit.period_days*(wall-hold)/duration
                    alpha = max(0.0, 1.0-(wall-hold)/.5)
                else:
                    elapsed, alpha = self.orbit.period_days, 0.0
                process.stdin.write(self.frame(elapsed, alpha).tobytes())
        finally:
            process.stdin.close()
        if process.wait():
            raise RuntimeError("FFmpeg encoding failed")
