#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Earth-centred comparison of an orbital Moon and a uniform-rate virtual Moon.

The presentation follows this repository's author-2d-animation renderer.
Source coordinates and analyst-supplied dynamics are documented in dremt-derived-model/docs/ANIMATION.md
and orbit_comparison.py. Latitude stays in the data layer for a later display.
The Earth-centred separation is a positional angle, not a body-attitude angle.
The left pendulum panel rotates and translates the same samples at the same scale.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
import shutil
import subprocess

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from orbit_comparison import EphemerisComparisonTrack, build_perigee_comparison_track
from pendulum_view import pendulum_view
from project_paths import ANIMATION_OUTPUT_DIR, project_path

BG, PANEL, PANEL_EDGE = "#07111f", "#101e31", "#29415f"
WHITE, MUTED, DIM = "#f2f5f9", "#a7b2c2", "#617086"
EARTH_BLUE, EARTH_LAND = "#3a8dde", "#72b879"
MOON, MOON_DARK = "#b8bdc7", "#858c98"
CYAN, GREEN, YELLOW = "#52d6e8", "#63d28d", "#ffe06b"
MAGENTA = "#ff63df"


def font(size, bold=False):
    suffix = "-Bold" if bold else ""
    for root, name in (("/usr/share/fonts/truetype/dejavu", f"DejaVuSans{suffix}.ttf"),
                       ("/usr/share/fonts/truetype/liberation2", f"LiberationSans{suffix or '-Regular'}.ttf")):
        path = Path(root) / name
        if path.exists():
            return ImageFont.truetype(str(path), size)
    return ImageFont.load_default(size=size)


def point(origin, radius, angle):
    return origin[0] + radius * math.cos(angle), origin[1] - radius * math.sin(angle)


def line(draw, start, end, color, width, dash=None):
    if dash is None:
        draw.line((start, end), fill=color, width=width)
        return
    distance = math.dist(start, end)
    if distance == 0:
        return
    dx, dy = (end[0] - start[0]) / distance, (end[1] - start[1]) / distance
    cursor = 0.0
    while cursor < distance:
        stop = min(cursor + dash, distance)
        draw.line((start[0] + dx * cursor, start[1] + dy * cursor,
                   start[0] + dx * stop, start[1] + dy * stop), fill=color, width=width)
        cursor += dash * 1.75


def arrow(draw, start, end, color, width, head):
    line(draw, start, end, color, width)
    angle = math.atan2(start[1] - end[1], end[0] - start[0])
    draw.polygon((end, point(end, head, angle + math.pi * .82),
                  point(end, head, angle - math.pi * .82)), fill=color)


def dotted_line(draw, start, end, color, radius, spacing):
    distance = math.dist(start, end)
    count = max(1, math.ceil(distance / spacing))
    for index in range(count + 1):
        fraction = index / count
        x = start[0] + fraction * (end[0] - start[0])
        y = start[1] + fraction * (end[1] - start[1])
        draw.ellipse((x-radius, y-radius, x+radius, y+radius), fill=color)


def dotted_circle(draw, centre, radius, color, dot_radius, spacing):
    count = max(12, math.ceil(2 * math.pi * radius / spacing))
    for index in range(count):
        x, y = point(centre, radius, 2 * math.pi * index / count)
        draw.ellipse((x-dot_radius, y-dot_radius, x+dot_radius, y+dot_radius), fill=color)


@dataclass(frozen=True)
class Layout:
    width: int
    height: int

    def p(self, x, y):
        return round(x * self.width / 1920), round(y * self.height / 1080)

    def n(self, value):
        return max(1, round(value * min(self.width / 1920, self.height / 1080)))


@dataclass(frozen=True)
class Presentation:
    subtitle: str = "Ephemeris Moon compared with a virtual Moon advancing at a constant mean angular rate"
    interval: str | None = None
    moon_legend: str = "Real Moon · JPL Horizons"
    footer: str = "JPL ephemeris · 0° is the starting perigee direction · Body sizes enlarged"


class Renderer:
    def __init__(self, track, width=1920, height=1080, presentation=Presentation()):
        self.track, self.layout = track, Layout(width, height)
        self.presentation = presentation
        L = self.layout
        self.title = font(L.n(43), True)
        self.subtitle = font(L.n(22))
        self.heading = font(L.n(23), True)
        self.label = font(L.n(18))
        self.value = font(L.n(20))
        self.small = font(L.n(16))
        self.big = font(L.n(46), True)
        self.earth = L.p(930, 555)
        self.pendulum_virtual = L.p(262, 794)
        self.moon_radius = L.n(32)
        # Body radii retain their physical ratio, with a separate distance scale.
        self.earth_radius = self.moon_radius * (6371.0 / 1737.4)
        self.earth_icon = self.make_earth_icon()
        radial_extent = np.max(np.linalg.norm(track.projected_orbit_km, axis=1))
        self.km_to_px = L.n(393) / radial_extent
        self.path = [self.position(p) for p in track.projected_orbit_km]
        self.static = self.background()

    def position(self, projected):
        return (self.earth[0] + self.km_to_px * projected[0],
                self.earth[1] - self.km_to_px * projected[1])

    def pendulum_position(self, projected):
        # Exactly the same km/pixel scale as the main orbital view.
        return (self.pendulum_virtual[0] + self.km_to_px*projected[0],
                self.pendulum_virtual[1] - self.km_to_px*projected[1])

    def make_earth_icon(self):
        r = self.earth_radius
        size = math.ceil(2*r) + 8
        icon = Image.new("RGBA", (size, size), (0, 0, 0, 0))
        draw = ImageDraw.Draw(icon)
        x = y = size/2
        draw.ellipse((x-r, y-r, x+r, y+r), fill=EARTH_BLUE, outline="#a9d8ff", width=self.layout.n(2))
        draw.ellipse((x-r*.55, y-r*.35, x-r*.05, y-r*.03), fill=EARTH_LAND)
        draw.ellipse((x+r*.05, y-r*.05, x+r*.52, y+r*.50), fill=EARTH_LAND)
        return icon

    def draw_earth(self, image, centre, angle=0.):
        icon = self.earth_icon
        if angle:
            icon = icon.rotate(math.degrees(angle) % 360, resample=Image.Resampling.BICUBIC)
        image.paste(icon, (round(centre[0]-icon.width/2), round(centre[1]-icon.height/2)), icon)

    def draw_real_moon(self, draw, pos, marker_angle):
        x, y = pos
        r = self.moon_radius
        draw.ellipse((x-r, y-r, x+r, y+r), fill=MOON, outline=WHITE, width=self.layout.n(2))
        for radial, theta, radius in ((.38, 2.4, 6), (.48, -.45, 4), (.51, .85, 3)):
            cx, cy = point(pos, r*radial, marker_angle + theta)
            rr = self.layout.n(radius)
            draw.ellipse((cx-rr, cy-rr, cx+rr, cy+rr), fill=MOON_DARK)

    def right_text(self, draw, x, y, text, style, color=WHITE):
        box = draw.textbbox((0, 0), text, font=style)
        px, py = self.layout.p(x, y)
        draw.text((px - (box[2] - box[0]), py), text, font=style, fill=color)

    def centered_text(self, draw, x, y, text, style, color=WHITE):
        box = draw.textbbox((0, 0), text, font=style)
        px, py = self.layout.p(x, y)
        draw.text((px - (box[2] - box[0]) / 2, py), text, font=style, fill=color)

    def background(self):
        L = self.layout
        image = Image.new("RGB", (L.width, L.height), BG)
        draw = ImageDraw.Draw(image)
        draw.text(L.p(64, 34), "PENDULAR LIBRATION MODEL", font=self.title, fill=WHITE)
        draw.text(L.p(66, 91), self.presentation.subtitle,
                  font=self.subtitle, fill=MUTED)
        draw.line(self.path, fill="#43526a", width=L.n(2))
        # Mark actual radial extrema, without placing apogee at an imposed 180°.
        for event in self.track.apses:
            # A closed ideal ellipse has one shared starting/ending perigee.
            if event["elapsed_days"] == self.track.duration_days:
                continue
            sample = self.track.sample(event["elapsed_days"])
            apsis_position = self.position(sample.projected_position_km)
            ax, ay = apsis_position
            rr = L.n(4)
            draw.ellipse((ax-rr, ay-rr, ax+rr, ay+rr), fill=DIM)
            if event["elapsed_days"] == 0:
                closed = self.track.apses[-1]["elapsed_days"] == self.track.duration_days
                label = "PERIGEE · 0° / 360°" if closed else "START PERIGEE · 0°"
                draw.text((ax-L.n(218), ay+L.n(44)), label, font=self.small, fill=MUTED)
                continue
            elif event["kind"] == "apogee":
                label = "APOGEE"
            else:
                label = "NEXT PERIGEE"
            direction = sample.real_angle_rad
            label_position = point(apsis_position, -L.n(110), direction)
            box = draw.textbbox((0, 0), label, font=self.small)
            draw.text((label_position[0]-(box[2]-box[0])/2, label_position[1]),
                      label, font=self.small, fill=MUTED)
        self.draw_earth(image, self.earth)
        self.centered_text(draw, 930, 690, "Earth", self.label, EARTH_BLUE)
        draw.rounded_rectangle((*L.p(32, 165), *L.p(490, 952)), radius=L.n(18),
                               fill=PANEL, outline=PANEL_EDGE, width=L.n(2))
        self.centered_text(draw, 262, 202, "PENDULUM VIEW", self.heading, CYAN)
        self.centered_text(draw, 262, 241, "Virtual Moon held vertical and fixed", self.small, MUTED)
        self.centered_text(draw, 262, 871, "Earth graphic marks the fixed frame", self.small, MUTED)
        self.centered_text(draw, 262, 898, "View rotates and pans · Same scale", self.small, MUTED)
        self.centered_text(draw, 262, 925, "Distance variations are preserved", self.small, MUTED)
        draw.rounded_rectangle((*L.p(1374, 165), *L.p(1870, 952)), radius=L.n(18),
                               fill=PANEL, outline=PANEL_EDGE, width=L.n(2))
        draw.text(L.p(1420, 202), "MOON & MEAN MOTION", font=self.heading, fill=CYAN)
        dates = self.presentation.interval or f"{self.track.start_epoch[:10]} → {self.track.end_epoch[:10]} (TDB)"
        draw.text(L.p(1420, 240), dates, font=self.small, fill=MUTED)
        draw.line((*L.p(1420, 272), *L.p(1823, 272)), fill=PANEL_EDGE, width=L.n(1))
        # Both rates are azimuth rates in the fixed camera plane.
        for i, label in enumerate(("Elapsed time", "Geometric angle", "Earth–Moon distance", "Orbital speed",
                                   "Virtual Moon angular rate", "Real Moon angular rate")):
            draw.text(L.p(1420, 300 + 56*i), label, font=self.label, fill=MUTED)
        draw.text(L.p(1420, 641), "ANGLE AT EARTH’S COM", font=self.heading, fill=YELLOW)
        self.centered_text(draw, 1622, 755, "Virtual Moon bearing − real Moon bearing", self.small, MUTED)
        draw.line((*L.p(1420, 802), *L.p(1823, 802)), fill=PANEL_EDGE, width=L.n(1))
        legend = ((CYAN, self.presentation.moon_legend, False),
                  (MAGENTA, "Virtual Moon · mean rate", True),
                  (YELLOW, "Earth-centred separation", False),
                  (GREEN, "Real Moon velocity", False))
        for i, (color, label, dashed) in enumerate(legend):
            y = 835 + 28*i
            if dashed:
                dotted_line(draw, L.p(1422, y), L.p(1470, y), color, L.n(2), L.n(10))
            else:
                line(draw, L.p(1422, y), L.p(1470, y), color, L.n(4))
            draw.text(L.p(1485, y-11), label, font=self.small, fill=MUTED)
        draw.text(L.p(65, 1013), self.presentation.footer,
                  font=self.small, fill=DIM)
        phase_note = ("Initially aligned" if abs(self.track.virtual_phase_rad) < 1e-12
                      else f"Initial offset: {math.degrees(self.track.virtual_phase_rad):+.2f}°")
        period = 2*math.pi/self.track.mean_angular_rate_rad_s/86400
        draw.text(L.p(65, 1040), f"Perigee start → full 360° revolution · Virtual mean period: {period:.6f} d · {phase_note}",
                  font=self.small, fill=DIM)
        return image

    def frame(self, elapsed_days):
        state = self.track.sample(elapsed_days)
        image = self.static.copy()
        draw, L = ImageDraw.Draw(image), self.layout
        pos = self.position(state.projected_position_km)
        virtual_pos = self.position(state.projected_virtual_position_km)
        # Both radial lines originate at Earth's COM. The comparison arc is
        # centred there too; it makes no assumption about a physical pivot.
        line(draw, pos, self.earth, CYAN, L.n(3), L.n(8))
        dotted_line(draw, self.earth, virtual_pos, MAGENTA, L.n(2), L.n(11))
        r = self.moon_radius
        self.draw_real_moon(draw, pos, state.body_marker_angle_rad)
        # This is the only arrow originating at the real Moon: the true
        # ephemeris velocity projected into the camera plane, at a fixed scale.
        velocity_end = (pos[0] + L.n(95)*state.projected_velocity_km_s[0],
                        pos[1] - L.n(95)*state.projected_velocity_km_s[1])
        arrow(draw, pos, velocity_end, GREEN, L.n(4), L.n(12))
        dotted_circle(draw, virtual_pos, r, MAGENTA, L.n(2), L.n(9))
        correction = state.angular_separation_rad
        if abs(correction) > math.radians(.1):
            count = max(5, math.ceil(abs(correction)*80))
            arc = [point(self.earth, L.n(174), state.real_angle_rad + correction*i/count)
                   for i in range(count+1)]
            draw.line(arc, fill=YELLOW, width=L.n(4))
        ex, ey = self.earth
        draw.ellipse((ex-L.n(3), ey-L.n(3), ex+L.n(3), ey+L.n(3)), fill=WHITE)
        rate = math.degrees(state.mean_angular_rate_rad_s)*86400
        real_rate = math.degrees(state.real_angular_rate_rad_s)*86400
        values = (f"{elapsed_days:,.2f} d", f"{math.degrees(state.geometric_angle_rad):.2f}°", f"{state.distance_km:,.0f} km",
                  f"{state.speed_km_s:.3f} km/s", f"{rate:.3f}°/d", f"{real_rate:.3f}°/d")
        for i, value in enumerate(values):
            self.right_text(draw, 1823, 325 + 56*i, value, self.value)
        self.centered_text(draw, 1622, 687, f"{math.degrees(correction):+.2f}°", self.big, YELLOW)
        self.draw_pendulum_view(image, state)
        return image

    def draw_pendulum_view(self, image, state):
        view = pendulum_view(state)
        L = self.layout
        earth = self.pendulum_position(view.earth_km)
        real = self.pendulum_position(view.real_moon_km)
        virtual = self.pendulum_position(view.virtual_moon_km)
        # Earth's identical graphic is fixed in main-view axes and turns here
        # with the coordinate transformation, not with a physical spin model.
        self.draw_earth(image, earth, view.earth_graphic_angle_rad)
        draw = ImageDraw.Draw(image)
        line(draw, earth, real, CYAN, L.n(3), L.n(8))
        dotted_line(draw, earth, virtual, MAGENTA, L.n(2), L.n(11))
        self.draw_real_moon(draw, real, view.real_marker_angle_rad)
        dotted_circle(draw, virtual, self.moon_radius, MAGENTA, L.n(2), L.n(9))
        correction = state.angular_separation_rad
        if abs(correction) > math.radians(.1):
            first_angle = state.real_angle_rad + view.earth_graphic_angle_rad
            count = max(5, math.ceil(abs(correction)*80))
            arc = [point(earth, L.n(174), first_angle+correction*i/count) for i in range(count+1)]
            draw.line(arc, fill=YELLOW, width=L.n(4))
        ex, ey = earth
        draw.ellipse((ex-L.n(3), ey-L.n(3), ex+L.n(3), ey+L.n(3)), fill=WHITE)

    def render(self, output, fps, duration, hold, ffmpeg):
        frames = round((duration + 2*hold)*fps)
        cmd = [ffmpeg, "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24",
               "-s", f"{self.layout.width}x{self.layout.height}", "-r", str(fps), "-i", "-",
               "-an", "-c:v", "libx264", "-preset", "medium", "-crf", "17", "-pix_fmt", "yuv420p",
               "-movflags", "+faststart", str(output)]
        # Stream frames; no large temporary image sequence is created.
        process = subprocess.Popen(cmd, stdin=subprocess.PIPE)
        assert process.stdin is not None
        try:
            for index in range(frames):
                playback = index/fps
                fraction = max(0.0, min(1.0, (playback-hold)/duration))
                process.stdin.write(self.frame(fraction*self.track.duration_days).tobytes())
                if index % max(1, fps*3) == 0:
                    print(f"Rendering {100*index/frames:.0f}%", flush=True)
            process.stdin.close()
            if process.wait() != 0:
                raise RuntimeError("FFmpeg encoding failed")
        except BaseException:
            # A failed encoder can also raise while closing its input pipe.
            # Preserve the original error and always reap the child process.
            try:
                process.stdin.close()
            except OSError:
                pass
            if process.poll() is None:
                process.terminate()
            process.wait()
            raise
        return frames


def find_ffmpeg(explicit=None):
    if explicit:
        return explicit
    executable = shutil.which("ffmpeg")
    if executable:
        return executable
    try:
        import imageio_ffmpeg
    except ImportError as exc:
        raise RuntimeError("Install dremt-derived-model/requirements/requirements-animation.txt or pass --ffmpeg PATH") from exc
    return imageio_ffmpeg.get_ffmpeg_exe()


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--resolution", default="1920,1080", help="WIDTH,HEIGHT, both even")
    p.add_argument("--fps", type=int, default=60)
    p.add_argument("--duration", type=float, default=30, help="playback seconds for one full revolution starting at perigee")
    p.add_argument("--hold", type=float, default=1.5, help="seconds at each endpoint")
    p.add_argument("--virtual-phase-deg", type=float, default=0, help="initial virtual bearing minus real bearing")
    p.add_argument("--orbit", choices=("ephemeris", "keplerian"), default="ephemeris",
                   help="cached actual orbit or the original animation's ideal two-body ellipse")
    p.add_argument("--frame-day", type=float, help="save one PNG at elapsed model days instead of video")
    p.add_argument("--output", type=Path, help="MP4 or PNG path relative to the clone root")
    p.add_argument("--reference", type=Path, help="clone-root-relative cached Horizons JSON covering a full orbit after perigee")
    p.add_argument("--ffmpeg", help="optional explicit encoder executable")
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        width, height = [int(n) for n in args.resolution.split(",")]
    except (ValueError, TypeError):
        raise SystemExit("resolution must be WIDTH,HEIGHT")
    if width < 640 or height < 360 or width % 2 or height % 2:
        raise SystemExit("resolution must be at least 640×360 with even dimensions")
    if args.fps <= 0 or args.duration <= 0 or args.hold < 0:
        raise SystemExit("fps and duration must be positive; hold must be nonnegative")
    if args.orbit == "keplerian":
        if args.reference:
            raise SystemExit("--reference applies only to the ephemeris orbit")
        from keplerian_comparison import KeplerianComparisonTrack
        print("Constructing the original animation's ideal Keplerian orbit", flush=True)
        track = KeplerianComparisonTrack(virtual_phase_rad=math.radians(args.virtual_phase_deg))
        presentation = Presentation(
            subtitle="Ideal Keplerian orbit compared with a virtual Moon advancing at the same mean angular rate",
            interval=f"One ideal orbit · {track.duration_days:.6f} days",
            moon_legend="Moon · ideal Keplerian orbit",
            footer="Ideal two-body ellipse · e = 0.0549 · Uniform orientation change · Body sizes enlarged",
        )
        default_name = "pendular_libration_keplerian_1080p.mp4"
    else:
        print("Loading the cached JPL ephemeris and constructing the mean-rate comparison", flush=True)
        track = build_perigee_comparison_track(args.reference, virtual_phase_rad=math.radians(args.virtual_phase_deg))
        presentation = Presentation()
        default_name = "dremt_pendulum_left_1080p.mp4"
    renderer = Renderer(track, width, height, presentation)
    output = project_path(args.output) if args.output else ANIMATION_OUTPUT_DIR / default_name
    output.parent.mkdir(parents=True, exist_ok=True)
    if args.frame_day is not None:
        if output.suffix.lower() != ".png":
            raise SystemExit("--frame-day requires an output filename ending in .png")
        renderer.frame(args.frame_day).save(output)
    else:
        if output.suffix.lower() != ".mp4":
            raise SystemExit("video output must end in .mp4")
        frames = renderer.render(output, args.fps, args.duration, args.hold, find_ffmpeg(args.ffmpeg))
        samples = [track.sample(float(day)) for day in np.linspace(0, track.duration_days, 113)]
        metadata = {
            "description": "JPL ephemeris Moon versus uniform-rate virtual Moon; angle at Earth's COM",
            "model": args.orbit, "title": "Pendular Libration Model",
            "interval_days": track.duration_days,
            "selection": "First perigee with a complete 360-degree projected revolution remaining in the cache",
            "end_definition": "First return to initial perigee bearing after a full revolution; not forced to be a later perigee",
            "geometric_zero": "Starting perigee direction in fixed axes; later apsidal drift retained",
            "apses": list(track.apses),
            "resolution": [width, height], "fps": args.fps,
            "frames": frames, "playback_seconds": frames/args.fps,
            "real_moon_source": "Cached NASA/JPL Horizons state and orientation samples, interpolated for playback",
            "virtual_mean_period_days": 2*math.pi/track.mean_angular_rate_rad_s/86400,
            "virtual_mean_angular_rate_deg_per_day": math.degrees(track.mean_angular_rate_rad_s)*86400,
            "virtual_initial_phase_deg": args.virtual_phase_deg,
            "angle_origin": "Earth's centre of mass",
            "angle_sign": "Virtual Moon bearing minus real Moon bearing, counterclockwise positive",
            "panel_order": ["pendulum", "top_down_orbit", "dashboard"],
            "pendulum_view": {
                "transform": "Rotate by -90 degrees minus virtual bearing, then translate the virtual Moon to a fixed anchor",
                "virtual_moon_anchor_px": list(renderer.pendulum_virtual),
                "main_km_to_px": renderer.km_to_px, "pendulum_km_to_px": renderer.km_to_px,
                "earth_radius_px_both_views": renderer.earth_radius,
                "moon_radius_px_both_views": renderer.moon_radius,
                "earth_graphic": "Marks fixed main-view inertial axes; apparent turning and panning are changes of view",
            },
            "assumptions": ["Virtual Moon's projected bearing advances uniformly at the supplied mean sidereal orbital rate",
                            "Virtual Moon shares the real Moon's instantaneous projected radius and hidden out-of-plane height",
                            "Earth-centred angular comparison; no physical suspension or pendulum force law inferred"],
            "earth_centred_angular_separation_range_deg": [min(math.degrees(s.angular_separation_rad) for s in samples),
                                                          max(math.degrees(s.angular_separation_rad) for s in samples)],
            "latitude_range_deg_not_displayed": [min(s.subearth_latitude_deg for s in samples),
                                                 max(s.subearth_latitude_deg for s in samples)],
        }
        if args.orbit == "ephemeris":
            metadata.update(ephemeris_start=track.start_epoch, ephemeris_end=track.end_epoch,
                            ephemeris_interval_days=track.duration_days)
        else:
            metadata.update(
                description="Ideal Keplerian Moon versus uniform-rate virtual Moon; angle at Earth's COM",
                selection="Original ideal-ellipse inputs; begin at perigee and show exactly one orbit",
                end_definition="Next perigee at 360 degrees after one ideal orbital period",
                geometric_zero="Perigee direction in fixed axes; no apsidal drift",
                real_moon_source="Analytic LunarOrbit from author-2d-animation/src/moon_libration/orbit.py",
                semi_major_axis_km=track.orbit.semi_major_axis_km, eccentricity=track.orbit.eccentricity,
                supplied_uniform_orientation_rate_rad_s=track.mean_angular_rate_rad_s,
                supplied_initial_orientation_rad=math.pi,
                effective_mu_km3_s2=track.mean_angular_rate_rad_s**2*track.orbit.semi_major_axis_km**3,
                relative_angle="M − f; the original animation displays its negative",
                assumptions=["Fixed planar two-body Keplerian ellipse using the original animation's a, e and T",
                             "Nonzero uniform orientation rate C = 2π/T, initially facing Earth at perigee",
                             "Virtual Moon shares instantaneous radius and advances uniformly at C",
                             "Coordinate pendulum view; no physical suspension, restoring force or energy transfer supplied"],
            )
        output.with_suffix(".json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2)+"\n")
    print(f"Saved {output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
