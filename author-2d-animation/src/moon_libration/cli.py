"""Command-line entry point for rendering the animation."""

from __future__ import annotations

import argparse
from pathlib import Path

from moon_libration.scene import VideoRenderer


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Render the Moon longitudinal-libration animation")
    p.add_argument("--start", choices=("perigee", "apogee"), default="perigee")
    p.add_argument("--duration", type=float, default=30.0, help="seconds for one orbit")
    p.add_argument("--hold", type=float, default=1.5, help="still-frame seconds at each end")
    p.add_argument("--fps", type=int, default=60)
    p.add_argument("--resolution", default="1920,1080", help="WIDTH,HEIGHT")
    p.add_argument("--output", type=Path, default=Path("author-2d-animation/outputs/moon_libration_1080p.mp4"))
    p.add_argument("--ffmpeg", default="ffmpeg", help="FFmpeg executable")
    return p


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.duration <= 0 or args.hold < 0 or args.fps <= 0:
        raise SystemExit("duration/fps must be positive and hold must be nonnegative")

    project = Path(__file__).resolve().parents[3]
    args.output = (project / args.output).resolve() if not args.output.is_absolute() else args.output
    args.output.parent.mkdir(parents=True, exist_ok=True)
    try:
        width, height = (int(part) for part in args.resolution.split(","))
    except ValueError as exc:
        raise SystemExit("resolution must be WIDTH,HEIGHT") from exc
    VideoRenderer(width, height, args.fps, args.start).render(
        args.output, args.duration, args.hold, args.ffmpeg
    )
    print(f"Rendered {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
