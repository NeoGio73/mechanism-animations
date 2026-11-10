import sys
from pathlib import Path

from .steps import STEPS, Reaction
from . import card, gas_phase, meshes, orbital_steps, video

STEPS.update({"orbital_sps": orbital_steps.orbital_sps,
              "orbital_track": orbital_steps.orbital_track,
              "orbital_cubes": orbital_steps.orbital_cubes,
              "meshes": meshes.build_meshes,
              "video_preview": video.video_preview,
              "video": video.video,
              "video_composite": video.video_composite,
              "gas_phase": gas_phase.gas_phase,
              "card": card.card})


def main(argv: list[str]) -> None:
    if len(argv) < 2:
        print("usage: python -m pipeline <reaction_dir> <step> [<step> ...]\n"
              f"steps: {', '.join(STEPS)}")
        sys.exit(1)
    rx = Reaction.load(Path(argv[0]))
    for step in argv[1:]:
        print(f"=== {rx.cfg['id']}: {step}", flush=True)
        STEPS[step](rx)


if __name__ == "__main__":
    main(sys.argv[1:])
