"""Stub: run a DRAMA ARES collision-probability analysis via pydrama.

Not exercised against real data yet -- ARES answers a different question
(collision probability against the debris environment) than OREM's re-entry
predictions, so it isn't part of the re-entry validation goal this folder
was set up for. Scaffolded here so the tool is wired up and importable if a
future task needs it.

See drama.ares.get_basic_config() for the full default config (mirrors
run_oscar_reentry.py's use of drama.oscar.get_basic_config()), and
drama.ares.run.__doc__ for the config schema.
"""

from drama import ares


def main():
    config = ares.get_basic_config()
    # TODO: override config fields for a real object/scenario before use.
    results = ares.run(config=config, parallel=False)
    print(results)


if __name__ == "__main__":
    main()
