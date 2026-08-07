"""Stub: run a DRAMA SARA re-entry risk (casualty area) analysis via pydrama.

Not exercised against real data yet. SARA answers a different question
(ground-casualty risk from surviving re-entry fragments) than OREM's
re-entry *timing* predictions -- it's a natural next step once OSCAR-based
re-entry-date validation (run_oscar_reentry.py) is established, since SARA
needs a re-entry epoch/state as one of its inputs. Scaffolded here so the
tool is wired up and importable.

Note: SARA-RISK requires Java 17+ (installed at E:\\Java\\jdk-17 as part of
this DRAMA setup). See drama.sara.run.__doc__ for the full config schema.
"""

from drama import sara


def main():
    config = sara.get_basic_config()
    # TODO: override config fields -- e.g. feed in an OSCAR final_state
    # (see run_oscar_reentry.py's output) as the re-entry epoch/orbit.
    results = sara.run(config=config, parallel=False)
    print(results)


if __name__ == "__main__":
    main()
