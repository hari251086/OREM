"""Shared TLE line pairing for the scratch_rpe analysis scripts.

`tle_pairs(lines)` returns every (line1, line2) pair in file order, plus the name line ('0 ...') that precedes it if
there is one. It works for 2-line and 3-line (name + lines 1/2) files alike: it scans line by line and accepts a pair
only when line 1 is directly followed by a line 2 of the same catalogue number. It replaces the older pattern
`for k in range(0, len(lines) - 1, 2)`, which stays in step only for 2-line files and silently read every second
record of 3-line files (issue #56, THEORY section 19)."""


def tle_pairs(lines):
    out = []
    k = 0
    while k < len(lines) - 1:
        l1, l2 = lines[k].rstrip('\r\n'), lines[k + 1].rstrip('\r\n')
        if l1.startswith('1 ') and l2.startswith('2 ') and l1[2:7] == l2[2:7]:
            name = lines[k - 1].rstrip('\r\n') if k > 0 and lines[k - 1].startswith('0 ') else None
            out.append((name, l1, l2))
            k += 2
        else:
            k += 1
    return out
