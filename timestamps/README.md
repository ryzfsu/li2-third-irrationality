# Timestamps

[OpenTimestamps](https://opentimestamps.org) proofs for the assets of release `v1`.
Each proof commits to the SHA-256 of the file below and is anchored in Bitcoin block
970653 (2026-10-09 15:43:29 UTC).

| file | SHA-256 |
|---|---|
| `li2-third-irrationality-v1.pdf` | `3a9affcb4ed49d874aa9305a64d70dd4f098c801461bc72e27be99b6042b7938` |
| `li2-third-irrationality-v1.tar.gz` | `debfb3316174bc3e80fe2f48935102d18ff782a906dd610da0cf9941b6acfb06` |

To verify, download the release asset into this directory and run

    pip install opentimestamps-client
    ots verify timestamps/<file>.ots
