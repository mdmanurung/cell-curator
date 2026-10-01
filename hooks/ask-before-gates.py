#!/usr/bin/env python3
"""Claude Code PreToolUse hook: ask before a command that records human approval.

`assumptions approve` passes Gate A and `write-back apply` writes labels into the
user's data. An agent must not run either without a person saying so.
"""

import json
import re
import sys

GATED = re.compile(r"\bcell-curator\b.*\b(assumptions\s+approve|write-back\s+apply)\b")

event = json.load(sys.stdin)
command = str(event.get("tool_input", {}).get("command", ""))
if GATED.search(command):
    json.dump(
        {
            "hookSpecificOutput": {
                "hookEventName": "PreToolUse",
                "permissionDecision": "ask",
                "permissionDecisionReason": (
                    "cell-curator: this records a human approval; a person must confirm it."
                ),
            }
        },
        sys.stdout,
    )
