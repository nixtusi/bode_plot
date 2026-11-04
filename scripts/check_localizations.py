#!/usr/bin/env python3
"""Check the string catalog and, optionally, an Xcode-built app bundle.

    python3 scripts/check_localizations.py
    python3 scripts/check_localizations.py --app /path/to/BodePlotApp.app

Only active catalog entries are checked. Xcode can retain obsolete entries with
extractionState=stale; those do not represent strings the app currently uses.
"""

import argparse
import json
import plistlib
import re
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
CATALOG = ROOT / "BodePlotApp" / "Localizable.xcstrings"
MINIMUM_LANGUAGES = 50  # Includes the source language.

# Foundation format specifiers, including positional variants such as %1$@.
FORMAT = re.compile(r"%(?:(\d+)\$)?[-+ #0']*(?:\d+)?(?:\.\d+)?(?:hh|h|ll|l|q|z|t|j|L)?([@diuoxXfFeEgGaAcCsSp])")


def placeholders(value):
    """Return argument positions and types, allowing translations to reorder them."""
    matches = []
    next_position = 1
    index = 0
    while index < len(value):
        if value[index] != "%":
            index += 1
            continue
        if value[index:index + 2] == "%%":
            index += 2
            continue
        match = FORMAT.match(value, index)
        if not match:
            # A literal percent sign is valid in a string with no substitutions.
            index += 1
            continue
        position = int(match.group(1)) if match.group(1) else next_position
        if not match.group(1):
            next_position += 1
        matches.append((position, match.group(2)))
        index = match.end()
    return Counter(matches)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--app", type=Path, help="also validate Localizable.strings in a built .app")
    args = parser.parse_args()

    with CATALOG.open(encoding="utf-8") as file:
        catalog = json.load(file)
    source = catalog["sourceLanguage"]
    active = {
        key: entry for key, entry in catalog["strings"].items()
        if entry.get("extractionState") != "stale" and entry.get("shouldTranslate") is not False
    }
    languages = {source}
    for entry in active.values():
        languages.update(entry.get("localizations", {}))

    errors = []
    if len(languages) < MINIMUM_LANGUAGES:
        errors.append(f"Expected at least {MINIMUM_LANGUAGES} languages including {source}; found {len(languages)}")
    for key, entry in active.items():
        localizations = entry.get("localizations", {})
        for language in sorted(languages - {source}):
            unit = localizations.get(language, {}).get("stringUnit", {})
            value = unit.get("value")
            if unit.get("state") != "translated" or not isinstance(value, str) or not value.strip():
                errors.append(f"{language}: missing translated value for {key!r}")
            elif placeholders(value) != placeholders(key):
                errors.append(f"{language}: placeholders do not match for {key!r}: {value!r}")

    if args.app:
        for language in sorted(languages):
            path = args.app / f"{language}.lproj" / "Localizable.strings"
            try:
                with path.open("rb") as file:
                    bundled = plistlib.load(file)
            except (OSError, plistlib.InvalidFileException) as error:
                errors.append(f"{language}: cannot read {path}: {error}")
                continue
            for key, entry in active.items():
                expected = key if language == source else entry.get("localizations", {}).get(language, {}).get("stringUnit", {}).get("value")
                # An explicit source-language value can differ from the key.
                if language == source:
                    expected = entry.get("localizations", {}).get(source, {}).get("stringUnit", {}).get("value", key)
                if expected is not None and bundled.get(key) != expected:
                    errors.append(f"{language}: built value missing or differs for {key!r}")

    if errors:
        print("Localization checks failed:", file=sys.stderr)
        for error in errors:
            print(f"  {error}", file=sys.stderr)
        return 1
    print(f"Validated {len(active)} active keys in {len(languages)} languages" + (" and the app bundle" if args.app else "") + ".")
    return 0


if __name__ == "__main__":
    sys.exit(main())
