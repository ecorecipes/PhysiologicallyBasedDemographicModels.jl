#!/usr/bin/env python3

import argparse
import base64
import json
import os
import pathlib
import sys
import urllib.request


def b64(path: pathlib.Path) -> str:
    return base64.b64encode(path.read_bytes()).decode("utf-8")


def compare(model: str, generated: pathlib.Path, reference: pathlib.Path) -> dict:
    prompt = (
        "Compare the first image (generated vignette figure) against the second image "
        "(reference paper figure). Focus on curve shapes, peak timing, ranking of scenarios, "
        "axes/labels, and major mismatches. Return strict JSON with keys "
        "overall_match, similarities, differences, and verdict."
    )
    payload = {
        "model": model,
        "stream": False,
        "messages": [
            {
                "role": "user",
                "content": prompt,
                "images": [b64(generated), b64(reference)],
            }
        ],
        "format": {
            "type": "object",
            "properties": {
                "overall_match": {"type": "string"},
                "similarities": {"type": "array", "items": {"type": "string"}},
                "differences": {"type": "array", "items": {"type": "string"}},
                "verdict": {"type": "string"},
            },
            "required": ["overall_match", "similarities", "differences", "verdict"],
        },
    }
    req = urllib.request.Request(
        os.environ.get("OLLAMA_API", "http://127.0.0.1:11434/api/chat"),
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req) as response:
        raw = json.load(response)
    content = raw["message"]["content"]
    return json.loads(content)


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare generated and reference figures with a local VLM.")
    parser.add_argument("generated")
    parser.add_argument("reference")
    parser.add_argument("--model", default=os.environ.get("PBDM_VLM_MODEL", "llama3.2-vision:latest"))
    args = parser.parse_args()

    generated = pathlib.Path(args.generated)
    reference = pathlib.Path(args.reference)
    if not generated.is_file():
        raise SystemExit(f"Missing generated figure: {generated}")
    if not reference.is_file():
        raise SystemExit(f"Missing reference figure: {reference}")

    result = compare(args.model, generated, reference)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
