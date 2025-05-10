#!/usr/bin/env python3
"""
Generate a biotechnological hypothesis for a BGC using OpenAI’s Chat Completion API via HTTP (no external Python packages required).
"""

import argparse
import json
import os
import sys
import urllib.request

PROMPT_TEMPLATE = """
Given the following information about a biosynthetic gene cluster (BGC), formulate a hypothesis about its biotechnological potential, including possible bioactivities or applications. Base your hypothesis on known structure-activity relationships, ecological roles, and taxonomic context.

- Domain composition: {domains}
- BGC class: {bgc_class}
- Closest MIBiG cluster: {mibig_hit}
- Taxonomic affiliation: {taxonomy}
- Environmental context: {environment}

Be concise but informative. Mention potential compound class and application.
"""


def generate_prompt(domains, bgc_class, mibig_hit, taxonomy, environment):
    return PROMPT_TEMPLATE.format(
        domains=domains,
        bgc_class=bgc_class,
        mibig_hit=mibig_hit,
        taxonomy=taxonomy,
        environment=environment,
    )


def load_descriptors(json_path):
    try:
        with open(json_path, "r") as f:
            return json.load(f)
    except (IOError, json.JSONDecodeError) as e:
        sys.exit(f"Error reading JSON file {json_path}: {e}")


def call_openai_api(prompt, api_key, model="gpt-4", temperature=0.5):
    url = "https://api.openai.com/v1/chat/completions"
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {api_key}",
    }
    body = {
        "model": model,
        "messages": [{"role": "user", "content": prompt}],
        "temperature": temperature,
    }
    data = json.dumps(body).encode("utf-8")
    req = urllib.request.Request(url, data=data, headers=headers, method="POST")
    try:
        with urllib.request.urlopen(req) as resp:
            result = json.load(resp)
            return result
    except urllib.error.HTTPError as e:
        error_body = e.read().decode()
        sys.exit(f"API request failed ({e.code}): {error_body}")
    except Exception as e:
        sys.exit(f"API connection error: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate biotechnological hypothesis for a BGC using HTTP to OpenAI API"
    )
    parser.add_argument(
        "--json",
        required=True,
        help="Path to a JSON file with BGC descriptors",
    )
    parser.add_argument(
        "--openai-api-key",
        dest="openai_api_key",
        help="Your OpenAI API key (or set the OPENAI_API_KEY env var).",
    )
    parser.add_argument(
        "--model",
        default="gpt-4",
        help="OpenAI model to use (e.g., gpt-4 or gpt-3.5-turbo)",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Run basic self-tests and exit",
    )
    args = parser.parse_args()

    if args.test:
        # Basic tests for generate_prompt
        sample = {
            "domains": ["PKS", "NRPS"],
            "class": "hybrid",
            "mibig_hit": "BGC0000001",
            "taxonomy": "Streptomyces",
            "environment": "marine sediment"
        }
        prompt = generate_prompt(sample["domains"], sample["class"], sample["mibig_hit"], sample["taxonomy"], sample["environment"])
        assert "Domain composition: ['PKS', 'NRPS']" in prompt
        assert "BGC class: hybrid" in prompt
        print("All tests passed.")
        sys.exit(0)

    data = load_descriptors(args.json)
    domains     = data.get("domains", "N/A")
    bgc_class   = data.get("class", "N/A")
    mibig_hit   = data.get("mibig_hit", "N/A")
    taxonomy    = data.get("taxonomy", "N/A")
    environment = data.get("environment", "N/A")

    api_key = args.openai_api_key or os.getenv("OPENAI_API_KEY")
    if not api_key:
        sys.exit("Error: No OpenAI API key provided. Use --openai-api-key or set $OPENAI_API_KEY.")

    prompt = generate_prompt(domains, bgc_class, mibig_hit, taxonomy, environment)
    result = call_openai_api(prompt, api_key, model=args.model)

    # Extract the assistant message
    try:
        hypothesis = result["choices"][0]["message"]["content"].strip()
    except (KeyError, IndexError) as e:
        sys.exit(f"Unexpected API response format: {e}\n{json.dumps(result)}")

    print("\n--- Hypothesis ---\n")
    print(hypothesis)
