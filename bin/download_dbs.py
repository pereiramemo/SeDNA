#!/usr/bin/env python3

"""
SeDNA database pre-fetch: populate ~/.sedna/db before a pipeline run.

Replaces the former install/get_databases.sh. Two databases are bind-mounted
into the module containers from --db_dir (params.db_dir in nextflow.config):

    <db_dir>/antismash/        module 1  (download-antismash-databases, in the module-1 image)
    <db_dir>/deepbgc/data/      module 2  (deepbgc download, in the module-2 image)

The antiSMASH and deepBGC downloaders exist only inside their images, so they
are run there with --db_dir mounted at the same path, by Docker or Apptainer
(--engine). Build or pull the images first (docker/dockerbuild_commands.sh or
bin/download_imgs.py). Fetching from inside the pipeline instead would make every
per-sample task race for the same download, and compute nodes often have no
outbound network at all, so run this once on a login node.

Every database is skipped when already present, so re-running is cheap.
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(__file__))
from utils import run, check_tools, check_file
# Reused so the SIF filename is derived in exactly one place
from download_imgs import APPTAINER_CACHE_DIR, IMAGE_TEMPLATE, cache_filename

DB_DIR_DEFAULT = os.path.join(os.path.expanduser("~"), ".sedna", "db")

ANTISMASH_MODULE = 1
DEEPBGC_MODULE = 2

###############################################################################
# 2. Define utility functions
###############################################################################

###############################################################################
# 2.1 Parse command-line arguments
###############################################################################

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"{os.path.basename(__file__)}: pre-fetch SeDNA databases",
        add_help=False)

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--db_dir", dest="db_dir", default=DB_DIR_DEFAULT,
        help=f"database root; must match params.db_dir (default: {DB_DIR_DEFAULT})")

    parser.add_argument("--databases", dest="databases", nargs="+",
        default=["antismash", "deepbgc"], choices=["antismash", "deepbgc"],
        help="databases to fetch (default: antismash deepbgc)")

    parser.add_argument("--engine", dest="engine", default="docker",
        choices=["docker", "apptainer"],
        help="container runtime used for the antiSMASH and deepBGC downloaders (default: docker)")

    parser.add_argument("--cache_dir", dest="cache_dir", default=APPTAINER_CACHE_DIR,
        help=f"Apptainer image cache, --engine apptainer only (default: {APPTAINER_CACHE_DIR})")

    return parser.parse_args()

###############################################################################
# 2.2 Run a command inside a module image, with db_dir mounted
###############################################################################

def run_in_image(module: int, cmd: list, db_dir: str, engine: str, cache_dir: str,
                 env: dict = None) -> None:
    img = IMAGE_TEMPLATE.format(n=module)
    env = env or {}

    if engine == "docker":
        check_tools(["docker"])
        # Run as the host user so the databases are not left owned by root;
        # HOME=/tmp gives the tools a writable home under that uid
        full_cmd = ["docker", "run", "--rm", "--volume", f"{db_dir}:{db_dir}",
                    "--user", f"{os.getuid()}:{os.getgid()}", "--env", "HOME=/tmp"]
        for key, value in env.items():
            full_cmd += ["--env", f"{key}={value}"]
        full_cmd += [img] + cmd
    else:
        check_tools(["apptainer"])
        sif = os.path.join(cache_dir, cache_filename(img))
        check_file(sif, f"Apptainer image {sif} (run bin/download_imgs.py first)")
        full_cmd = ["apptainer", "exec", "--bind", db_dir]
        for key, value in env.items():
            full_cmd += ["--env", f"{key}={value}"]
        full_cmd += [sif] + cmd

    run(full_cmd)

###############################################################################
# 2.4 antiSMASH databases
###############################################################################

def fetch_antismash(db_dir: str, engine: str, cache_dir: str) -> None:
    antismash_dir = os.path.join(db_dir, "antismash")
    # download-antismash-databases is itself idempotent, but needs the image;
    # a populated pfam/ subdir is enough to skip it entirely
    if os.path.isdir(os.path.join(antismash_dir, "pfam")):
        print(f"{antismash_dir} already present; skipping")
        return

    os.makedirs(antismash_dir, exist_ok=True)
    print(f"Downloading antiSMASH databases into {antismash_dir} ...")
    try:
        run_in_image(ANTISMASH_MODULE,
                     ["download-antismash-databases", "--database-dir", antismash_dir],
                     db_dir, engine, cache_dir)
    except subprocess.CalledProcessError:
        print("antiSMASH database download failed", file=sys.stderr)
        sys.exit(1)

###############################################################################
# 2.5 deepBGC models and Pfam copy
###############################################################################

def fetch_deepbgc(db_dir: str, engine: str, cache_dir: str) -> None:
    # deepBGC's own default is ~/.local/share/deepbgc/data, and it stores files under
    # <DEEPBGC_DOWNLOADS_DIR>/<release>/; matching that keeps --db_dir ~/.local/share valid
    deepbgc_dir = os.path.join(db_dir, "deepbgc", "data")
    if os.path.isdir(deepbgc_dir) and os.listdir(deepbgc_dir):
        print(f"{deepbgc_dir} already present; skipping")
        return

    os.makedirs(deepbgc_dir, exist_ok=True)
    print(f"Downloading deepBGC models into {deepbgc_dir} ...")
    try:
        run_in_image(DEEPBGC_MODULE, ["deepbgc", "download"], db_dir, engine, cache_dir,
                     env={"DEEPBGC_DOWNLOADS_DIR": deepbgc_dir})
    except subprocess.CalledProcessError:
        print("deepBGC download failed", file=sys.stderr)
        sys.exit(1)

###############################################################################
# 3. Define the main function
###############################################################################

def main() -> None:

    args = parse_args()
    db_dir = os.path.abspath(args.db_dir)
    os.makedirs(db_dir, exist_ok=True)

    if "antismash" in args.databases:
        fetch_antismash(db_dir, args.engine, args.cache_dir)

    if "deepbgc" in args.databases:
        fetch_deepbgc(db_dir, args.engine, args.cache_dir)

    print(f"{os.path.basename(__file__)} exited successfully")
    sys.exit(0)

###############################################################################
# 4. Run the main function
###############################################################################

if __name__ == "__main__":
    main()
