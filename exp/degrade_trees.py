import os
import argparse
import shutil

from lib.util import check_make_dir
from lib.log_config import setup_logging
from lib.log_config import append_date

setup_logging()
import logging


def is_graph_file(f: str) -> bool:
    # f is a filename / basename
    return f.split(".")[-1] == "txt"


def main():
    # default number of degradations
    n: int = 10

    # args: setname (assumed to only contain trees)
    parser = argparse.ArgumentParser(
        prog="aspi_degrade_trees",
        description="Degrade pseudotrees incrementally.",
        usage="Pass path to -s (required) for src (existing set of trees). Pass an int to -n (optional) for number of degradations, default is 10.",
    )
    parser.add_argument(
        "-s", "--src_set", type=str, help="Path to the src set.", required=True
    )
    parser.add_argument(
        "-n", "--degradations", type=int, help="Number of degradations."
    )

    args = parser.parse_args()
    if args.degradations:
        n = args.degradations
    src_path: str = os.path.abspath(args.src_set)
    src_sn: str = os.path.basename(src_path)
    print(src_path)
    print(src_sn)

    # create target directory target
    temp_path: str = os.path.join(os.path.dirname(src_path), f"{src_sn}deg{str(n)}")
    dst_path: str = check_make_dir(temp_path, 0, delim="_")
    print(dst_path)
    dst_sn: str = os.path.basename(dst_path)

    # for every .txt (graph) file in set
    for fn in os.listdir(src_path):
        src_fp = os.path.join(src_path, fn)
        if not os.path.isfile(src_fp):
            continue
        if not is_graph_file(fn):
            continue

        print(fn)
        dst_fn: str = fn.replace(src_sn, dst_sn)
        dst_fp: str = os.path.join(dst_path, dst_fn)
        shutil.copyfile(src_fp, dst_fp)

        if os.path.isfile(dst_fp):
            print(f"copied {fn} to {dst_fp}")

        # degrade the graph (cumulatively) n times

        # save result of each degradation to target


if __name__ == "__main__":
    main()
