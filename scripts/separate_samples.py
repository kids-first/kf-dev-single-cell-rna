"""Take a tar file with fastqs return fastq1s, fastq2s, and sample names."""
import sys
import os
import pathlib
import shutil
import argparse


def parse_args(args):
    """Get arguments."""
    # required args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        "--paired",
        help="Flag for paired end. Will attempt to find paried reads based on \
        fastq names",
        action="store_true",
        required=False,
    )
    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument(
        "--dir", help="Input tarball containing a directory with fastqs", required=True
    )

    # parse and return arguments
    args = parser.parse_args()
    paired = args.paired
    in_dir = args.dir

    return in_dir, paired


def make_sample_name(file):
    "Extract sample name from fastq file"
    parts = file.name.split(".")
    parts = parts[: len(parts) - 2]
    sample_name = ".".join(parts)
    return sample_name


def main(args):
    """Main."""
    # parse arguments
    in_dir, paired = parse_args(args)

    # make output directories
    fastq1 = "./fastq1"
    fastq2 = "./fastq2"
    os.makedirs(fastq1, exist_ok=True)
    os.makedirs(fastq2, exist_ok=True)

    # open tarfile, find sample names and extract to correct folder
    sample_names = []
    sample_names_1 = []
    sample_names_2 = []
    for entry in pathlib.Path(in_dir).iterdir():
        if entry.is_file():
            # get the sample name (without .fastq.gz)
            sample_name = make_sample_name(entry)
            # remove directory from path to extract file to output directory
            basename = os.path.basename(entry.name)
            move_path = fastq1 + "/" + basename
            if not paired:
                shutil.move(entry, move_path)
                sample_names.append(sample_name)
            else:
                # figure out if the file is read 1 or read 2
                sample_split = sample_name.split(".")
                read_number = sample_split[-1]
                sample_split.pop()
                sample_name = ".".join(sample_split)
                # extract to the correct folder
                if read_number not in ["R1", "R2"]:
                    raise ValueError(
                        sample_name + " has an invalid read number: " + read_number
                    )
                elif read_number == "R1":
                    shutil.move(entry, move_path)
                    sample_names_1.append(sample_name)
                elif read_number == "R2":
                    move_path = fastq2 + "/" + basename
                    shutil.move(entry, move_path)
                    sample_names_2.append(sample_name)

    if paired:
        # find unmatched samples
        unmatched = list(set(sample_names_1) - set(sample_names_2)) + list(
            set(sample_names_2) - set(sample_names_1)
        )
        if len(unmatched) != 0:
            raise ValueError(
                "The following samples are missing a read pair " + ",".join(unmatched)
            )
        sample_names = sample_names_1

    print("\n".join(sorted(sample_names)))


if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv)
