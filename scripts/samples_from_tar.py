'''Take a tar file with fastqs return fastq1s, fastq2s, and sample names.'''
import sys
import os
import tarfile
import argparse

def parse_args(args):
    '''Get arguments.'''
    #required args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p", "--paired",
        help='Flag for paired end. Will attempt to find paried reads based on \
        fastq names',
        action='store_true',
        required=False
        )
    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument(
        "--tar",
        help="Input tarball containing a directory with fastqs",
        required=True
        )

    #parse and return arguments
    args = parser.parse_args()
    paired = args.paired
    in_tar = args.tar

    return in_tar, paired

def main(args):
    '''Main.'''
    #parse arguments
    in_tar, paired = parse_args(args)

    #make output directories
    fastq1 = "./fastq1"
    fastq2 = "./fastq2"
    os.makedirs(fastq1, exist_ok=True)
    os.makedirs(fastq2, exist_ok=True)

    #open tarfile, find sample names and extract to correct folder
    sample_names = []
    sample_names_1 = []
    sample_names_2 = []
    with tarfile.open(in_tar, "r:gz") as tar:
        for entry in tar:
            if entry.isfile():
                #get the sample name (without .fastq.gz)
                sample_name = os.path.basename(entry.name.split('.')[0])
                #remove directory from path to extract file to output directory
                entry.name = os.path.basename(entry.name)
                if not paired:
                    tar.extract(entry, fastq1)
                    sample_names.append(sample_name)
                else:
                    #figure out if the file is read 1 or read 2
                    sample_split = sample_name.split("_")
                    read_number = sample_split[-1]
                    sample_split.pop()
                    sample_name = "_".join(sample_split)
                    #extract to the correct folder
                    if read_number not in ['1', '2']:
                        raise ValueError(sample_name + " has an invalid read number: " + read_number)
                    elif read_number == '1':
                        tar.extract(entry, fastq1)
                        sample_names_1.append(sample_name)
                    elif read_number == '2':
                        tar.extract(entry, fastq2)
                        sample_names_2.append(sample_name)

    if paired:
        #find unmatched samples
        unmatched = list(set(sample_names_1) - set(sample_names_2)) + \
            list(set(sample_names_2) - set(sample_names_1))
        if len(unmatched) != 0:
            raise ValueError("The following samples are missing a read pair " + ",".join(unmatched))
        sample_names = sample_names_1

    print(sample_names)

if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv)
