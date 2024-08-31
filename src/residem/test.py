import argparse

def create_parser():
    parser = argparse.ArgumentParser(description="Process some files.")
    parser.add_argument('-i', '--input', required=True, help='Specifies the input file.')
    parser.add_argument('-o', '--output', required=True, help='Specifies the output file.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Run in verbose mode.')
    return parser

if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    # Insert your file processing logic here
