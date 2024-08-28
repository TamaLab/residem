import urllib.request


def get_pdb_mtz(pdb_id):
   url_mtz =  f"https://edmaps.rcsb.org/coefficients/%s.mtz"%pdb_id
   url_pdb = f"https://files.rcsb.org/download/%s.pdb"%pdb_id
   urllib.request.urlretrieve(url_mtz, f'{pdb_id}.mtz')
   urllib.request.urlretrieve(url_pdb, f'{pdb_id}.pdb')



def main():
    import argparse

    parser = argparse.ArgumentParser(prog='residem_pdb_mtz', formatter_class=argparse.RawTextHelpFormatter,
                                     description=f""" This is script downloads pdb and corresponding mtz file from pdb.
    Example: \n\n
    $ residem_pdb_mtz -r 5b6v
    
    """)


    parser.add_argument("-r", "--ref", type=str, help="reference pdb file to download")
    args = parser.parse_args()
    if len(args.ref) == 4:
        get_pdb_mtz(args.ref)
    else:
        print("Enter a Valid pdb-id, example: \n residem_pdb_mtz -r 5b6v")





if __name__ == "__main__":
    main()