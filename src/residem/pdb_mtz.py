import urllib.request
import subprocess


def run_result(template):
    try:
        result = subprocess.run(template, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True, shell=True, check=True)
        return result.returncode
    except subprocess.CalledProcessError as e:
        print(f"Background error: {e.stderr}")
        return e.returncode


def get_pdb_mtz(pdb_id):
   url_mtz =  f"https://files.rcsb.org/download/%s-sf.cif"%pdb_id
   url_pdb = f"https://files.rcsb.org/download/%s.pdb"%pdb_id
   urllib.request.urlretrieve(url_mtz, f'{pdb_id}.mtz')
   urllib.request.urlretrieve(url_pdb, f'{pdb_id}.pdb')
   run_result(f"sf_convert -o mtz -sf %s-sf.cif -out %s.mtz -pdb %s.pdb"%(pdb_id,pdb_id,pdb_id))
   run_result("rm sf_format_guess.text")
   run_result("rm sf_information.cif")




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