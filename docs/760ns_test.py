import os
import subprocess
import tqdm
import shutil
import time


global working_dir
working_dir = os.getcwd()

os.makedirs("SVD/760ns", exist_ok=True)
shutil.copyfile(f'{working_dir}/05_760ns/Data_folder_0/map_dump_default/chain_A_U_csv/map_dump_full_negative_chain_A_U.csv',
                f'{working_dir}/SVD/760ns/05_760ns_map_dump_full_negative_chain_A_U.csv')

os.chdir("SVD/760ns")
template = f"residem_svd_csv -r {working_dir}/dark/5b6v.pdb -f 05_760ns_map_dump_full_positive_negative_chain_A_U.csv " \
    f"-n 93,94,111,174,175,176,177,178,179,180,181,182,211,212,213,214,215,216,217,218,219,300 -t density"


result = subprocess.Popen(template, shell=True)
os.chdir(working_dir)


# residem_svd_csv -r 5b6v.pdb -f map_dump_full_negative_chain_A_U.csv -n 93,94,111,174,175,176,177,178,179,180,181,182,211,212,213,214,215,216,217,218,219,300 -t density
