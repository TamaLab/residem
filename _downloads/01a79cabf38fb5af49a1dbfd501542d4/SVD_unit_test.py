import os
import subprocess
import tqdm
import shutil
import time

br_reference = {"dark": "5b6v", "SVD": None, "SVD/all": None,
                "SVD/positive": None, "SVD/negative": None, "SVD/negative_172_614": None, "SVD/positive_172_614": None}

bR_dict = {"01_16ns": "5b6w", "02_40ns": "5h2h", "03_110ns": "5h2i", "04_290ns": "5h2j",
           "05_760ns": "5b6x", "06_2us": "5h2k", "07_5.25us": "5h2l", "08_13.8us": "5h2m", "09_36.2us": "5b6y",
           "10_95.2": "5h2n", "11_250us": "5h2o", "12_657us": "5h2p", "13_1.725ms": "5b6z"}

# bR_dict = {"01_16ns": "5b6w",
#            "05_760ns": "5b6x"}

global working_dir
working_dir = os.getcwd()
print(working_dir)


def cp_all_svd_files(x):
    positive_file = f'{working_dir}/{x}/Data_folder_0/map_dump_default/chain_A_U_csv/map_dump_full_positive_chain_A_U.csv'
    positive_svd = f'{working_dir}/SVD/positive/{x}_map_dump_full_positive_chain_A_U.csv'
    positive_svd_172_614 = f'{working_dir}/SVD/positive_172_614/{x}_map_dump_full_positive_chain_A_U.csv'

    negative_file = f'{working_dir}/{x}/Data_folder_0/map_dump_default/chain_A_U_csv/map_dump_full_negative_chain_A_U.csv'
    negative_svd = f'{working_dir}/SVD/negative/{x}_map_dump_full_negative_chain_A_U.csv'
    negative_svd_172_614 = f'{working_dir}/SVD/negative_172_614/{x}_map_dump_full_negative_chain_A_U.csv'

    all_file = f'{working_dir}/{x}/Data_folder_0/map_dump_default/chain_A_U_csv/map_dump_full_positive_negative_chain_A_U.csv'
    all_svd = f'{working_dir}/SVD/all/{x}_map_dump_full_positive_negative_chain_A_U.csv'

    shutil.copyfile(positive_file, positive_svd)
    shutil.copyfile(negative_file, negative_svd)
    shutil.copyfile(positive_file, positive_svd_172_614)
    shutil.copyfile(negative_file, negative_svd_172_614)
    shutil.copyfile(all_file, all_svd)
    return None


length = len(list(bR_dict.keys()))


def run_svd(x, y, all="all"):
    os.chdir(x)
    f = open("time_and_file_name.txt", "w")
    f.write("{")
    for i, x in enumerate(list(bR_dict.keys())):
        if i+1 != length:
            f.write(f"\"{x[3:]}\":\"{x}_{y}\",")
        else:
            f.write(f"\"{x[3:]}\":\"{x}_{y}\"")

    f.write("}")
    f.close()
    if all == "all":
        template = f"residem_svd_csv -r {working_dir}/dark/5b6v.pdb --dict time_and_file_name.txt " \
            f"-n 172-614 -t density"
    else:
        template = f"residem_svd_csv -r {working_dir}/dark/5b6v.pdb --dict time_and_file_name.txt " \
            f"-n 93,94,111,174,175,176,177,178,179,180,181,182,211,212,213,214,215,216,217,218,219,300 -t density"

    result = subprocess.Popen(template, shell=True)
    os.chdir(working_dir)
    # return #result.returncode


[cp_all_svd_files(x) for x in list(bR_dict.keys())]

map_list = ["map_dump_full_positive_negative_chain_A_U.csv", "map_dump_full_positive_chain_A_U.csv",
            "map_dump_full_negative_chain_A_U.csv"]
for x in list(br_reference.keys())[2:]:
    if x == "SVD/all":
        run_svd(x, map_list[0], "no")
    elif x == "SVD/positive":
        run_svd(x, map_list[1], "no")
    elif x == "SVD/negative":
        run_svd(x, map_list[2], "no")
    elif x == "SVD/positive_172_614":
        run_svd(x, map_list[1], "all")
    elif x == "SVD/negative_172_614":
        run_svd(x, map_list[2], "all")


os.makedirs("SVD/760ns", exist_ok=True)
shutil.copyfile(f'{working_dir}/05_760ns/Data_folder_0/map_dump_default/chain_A_U_csv/map_dump_full_negative_chain_A_U.csv',
                f'{working_dir}/SVD/760ns/05_760ns_map_dump_full_negative_chain_A_U.csv')

os.chdir("SVD/760ns")
template = f"residem_svd_csv -r {working_dir}/dark/5b6v.pdb -f 05_760ns_map_dump_full_negative_chain_A_U.csv " \
    f"-n 93,94,111,174,175,176,177,178,179,180,181,182,211,212,213,214,215,216,217,218,219,300 -t density"

result = subprocess.Popen(template, shell=True)
os.chdir(working_dir)
exit()
