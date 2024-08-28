import os
import subprocess
import tqdm
import shutil
import time
import subprocess
import concurrent.futures


def run_result(template):
    try:
        result = subprocess.run(template, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True, shell=True, check=True)
        return result.returncode
    except subprocess.CalledProcessError as e:
        print(f"Background error: {e.stderr}")
        return e.returncode


br_reference = {"dark": "5b6v", "SVD": None, "SVD/all": None,
                "SVD/positive": None, "SVD/negative": None, "SVD/positive_172_614": None, "SVD/negative_172_614": None, "SVD/760ns": None}

bR_dict = {"01_16ns": "5b6w", "02_40ns": "5h2h", "03_110ns": "5h2i", "04_290ns": "5h2j", "05_760ns": "5b6x", "06_2us": "5h2k", "07_5.25us": "5h2l", "08_13.8us": "5h2m",
           "09_36.2us": "5b6y", "10_95.2": "5h2n", "11_250us": "5h2o", "12_657us": "5h2p", "13_1.725ms": "5b6z"}

# bR_dict = {"01_16ns": "5b6w", "02_40ns": "5h2h",
#            "05_760ns": "5b6x"}
global working_dir
working_dir = os.getcwd()
# print(working_dir)


def mkdir_bR(x, y):
    os.makedirs(x, exist_ok=True)
    os.chdir(x)
    if y is not None:
        template = f"residem_pdb_mtz -r {y}"
        # result = run_result(template)
        p1 = subprocess.Popen(template, shell=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    os.chdir(working_dir)
    return


def run_residem(x):
    os.chdir(os.path.join(working_dir, x))
    pathdir_dark_pdb = os.path.join(working_dir, 'dark/5b6v.pdb')
    pathdir_dark_mtz = os.path.join(working_dir, 'dark/5b6v.mtz')
    template = f"residem -r %s -m %s -t %s.mtz -v max" % (
        pathdir_dark_pdb, pathdir_dark_mtz, bR_dict[x])
    result = run_result(template)
    # result = subprocess.Popen(template, shell=True)
    os.chdir(working_dir)
    return result


[mkdir_bR(x, y) for x, y in bR_dict.items()]
[mkdir_bR(x, y) for x, y in br_reference.items()]
# time.sleep(20)


for x in list(bR_dict.keys()):
    run_residem(x)
