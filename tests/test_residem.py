import subprocess
import os
import pytest




def run_residem():
    pdb_file = "tests/5b6v.pdb"
    mtz_file1 = "tests/5b6v.mtz"
    mtz_file2 = "tests/5b6x.mtz"

    command = ["residem", "-r", pdb_file, "-m", mtz_file1, "-t", mtz_file2]
    
    # Run the command
    result = subprocess.run(command, capture_output=True, text=True)

    # Check the return code
    assert result.returncode == 0, f"Command failed with return code {result.returncode}"

    return result


def check_output_files():
    # Expected files
    expected_files = ["Data_folder_0/F_obs_minus_F_obs.ccp4"]

    for filename in expected_files:
        assert os.path.isfile(filename), f"Expected file {filename} not found"


def test_residem_pdb_mtz():
    # Run the commands
    run_residem()

    # Check for expected output files
    check_output_files()
