import subprocess
import os
import pytest


def run_residem_pdb_mtz(pdb):
    # Command to be tested
    command = ["residem_pdb_mtz", "-r", f"{pdb}"]

    # Run the command
    result = subprocess.run(command, capture_output=True, text=True)

    # Check the return code
    assert result.returncode == 0, f"Command failed with return code%s" % result.returncode

    # Check for errors in stderr, but don't fail the test if there are warnings
    if result.stderr:
        print(f"Stderr output: {result.stderr}")

    return result


def run_residem():
    command = ["residem", "-r", "5b6v.pdb", "-m", "5b6v.mtz", "-t", "5b6x.mtz"]
    # Run the command
    result = subprocess.run(command, capture_output=True, text=True)

    # Check the return code
    assert result.returncode == 0, f"Command failed with return code%s" % result.returncode

    # Check for errors in stderr, but don't fail the test if there are warnings
    if result.stderr:
        print(f"Stderr output: {result.stderr}")

    return result


def check_output_files():
    # Expected files
    expected_files = ["5b6v.pdb", "5b6v.mtz",
                      "Data_folder_0/F_obs_minus_F_obs.ccp4"]

    for filename in expected_files:
        assert os.path.isfile(filename), f"Expected file {filename} not found"


def test_residem_pdb_mtz():
    # Run the commands
    run_residem_pdb_mtz("5b6v")
    run_residem_pdb_mtz("5b6x")
    run_residem()

    # Check for expected output files
    check_output_files()
