import sys
sys.path.append('path_to_psspy_module')
import psspy

import psspy
import psse34  # Change this based on your PSSE version

def initialize_psse():
    # Initialize PSSE
    psspy.psseinit(50)  # Initialize PSSE with 50 buses
    print("PSSE initialized.")

def load_case(case_file):
    # Load a power system case
    ierr = psspy.read(0, case_file)
    if ierr != 0:
        print(f"Error loading case: {ierr}")
    else:
        print(f"Case loaded from {case_file}.")

def run_power_flow():
    # Run a power flow solution
    ierr = psspy.fdns([0, 0, 0, 1, 0, 0, 0, 0, 0])
    if ierr != 0:
        print(f"Power flow solution failed: {ierr}")
    else:
        print("Power flow solution completed.")

def save_results(output_file):
    # Save the results
    ierr = psspy.save(output_file)
    if ierr != 0:
        print(f"Error saving results: {ierr}")
    else:
        print(f"Results saved to {output_file}.")

if __name__ == "__main__":
    initialize_psse()
    load_case("your_case_file.sav")  # Replace with your actual case file
    run_power_flow()
    save_results("results_file.sav")  # Replace with your desired results file
