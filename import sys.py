import sys

# Add the PSSE Python directory to the Python path
sys.path.append('C:\\Program Files (x86)\\PTI\\PSSE34\\Python\\Lib\\site-packages')

import psspy

def initialize_psse():
    psspy.psseinit(50)
    print("PSSE initialized.")

def load_case(case_file):
    ierr = psspy.read(0, case_file)
    if ierr != 0:
        print(f"Error loading case: {ierr}")
    else:
        print(f"Case loaded from {case_file}.")

if __name__ == "__main__":
    initialize_psse()
    load_case("your_case_file.sav")
