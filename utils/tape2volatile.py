import os
import subprocess


def run_command(command):
    try:
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()
        output = output.decode().strip().split('\n')
        error = error.decode().strip()
        return output if process.returncode == 0 else error
    except Exception as e:
        return str(e)


path = "/mss/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4"

# List all files in the specified path
try:
    files = os.listdir(path)
    for file in files:
        # Build the command for each file
        full_file_path = os.path.join(path, file)
        command = f"jcache get {full_file_path}"
        print(f"Running command: {command}")

        # Run the command
        result = run_command(command)

        # Print the result
#        if isinstance(result, list):
#            print("\n".join(result))
#        else:
#            print(result)
except FileNotFoundError:
    print(f"Path not found: {path}")
except Exception as e:
    print(f"An error occurred: {e}")
