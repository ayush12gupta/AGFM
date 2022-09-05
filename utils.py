import subprocess, sys

def execute(cmd):
    subprocess.check_call(cmd, shell=True, stdout=sys.stdout, stderr=subprocess.STDOUT)