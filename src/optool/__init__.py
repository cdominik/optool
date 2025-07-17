import os
import pathlib
import shutil
import subprocess
from .optool import *


MODULE_PATH = pathlib.Path(__file__).parent.resolve()
OPTOOL = shutil.which(f"{MODULE_PATH}/optool")


def make_optool(options=None):
    os.chdir(MODULE_PATH)
    cmd = ["make",]
    if options is not None:
        cmd.append(options)
    subprocess.Popen(cmd).wait()


if OPTOOL is None:
    q = input(
        "Cannot find 'optool' executable, shall I try to compile it using "
        "default options (make clean; make multi)? (y/n)"
    )
    if q.lower() == "y":
        make_optool(options="clean")
        make_optool(options="multi")
