from B_reactor import Reactor
import time
import subprocess
import os

banner = r"""
          +=============================================================+
          | ██████╗  ██████╗  ██████╗ ███████╗████████╗███████╗██████╗  |
          | ██╔══██╗██╔═══██╗██╔═══██╗██╔════╝╚══██╔══╝██╔════╝██╔══██╗ |
          | ██████╔╝██║   ██║██║   ██║███████╗   ██║   █████╗  ██████╔╝ |
          | ██╔══██╗██║   ██║██║   ██║╚════██║   ██║   ██╔══╝  ██╔══██╗ |
          | ██║  ██║╚██████╔╝╚██████╔╝███████║   ██║   ███████╗██║  ██║ |
          | ╚═╝  ╚═╝ ╚═════╝  ╚═════╝ ╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝ |
          +=============================================================+

"""
name      = "          | Robust Object-Oriented Solver of Transport Equations in a Reactor\n"
copyright = "Copyright | 2020-2023, Advanced Nuclear Systems Group, Paul Scherrer Institute\n"
author    = "Author    | Konstantin Mikityuk, konstantin.mikityuk@psi.ch\n"
date      = "Date      | " + time.ctime() + "\n"
try:
    if os.environ.get('ROOSTER_COMMIT_ID') is not None:
        commit_id = os.environ.get('ROOSTER_COMMIT_ID')
    else:
        commit_id = subprocess.check_output(["git","rev-parse","HEAD"]).strip().decode("utf-8")
except:
    commit_id = "unkown git hash, please install git or run in the repo."

gitHash   = "Git HASH  | " + commit_id +"\n"

rooster = banner+name+copyright+author+gitHash+date

print(rooster)

# create and solve
r = Reactor()
