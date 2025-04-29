from inftools.tistools.concatenate import trjcat
from inftools.tistools.trajtxt_conv import trajtxt_conv
from infretis.asyncrunner import aiorunner, future_list
import os

"""Merge multiple paths (trjcat) simultaneously.

This script can be adapted to also calculate CVs etc.

"""

TOP = "path_to_conf.gro"
SUBC = 5
WORKERS = 20

def wrapper(dic):
	i = dic["ps"]
	if not os.path.exists(f"{i}/traj_xtc.txt"):
		trajtxt_conv(i=f"{i}/traj.txt", o=f"{i}/traj_xtc.txt", r=SUBC)

	if not os.path.exists(f"{i}/accepted/{i}.xtc"):
		trjcat(out=f"{i}/accepted/{i}.xtc", traj=f"{i}/traj_xtc.txt",
		       traj_format="xtc", topology=TOP)


runner = aiorunner({}, n_workers=WORKERS)
runner.set_task(wrapper)
runner.start()
futures = future_list()

# here we only merge paths 101 to 1001, dummy variables
start, end = 101, 1001
for i in range(start, end):
	futures.add(runner.submit_work({"ps":i}))
    print("submitted:", i)


for i in range(start, end):
	future = futures.as_completed()
    print("Done:", i, future)


runner.stop()
