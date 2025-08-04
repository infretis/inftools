from inftools.tistools.concatenate import trjcat
from inftools.tistools.trajtxt_conv import trajtxt_conv
from infretis.asyncrunner import aiorunner, future_list
import os

"""Merge multiple paths (trjcat) simultaneously.

This script can be adapted to also calculate CVs etc.

To run, put this script into load/.
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

n_loops = end - start
while loop_cnt < n_loops + n_workers:
    if loop_cnt < n_workers:
        print("submit", loop_cnt + start)
        futures.add(runner.submit_work({"ps": loop_cnt + start}))
    else:
        future = futures.as_completed()
        print("completed", future, "keep going", loop_cnt <= n_loops, loop_cnt, n_loops)
        if loop_cnt <= n_loops:
            print("submit", loop_cnt + start)
            futures.add(runner.submit_work({"ps":loop_cnt + start}))
    loop_cnt += 1

runner.stop()
