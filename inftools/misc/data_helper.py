def data_reader(inp):
    """Read infretis_data.txt type of file."""
    # collect all data into list of path dicts.
    paths = []

    # in case we get .gz data file..
    if inp[-3:] == ".gz":
        import gzip
        oopen = gzip.open
        readmode = "rt"
    else:
        oopen = open
        readmode = "r"

    with oopen(inp, readmode) as read:
        ensl = 0

        # read line
        for line in read:
            if line.startswith("#"):
                continue

            rip = line.rstrip().split()
            # get num ensembles from first-non comment line
            if not ensl:
                ensl = int(len(rip[3:])/2)

            # get line data
            pn, len0, max_op = rip[:3]
            path = {"pn": pn, "len": len0, "max_op": max_op}
            path["cols"] = {}
            f0l, w0l = rip[3:ensl+3], rip[3+ensl:2*ensl + 3]

            # skip if no weights
            if set(f0l) == set(w0l) == set(("----",)):
                continue

            # store only the weights based on col
            for col, (f0, w0) in enumerate(zip(f0l, w0l)):
                if '----' in (f0, w0):
                    continue
                path["cols"][col] = [f0, w0]

            # append to paths list
            paths.append(path)
        return paths


