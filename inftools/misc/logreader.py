def log_reader1(inp, skip):
    """Computes acceptance probability per ensemble.

    Only shooting moves for now."""
    enss = {}
    accs = 0
    with open(inp, "r") as read:
        for idx, line in enumerate(read):
            if "shooted" in line and "sh sh" not in line:
                rip = line.rstrip().split()
                ens = int(rip[5])
                if ens not in enss:
                    enss[ens] = []
                status = read.readline().split()
                if status[3] == "ACC":
                    accs += 1
                if accs <= skip:
                    continue
                enss[ens].append(status[3] == "ACC")
                # if status[3] != "ACC" and ens not in (0, 1):
                #     print(ens, status[3])

    keys = list(sorted(enss.keys()))
    pacc = [100*sum(enss[key])/len(enss[key]) for key in keys]
    return keys, pacc, [len(enss[key]) for key in keys]
