def assign_code(log):
    paths = {}
    codes = {}
    inits = []
    inits0 = []
    init_b = False
    shootings = []
    shooting = 0

    with open(log) as read:
        for line in read:
            rip = line.rstrip().split()
            if "stored" in line:
                init_b = True

            if "shooted" in line:
                shooting += 1

            if "|" in line and rip[-4] == "|":
                pn = int(rip[1][1:])
                if init_b:
                    inits0.append(pn)
                if pn not in paths:
                    code = f"{rip[-1]}_{rip[-2]}_{rip[-3]}"
                    paths[pn] = code
                    shootings.append(shooting)
                    if code not in codes:
                        codes[code] = [pn]
                    else:
                        codes[code].append(pn)
            if "submit" in line and init_b:
                inits.append(inits0)
                inits0 = []
                init_b = False

    return paths, codes, inits, shootings


def add_move(replica, ens_arr, old, new, ens_num):
    replica_new = []
    ens_new = []

    for idx, (rep, ens)in enumerate(zip(replica, ens_arr)):
        last = rep.split("+")[-1]
        # last_ens = ens.split("+")[-1]
        if old == last:
            replica_new.append(rep + "+" + new)
            ens_new.append(ens + "+" + ens_num)
        else:
            replica_new.append(rep)
            ens_new.append(ens)
    return replica_new, ens_new


def read_log(log, pn):
    # what we want:
    # p0 -> p10 -> p16 -> p20 -> p23 -> p31 -> p34 -> p50
    # e0 ->  e1 ->  e0 ->  e1 ->  e1 ->  e0 ->  e0 ->  e1

    paths, codes, inits, shootings = assign_code(log)

    istep = -1
    replicas = []
    ensembles = []
    with open(log) as read:
        for line in read:
            rip = line.rstrip().split()
            if "stored" in line:
                if istep > -1:
                    replicas.append(replicas0)
                    ensembles.append(ens0)
                istep += 1
                replicas0 = [f"{i}" for i in inits[istep]]
                ens0 = [str(i) for i in range(len(replicas0))]

            if "shooted" not in line:
                continue

            # zero swap
            if "shooted" in line and len(rip) == 15:
                pno1, pno2 = rip[-5], rip[-4]
                pnn1, pnn2 = rip[-2], rip[-1]

                replicas0, ens0 = add_move(replicas0, ens0, pno1, pnn2, str(1))
                replicas0, ens0 = add_move(replicas0, ens0, pno2, pnn1, str(0))

            # shoot
            elif "shooted" in line:
                pno, pnn = rip[-3], rip[-1]
                ens = str(int(rip[-6]))
                if pno == pnn:
                    continue

                replicas0, ens0 = add_move(replicas0, ens0, pno, pnn, ens)
    replicas.append(replicas0)
    ensembles.append(ens0)

    follow_pn = []
    follow_en = []
    follow_size = []
    code = paths[pn]
    for idx, replica, ensemble in zip(range(len(replicas))[::-1], replicas[::-1], ensembles[::-1]):
        for rep, ens in zip(replica, ensemble):
            rep_codesl = [paths[int(pn0)] for pn0 in rep.split("+")]
            rep_codes = "+".join(rep_codesl)
            if code in rep_codes:
                follow_pn.append(rep)
                follow_en.append(ens)
                code = paths[int(rep.split("+")[0])]
                break
        follow_size.append((len(replica), len(rep_codesl)))

    pns = ("+".join(follow_pn[::-1])).split("+")
    ens = ("+".join(follow_en[::-1])).split("+")

    paths = list(paths.keys())
    shootings0 = {pn: sht for pn, sht in zip(paths, shootings)}
    shootings1 = []
    for pn in pns:
        shootings1.append(shootings0[int(pn)])

    return ens, pns, follow_size[::-1], shootings1
