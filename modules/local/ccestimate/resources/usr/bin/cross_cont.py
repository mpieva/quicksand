#! /usr/bin/env pypy3

import sys

class Event:
    def __init__(self, rg, source, count, cont):
        self.rg = rg
        self.source = source
        self.count = count
        self.cont = cont

class RG:
    def __init__(self, rg, p7, p5, count):
        self.rg = rg
        self.p7 = p7
        self.p5 = p5
        self.count = count
        self.cont = 0

    def print_summary(self):
        perc = (self.cont/self.count)*100
        print(self.rg, round(self.cont,1), "0" if perc==0 else f"{perc:.6f}",
            sep='\t',
            file=sys.stdout
        )

def main(rgs):
    all_rg = [x for x in rgs if x.rg not in ['unexpected', 'unknown']]

    # generate lookup tables
    p7_p5_rg_dict = {(rg.p7,rg.p5):rg.rg for rg in all_rg}
    p7_p5_count_dict = {(rg.p7,rg.p5):rg.count for rg in rgs}

    # get the corner combinations
    p5_col_dict = {}
    p7_row_dict = {}

    for rg in all_rg:
        p5_col_dict[rg.p7] = set([x.p5 for x in rgs if rg.p7 == x.p7])
        p7_row_dict[rg.p5] = set([x.p7 for x in rgs if rg.p5 == x.p5])

    # walk through all possible index combinations
    res = []
    top = []

    for rg in all_rg:
        for op7 in p7_row_dict[rg.p5]:
             if op7==rg.p7:
                 continue
             for op5 in p5_col_dict[rg.p7]:
                 if op5 == rg.p5:
                     continue

                 corner1 = p7_p5_count_dict[(rg.p7,op5)]
                 corner2 = p7_p5_count_dict[(op7,rg.p5)]
                 min_corner = corner1 if corner1 < corner2 else corner2
                 cont = (min_corner/rg.count)**2 * rg.count
                 if cont > 0.5:
                     rg.cont += cont
                 source = p7_p5_rg_dict.get((op7,op5), f'{op7}/{op5}')
                 top.append(Event(rg.rg, source, rg.count, cont))
        res.append(rg)

    res.sort(key=lambda x: x.rg)
    top.sort(key=lambda x: x.cont, reverse=True)

    print('#### TOP 10 EVENTS #####' , file=sys.stdout)
    print('# Source','Into_RG','CC_reads','Total_RG','CC_percent', sep='\t', file=sys.stdout)
    for n,ev in enumerate(top):
        if n==10:
            break
        print(
            f"# {ev.source}",ev.rg, f"{ev.cont:.2f}",
            ev.count, f"{(ev.cont/ev.count)*100:.6f}",
            sep='\t',
            file=sys.stdout
        )
    print('#', file=sys.stdout)

    #print each RG summary
    print('#### SUMMARY #####' , file=sys.stdout)
    print('#RG',"CC_readsum","CC_percent", sep='\t', file=sys.stdout)
    for x in res:
        x.print_summary()


if __name__ == '__main__':
    #import the dataframe
    rgs = []
    with open(sys.argv[1]) as infile:
        for line in infile.read().splitlines():
            if not line.startswith('#'):
                count, _, p7, _, p5, rg = line.split('\t',5)
                try:
                    int(p7)
                    int(p5)
                except ValueError:
                    continue
                rgs.append(RG(rg, int(p7), int(p5), int(count)))

    main(rgs)
