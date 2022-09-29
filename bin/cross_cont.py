#! /usr/bin/env python3

import pandas as pd
import sys


class RG:
    def __init__(self,rg,p7,p5,count):
        self.rg = rg
        self.p7 = p7
        self.p5 = p5
        self.count = count
        self.cont = 0

    def print_summary(self):
        print(self.rg, round(self.cont,1), round(self.cont/self.count*100,4),
            sep='\t',
            file=sys.stdout
        )


def main(df):

    all_rg = df[df.RG.isin(['unexpected', 'unknown'])==False]

    #generate look-up tables
    p7_p5_rg_dict = {(p7,p5):rg for p7,p5,rg in zip(all_rg.p7ind, all_rg.p5ind, all_rg.RG)}
    p7_p5_count_dict = {(p7,p5):n for p7,p5,n in zip(df.p7ind, df.p5ind, df['#seqs'])}
    
    # get all corner combinations 
    p5_col_dict = {p7:set(x.p5ind) for p7,x in df.groupby('p7ind')}
    p7_row_dict = {p5:set(x.p7ind) for p5,x in df.groupby('p5ind')}

    # all rows and columns that have a shared p5 or p7 with the lib under observations
    # e.g. lib (1190, 1275)
    #p5ind     1233        1274        1275  1204   1381  
    #p7ind                                                                                  
    #1056       NaN         NaN         1.0   NaN    NaN  
    #1151       NaN         NaN       330.0   NaN    NaN 
    #1152       NaN         NaN        29.0   NaN    NaN 
    #1190      11.0      7703.0  18645395.0   1.0   99.0  
    #1191       NaN         NaN       155.0   NaN    NaN  
    #31         NaN         NaN       277.0   NaN    NaN  .
    #408        NaN         NaN         8.0   NaN    NaN

    # walk through all possible index combinations
    res = []
    top = []

    for rg,p7,p5,count in zip(all_rg.RG, all_rg.p7ind, all_rg.p5ind, all_rg['#seqs']):
         tmp = RG(rg,p7,p5,count)

         # now walk through the corner libraries
         
         for op7 in p7_row_dict[p5]:
             if op7==p7:
                 continue
             for op5 in p5_col_dict[p7]:
                 if op5 == p5:
                     continue
                 
                 corner1 = p7_p5_count_dict[(p7,op5)]
                 corner2 = p7_p5_count_dict[(op7,p5)]
                 min_corner = corner1 if corner1 < corner2 else corner2
                 cont = (min_corner/count)**2 * count
                 if cont > 0.5:
                     tmp.cont += cont
                 source =  p7_p5_rg_dict.get((op7,op5), f'{op7}/{op5}')
                 top.append((rg,source,cont,count))
         res.append(tmp)

    #sort the results
    res.sort(key=lambda x: x.rg)
    top.sort(key=lambda x: x[2], reverse=True)

    # print top10 contamination events
    print('#### TOP 10 EVENTS #####' , file=sys.stdout)
    print('# Source','Into_RG','CC_reads','Total_RG','CC_percent', sep='\t', file=sys.stdout) 
    for n,(rg,src,cont,count) in enumerate(top):
        if n==10:
            break
        print(
            f"# {src}",rg, f"{cont:.2f}",
            count, f"{(cont/count)*100:.6f}", 
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
    df = pd.read_csv(sys.argv[1], sep='\t', low_memory=False)
    df = df[df['#seqs'].astype(str).str.contains('^\d')]
    
    #filter phix-indices or missing indices
    df = df[(df['p5ind'].str.contains('^[0-9]+$'))&(df['p7ind'].str.contains('^[0-9]+$'))]
    df.convert_dtypes()
    main(df)


