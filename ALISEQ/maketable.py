"""
obtained from:
https://stackoverflow.com/questions/11347505/what-are-some-approaches-to-outputting-a-python-data-structure-to-restructuredte

example:

print make_table( [  ['Name', 'Favorite Food', 'Favorite Subject'],
                     ['Joe', 'Hamburgrs', 'I like things with really long names'],
                     ['Jill', 'Salads', 'American Idol'],
                     ['Sally', 'Tofu', 'Math']])

    ===== ============= ==================================== 
    Name  Favorite Food Favorite Subject                     
    ===== ============= ==================================== 
    Joe   Hamburgrs     I like things with really long names 
    ----- ------------- ------------------------------------ 
    Jill  Salads        American Idol                        
    ----- ------------- ------------------------------------ 
    Sally Tofu          Math                                 
    ===== ============= ====================================
"""

import pandas as pd


def make_table_from_csv(csv, sep=',', header=0):
    df = pd.read_csv(csv, sep=sep, header=header)
    return make_table_from_df(df)


def make_table_from_df(df):
    try:
        df = df.to_frame().astype(str)   # just in case it's actually a series
        grid = list([('attribute', 'value')]+[[str(k),str(v)] for k,v in df.itertuples()])
    except AttributeError:
        df = df.astype(str)
        grid = [list(df.columns),] + df.values.tolist()
    return make_table(grid)


def make_table(grid):
    max_cols = [max(out) for out in map(list, zip(*[[len(item) for item in row] for row in grid]))]
    rst = table_div(max_cols, 1)

    for i, row in enumerate(grid):
        header_flag = False
        if i == 0 or i == len(grid)-1: header_flag = True
        rst += normalize_row(row,max_cols)
        rst += table_div(max_cols, header_flag )
    return rst

def table_div(max_cols, header_flag=1):
    out = ""
    if header_flag == 1:
        style = "="
    else:
        style = "-"

    for max_col in max_cols:
        out += max_col * style + " "

    out += "\n"
    return out


def normalize_row(row, max_cols):
    r = ""
    for i, max_col in enumerate(max_cols):
        r += row[i] + (max_col  - len(row[i]) + 1) * " "

    return r + "\n"





