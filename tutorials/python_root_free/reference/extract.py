"""Extract Figure 4 central TMD curves from pdftocairo's page-9 SVG.

Input: pdftocairo -svg -f 9 -l 9 1707.00913.pdf figures/paper_page9.svg
The values are published plot coordinates, NOT author-supplied numerical tables.
Tick positions below were identified in the vector paths and checked visually.
"""
from pathlib import Path
from xml.etree import ElementTree as ET
import re
import csv
import argparse

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--page9',type=Path,required=True)
p.add_argument('--page10',type=Path,required=True)
p.add_argument('--output',type=Path,default=Path('reference'))
args=p.parse_args()
args.output.mkdir(parents=True,exist_ok=True)
paths=ET.parse(args.page9).findall('.//{http://www.w3.org/2000/svg}path')
specs=[(.13,1,537,137.812886,448.365639,228.560570,2992.120775,False),
       (.13,-1,559,137.812886,448.365639,228.560570,2992.120775,True),
       (.37,1,629,137.809817,448.362571,223.402835,2063.358450,False),
       (.37,-1,647,137.809817,448.362571,223.402835,2063.358450,True)]
with (args.output/'paper_figure4_vector.csv').open('w') as f:
    w=csv.writer(f,lineterminator='\n'); w.writerow(['z','charge','jt_GeV','A_UT_paper'])
    for z,charge,index,x01,x1,y0,yscale,endpoints in specs:
        path=paths[index]
        assert path.get('stroke') == ('rgb(0%, 0%, 0%)' if endpoints else 'rgb(100%, 0%, 0%)')
        tokens=path.get('d').split()
        assert not set(tokens)&{'C','Z'}
        xy=[(float(tokens[i+1]),float(tokens[i+2])) for i in range(0,len(tokens),3)]
        if endpoints: xy=xy[1::2]
        assert len(xy)==15 and abs(xy[0][0]-x01)<1e-5
        for x,y in xy:
            w.writerow([z,charge,10**((x-x01)/(x1-x01)-1),(y-y0)/yscale])
print('Saved vector-extracted central curves (15 published sampling points per curve).')

paths=ET.parse(args.page10).findall('.//{http://www.w3.org/2000/svg}path')
with (args.output/'paper_figure5_vector.csv').open('w') as f:
    w=csv.writer(f,lineterminator='\n'); w.writerow(['evolution','jt_GeV','cross_section_paper'])
    for evolution,index in [(1,350),(0,352)]:
        tokens=paths[index].get('d').split()
        xy=[(float(tokens[i+1]),float(tokens[i+2])) for i in range(0,len(tokens),3)][1::2]
        assert len(xy)==15
        for x,y in xy:
            j=.2+(x-138.915826)/(190.915946-138.915826)*.2
            value=10**(-13+(y-78.954496)/(378.166573-78.954496)*8)
            w.writerow([evolution,j,value])
