# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 16:42:58 2022

@author: yk5
"""

import os
from os import path

def line_process(line):
    line_=line.strip('\n').strip(' ')
    if len(line_)==16:
        serial_l=[line[6:11].strip(' '),line[11:].strip('\n').strip(' ')]
    if len(line_)==21:
        serial_l=[line[6:11].strip(' '),line[11:16].strip(' '),line[16:].strip('\n').strip(' ')]
    if len(line_)==26:
        serial_l=[line[6:11].strip(' '),line[11:16].strip(' '),line[16:21].strip(' '),line[21:].strip('\n').strip(' ')]
    if len(line_)==31:
        serial_l=[line[6:11].strip(' '),line[11:16].strip(' '),line[16:21].strip(' '),line[21:26].strip(' '),line[26:].strip('\n').strip(' ')]
    return serial_l

output_dir=r"x_lig"
filenames=os.listdir(r"covint_complex")
print(filenames)
for file in filenames:
    report_file = path.join(r'covint_report',file[:-4]+'.txt')
    with open (report_file,'r') as r:
        lines_r=r.readlines()
        lig_names=[]
        lig_position=[]
        gap_i=[]
        for line_r in lines_r:
            if 'SMALLMOLECULE' in line_r:
                lig_names.append(line_r[:3])
            if 'POLYMER' in line_r:
                lig_names.append(line_r[:3])
        input_file = path.join(r'covint_complex',file)
        with open (input_file,'r') as t:
            lines=t.readlines()
            for line in lines:
                if line[:6]=='ANISOU':
                    lines.remove(line)
            serial_atom=[]
            for line in lines:
                line__=line.split(' ')
                line_=[x.strip() for x in line__ if len(x.strip()) > 0]
                if line_[0]=='ATOM':
                    serial_atom.append(line_[1])
            lig_serial_line=[]
            lig_serial=[]
            serial_l=[]
            for line in lines:
                x=0
                y=0                
                if 'CONECT' == line[:6]:    
                    serial_l.append(line_process(line))
                    for line_i in line_process(line):
                        if line_i in serial_atom:
                            x+=1
                        if line_i not in serial_atom:
                            y+=1
                            lig_serial_=line_i
                    if x*y!=0:
                        lig_serial.append(lig_serial_)
                        lig_serial_line.append(line)
            lig_res_all=[]
            for lig_serial_ in lig_serial:
                for line in lines:
                    if line[:6]=='HETATM':
                        if line[6:11].strip(' ')==lig_serial_:                      
                            if line[17:20] in lig_names:
                                line_idx=lines.index(line)
                                res_idx=line[21:26]
                                line_idxa=line_idx
                                line_res=[]
                                while True:
                                    linea=lines[line_idxa]
                                    if linea[:6]=='HETATM':
                                        if linea[21:26]==res_idx:
                                            line_res.insert(0,linea)
                                        else:
                                            break
                                    else:
                                        break
                                    line_idxa-=1
                                line_idxa=line_idx
                                while True:
                                    line_idxa+=1
                                    linea=lines[line_idxa]
                                    if linea[:6]=='HETATM':
                                        if linea[21:26]==res_idx:
                                            line_res.append(linea)
                                        else:
                                            break
                                    else:
                                        break
                try:
                    lig_res_all.append(line_res)
                except:
                    continue
            sorted_list = sorted(lig_res_all, key=lambda x: len(x), reverse=True)
            try:
                line_lig=sorted_list[0]
            except:
                continue
            lig_index_l=[]
            for line in line_lig:
                lig_index_l.append(line[6:11].strip())
            for line in line_lig:
                lig_index=line[6:11].strip()
                for k in serial_l:       
                    if lig_index in k:
                        for i in lig_serial_line:
                            serial_l_=line_process(i)
                            if lig_index in serial_l_:                                
                                line_con=i
                                for j in serial_l_:
                                    if j not in lig_index_l:
                                        res_serial=j                                    
            for line in lines:
                if line[:4]=='ATOM':
                    if line[6:11].strip(' ')==res_serial:
                        print (1)
                        line_idx=lines.index(line)
                        res_idx=line[21:26]
                        line_idxa=line_idx
                        line_res=[]
                        while True:
                            linea=lines[line_idxa]
                            if linea[:4]=='ATOM':
                                print (line_res)
                                if linea[21:26]==res_idx:
                                    line_res.insert(0,linea)                                  
                                else:
                                    break
                            else:
                                break
                            line_idxa-=1
                        line_idxa=line_idx
                        while True:
                            line_idxa+=1
                            linea=lines[line_idxa]
                            if linea[:4]=='ATOM':
                                if linea[21:26]==res_idx:
                                    line_res.append(linea)
                                else:
                                    break
                            else:
                                break
            conect_idx=[]
            line_conect=[]
            output_file=path.join(output_dir,file)
            for j in line_lig:
                if j[-4]!='H':
                    conect_idx.append(j[6:11].strip())
            for k in lines:
                if k[:6]=='CONECT':
                    for k_ in line_process(k):
                        if k_ in conect_idx:
                            line_conect.append(k)
                            break
            with open (output_file,'w+') as g:
                for i in line_res:
                    if i[-4]!='H':
                        g.write(i)
                for j in line_lig:
                    if j[-4]!='H':
                        g.write(j)
                for k in line_conect:
                    g.write(k)
                # g.write(line_con)

                    
            


                    
            