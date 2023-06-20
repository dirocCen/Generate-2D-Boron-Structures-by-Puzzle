from ase.io.vasp import read_vasp
import numpy as np
from pdb import set_trace
from sagar.io.vasp import read_vasp
from pyvaspflow.utils import write_poscar
from itertools import combinations
from sagar.crystal.structure import symbol2number as s2n
from sagar.crystal.structure import Cell
import numpy as np
import os, shutil, time, random, string


def main():

    b1 = 0.122
    b2 = 0.123
    res = []
    for ii in range(11,21):
        for jj in range(100,201):
            tmp1 = ii/jj
            if b1 < tmp1 < b2:
                res.append(tuple([ii, jj, tmp1]))
    print(res)
    # set_trace()


def read_vacc_from_poscar():

    path0 = 'poscar/poscar-114/poscar-114-2'
    files = os.listdir(path0)
    vacr = np.array([])
    for ii in range(len(files)):
        path1 = os.path.join(path0, 'POSCAR' + str(ii))
        atoms = read_vasp(path1)
        tmp1 = len(atoms.numbers)
        tmp2 = 114 - tmp1
        tmp3 = tmp2 / 114
        vacr = np.hstack((vacr, tmp2))

    vacr2 = np.unique(vacr)
    print(vacr2)
    itp = np.where(vacr==14)[0]
    print(itp)
    # for ii in itp:


    set_trace()


def num2poscar(origin_poscar, seq, insert_ele, folder, idx):
    seq = seq.astype('int')
    orig_cell = read_vasp(origin_poscar)
    new_atoms = orig_cell.atoms
    new_atoms[seq] = s2n(insert_ele)
    new_cell = Cell(orig_cell.lattice, orig_cell.positions, new_atoms)
    if np.linalg.det(new_cell.lattice) < 0:
        new_cell.lattice[2] = -new_cell.lattice[2]
    write_poscar(new_cell, idx=idx, folder=folder)



def main_run(origin_poscar_path, seq_path, seq_line, save_path, idx):
    f = open(seq_path, 'r')
    # idx = seq_line
    data = f.readlines()
    line = data[seq_line-1]
    seq = line.strip()   # 删除尾部空格
    seq = seq.split()    # 以空格为分界线，分成数组
    seq = np.array(seq)  #
    # set_trace()
    num2poscar(origin_poscar_path, seq, 'Vac', save_path, idx)


def read_vacc_from_mag():

    atoms_num = 114
    vac_num = 14

    cor_path0 = 'mag_cor/poscar_%s' % atoms_num
    files = os.listdir(cor_path0)
    for file in files:
        num1 = file.split('_')[1]
        num2 = num1.split('-')[1]
        origin_poscar = 'primitive/extend_cell-%s/POSCAR-v-%s-%s' % (atoms_num, atoms_num, num2)
        save_path = 'vacc_0.1225_poscar/poscar_%s/poscar_%s_%s' % (atoms_num, atoms_num, num2)
        seq_path = 'seq/seq_v-%s/seq_%s' % (atoms_num, num2)
        cor_path1 = os.path.join(cor_path0, file)
        mag = np.load(cor_path1)['Mag']
        itp = np.where(np.sum(mag, axis=1)==(atoms_num-2*vac_num))[0]
        # set_trace()
        if len(itp) != 0:
            try:
                os.makedirs(save_path)
            except Exception as e:
                pass
            for ii, seq_line in enumerate(itp):
                # set_trace()
                seq_line = seq_line + 1
                main_run(origin_poscar, seq_path, seq_line, save_path, ii)



if __name__ == '__main__':
    read_vacc_from_mag()


