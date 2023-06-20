# -*- coding:utf-8 -*-
from multiprocessing.pool import Pool
from pdb import set_trace
from sagar.io.vasp import read_vasp
from pyvaspflow.utils import write_poscar
from itertools import combinations
from sagar.crystal.structure import symbol2number as s2n
from sagar.crystal.structure import Cell
from sagar.crystal.derive import ConfigurationGenerator, cells_nonredundant
from sagar.toolkit.mathtool import refine_positions
from sagar.crystal.utils import non_dup_hnfs, snf
from sagar.crystal.derive import PermutationGroup as PG
from sagar.toolkit.mathtool import is_int_np_array
from sagar.crystal.derive import PermutationGroup
import multiprocessing as mp
import numpy as np
import os, shutil, time, random, string
from scipy import io



def sortrows(mat):
    n = np.shape(mat)[1]
    # import pdb; pdb.set_trace()
    idx = np.lexsort(tuple(mat[:, i] for i in range(n - 1, -1, -1)))
    return mat[idx]


def num2poscar(origin_poscar, seq, insert_ele, folder, idx):
    seq = seq.astype('int')
    orig_cell = read_vasp(origin_poscar)
    new_atoms = orig_cell.atoms
    new_atoms[seq] = s2n(insert_ele)
    new_cell = Cell(orig_cell.lattice, orig_cell.positions, new_atoms)
    if np.linalg.det(new_cell.lattice) < 0:
        new_cell.lattice[2] = -new_cell.lattice[2]
    write_poscar(new_cell, idx=idx, folder=folder)


def main_import_pool(w_ind):
    origin_poscar_list = []
    seq_list = []
    seq_line = []
    idx = []
    save_path = '/home/cenyj/study_document/pycode/ACM/FCC/acm_1_20Convexposcar'
    os.makedirs(save_path)
    for ii,line in enumerate(w_ind):
        if ii>0:
            idx.append(ii+1)
            path1 = '/home/cenyj/study_document/pycode/ACM/FCC/seq/poscar-v-{0}'.format(line[0])   # seq的一级目录
            path2 = path1 + '/poscar-v-{0}-{1}'.format(line[0], line[1])    # seq二级目录
            path3 = path2 + '/all_seq_{0}'.format(line[2])
            seq_list.append(path3)
            path4 = '/home/cenyj/study_document/pycode/ACM/FCC/primitive/extend_cell-{0}/cell_v{1}_id{2}.vasp'.format(line[0], line[0], line[1])
            origin_poscar_list.append(path4)
            seq_line.append(line[3])
    return origin_poscar_list, seq_list, seq_line, save_path, idx


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


def same_atoms_to_one_dir():
    atoms_num = 105
    files = os.listdir('seq/seq_v-%s' % atoms_num)
    ind = 0
    save_path = 'poscar/poscar_multiple_7/poscar-%s' % atoms_num
    try:
        os.makedirs(save_path)
    except Exception as e:
        print(e)
    for line in files:
        ii = line.split('_')[1]
        # set_trace()
        origin_poscar = 'primitive/extend_cell-%s/POSCAR-v-%s-%s' % (atoms_num, atoms_num, ii)
        seq = 'seq/seq_v-%s/seq_%s' % (atoms_num, ii)
        Mag = np.load('mag_cor/poscar_%s/poscar_%s-%s_Mag_cor.npz' % (atoms_num, atoms_num, ii))['Mag']
        num = len(Mag)
        # set_trace()
        for seq_line in range(1, num + 1):
            idx = ind
            main_run(origin_poscar, seq, seq_line, save_path, idx)
            ind = ind + 1
            # set_trace()


def multi_main():
    atoms_num = 79
    files = os.listdir('seq/seq_v-%s' % atoms_num)
    for line in files:
        ii = line.split('_')[1]
        # set_trace()
        origin_poscar = 'primitive/extend_cell-%s/POSCAR-v-%s-%s' % (atoms_num, atoms_num, ii)
        seq = 'seq/seq_v-%s/seq_%s' % (atoms_num, ii)
        save_path = 'poscar/poscar-%s/poscar-%s-%s' % (atoms_num, atoms_num, ii)
        try:
            os.makedirs(save_path)
        except Exception as e:
            print(e)
        Mag = np.load('mag_cor/poscar_%s/poscar_%s-%s_Mag_cor.npz' % (atoms_num, atoms_num, ii))['Mag']
        num = len(Mag)
        # set_trace()
        for seq_line in range(1, num + 1):
            idx = seq_line - 1
            main_run(origin_poscar, seq, seq_line, save_path, idx)


def single_main():
    num = 144
    name = '%s-test' % num
    origin_poscar = 'primitive/POSCAR-%s' % name
    seq = 'seq/seq_%s' % name
    save_path = 'poscar-test/poscar-%s' % name
    try:
        os.makedirs(save_path)
    except Exception as e:
        print(e)
    mag_path = 'mag_cor/'

    Mag = np.load('mag_cor/poscar_%s_Mag_cor.npz' % name)['Mag']
    num = len(Mag)
    # set_trace()
    for seq_line in range(1, num + 1):
        idx = seq_line - 1
        main_run(origin_poscar, seq, seq_line, save_path, idx)
        # set_trace()


if __name__ == '__main__':
    same_atoms_to_one_dir()


