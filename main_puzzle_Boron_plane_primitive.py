import multiprocessing
from pdb import set_trace
import pretty_errors
import numpy as np
from scipy import io
from numpy.matlib import repmat
from ase.io.vasp import read_vasp, write_vasp
from collections import Counter
import os
import time
from scipy import io
# import matplotlib.pyplot as plt
import itertools
import time
from sagar.io.vasp import read_vasp as read
from sagar.toolkit.mathtool import refine_positions



def expand_cell(basis, bb, n1, n2, n3):
    b = np.dot(bb, basis)
    all = []
    for ii in range(n1):
        for jj in range(n2):
            for kk in range(n3):
                tb = b + repmat(np.dot([ii, jj, kk], basis), len(b), 1)
                if len(all) == 0:
                    all = tb
                else:
                    all = np.vstack((all, tb))
    basis = np.array([n1 * basis[0], n2 * basis[1], n3 * basis[2]])
    bb = np.dot(all, np.linalg.inv(basis))
    b = all

    return bb, basis


def sortrows(a):
    return a[np.lexsort(np.rot90(a))]


def findByRow(mat, row):
    return np.where((mat == row).all(1))[0]


def join_B_A(A, B):
    tmp = []
    # set_trace()
    B = np.unique(B)
    for ii, ele in enumerate(B):
        if np.sum(np.where(A == ele, 0, 0)) > 0:
            tmp.append(ii)
    B = np.delete(B, tmp)
    if len(B) != -1:
        r = np.hstack((A, B))
    else:
        r = np.copy(A)
    return r


def first_distance(basis, bb):
    b = np.dot(bb, basis)
    tmp = np.unique(np.round(np.sqrt(np.sum((b - repmat(b[0], len(b), 1)) ** 2, axis=1)), 3))
    return tmp[1]


def get_perms(cell,str_type='crystal',symprec=1e-3):
    latt = cell.lattice
    pos = cell.positions
    pos = np.dot(pos,latt)
    if str_type == "crystal":
        symm = cell.get_symmetry(symprec=symprec)
        trans,rots = symm['translations'],symm['rotations']
        perms = np.zeros((np.shape(trans)[0],len(cell.atoms)))
        origin_positions = refine_positions(cell.positions)
        for ix, rot in enumerate(rots):
            for iy,o_pos in enumerate(origin_positions):
                new_pos = np.dot(rot,o_pos.T) + trans[ix]
                new_pos = np.mod(new_pos,1)
                new_pos = refine_positions(new_pos)
                idx = np.argmin(np.linalg.norm(new_pos-origin_positions,axis=1))
                perms[ix,iy] = idx
        perms_table = np.unique(perms,axis=0)
    else:
        pass
        # mol = Molecule(pos,cell.atoms)
        # perms_table = mol.get_symmetry_permutation(symprec)
    return perms_table


def find_order(NNindex):
    order = np.array([0])
    tmp = np.array([0])
    tmp1 = join_B_A(tmp, NNindex[0])
    tmp2 = np.array(range(len(NNindex)))
    tmp2 = tmp2.reshape(tmp2.shape[0], 1)
    tmp2 = np.hstack((tmp2, NNindex))
    ii = 0
    while len(order) < len(NNindex):
        tmp3 = np.delete(tmp2, order, axis=0)
        tmp4 = np.copy(tmp3)
        # tmp2 = np.delete(tmp2, order[ii] - ii, axis=0)
        tmp1 = join_B_A(tmp1, tmp2[order[ii]])
        # set_trace()
        for jj in tmp1:
            tmp4[tmp4 == jj] = -1
        tmp4 = np.sum(np.where(tmp4 == -1, 1, 0), 1)
        # set_trace()
        tmp4 = tmp3[np.where(tmp4==np.max(tmp4))[0][0]][0]
        # set_trace()
        order = np.hstack((order, tmp4))

        ii += 1
    return order


def NNfind_3(basis, a):
    b = np.dot(a, basis)  # direct to xyz
    # b = a   # for Car

    NN = 0.68 * 1.25 * np.array([[1, 1.732050807, 0],[2,  0, 0],[1, -1.732050807, 0],
       [-1, -1.732050807, 0],[-2, 0, 0],[-1, 1.732050807, 0],[2, 3.464101615, 0],
       [3, 1.732050807, 0],[4, 0, 0],[3, -1.732050807, 0],[2, -3.464101615, 0],
       [0, -3.464101615, 0],[-2, -3.464101615, 0],[-3, -1.732050807, 0], [-4, 0.000000000, 0],
       [-3, 1.732050807, 0],[-2, 3.464101615, 0],[0, 3.464101615, 0]])
    NNindex = np.zeros((len(b), len(NN)), dtype=np.int16)
    all = []
    for ii in range(-5, 6):
        for jj in range(-5, 6):
            for kk in range(-5, 6):
                tb = b + repmat(np.dot([ii, jj, kk], basis), len(b), 1)
                if len(all) == 0:
                    all = tb
                else:
                    all = np.vstack((all, tb))

    tmp = np.array(range(len(b)))
    xyz = np.hstack((all, repmat(tmp.reshape(tmp.shape[0], 1), int(len(all)/len(b)), 1)))
    for ii in range(len(b)):
        for jj in range(len(NN)):
            # set_trace()
            tp = b[ii] + NN[jj]
            tmp1 = (repmat(tp, len(all), 1) - all) ** 2
            # tmp1 = () ** 2
            dtp = tmp1.sum(axis=1)
            q = dtp.min()
            qq = np.argmin(dtp)
            NNindex[ii, jj] = xyz[qq, 3]
            # set_trace()
            if q > 1:
                print('error')
    return NNindex


def prep_primitive_path():
    primitive = []
    path0 = 'primitive/extend_cell-'
    w1_2 = []
    for ii in range(2, 31):
        path1 = path0 + str(ii)
        files = os.listdir(path1)
        for jj, file in enumerate(files):
            path2 = os.path.join(path1, 'POSCAR-v-' + str(ii) + '-' + str(jj))
            primitive.append(path2)
            w1_2.append([ii, jj])
            # set_trace()
    return primitive, w1_2


def transform_ele2all_nonsymmetry(ele, perms):
    ele = ele[perms]
    ele = np.unique(ele, axis=0)
    return ele


def puzzle(Mag, cor, ele, ele_pos, idx, NNindex):
    pos = np.hstack((idx, NNindex[idx]))
    # ele = np.array([[1, 1, -1, 1, -1, 1], [1, -1, 1, 1, -1, 1]])
    # pos = np.array([13, 2, 5, 13, 2, 1])
    tmp1, tmp2 = np.unique(pos, return_counts=True)
    if len(tmp1)<len(pos):
        tmp3 = tmp1[np.where(tmp2 > 1)[0]]
        for ii, ind in enumerate(tmp3):
            if ii==0:
                tmp4 = ele[:, np.where(pos == ind)[0]]
                tmp5 = np.abs(np.sum(tmp4, axis=1))==len(tmp4[0])
            else:
                tmp4 = ele[:, np.where(pos == ind)[0]]
                tmp4 = np.abs(np.sum(tmp4, axis=1))==len(tmp4[0])
                tmp5 = (tmp4 & tmp5)
        ele = ele[tmp5]
        ele_pos = ele_pos[tmp5]

    puz = Mag[:, pos]
    puz1 = np.repeat(puz, len(ele), axis=0)
    ele1 = np.tile(ele, (len(puz), 1))
    r1 = (puz1 + ele1).all(1)

    tmp1 = np.where(r1 == True)[0]
    if tmp1.size == 0:
        Mag = np.array([])
        cor = np.array([])
    else:
        tmp2 = tmp1 // len(ele)
        tmp3 = np.mod(tmp1, len(ele))
        tmp4 = Counter(tmp2)
        # set_trace()
        Mag = np.repeat(Mag[list(tmp4.keys())], tuple(tmp4.values()), axis=0)
        cor = np.repeat(cor[list(tmp4.keys())], tuple(tmp4.values()), axis=0)

        tmp_ele = ele[tmp3]
        Mag[:, pos] = tmp_ele

        tmp_ele_pos = ele_pos[tmp3].astype(np.int32)
        # set_trace()
        cor[range(len(cor)), tmp_ele_pos] += 1
    return Mag, cor


def exceed_boundary(Mag1, cor1, NNindex, all_ele, ele_pos, order, itp,
                    savelog_path, cell_perms):
    ind1 = 20
    ind2 = len(Mag1) // ind1  #
    # cell_perms = np.loadtxt(read_perms_path).astype(np.int16) - 1


    Mag = np.array([])
    cor = np.array([])
    for jj in range(ind1):
        # set_trace()
        with open(savelog_path, 'a') as f:
            f.write('\n' + str(jj) + '/' + str(ind1) + ': ')
        if jj != ind1-1:
            tmp_M = Mag1[jj*ind2:(jj+1)*ind2]
            tmp_cor = cor1[jj*ind2:(jj+1)*ind2]
            for kk, idx in enumerate(order[itp+1:]):
                tmp_M, tmp_cor = puzzle(tmp_M, tmp_cor, all_ele, ele_pos, idx, NNindex)
                with open(savelog_path, 'a') as f:
                    f.write(str(itp + kk + 1) + ' ')
                if tmp_M.size==0:
                    with open(savelog_path, 'a') as f:
                        f.write('\n' + 'when ind==' + str(jj) + ' there are no correspongding strucrt' + '\n')
                    continue
            if tmp_M.size !=0:
                for ii in range(len(tmp_M)):
                    tmp = sortrows(tmp_M[ii][cell_perms])[0]
                    tmp_M[ii] = np.copy(tmp)
                tmp_M, itp1 = np.unique(tmp_M, return_index=True, axis=0)
                tmp_cor = tmp_cor[itp1] / len(order)
                if Mag.size==0:
                    Mag = tmp_M
                    cor = tmp_cor
                else:
                    Mag = np.vstack((Mag, tmp_M))
                    cor = np.vstack((cor, tmp_cor))
        else:
            tmp_M = Mag1[jj*ind2:]
            tmp_cor = cor1[jj*ind2:]
            for kk, idx in enumerate(order[itp+1:]):
                tmp_M, tmp_cor = puzzle(tmp_M, tmp_cor, all_ele, ele_pos, idx, NNindex)
                with open(savelog_path, 'a') as f:
                    f.write(str(itp + kk + 1) + ' ')
                if tmp_M.size==0:
                    with open(savelog_path, 'a') as f:
                        f.write('\n' + 'when ind==' + str(jj) + ' there are no correspongding strucrt' + '\n')
                    continue
            if tmp_M.size !=0:
                for ii in range(len(tmp_M)):
                    tmp = sortrows(tmp_M[ii][cell_perms])[0]
                    tmp_M[ii] = np.copy(tmp)
                tmp_M, itp1 = np.unique(tmp_M, return_index=True, axis=0)
                tmp_cor = tmp_cor[itp1] / len(order)
                if Mag.size==0:
                    Mag = tmp_M
                    cor = tmp_cor
                else:
                    Mag = np.vstack((Mag, tmp_M))
                    cor = np.vstack((cor, tmp_cor))
    if Mag.size != 0:
        for ii in range(len(Mag)):
            tmp = sortrows(Mag[ii][cell_perms])[0]
            Mag[ii] = np.copy(tmp)
        Mag, itp1 = np.unique(Mag, return_index=True, axis=0)
        cor = cor[itp1]
    return Mag, cor


def prep_path():
    name = '144-test'
    primitive_path = 'primitive/POSCAR-' + name
    savelog_path = 'log/main_%s.log' % name
    write_seq_path = 'seq/seq_%s' % name
    save_Magcor_path = 'mag_cor/poscar_%s_Mag_cor.npz'

    return primitive_path, savelog_path, write_seq_path, save_Magcor_path


def main():
    data = io.loadmat('19atoms_boron_envir_notINatoms.mat')
    SSS = data['cluster'].astype(np.int8) * 2 - 1
    perms = data['perms_1'] - 1
    ele = SSS
    primitive_path, savelog_path, write_seq_path, save_Magcor_path = prep_path()

    t0 = time.time()

    atoms = read(primitive_path)
    cell_perms = get_perms(atoms, symprec=1e-3).astype(np.int16)
    # set_trace()
    # os.mkdir(write_seq_path)

    for kk, line in enumerate(ele):
        tmp = transform_ele2all_nonsymmetry(line, perms)
        if kk == 0:
            all_ele = tmp
            ele_pos = np.ones((1, len(tmp))) * findByRow(SSS, line)
        else:
            all_ele = np.vstack((all_ele, tmp))
            ele_pos = np.hstack((ele_pos, np.ones((1, len(tmp))) * findByRow(SSS, line)))
    ele_pos = ele_pos[0]

    atoms = read_vasp(primitive_path)
    bb = atoms.get_scaled_positions()
    basis = atoms.cell

    NNindex = NNfind_3(basis, bb)
    # set_trace()
    order = find_order(NNindex)

    Mag1 = np.zeros((1, len(order))).astype(np.int8)
    cor1 = np.zeros((1, len(SSS))).astype(np.int8)
    for ii,idx in enumerate(order):
        Mag1, cor1 = puzzle(Mag1, cor1, all_ele, ele_pos, idx, NNindex)
        # print(ii)
        with open(savelog_path, 'a') as f:
            f.write(str(ii) + '\n')
        if len(Mag1)>10**6:
            with open(savelog_path, 'a') as f:
                f.write('now start pruning and parallel' + '\n')
            Mag1, cor1 = exceed_boundary(Mag1, cor1, NNindex, all_ele, ele_pos, order, ii,
                                         savelog_path, cell_perms)
            break
        # set_trace()
        if Mag1.size == 0:
            print('no corresponding struct for %s')
            break
    # set_trace()
    if Mag1.size != 0:
        for ii in range(len(Mag1)):
            tmp = sortrows(Mag1[ii][cell_perms])[0]
            Mag1[ii] = np.copy(tmp)
        Mag1, itp1 = np.unique(Mag1, return_index=True, axis=0)
        cor1 = cor1[itp1] / len(order)
        # print(len(Mag1))
        with open(savelog_path, 'a') as f:
            f.write('\ntotal struct number:' + str(len(Mag1)))
        # set_trace()
        for tmp in Mag1:
            tmp1 = list(np.where(tmp==-1)[0])
            # set_trace()
            with open(write_seq_path, "a") as f:
                for ii in tmp1:
                    f.write(str(ii) + ' ')
                f.write('\n')
        np.savez(save_Magcor_path, Mag=Mag1, cor=cor1)
            # set_trace()
    t1 = time.time()
    with open(savelog_path, 'a') as f:
        f.write('\ntotal spend time:' + str(t1 - t0) + '\n')


if __name__ == '__main__':
    main()



