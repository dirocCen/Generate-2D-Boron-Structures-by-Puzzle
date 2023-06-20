import numpy as np
from pdb import set_trace
import os
from scipy import io


def sortrows(a):
    return a[np.lexsort(np.rot90(a))]


def findByRow(mat, row, prec=1e-4):
    mat = np.round(mat, 4)
    row = np.round(row, 4)
    return np.where((mat == row).all(1))[0]


def read_cor():
    path0 = 'mag_cor'
    # res = np.unique(data, axis=1)
    cor = np.array([])
    for ii in range(91):
        path1 = os.path.join(path0, 'poscar_%s' % ii)
        if os.path.exists(path1):
            for jj,ind1 in enumerate(os.listdir(path1)):
                path2 = os.path.join(path1, ind1)
                # set_trace()
                if cor.size==0:
                    cor = np.load(path2)['cor']
                else:
                    cor = np.vstack((cor, np.load(path2)['cor']))
                    # set_trace()
                    cor = np.unique(cor, axis=0)
    io.savemat('B_1_90_cor.mat', {'cor': cor})  # 保存成mat格式
    set_trace()


def get_poscar_from_cor():
    M_convex = io.loadmat('B_1_80_M_convex.mat')['M_convex']
    path0 = 'mag_cor'
    pos = []
    for line in M_convex:
        itp = np.array([])
        for ii in range(81):
            path1 = os.path.join(path0, 'poscar_%s' % ii)
            if os.path.exists(path1):
                for jj, ind1 in enumerate(os.listdir(path1)):
                    path2 = os.path.join(path1, ind1)
                    cor = np.load(path2)['cor']
                    itp = findByRow(cor, line)
                    # set_trace()
                    if itp.size!=0:
                        pos.append({path2: itp})
                        # set_trace()
                        break
            if itp.size != 0:
                break
                    # set_trace()
    set_trace()




if __name__ == '__main__':
    get_poscar_from_cor()
