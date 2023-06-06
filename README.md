# Generate-2D-Boron-Structures-by-Puzzle


## 第一步提前准备扩胞后的超胞作为拼图的背景

可使用http://sagar.compphys.cn/sagar 网站上的 Specific Volume Supercell Generator 功能。  


## 获取拼图元素
在该例中使用的"19atoms_boron_envir_notINatoms.mat"是符合B规则，六角晶格19个原子的拼图元素，已经提前生成好。  

其中的"cluster"是拼图元素，"perms_1"是该碎片的对称操作。  

如果是要算自己的例子，碎片的编号顺序可以自己定义。  


## 根据碎片和背景超胞生成结构

```
python main_puzzle_Boron_plane_primitive.py
```  
功能：  
根据选择的拼图碎片和原胞背景，拼出所有可能的结构，其中任意的碎片都可以用任意次，且满足边界周期条件。



程序说明：  
可以直接从main.py函数开始看起，其中的几个关键的函数：  

prep_path()：准备好各种路径。  

transform_ele2all_nonsymmetry（）：原先的cluster碎片中每个碎片都是对称不等价的，但实际上拼图的时候碎片是有方向性的，所以在具体拼的时候需要拿回所有对称等价的碎片。  

NNfind_3(basis, bb)：生成这个poscar中每个原子的近邻关系NNindex。  

find_order（NNindex）：根据每个原子的近邻关系确定一个拼图的顺序，尽可能地减少可能的中间结构的分支。  

exceed_boundary（）：在拼的过程中如果中间结构太多了内存会不够用，所以设定当结构数len(Mag1)>10^6时，开始分类计算深度优先，时间换空间。  



参数说明：  
Mag1存储的是每个位置的元素信息atom_type。  

cor1存储的是使用的碎片信息。  




