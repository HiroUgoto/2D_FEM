# Input File Format - mesh.in
---

Input file `mesh.in` consists of three definition blocks, namely **header block**, **node block**, **element block**, and **material block**. The detail formats are described as follows.

---
## Header block

Header block defines a number of nodes and elements, and degrees of freedom of the FEM model. The header block must be writtern in a single line.

```
nnode nelem nmaterial dof
```

**nnode** : Number of nodes  
**nelem** : Number of elements  
**nmaterial** : Number of materials  
**dof** : Degrees of freedom
- 2D-SH (anti-plane) problem : dof = 1
- 2D-PSV (in-plane) problem : dof = 2
- 2D-SH,PSV problem : dof = 3


### Example

```
121 100 2 2
```
FEM model consists 121 nodes and 100 elements. 2 materials are given.
The problem targets 2D in-plane deformation.

---
## Node block

Node block defines locations of nodes and its degree of freedom (constrains).

```
id x y dof0 [dof1] [dof2]
```

**id** : Node ID. Prefer to define it by sequential order.  
**x** : x coordinate [m]  
**y** : y coordinate [m]  
**dof0** : 1st degree of freedom [fix: 0, no fix: 1]  
**dof1** : 2nd degree of freedom (if dof >= 2)  
**dof2** : 3rd degree of freedom (if dof =3)

### Example

```
0 0.0 0.0 0 0
1 1.0 0.0 1 1
2 2.0 0.0 0 0
```
3 Nodes are defined at (x,y) = (0.0,0.0), (1.0,0.0), (2.0,0.0). Node ID 0 and 2 are fixed both x and y components.  

---
## Element block

Element block defines element types, material constants, and belonging nodes.

```
id style material_id node_id
```

**id** : Element ID. Prefer to define it by sequential order.  
**style** : Element style as listed below.
- 2d4solid : 4-node isoparametric solid element (2D)
- 2d8solid : 8-node serendipity solid element (2D)
- 2d9solid : 9-node isoparametric solid element (2D)
- 1d2line : 2-node line element (1D)
- 1d3line : 3-node line element (1D)
- 1d2input : 2-node line input boundary element
- 1d3input : 3-node line input boundary element  

**material_id** : Material ID defined in material block.  
**node_id** : Node ID list in the order corresponding to the element style.

### Example

```
0 2d4solid 0 0 1 3 2
```

2D 4-node isoparametric element is defined.
Material parameters are given in ID 0.
The element consists of node ID 0, 1, 3, 2 in the order.

```
1 1d2input 1 0 1
```

Input wave boundary is defined on the boundary consisting of node ID 0 and 1.

---
## Material block

Material block defines material types and its physical parameters.

```
id style param
```

**id** : Material ID. Prefer to define it by sequential order.  
**style** : Material style as listed below.  
**param** : List of physical parameters. The order is described as below.  
- vs_vp_rho : Linear elastic material. Parameters are given in the order of  
  + S-wave velocity [m/s], P-wave velocity [m/s], and density [kg/m^3]
- nu_vp_rho : Linear elastic material. Parameters are given in the order of
  + Poisson's ratio, P-wave velocity [m/s], and density [kg/m^3]

### Example

```
0 vs_vp_rho 250. 1500. 1750.
```

Linear elastic material is defined. The parameters are defined in the order of S-wave velocity (250m/s), P-wave velocity (1500m/s), and density (1750kg/m^3).
