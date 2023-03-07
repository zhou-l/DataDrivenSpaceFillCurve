# DataDrivenSpaceFillCurve
Research codes of "Data-Driven Space-Filling Curves" (https://arxiv.org/abs/2009.06309).</br>
Please cite the paper if you use the codes:</br>

**L. Zhou, C. R. Johnson and D. Weiskopf, "Data-Driven Space-Filling Curves," in _IEEE Transactions on Visualization and Computer Graphics_, doi: 10.1109/TVCG.2020.3030473.**

![Linearizations of a synthetic volume data of a sphere](/sfc_reprod_full.png "Linearizations (top) of a synthetic volume data of a sphere.")

# File name convention
Important: main functions are named as "figGenXXX.m". </br>
For example, </br>
**figGenSFC2DCase.m**---2D regular grids;</br> 
**figGenSFC3DCase.m**---3D regular grids; </br>
**figGenSFC2DQuadTree.m**---2D multiscale; </br>
**figGenSFC3DOctree.m**---3D multiscale. </br>

# Paper Figures Reproducibility
Paper reproducibility: </br>
Figures 6, 7, 8, 10, 12 of the paper can be readily reproduced with scripts named **paperFigX.m** where X is the number of the figure. </br>
**paperFig6.m**---Figure 6;</br>
**paperFig7.m**---Figure 7;</br>
**paperFig8.m**---Figure 8;</br>
**paperFig10.m**---Figure 10;</br>
**paperFig12.m**---Figure 12;</br>
**paperFigsAll.m** can generate results of all figures above. 
