%convert the octree level to the SFC level
%The octree always use 0 to denote coarest level of a number of 1*1*1
%block, whereas the sfcLvl starts from 1-finest to sfcMaxLvl-coarsest leaf
%node.
function sfcLvl = convertOctLvlToSFCLvl(octLvl, octMaxLvl)
   sfcLvl = octMaxLvl - octLvl + 1;
end