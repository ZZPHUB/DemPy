from PtcDat import *
from DatType import *


def GenCube(Gemo:Gemo_t,Len:float,Dx:float,CPos:Vec_t,Type:int,Mat:Material_t):
    for i in range(int(Len/Dx)):
        for j in range(int(Len/Dx)):
            for k in range(int(Len/Dx)):
                Ptc = Ptc_t()
                Ptc.SetMaterial(Mat)
                Ptc.Pos = CPos - Vec_t(Len/2.0,Len/2.0,Len/2.0) + Vec_t(i*Dx,j*Dx,k*Dx)
                Ptc.Type = Type
                Ptc.Id = Gemo.size() + i*int(Len/Dx)*int(Len/Dx) + j*int(Len/Dx) + k
                Gemo.AddPtc(Ptc)