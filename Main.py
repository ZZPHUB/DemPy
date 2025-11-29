from SimDem import SimDem_t
from PtcDat import *
from DatType import *
from PtcGen import *

if __name__ == "__main__":


    Gemo = Gemo_t()

    GenCube(Gemo,1.0,0.1,Vec_t(-0.6,0.0,0.0),1,Material_t.Rubber)
    GenCube(Gemo,1.0,0.1,Vec_t(0.6,0.0,0.0),2,Material_t.Rubber)

    for i in range(Gemo.size()):
        if Gemo.PtcList[i].Type == 1 :
            Gemo.PtcList[i].Vel = Vec_t(2.0,0.0,0.0)
        else:
            Gemo.PtcList[i].Vel = Vec_t(-2.0,0.0,0.0)

    Sim = SimDem_t()
    Sim.TotalStep = 2000
    Sim.OutStep = 1

    Sim.AddGemo(Gemo)

    Sim.SetPtcR(0.05)
    Sim.CalDt()

    print(Sim.Dt)

    Sim.Dt = 0.0005

    Sim.Sim()