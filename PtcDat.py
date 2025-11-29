from DatType import Vec_t, Mat_t
from copy import copy,deepcopy
from enum import Enum

class OverLapT_t:
    def __init__(self,Id:int,OverLap:Vec_t):
        self.Id = Id
        self.OverLap = OverLap
    
class Material_t(Enum):
    Ice = 1
    Aluminum = 2
    Rubber = 3
    Rigid = 4

class Ptc_t:
    def __init__(self):
        self.Pos:Vec_t = Vec_t()
        self.Vel:Vec_t = Vec_t()
        self.Omega:Vec_t = Vec_t()
        self.Acc:Vec_t = Vec_t()
        self.Alapha:Vec_t = Vec_t()
        #self.Raidus:float = 0.0
        #self.Mass:float = 0.0
        self.Id:int = 0
        self.OverLapList:list[OverLapT_t] = []
        self.Type:int = 0 
    
    def SetMaterial(self,Mat:Material_t):
        match Mat:
            case Material_t.Ice:
                self.E:float = 10.0e0
                self.G:float = 3.5e9 
                self.Poisson:float = 0.3 
                self.Rho:float = 910.0
            case Material_t.Aluminum:
                self.E:float = 17.0e9
                self.G:float = 7.0e9
                self.Poisson:float = 0.42
                self.Rho:float = 2700.0
            case Material_t.Rubber:
                self.E:float = 0.00784e9
                self.G:float = 2.9e9 
                self.Poisson:float = 0.47 
                self.Rho:float = 920.0
            case Material_t.Rigid:
                self.E:float = 0.0
                self.G:float = 0.0
                self.Poisson:float = 0.0
                self.Rho:float = 0.0


class Gemo_t:
    def __init__(self):
        self.PtcList:list[Ptc_t] = []
    
    def __size__(self):
        return len(self.PtcList)
    
    def size(self):
        return self.__size__()
    
    def AddPtc(self,Ptc:Ptc_t):
        self.PtcList.append(deepcopy(Ptc))
    

if __name__ == "__main__":
    Ptc = Ptc_t(Vec_t(0,0,0),Vec_t(0,0,0),Vec_t(0,0,0),Vec_t(0,0,0),Vec_t(0,0,0),1.0,1.0,0)
    Gemo = Gemo_t()
    Gemo.AddPtc(Ptc)
    print(Gemo.PtcList[0].Acc)
    Ptc.Acc = Vec_t(0,0,1)
    print(Gemo.PtcList[0].Acc)
    Gemo.AddPtc(Ptc)
    print(Gemo.PtcList[1].Acc)
    print(Gemo.size())
