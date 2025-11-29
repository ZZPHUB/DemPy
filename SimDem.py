from PtcDat import *
import vtk as Vtk
from math import *
class SimDem_t:
    def __init__(self):
        self.CurrentStep = int(0)
        self.TotalStep = int(0)
        self.OutStep = int(0)
        #self.Dt = float(0.0)

        self.Gravity:Vec_t = Vec_t(0.0,0.0,0.0)

        self.Restitution:float = 0.7
        self.Cn = - log(self.Restitution)/sqrt(pi**2 + log(self.Restitution)**2)
        self.Friction:float = 0.3
    
    def SetPtcR(self,R:float):
        self.Radius:float = R
    
    def CalDt(self):
        EMax = 0.0
        Rho = 0.0
        for i in range(self.Gemo.size()):
            PtcI = self.Gemo.PtcList[i]
            EMax = max(PtcI.E,EMax)
            Rho = max(PtcI.Rho,Rho)
        self.Dt = 0.5 * pi/50.0 * sqrt(sqrt(4.0/3.0*EMax*sqrt(self.Radius))/(4.0/3.0 * pi * self.Radius**3 * Rho)) 
    
    def AddGemo(self,Gemo:Gemo_t):
        self.Gemo:Gemo_t = Gemo
    
    def CalContact(self):
        for i in range(self.Gemo.size()):
            PtcI = self.Gemo.PtcList[i]
            PtcI.Acc = Vec_t()
            PtcI.Alapha = Vec_t()

        for i in range(self.Gemo.size()):
            PtcI = self.Gemo.PtcList[i]
            for j in range(i+1,self.Gemo.size()):
                PtcJ = self.Gemo.PtcList[j]

                VecIJ = PtcJ.Pos - PtcI.Pos
                Dx = sqrt(VecIJ * VecIJ)

                if(abs(Dx) > 1e-6):
                    VecIJ = VecIJ / Dx

                NormOL = self.Radius*2.0 - Dx

                if (NormOL < 0.0 and PtcI.Type != PtcJ.Type):
                    continue
                else:
                    if(NormOL < 0.0 and PtcI.Type == PtcJ.Type):
                        if NormOL < -self.Radius/2.0:
                            continue
                        else:
                            VecIJ = VecIJ * (-1.0)
                            NormOL = -1.0 * NormOL
                    
                    VelIJNorm = VecIJ * ((PtcJ.Vel - PtcI.Vel) * VecIJ)

                    EStar = 1.0 /((1.0 - PtcI.Poisson**2)/PtcI.E + (1.0 - PtcJ.Poisson**2)/PtcJ.E)
                    RStar = self.Radius/2.0 
                    Mass = 4.0/3.0 * self.Radius**3 * (PtcI.Rho * PtcJ.Rho) / (PtcI.Rho + PtcJ.Rho)

                    KNorm = 4.0/3.0*EStar*sqrt(RStar)
                    Gamma = 1e-5 * sqrt(6.0 * Mass * EStar * sqrt(RStar))

                    ForceNorm = VecIJ * KNorm*sqrt(NormOL**3)
                    ForceDampNorm = VelIJNorm * (-Gamma*NormOL**(0.25))

                    #ForceNorm = VecIJ * (4.0/3.0)*EStar*sqrt(RStar)*sqrt(NormOL**3)
                    #ForceDampNorm = VelIJNorm * (-self.Cn*sqrt(8.0*Mass*EStar*sqrt(RStar*NormOL)))

                    #if(PtcI.Type == PtcJ.Type):
                    #    ForceDampNorm = ForceDampNorm * (-1.0)

                    Force = ForceNorm + ForceDampNorm

                    PtcI.Acc = PtcI.Acc - Force/Mass
                    PtcJ.Acc = PtcJ.Acc + Force/Mass

    def CalUpdate(self):
        for i in range(self.Gemo.size()):
            Ptc = self.Gemo.PtcList[i]
            Ptc.Vel = Ptc.Vel + Ptc.Acc * self.Dt
            Ptc.Pos = Ptc.Pos + Ptc.Vel * self.Dt

    def SimStep(self):
        self.CalContact()
        self.CalUpdate()
    def Sim(self):
        for i in range(self.TotalStep):
            self.SimStep()
            self.CurrentStep = i 
            if self.CurrentStep % self.OutStep == 0:
                self.Output()

    def Output(self):
        Ptc = Vtk.vtkPoints()
        Ptc.SetNumberOfPoints(self.Gemo.size())
        Poly = Vtk.vtkPolyData()
        Poly.SetPoints(Ptc)

        Type = Vtk.vtkIntArray()
        Type.SetNumberOfComponents(1)
        Type.SetName("Type")
        Type.SetNumberOfTuples(self.Gemo.size())

        Vel = Vtk.vtkFloatArray()
        Vel.SetNumberOfComponents(3)
        Vel.SetName("Vel")
        Vel.SetNumberOfTuples(self.Gemo.size())

        Omega = Vtk.vtkFloatArray()
        Omega.SetNumberOfComponents(3)
        Omega.SetName("Omega")
        Omega.SetNumberOfTuples(self.Gemo.size())

        Acc = Vtk.vtkFloatArray()
        Acc.SetNumberOfComponents(3)
        Acc.SetName("Acc")
        Acc.SetNumberOfTuples(self.Gemo.size())

        Alapha = Vtk.vtkFloatArray()
        Alapha.SetNumberOfComponents(3)
        Alapha.SetName("Alapha")
        Alapha.SetNumberOfTuples(self.Gemo.size())

        Id = Vtk.vtkIntArray()
        Id.SetNumberOfComponents(1)
        Id.SetName("Id")
        Id.SetNumberOfTuples(self.Gemo.size())

        for i,P in enumerate(self.Gemo.PtcList):
            Ptc.SetPoint(i,P.Pos.x,P.Pos.y,P.Pos.z)
            Type.SetTuple1(i,P.Type)
            Vel.SetTuple3(i,P.Vel.x,P.Vel.y,P.Vel.z)
            Omega.SetTuple3(i,P.Omega.x,P.Omega.y,P.Omega.z)
            Acc.SetTuple3(i,P.Acc.x,P.Acc.y,P.Acc.z)
            Alapha.SetTuple3(i,P.Alapha.x,P.Alapha.y,P.Alapha.z)
            Id.SetTuple1(i,P.Id)
        Ptc.Modified()
        Poly.GetPointData().AddArray(Type)
        Poly.GetPointData().AddArray(Vel)
        Poly.GetPointData().AddArray(Omega)
        Poly.GetPointData().AddArray(Acc)
        Poly.GetPointData().AddArray(Alapha)
        Poly.GetPointData().AddArray(Id)

        #Poly.GetPointData().SetScalars(Id)
        Poly.GetPointData().SetScalars(Type)
        Poly.GetPointData().SetVectors(Vel)
        #Poly.GetPointData().SetVectors(Omega)
        #Poly.GetPointData().SetVectors(Acc)
        #Poly.GetPointData().SetVectors(Alapha)

        FName = f"./OutFile/PtcOut{self.CurrentStep}.vtp"
        Wri = Vtk.vtkXMLPolyDataWriter()
        Wri.SetFileName(FName)
        Wri.SetDataModeToAscii()
        Wri.SetInputData(Poly)
        if Wri.Write():
                print(f"Write {FName} success")
        else:
            print(f"Write {FName} failed")

        
