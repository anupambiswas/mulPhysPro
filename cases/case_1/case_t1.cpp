#include<iostream>
#include<cmath>
#include<cstdlib>
#include<fstream>
#include"../../headers/flowDomain_2.h"

using namespace std;

const int HOR=0,VER=1,STA=0,END=1;
const int WALL=0,VINLET=1,POUTLET=2,IFACE=3;

class fDomain:public flowDomain_2
{
public:
  fDomain(std::string name);
  friend void solvePSync(fDomain **domList,int numDom);
  void vmStore();
};

fDomain::fDomain(std::string name):flowDomain_2(name)
{
}

void solvePSync(fDomain **domList,int numDom)
{
  int i;
  double maxDif;

  do
    {
      maxDif=-1.0e40;
      for(i=0;i<numDom;i++)
	domList[i]->getPCorrectionFromNeibs();
      for(i=0;i<numDom;i++)
	domList[i]->solvePSingle(0.666,maxDif);
    }
  while(maxDif>domList[0]->iTolP);
}

void fDomain::vmStore()
{
  int i,j;
  double vm;
  std::ofstream F((domainName+"vm").c_str(),ios::out);
  for(i=2;i<nxp2[0];i++)
    for(j=2;j<nxp2[1];j++)
      {
	vm=pow(u[i][j]*u[i][j]+v[i][j]*v[i][j],0.5);
	F<<(i-1.5)*dx[0]+xOri<<" "<<(j-1.5)*dx[1]+yOri<<" "<<vm<<std::endl;
      }
  F.close();
}

int main()
{
  remove("res_u");
  remove("res_v");
  remove("res_p");

  int NH,NR1,NR2,Nfw,NG,Nfh;
  double L;
  double dx,dt;

  NH=60;
  NR1=60;
  NR2=60;
  Nfw=10;
  NG=10;
  Nfh=24;
  dx=1.0/NH;
  dt=1.0e-3;
  L=(NR1+NR2+3*NG+4*Nfw)*dx;

  int NfinH=NH-Nfh;

  fDomain R1("R1");
  R1.setGridProp(NR1,NH,dx,dx,dt);
  R1.setFlowProp(1.0,0.71,0);
  R1.setIterTol(1.0e-8,1.0e-8,1.0e-8);
  R1.setSoluTol(1.0e-8,1.0e-8,1.0e-8);
  R1.setRelPar(0.3,0.3,0.7);
  R1.setOrigin(0.0,0.0);

  fDomain R2("R2");
  R2.setGridProp(NR2,NH,dx,dx,dt);
  R2.setFlowProp(1.0,0.71,0);
  R2.setIterTol(1.0e-8,1.0e-8,1.0e-8);
  R2.setSoluTol(1.0e-8,1.0e-8,1.0e-8);
  R2.setRelPar(0.3,0.3,0.7);
  R2.setOrigin(L-dx*NR2,0.0);

  fDomain F1("F1");
  F1.setGridProp(Nfw,Nfh,dx,dx,dt);
  F1.setFlowProp(1.0,0.71,0);
  F1.setIterTol(1.0e-8,1.0e-8,1.0e-8);
  F1.setSoluTol(1.0e-8,1.0e-8,1.0e-8);
  F1.setRelPar(0.3,0.3,0.7);
  F1.setOrigin(NR1*dx,NfinH*dx);

  fDomain F2("F2");
  F2.setGridProp(Nfw,Nfh,dx,dx,dt);
  F2.setFlowProp(1.0,0.71,0);
  F2.setIterTol(1.0e-8,1.0e-8,1.0e-8);
  F2.setSoluTol(1.0e-8,1.0e-8,1.0e-8);
  F2.setRelPar(0.3,0.3,0.7);
  F2.setOrigin((NR1+Nfw+NG)*dx,0.0);

  fDomain F3("F3");
  F3.setGridProp(Nfw,Nfh,dx,dx,dt);
  F3.setFlowProp(1.0,0.71,0);
  F3.setIterTol(1.0e-8,1.0e-8,1.0e-8);
  F3.setSoluTol(1.0e-8,1.0e-8,1.0e-8);
  F3.setRelPar(0.3,0.3,0.7);
  F3.setOrigin((NR1+2*Nfw+2*NG)*dx,NfinH*dx);

  fDomain F4("F4");
  F4.setGridProp(Nfw,Nfh,dx,dx,dt);
  F4.setFlowProp(1.0,0.71,0);
  F4.setIterTol(1.0e-8,1.0e-8,1.0e-8);
  F4.setSoluTol(1.0e-8,1.0e-8,1.0e-8);
  F4.setRelPar(0.3,0.3,0.7);
  F4.setOrigin((NR1+3*Nfw+3*NG)*dx,0.0);

  fDomain G12("G12");
  G12.setGridProp(NG,NH,dx,dx,dt);
  G12.setFlowProp(1.0,0.71,0);
  G12.setIterTol(1.0e-8,1.0e-8,1.0e-8);
  G12.setSoluTol(1.0e-8,1.0e-8,1.0e-8);
  G12.setRelPar(0.3,0.3,0.7);
  G12.setOrigin((NR1+Nfw)*dx,0.0);

  fDomain G23("G23");
  G23.setGridProp(NG,NH,dx,dx,dt);
  G23.setFlowProp(1.0,0.71,0);
  G23.setIterTol(1.0e-8,1.0e-8,1.0e-8);
  G23.setSoluTol(1.0e-8,1.0e-8,1.0e-8);
  G23.setRelPar(0.3,0.3,0.7);
  G23.setOrigin((NR1+2*Nfw+NG)*dx,0.0);

  fDomain G34("G34");
  G34.setGridProp(NG,NH,dx,dx,dt);
  G34.setFlowProp(1.0,0.71,0);
  G34.setIterTol(1.0e-8,1.0e-8,1.0e-8);
  G34.setSoluTol(1.0e-8,1.0e-8,1.0e-8);
  G34.setRelPar(0.3,0.3,0.7);
  G34.setOrigin((NR1+3*Nfw+2*NG)*dx,0.0);

  fDomain **domList=new fDomain* [9];
  domList[0]=&R1;
  domList[1]=&F1;
  domList[2]=&G12;
  domList[3]=&F2;
  domList[4]=&G23;
  domList[5]=&F3;
  domList[6]=&G34;
  domList[7]=&F4;
  domList[8]=&R2;

  int i,j;
  double **velW,**velI;
  double *pres;

  velW=new double* [NR1]; // Choose biggest size
  for(i=0;i<NR1;i++)
    {
      velW[i]=new double [2];
      velW[i][0]=velW[i][1]=0.0;
    }

  velI=new double* [NH];
  for(i=0;i<NH;i++)
    {
      velI[i]=new double [2];
      velI[i][0]=1.0;
      velI[i][1]=0.0;
    }

  pres=new double [NH];
  for(i=0;i<NH;i++)
    pres[i]=0.0;


  R1.defineBoundary(HOR,STA,2,NR1+2,WALL,velW);
  R1.defineBoundary(HOR,END,2,NR1+2,WALL,velW);
  R1.defineBoundary(VER,STA,2,NH+2,VINLET,velI);
  R1.defineBoundary(VER,END,2,NfinH+2,WALL,velW);
  R1.defineBoundary(VER,END,NfinH+2,NH+2,IFACE,NULL,NULL,NULL,0.0,&F1,2);

  F1.defineBoundary(HOR,STA,2,Nfw+2,WALL,velW);
  F1.defineBoundary(HOR,END,2,Nfw+2,WALL,velW);
  F1.defineBoundary(VER,STA,2,Nfh+2,IFACE,NULL,NULL,NULL,0.0,&R1,NfinH+2);
  F1.defineBoundary(VER,END,2,Nfh+2,IFACE,NULL,NULL,NULL,0.0,&G12,NfinH+2);

  G12.defineBoundary(HOR,STA,2,NG+2,WALL,velW);
  G12.defineBoundary(HOR,END,2,NG+2,WALL,velW);
  G12.defineBoundary(VER,STA,2,NfinH+2,WALL,velW);
  G12.defineBoundary(VER,STA,NfinH+2,NH+2,IFACE,NULL,NULL,NULL,0.0,&F1,2);
  G12.defineBoundary(VER,END,2,Nfh+2,IFACE,NULL,NULL,NULL,0.0,&F2,2);
  G12.defineBoundary(VER,END,Nfh+2,NH+2,WALL,velW);

  F2.defineBoundary(HOR,STA,2,Nfw+2,WALL,velW);
  F2.defineBoundary(HOR,END,2,Nfw+2,WALL,velW);
  F2.defineBoundary(VER,STA,2,Nfh+2,IFACE,NULL,NULL,NULL,0.0,&G12,2);
  F2.defineBoundary(VER,END,2,Nfh+2,IFACE,NULL,NULL,NULL,0.0,&G23,2);

  G23.defineBoundary(HOR,STA,2,NG+2,WALL,velW);
  G23.defineBoundary(HOR,END,2,NG+2,WALL,velW);
  G23.defineBoundary(VER,STA,2,Nfh+2,IFACE,NULL,NULL,NULL,0.0,&F2,2);
  G23.defineBoundary(VER,STA,Nfh+2,NH+2,WALL,velW);
  G23.defineBoundary(VER,END,2,NfinH+2,WALL,velW);
  G23.defineBoundary(VER,END,NfinH+2,NH+2,IFACE,NULL,NULL,NULL,0.0,&F3,2);

  F3.defineBoundary(HOR,STA,2,Nfw+2,WALL,velW);
  F3.defineBoundary(HOR,END,2,Nfw+2,WALL,velW);
  F3.defineBoundary(VER,STA,2,Nfh+2,IFACE,NULL,NULL,NULL,0.0,&G23,NfinH+2);
  F3.defineBoundary(VER,END,2,Nfh+2,IFACE,NULL,NULL,NULL,0.0,&G34,NfinH+2);

  G34.defineBoundary(HOR,STA,2,NG+2,WALL,velW);
  G34.defineBoundary(HOR,END,2,NG+2,WALL,velW);
  G34.defineBoundary(VER,STA,2,NfinH+2,WALL,velW);
  G34.defineBoundary(VER,STA,NfinH+2,NH+2,IFACE,NULL,NULL,NULL,0.0,&F3,2);
  G34.defineBoundary(VER,END,2,Nfh+2,IFACE,NULL,NULL,NULL,0.0,&F4,2);
  G34.defineBoundary(VER,END,Nfh+2,NH+2,WALL,velW);

  F4.defineBoundary(HOR,STA,2,Nfw+2,WALL,velW);
  F4.defineBoundary(HOR,END,2,Nfw+2,WALL,velW);
  F4.defineBoundary(VER,STA,2,Nfh+2,IFACE,NULL,NULL,NULL,0.0,&G34,2);
  F4.defineBoundary(VER,END,2,Nfh+2,IFACE,NULL,NULL,NULL,0.0,&R2,2);

  R2.defineBoundary(HOR,STA,2,NR2+2,WALL,velW);
  R2.defineBoundary(HOR,END,2,NR2+2,WALL,velW);
  R2.defineBoundary(VER,STA,2,Nfh+2,IFACE,NULL,NULL,NULL,0.0,&F4,2);
  R2.defineBoundary(VER,STA,Nfh+2,NH+2,WALL,velW);
  R2.defineBoundary(VER,END,2,NH+2,POUTLET,NULL,pres);

  int iter=0;
  bool resid;

  while(1)
    {
      R1.resetBoundary();
      R1.computeFlux();
      R1.computeCoefficients();

      F1.resetBoundary();
      F1.computeFlux();
      F1.computeCoefficients();

      G12.resetBoundary();
      G12.computeFlux();
      G12.computeCoefficients();

      F2.resetBoundary();
      F2.computeFlux();
      F2.computeCoefficients();

      G23.resetBoundary();
      G23.computeFlux();
      G23.computeCoefficients();

      F3.resetBoundary();
      F3.computeFlux();
      F3.computeCoefficients();

      G34.resetBoundary();
      G34.computeFlux();
      G34.computeCoefficients();

      F4.resetBoundary();
      F4.computeFlux();
      F4.computeCoefficients();

      R2.resetBoundary();
      R2.computeFlux();
      R2.computeCoefficients();

      R1.resetSolverCoef(0);
      R1.solveGS(0,0.666);
      R1.solveGS(1,0.666);

      F1.resetSolverCoef(0);
      F1.solveGS(0,0.666);
      F1.solveGS(1,0.666);

      G12.resetSolverCoef(0);
      G12.solveGS(0,0.666);
      G12.solveGS(1,0.666);

      F2.resetSolverCoef(0);
      F2.solveGS(0,0.666);
      F2.solveGS(1,0.666);

      G23.resetSolverCoef(0);
      G23.solveGS(0,0.666);
      G23.solveGS(1,0.666);

      F3.resetSolverCoef(0);
      F3.solveGS(0,0.666);
      F3.solveGS(1,0.666);

      G34.resetSolverCoef(0);
      G34.solveGS(0,0.666);
      G34.solveGS(1,0.666);

      F4.resetSolverCoef(0);
      F4.solveGS(0,0.666);
      F4.solveGS(1,0.666);

      R2.resetSolverCoef(0);
      R2.solveGS(0,0.666);
      R2.solveGS(1,0.666);

      R1.manageInterface();
      F1.manageInterface();
      G12.manageInterface();
      F2.manageInterface();
      G23.manageInterface();
      F3.manageInterface();
      G34.manageInterface();
      F4.manageInterface();
      R2.manageInterface();

      R1.computePCoefReq();
      F1.computePCoefReq();
      G12.computePCoefReq();
      F2.computePCoefReq();
      G23.computePCoefReq();
      F3.computePCoefReq();
      G34.computePCoefReq();
      F4.computePCoefReq();
      R2.computePCoefReq();

      R1.computePCoefficients();
      F1.computePCoefficients();
      G12.computePCoefficients();
      F2.computePCoefficients();
      G23.computePCoefficients();
      F3.computePCoefficients();
      G34.computePCoefficients();
      F4.computePCoefficients();
      R2.computePCoefficients();

      R1.resetSolverCoef(2);
      F1.resetSolverCoef(2);
      G12.resetSolverCoef(2);
      F2.resetSolverCoef(2);
      G23.resetSolverCoef(2);
      F3.resetSolverCoef(2);
      G34.resetSolverCoef(2);
      F4.resetSolverCoef(2);
      R2.resetSolverCoef(2);

      solvePSync(domList,9);//return(0);

      R1.renewVariables();
      F1.renewVariables();
      G12.renewVariables();
      F2.renewVariables();
      G23.renewVariables();
      F3.renewVariables();
      G34.renewVariables();
      F4.renewVariables();
      R2.renewVariables();

      R1.swapVarMems();
      F1.swapVarMems();
      G12.swapVarMems();
      F2.swapVarMems();
      G23.swapVarMems();
      F3.swapVarMems();
      G34.swapVarMems();
      F4.swapVarMems();
      R2.swapVarMems();

      resid=R1.computeResiduals() && F1.computeResiduals() && G12.computeResiduals() && F2.computeResiduals();
      resid=resid && G23.computeResiduals() && F3.computeResiduals() && G34.computeResiduals();
      resid=resid && F4.computeResiduals() && R2.computeResiduals();
      if(resid)break;

      R1.storeResiduals(iter);
      if(iter%200==0){cout<<iter<<" ";R1.displayResiduals();}
      if(iter%500==0)
	{
	  R1.vmStore();
	  F1.vmStore();
	  G12.vmStore();
	  F2.vmStore();
	  G23.vmStore();
	  F3.vmStore();
	  G34.vmStore();
	  F4.vmStore();
	  R2.vmStore();

	  R1.dataStore(0);R1.dataStore(1);R1.dataStore(2);
	  F1.dataStore(0);F1.dataStore(1);F1.dataStore(2);
	  G12.dataStore(0);G12.dataStore(1);G12.dataStore(2);
	  F2.dataStore(0);F2.dataStore(1);F2.dataStore(2);
	  G23.dataStore(0);G23.dataStore(1);G23.dataStore(2);
	  F3.dataStore(0);F3.dataStore(1);F3.dataStore(2);
	  G34.dataStore(0);G34.dataStore(1);G34.dataStore(2);
	  F4.dataStore(0);F4.dataStore(1);F4.dataStore(2);
	  R2.dataStore(0);R2.dataStore(1);R2.dataStore(2);
	}
      iter++;//return(0);///////////////////////////////
    }

  R1.resetBoundary();
  F1.resetBoundary();
  G12.resetBoundary();
  F2.resetBoundary();
  G23.resetBoundary();
  F3.resetBoundary();
  G34.resetBoundary();
  F4.resetBoundary();
  R2.resetBoundary();

  return(0);
}
