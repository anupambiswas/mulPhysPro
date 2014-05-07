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

int main()
{
  int chanWall=100,exteWall=80,inletSize=60,outletSize=80,ifcS;
  double dx=1.0/inletSize;

  //double *a=new double [0];return(0);

  ifcS=outletSize-inletSize;

  fDomain chan("Channel");
  chan.setGridProp(chanWall,inletSize,dx,dx,0);
  chan.setFlowProp(1.0,0);
  chan.setIterTol(1.0e-8,1.0e-8,1.0e-8);
  chan.setSoluTol(1.0e-8,1.0e-8,1.0e-8);
  chan.setRelPar(0.3,0.3,0.7);
  chan.setOrigin(0.0,0.0);

  fDomain exte("Extension");
  exte.setGridProp(exteWall,outletSize,dx,dx,0);
  exte.setFlowProp(1.0,0);
  exte.setIterTol(1.0e-8,1.0e-8,1.0e-8);
  exte.setSoluTol(1.0e-8,1.0e-8,1.0e-8);
  exte.setRelPar(0.3,0.3,0.7);
  exte.setOrigin(chanWall*dx,-ifcS*dx);

  fDomain **domList=new fDomain* [2];
  domList[0]=&chan;
  domList[1]=&exte;

  int i,j;
  double **velW,**velI;
  double *pres;

  velW=new double* [chanWall];
  for(i=0;i<chanWall;i++)
    {
      velW[i]=new double [2];
      velW[i][0]=velW[i][1]=0.0;
    }

  velI=new double* [inletSize];
  for(i=0;i<inletSize;i++)
    {
      velI[i]=new double [2];
      velI[i][0]=1.0;
      velI[i][1]=0.0;
    }

  pres=new double [outletSize];
  for(i=0;i<outletSize;i++)
    pres[i]=0.0;


  chan.defineBoundary(HOR,STA,2,chanWall+2,WALL,velW);
  chan.defineBoundary(HOR,END,2,chanWall+2,WALL,velW);
  chan.defineBoundary(VER,STA,2,inletSize+2,VINLET,velI);
  chan.defineBoundary(VER,END,2,inletSize+2,IFACE,NULL,NULL,NULL,0,&exte,ifcS+2);

  exte.defineBoundary(HOR,STA,2,exteWall+2,WALL,velW);
  exte.defineBoundary(HOR,END,2,exteWall+2,WALL,velW);
  exte.defineBoundary(VER,STA,2,ifcS+2,WALL,velW);
  exte.defineBoundary(VER,STA,ifcS+2,outletSize+2,IFACE,NULL,NULL,NULL,0,&chan,2);
  exte.defineBoundary(VER,END,2,outletSize+2,POUTLET,NULL,pres);

  int iter=0;
  while(1)
    {
      chan.resetBoundary();
      chan.computeFlux();
      chan.computeCoefficients();

      exte.resetBoundary();
      exte.computeFlux();
      exte.computeCoefficients();

      chan.resetSolverCoef(0);
      chan.solveGS(0,0.666);
      chan.solveGS(1,0.666);

      exte.resetSolverCoef(0);
      exte.solveGS(0,0.666);
      exte.solveGS(1,0.666);

      chan.manageInterface();
      exte.manageInterface();

      chan.computePCoefReq();
      exte.computePCoefReq();

      chan.computePCoefficients();
      chan.resetSolverCoef(2);

      exte.computePCoefficients();
      exte.resetSolverCoef(2);

      solvePSync(domList,2);//return(0);

      chan.renewVariables();
      exte.renewVariables();

      chan.swapVarMems();
      exte.swapVarMems();

      if(chan.computeResiduals() && exte.computeResiduals())break;

      chan.storeResiduals(iter);
      if(iter%200==0){cout<<iter<<" ";chan.displayResiduals();}
      if(iter%500==0)
	{
	  chan.dataStore(0);
	  chan.dataStore(1);
	  chan.dataStore(2);
	  exte.dataStore(0);
	  exte.dataStore(1);
	  exte.dataStore(2);
	}
      iter++;//return(0);///////////////////////////////
    }
  chan.resetBoundary();
  exte.resetBoundary();

  cout<<"hello\n";
  return(0);
}
