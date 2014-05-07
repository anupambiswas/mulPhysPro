#ifndef FLOW_DOMAIN_H
#define FLOW_DOMAIN_H
#include<iostream>
#include<cstdlib>
#include<string>
#include<cmath>
#include<fstream>
class flowDomain
{
  typedef struct boundaryZone
  {
    int bSize,bcType,aln,sid,sidFac;
    int **inter,**cell,**ifc,**first,**second;
    double **vel,*pr,*Te;
    double mulFac0,mulFac1;
    struct boundaryZone *next;
    boundaryZone(){vel=NULL;pr=Te=NULL;next=NULL;}
  }bZone;

 protected:

  std::string domainName;

  int nx[2],nxa[2],nxp1[2],nxp2[2];
  double dx,dy,dt,dxbdt;
  double cellAR,dxReInv,dxRePrInv,c0,c1,c2,c3;
  double relU,relV,relP;
  double pstab;
  int simType,numDat;

  std::string varName[4];
  int iStart[4],iEnd[4],jStart[4],jEnd[4];
  double iCor[4],jCor[4];

  double iTolU,iTolV,iTolP,iTolT;
  double sTolU,sTolV,sTolP,sTolT;

  double Re,Pr;
  double **u,**v,**p,**un,**vn,**pc,**T,**Tn;
  double **fxu,**fyu,**fxv,**fyv;
  double **mcu,**mcv,**mcp,**scw,**sce,**scs,**scn,**souu,**souv,**soup;
  double **sc[2][2],**source[4];
  double **mcof[4];
  double **mcT,**souT;
  double **uRes,**vRes;

  double errorU,errorV,errorP,errorT;
  double xOri,yOri;

  bZone *allBoun;

  void allocateMemoryFor(double ***var);
  void allocateMemoryForAll();

 public:
  void computeFlux();
  void computeCoefficients(int dec);
  void computePCoefficients();

  void resetBoundary();
  void resetSolverCoef(int varIde);
  virtual double fluxLim(double rat);
  void renewVariables();
  void swapVarMems();
  bool computeResiduals();
  void storeResiduals(int iter);
  void displayResiduals();
  void dataStore(int dec);
  void displayResidualsT();
  void copyVars();

  void solveGS(int dec,double rel);
  void initialize(int dec,double **val);
  double getT(const int& i,const int& j);

  std::string getFilename(int tstep,std::string filetype,int m=100000);
  void initialize(std::string file);
  void writeData(int tstep=0);

  flowDomain(std::string nam="flow_domain");

  void setGridProp(int nx,int ny,double dxx,double dyy,double dtt=0.0);
  void setFlowProp(double Rey,double Pra,int simtype=0);
  void setIterTol(double tolu,double tolv,double tolp,double tolT=1.0e-8);
  void setSoluTol(double tolu,double tolv,double tolp,double tolT=1.0e-8);
  void setRelPar(double relu,double relv,double relp);

  void setOrigin(double x0,double y0);

  void setPstab(double stab);

  void defineBoundary(int align,int side,int start,int end,int bctype,double **velSource=NULL,double *pSource=NULL,double *TSource=NULL,double rk=0);

  void salvage();

  void resetBoundaryT();
  void resetSolverCoefT();

  void dataStoreT();
  void swapVarMemsT();
  bool computeResidualsT();
  void setMulFacs(double mFac0=-1.0,double mFac1=2.0);
};
#endif
