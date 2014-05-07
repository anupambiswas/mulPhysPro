#ifndef FLOW_DOMAIN_2_H
#define FLOW_DOMAIN_2_H
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<cmath>
#include<string>

class flowDomain_2
{
  typedef struct boundaryZone
  {
    int bSize,bcType,aln,sid,sidFac;
    int **inter,**cell,**first,**second;
    double **vel,*pr,*Te;
    double mulFac0,mulFac1;
    flowDomain_2 *neibDom;
    int neibStart;
    struct boundaryZone *next;
    boundaryZone(){vel=NULL;pr=NULL;Te=NULL;next=NULL;}
  }bZone;

 protected:

  std::string domainName;

  int nx[2],nxa[2],nxp1[2],nxp2[2],nxp3[2];
  double dx[2],dt,dxbdt;
  double cellAR,dxReInv,dxRePrInv,c0,c1,c2,c3;
  double relU,relV,relP;
  double pstab;
  int simType,numDat;

  double iTolU,iTolV,iTolP,iTolT;
  double sTolU,sTolV,sTolP,sTolT;
  double xOri,yOri;

  double Re,Pr;
  double **u,**v,**p,**un,**vn,**pc,**T,**Tn;
  double **fx,**fy;
  double **mcu,**mcv,**mcp,**scw,**sce,**scs,**scn,**souu,**souv,**soup;
  double **mcuv,**extu,**extv;
  double **sc[2][2],**source[3];
  double **mcof[3];
  double **otu,**otv,**dAvgU,**dAvgV;
  double **uFace,**vFace;
  double **mcT,**souT;

  double errorU,errorV,errorP;

  bZone *allBoun;

  void allocateMemoryFor(double ***var);
  void allocateMemoryForAll();

 public:
  void computeFlux();
  void computeCoefficients();
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

  void solveGS(int dec,double rel);

  flowDomain_2(std::string nam="flow_domain");

  void setGridProp(int nx,int ny,double dxx,double dyy,double dtt=0.0);
  void setFlowProp(double Rey,double Pra,int simtype=0);
  void setIterTol(double tolu,double tolv,double tolp);
  void setSoluTol(double tolu,double tolv,double tolp);
  void setRelPar(double relu,double relv,double relp);

  void setPstab(double stab);

  void defineBoundary(int align,int side,int start,int end,int bctype,double **velSource=NULL,double *pSource=NULL,double *TSource=NULL,double rk=0.0,flowDomain_2 *neib=NULL,int nStart=0);

  void computePCoefReq();
  void getPCorrectionFromNeibs();
  void setOrigin(double x0,double y0);
  void solvePSingle(double rel,double &maxDif);
  void manageInterface();
 
  std::string getFilename(int tstep,std::string filetype,int m=100000);
  void readForPP(std::string file);
  void writeData(int tstep=0);

  void setMulFacs(double mFac0,double mFac1);
  double getT(const int& i,const int& j);
  void swapVarMemsT();
  void resetBoundaryT();
  void resetSolverCoefT();
  void computeCoefficientsT();
};

#endif
