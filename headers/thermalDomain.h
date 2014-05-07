#ifndef THERMAL_DOMAIN_H
#define THERMAL_DOMAIN_H
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<string>
#include<fstream>

class thermalDomain
{
  typedef struct boundaryZone
  {
    int bSize,bcType,aln,sid;
    double mulFac0,mulFac1;
    int **cell,**first;
    double *val;
    struct boundaryZone *next;
    boundaryZone(){cell=first=NULL;val=NULL;next=NULL;}
  }bZone;
 protected:
  bZone *allBoun;

  int nx[2],nxp1[2],nxa[2];
  double **T,**Tn,**scw,**sce,**scs,**scn,**mc,**sou;
  double **sc[2][2];
  std::string domainName;
  double dx[2],dt,dxbdt,cellAR;
  double Re,Pr;
  double cwe,csn,cmc;
  int simType;
  double TTol,rel;
  double xOri,yOri;

  void allocateMemory(double ***var);
  void allocateMemoryForAll();
  void del(double ***var);
public:
  thermalDomain(std::string name="thermal_domain");
  ~thermalDomain();
  void setGridProp(int Nx,int Ny,double dxx,double dyy,double dtt=0.001);
  void setThermProp(double Rey,double Pra,int simtype=0);
  void setSoluTol(double tol);
  void setRelPar(double relp);

  void defineBoundary(int align,int side,int start,int end,int bctype,double *valSource,double rk=0.0);
  void resetBoundary();
  void resetSolverCoef();
  void solveGS();
  void computeCoefficients();
  void resetCoefficients();
  bool computeError();
  void swapVarMems();

  double getT(const int& i,const int& j);
  void setOrigin(double x0,double y0);
  void initialize(double **val);

  virtual void computeSource();

  std::string getFilename(int tstep,std::string filename,int m=100000);
  void initialize(std::string fileName);
  void writeData(int tstep=0);
  void dataStore();
};
#endif
