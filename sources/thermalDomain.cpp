#include"../headers/thermalDomain.h"

void thermalDomain::allocateMemory(double ***var)
{
  int i,j;
  *var=new double* [nxa[0]];
  for(i=0;i<nxa[0];i++)
    (*var)[i]=new double [nxa[1]];
  for(i=0;i<nxa[0];i++)
    for(j=0;j<nxa[1];j++)
      (*var)[i][j]=0.0;
}

void thermalDomain::allocateMemoryForAll()
{
  allocateMemory(&T);
  allocateMemory(&Tn);
  allocateMemory(&scw);
  allocateMemory(&sce);
  allocateMemory(&scs);
  allocateMemory(&scn);
  allocateMemory(&mc);
  allocateMemory(&sou);

  sc[0][0]=scs;
  sc[0][1]=scn;
  sc[1][0]=scw;
  sc[1][1]=sce;
}

thermalDomain::thermalDomain(std::string name)
{
  domainName=name;
  allBoun=NULL;
  TTol=1.0e-8;
  rel=0.7;
  xOri=yOri=0.0;
}

void thermalDomain::setGridProp(int Nx,int Ny,double dxx,double dyy,double dtt)
{
  dx[0]=dxx;
  dx[1]=dyy;
  dt=dtt;
  cellAR=dx[0]/dx[1];
  dxbdt=dx[0]/dt;

  nx[0]=Nx;
  nx[1]=Ny;
  nxp1[0]=Nx+1;
  nxp1[1]=Ny+1;
  nxa[0]=Nx+2;
  nxa[1]=Ny+2;

  allocateMemoryForAll();
}

void thermalDomain::setThermProp(double Rey,double Pra,int simtype)
{
  Re=Rey;
  Pr=Pra;

  simType=simtype;
  if(!(simType==0 || simType==1))
    {
      std::cout<<"Invalid simulation type chosen. 0 for steady-state, 1 for unsteady. Terminating.\n";
      exit(0);
    }

  if(simType==0)
    {
      dxbdt=0.0;
      Re=Pr=1.0;
    }

  cwe=1.0/(dx[0]*Re*Pr);
  csn=cellAR*cellAR*cwe;
  cmc=-dxbdt-2.0*(1+cellAR*cellAR)*cwe;
  std::cout<<cwe<<" "<<csn<<" "<<cmc<<std::endl;
}

void thermalDomain::setSoluTol(double tol)
{
  TTol=tol;
}

void thermalDomain::setRelPar(double relp)
{
  rel=relp;
}

void thermalDomain::setOrigin(double x0,double y0)
{
  xOri=x0;
  yOri=y0;
}

void thermalDomain::del(double ***var)
{
  for(int i=0;i<nxa[0];i++)
    delete [] (*var)[i];
  delete [] (*var);
}

thermalDomain::~thermalDomain()
{
  del(&T);
  del(&Tn);
}

void thermalDomain::computeCoefficients()
{
  int i,j;
  for(i=1;i<nxp1[0];i++)
    for(j=1;j<nxp1[1];j++)
      {
	scw[i][j]=sce[i][j]=cwe;
	scs[i][j]=scn[i][j]=csn;
	mc[i][j]=cmc;
      }
}

void thermalDomain::solveGS()
{
  double maxErr,err,newval,inc;
  int i,j;
  do
    {
      maxErr=-1.0e40;
      for(i=1;i<nxp1[0];i++)
	for(j=1;j<nxp1[1];j++)
	  {
	    newval=(sou[i][j]-(scw[i][j]*Tn[i-1][j]+sce[i][j]*Tn[i+1][j]+scs[i][j]*Tn[i][j-1]+scn[i][j]*Tn[i][j+1]))/mc[i][j];
	    inc=newval-Tn[i][j];
	    err=fabs(inc);
	    Tn[i][j]+=rel*inc;
	    if(maxErr<err)maxErr=err;
	  }
    }
  while(maxErr>TTol);
}

void thermalDomain::defineBoundary(int align,int side,int start,int end,int bctype,double *valSource,double rk)
{
  bZone *bAdd=new bZone;
  int i;
  bAdd->next=allBoun;
  allBoun=bAdd;
  bAdd->bSize=end-start;
  bAdd->bcType=bctype;
  bAdd->val=valSource;
  bAdd->aln=align;
  bAdd->sid=side;
  switch(bctype)
    {
    case(0):
      bAdd->mulFac0=-1.0;
      bAdd->mulFac1=2.0;
      break;
    case(1):
      bAdd->mulFac0=1.0;
      bAdd->mulFac1=dx[!align]*(2*!side-1);
      break;
    case(2):
      bAdd->mulFac0=(1.0-rk)/(1.0+rk);
      bAdd->mulFac1=2.0*rk/(1.0+rk);
      break;
    default:
      std::cout<<"Invalid boundary type. Permissible types: 0 for Dirichlet, 1 for Neumann, 2 for interface. Terminating for now.\n";
      exit(0);
      break;
    }
  bAdd->cell=new int* [bAdd->bSize];
  bAdd->first=new int* [bAdd->bSize];
  for(i=0;i<bAdd->bSize;i++)
    {
      bAdd->cell[i]=new int [2];
      bAdd->first[i]=new int [2];
    }
  for(i=start;i<end;i++)
    {
      bAdd->cell[i-start][align]=bAdd->first[i-start][align]=i;
      bAdd->cell[i-start][!align]=side*nxp1[!align];
      bAdd->first[i-start][!align]=side*nxp1[!align]+(2*!side-1);
    }
}

void thermalDomain::resetBoundary()
{
  bZone *bList=allBoun;
  int **cel,**fir;
  int i;
  while(bList)
    {
      cel=bList->cell;
      fir=bList->first;
      for(i=0;i<bList->bSize;i++)
	T[cel[i][0]][cel[i][1]]=bList->mulFac0*T[fir[i][0]][fir[i][1]]+bList->mulFac1*bList->val[i];
      bList=bList->next;
    }
}

void thermalDomain::resetCoefficients()
{
  int i,j;
  for(i=1;i<nxp1[0];i++)
    {
      mc[i][1]=cmc;
      scs[i][1]=csn;
      sou[i][1]=0.0;
      mc[i][nx[1]]=cmc;
      scn[i][nx[1]]=csn;
      sou[i][nx[1]]=0.0;
    }
  for(j=1;j<nxp1[1];j++)
    {
      mc[1][j]=cmc;
      scw[1][j]=cwe;
      sou[1][j]=0.0;
      mc[nx[0]][j]=cmc;
      sce[nx[0]][j]=cwe;
      sou[nx[0]][j]=0.0;
    }
}

void thermalDomain::resetSolverCoef()
{
  bZone *bList=allBoun;
  int **cel,**fir;
  int i,is,js;
  while(bList)
    {
      cel=bList->cell;
      fir=bList->first;
      for(i=0;i<bList->bSize;i++)
	{
	  is=bList->first[i][0];
	  js=bList->first[i][1];
	  mc[is][js]+=sc[bList->aln][bList->sid][is][js]*bList->mulFac0;
	  sou[is][js]-=sc[bList->aln][bList->sid][is][js]*bList->mulFac1*bList->val[i];
	  sc[bList->aln][bList->sid][is][js]=0.0;
	}
      bList=bList->next;
    }
}

void thermalDomain::swapVarMems()
{
  double **temp;
  temp=T;
  T=Tn;
  Tn=temp;
}

void thermalDomain::computeSource()
{
  int i,j;
  for(i=1;i<nxp1[0];i++)
    for(j=1;j<nxp1[1];j++)
      sou[i][j]=-T[i][j]*dxbdt;//+other terms if needed (-qprime*dx[0]/(Re*Pr))
}

double thermalDomain::getT(const int& i,const int& j)
{
  return(T[i][j]);
}

bool thermalDomain::computeError()
{
  int i,j;
  double maxErr=-1.0e40,err;
  for(i=1;i<nxp1[0];i++)
    for(j=1;j<nxp1[1];j++)
      {
	err=fabs(Tn[i][j]-T[i][j]);
	if(maxErr<err)maxErr=err;
      }
  std::cout<<" "<<maxErr<<" ";
  if(maxErr<TTol)return(false);
  return(true);
}

void thermalDomain::initialize(double **val)
{
  int i,j;
  for(i=0;i<nxa[0];i++)
    for(j=0;j<nxa[1];j++)
      T[i][j]=val[i][j];
}

void thermalDomain::initialize(std::string fileName)
{
  int nXa[2];
  std::ifstream F(fileName.c_str(),std::ios::in|std::ios::binary);
  F.read((char*)&nXa[0],sizeof(int));
  F.read((char*)&nXa[1],sizeof(int));
  if(!(nXa[0]==nxa[0] && nXa[1]==nxa[1]))
    {
      std::cout<<"Initializing file size discrepancy. Terminating.\n";
      exit(0);
    }
  for(int i=0;i<nxa[0];i++)
    for(int j=0;j<nxa[1];j++)
      F.read((char*)&T[i][j],sizeof(double));
  F.close();
}

std::string thermalDomain::getFilename(int tstep,std::string filename,int m)
{
  int dig;
  while(m>0)
    {
      dig=tstep/m;
      tstep-=dig*m;
      m/=10;
      filename=filename+static_cast<char>(dig+48);
    }
  return(filename);
}

void thermalDomain::writeData(int tstep)
{
  std:: string filename=getFilename(tstep,domainName);
  std::ofstream F(filename.c_str(),std::ios::out|std::ios::binary);

  F.write((char*)&nxa[0],sizeof(int));
  F.write((char*)&nxa[1],sizeof(int));

  for(int i=0;i<nxa[0];i++)
    for(int j=0;j<nxa[1];j++)
      F.write((char*)&T[i][j],sizeof(double));
  F.close();
}

void thermalDomain::dataStore()
{
  int i,j;
  std::ofstream F((domainName+"T").c_str());
  for(i=1;i<nxp1[0];i++)
    for(j=1;j<nxp1[1];j++)
      F<<(i-0.5)*dx[0]+xOri<<" "<<(j-0.5)*dx[1]+yOri<<" "<<T[i][j]<<std::endl;
  F.close();
}
