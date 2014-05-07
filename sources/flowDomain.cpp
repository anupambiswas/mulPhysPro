#include"../headers/flowDomain.h"

void flowDomain::allocateMemoryFor(double ***var)
{
  int i,j;
  *var=new double* [nxa[0]];
  for(i=0;i<nxa[0];i++)
    (*var)[i]=new double [nxa[1]];
  for(i=0;i<nxa[0];i++)
    for(j=0;j<nxa[1];j++)
      (*var)[i][j]=0.0;
}

void flowDomain::allocateMemoryForAll()
{
  allocateMemoryFor(&u);
  allocateMemoryFor(&v);
  allocateMemoryFor(&p);

  allocateMemoryFor(&un);
  allocateMemoryFor(&vn);
  allocateMemoryFor(&pc);

  allocateMemoryFor(&fxu);
  allocateMemoryFor(&fyu);
  allocateMemoryFor(&fxv);
  allocateMemoryFor(&fyv);

  allocateMemoryFor(&mcu);
  allocateMemoryFor(&mcv);
  allocateMemoryFor(&mcp);
  allocateMemoryFor(&scw);
  allocateMemoryFor(&sce);
  allocateMemoryFor(&scs);
  allocateMemoryFor(&scn);
  allocateMemoryFor(&souu);
  allocateMemoryFor(&souv);
  allocateMemoryFor(&soup);

  allocateMemoryFor(&T);
  allocateMemoryFor(&Tn);

  allocateMemoryFor(&mcT);
  allocateMemoryFor(&souT);

  allocateMemoryFor(&uRes);
  allocateMemoryFor(&vRes);

  sc[0][0]=scs;
  sc[0][1]=scn;
  sc[1][0]=scw;
  sc[1][1]=sce;

  source[0]=souu;
  source[1]=souv;
  source[2]=soup;

  mcof[0]=mcu;
  mcof[1]=mcv;
  mcof[2]=mcp;

  //mcT=mcp;
  mcof[3]=mcT;
  //souT=soup;
  source[3]=souT;
  //Tn=pc;
}

flowDomain::flowDomain(std::string nam)
{
  domainName=nam;
  allBoun=NULL;
  pstab=0.1;
  numDat=3;

  // Automatic initialization for some simulation parameters
  iTolU=iTolV=iTolP=iTolT=1.0e-8;
  sTolU=sTolV=sTolP=sTolT=1.0e-8;
  relU=relV=0.3;
  relP=0.7;

  // Automatic initialization of geometric origin
  xOri=yOri=0.0;
}

void flowDomain::setGridProp(int Nx,int Ny,double dxx,double dyy,double dtt)
{
  dx=dxx;
  dy=dyy;
  dt=dtt;
  cellAR=dx/dy;
  dxbdt=dx/dt;

  nx[0]=Nx;
  nx[1]=Ny;
  nxa[0]=Nx+4;
  nxa[1]=Ny+4;
  nxp1[0]=Nx+1;
  nxp1[1]=Ny+1;
  nxp2[0]=Nx+2;
  nxp2[1]=Ny+2;

  allocateMemoryForAll();

  //xvel
  iStart[0]=3;
  iEnd[0]=nxp2[0];
  jStart[0]=2;
  jEnd[0]=nxp2[1];
  iCor[0]=2.0;
  jCor[0]=1.5;
  varName[0]="u";
  //yvel
  iStart[1]=2;
  iEnd[1]=nxp2[0];
  jStart[1]=3;
  jEnd[1]=nxp2[1];
  iCor[1]=1.5;
  jCor[1]=2.0;
  varName[1]="v";
  //pres
  iStart[2]=2;
  iEnd[2]=nxp2[0];
  jStart[2]=2;
  jEnd[2]=nxp2[1];
  iCor[2]=1.5;
  jCor[2]=1.5;
  varName[2]="p";
  //Temp
  iStart[3]=2;
  iEnd[3]=nxp2[0];
  jStart[3]=2;
  jEnd[3]=nxp2[1];
  iCor[3]=1.5;
  jCor[3]=1.5;
  varName[3]="T";

}

void flowDomain::setFlowProp(double Rey,double Pra,int simtype)
{
  Re=Rey;
  Pr=Pra;

  dxReInv=1.0/(dx*Re);
  c0=2.0*dxReInv*(1.0+cellAR*cellAR);
  c1=cellAR*cellAR*dxReInv;

  dxRePrInv=dxReInv/Pr;
  c2=2.0*dxRePrInv*(1.0+cellAR*cellAR);
  c3=cellAR*cellAR*dxRePrInv;

  simType=simtype;
  if(!(simtype==0 || simtype==1))
    {
      std::cout<<"Invalid flowtype. 0 for steady, 1 for unsteady. Terminating.\n";
      exit(0);
    }
  if(simType==0)
    dxbdt=0.0;
}


void flowDomain::setOrigin(double x0,double y0)
{
  xOri=x0;
  yOri=y0;
}

void flowDomain::setIterTol(double tolu,double tolv,double tolp,double tolT)
{
  iTolU=tolu;
  iTolV=tolv;
  iTolP=tolp;
  iTolT=tolT;

}

void flowDomain::setSoluTol(double tolu,double tolv,double tolp,double tolT)
{
  sTolU=tolu;
  sTolV=tolv;
  sTolP=tolp;
  sTolT=tolT;
}


void flowDomain::setRelPar(double relu,double relv,double relp)
{
  relU=relu;
  relV=relv;
  relP=relp;
}

void flowDomain::setPstab(double stab)
{
  pstab=stab;
}

void flowDomain::computeFlux()
{
  int i,j;

  for(i=3;i<=nxp2[0];i++)
    for(j=2;j<nxp2[1];j++)
      fxu[i][j]=0.5*(u[i-1][j]+u[i][j]);

  for(i=3;i<nxp2[0];i++)
    for(j=2;j<=nxp2[1];j++)
      fyu[i][j]=0.5*(v[i-1][j]+v[i][j]);

  for(i=2;i<=nxp2[0];i++)
    for(j=3;j<nxp2[1];j++)
      fxv[i][j]=0.5*(u[i][j]+u[i][j-1]);

  for(i=2;i<nxp2[0];i++)
    for(j=3;j<=nxp2[1];j++)
      fyv[i][j]=0.5*(v[i][j-1]+v[i][j]);
}

void flowDomain::computeCoefficients(int dec)
{
  int i,j,is,js,ibeg,jbeg;
  double Fw,Fe,Fs,Fn,WCT;
  int aw,ae,as,an;
  double veldif,rat,pfac,con0,con1,con2;
  double **fluxX,**fluxY,**var,**mc,**sou;
  int pdecI,pdecJ;
  switch(dec)
    {
    case(0):
      fluxX=fxu;
      fluxY=fyu;
      var=uRes;
      mc=mcu;
      sou=souu;
      ibeg=3;
      jbeg=2;
      pdecI=-1;
      pdecJ=0;
      pfac=1.0;
      con0=c0;
      con1=c1;
      con2=dxReInv;
      break;
    case(1):
      fluxX=fxv;
      fluxY=fyv;
      var=vRes;
      mc=mcv;
      sou=souv;
      ibeg=2;
      jbeg=3;
      pdecI=0;
      pdecJ=-1;
      pfac=cellAR;
      con0=c0;
      con1=c1;
      con2=dxReInv;
      break;
    case(3):
      fluxX=u;
      fluxY=v;
      var=T;
      mc=mcT;
      sou=souT;
      ibeg=2;
      jbeg=2;
      pdecI=0;
      pdecJ=0;
      pfac=0.0;
      con0=c2;
      con1=c3;
      con2=dxRePrInv;
      break;
    default:
      std::cout<<"\nSome error in option-selection (0 for u, 1 for v, 3 for Temp) of computeVelCoefficients(int).\nTerminating from computeCoefficients().\n\n";
      exit(0);
      break;
    }


  for(i=ibeg;i<nxp2[0];i++)
    for(j=jbeg;j<nxp2[1];j++)
      {
	Fw=fluxX[i][j];
	Fe=fluxX[i+1][j];
	Fs=fluxY[i][j];
	Fn=fluxY[i][j+1];

	aw=Fw>0.0;
	ae=Fe>0.0;
	as=Fs>0.0;
	an=Fn>0.0;

	mc[i][j]=Fe*ae-Fw*!aw+cellAR*(Fn*an-Fs*!as)+con0+dxbdt;
	scw[i][j]=-Fw*aw-con2;
	sce[i][j]=Fe*!ae-con2;
	scs[i][j]=-cellAR*Fs*as-con1;
	scn[i][j]=cellAR*Fn*!an-con1;

	WCT=0.0;
	//west face
	is=i-1+2*!aw;
	veldif=var[i][j]-var[i-1][j];
	rat=(var[is][j]-var[is-1][j])/veldif;
	WCT-=0.5*Fw*fluxLim(rat)*(2.0*aw-1)*veldif;

	//east face
	is=i+2*!ae;
	veldif=var[i+1][j]-var[i][j];
	rat=(var[is][j]-var[is-1][j])/veldif;
	WCT+=0.5*Fe*fluxLim(rat)*(2.0*ae-1)*veldif;

	//south face
	js=j-1+2*!as;
	veldif=var[i][j]-var[i][j-1];
	rat=(var[i][js]-var[i][js-1])/veldif;
	WCT-=cellAR*0.5*Fs*fluxLim(rat)*(2.0*as-1)*veldif;

	//north face
	js=j+2*!an;
	veldif=var[i][j+1]-var[i][j];
	rat=(var[i][js]-var[i][js-1])/veldif;
	WCT+=cellAR*0.5*Fn*fluxLim(rat)*(2.0*an-1)*veldif;

	sou[i][j]=pfac*(p[i+pdecI][j+pdecJ]-p[i][j])-WCT+dxbdt*var[i][j];
      }
}

void flowDomain::computePCoefficients()
{
  int i,j;
  for(i=2;i<nxp2[0];i++)
    for(j=2;j<nxp2[1];j++)
      {
	scw[i][j]=-1.0/mcu[i][j];
	sce[i][j]=-1.0/mcu[i+1][j];
	scs[i][j]=-cellAR/mcv[i][j];
	scn[i][j]=-cellAR/mcv[i][j+1];
	soup[i][j]=-(u[i+1][j]-u[i][j]+cellAR*(v[i][j+1]-v[i][j]));
      }

  for(i=2;i<nxp2[0];i++)
    for(j=2;j<nxp2[1];j++)
      mcp[i][j]=-(scw[i][j]+sce[i][j]+scs[i][j]+scn[i][j])+pstab;

}

void flowDomain::solveGS(int dec,double rel)
{
  int i,j;
  double newVal,absVal,maxDif,inc;
  double **mc,**var,**sou;
  double toler;
  int is,ie,js,je;
  switch(dec)
    {
    case(0):
      mc=mcu;
      var=un;
      sou=souu;
      toler=iTolU;
      is=3;
      ie=nxp2[0];
      js=2;
      je=nxp2[1];
      break;
    case(1):
      mc=mcv;
      var=vn;
      sou=souv;
      toler=iTolV;
      is=2;
      ie=nxp2[0];
      js=3;
      je=nxp2[1];
      break;
    case(2):
      mc=mcp;
      var=pc;
      sou=soup;
      toler=iTolP;
      is=2;
      ie=nxp2[0];
      js=2;
      je=nxp2[1];
      break;
    case(3):
      mc=mcT;
      var=Tn;
      sou=souT;
      toler=iTolT;
      is=2;
      ie=nxp2[0];
      js=2;
      je=nxp2[1];
      break;
    }

  do
    {
      maxDif=-1.0e40;
      for(i=is;i<ie;i++)
	for(j=js;j<je;j++)
	  {
	    newVal=(sou[i][j]-(scw[i][j]*var[i-1][j]+sce[i][j]*var[i+1][j]+scs[i][j]*var[i][j-1]+scn[i][j]*var[i][j+1]))/mc[i][j];
	    inc=rel*(newVal-var[i][j]);
	    absVal=fabs(inc);
	    var[i][j]+=inc;
	    if(maxDif>absVal)maxDif=absVal;
	  }
    }
  while(maxDif>toler);
}

void flowDomain::defineBoundary(int align,int side,int start,int end,int bctype,double **velSource,double *pSource,double *TSource,double rk)
{
  bZone *bAdd=new bZone;
  int i;
  bAdd->next=allBoun;
  allBoun=bAdd;
  bAdd->bSize=end-start;
  bAdd->bcType=bctype;
  bAdd->vel=velSource;
  bAdd->pr=pSource;
  bAdd->Te=TSource;
  bAdd->aln=align;
  bAdd->sid=side;
  bAdd->sidFac=2*!side-1;

  switch(bctype)
    {
    case(0)://WALL
      bAdd->mulFac0=(1-rk)/(1+rk);
      bAdd->mulFac1=2.0*rk/(1+rk);
      break;
    case(1)://VINLET
      bAdd->mulFac0=-1.0;
      bAdd->mulFac1=2.0;
      break;
    case(2)://POUTLET
      break;
    default:
      std::cout<<"\nInvalid BC type. 0-wall, 1-velocity inlet, 2-pressure outlet. Terminating.\n\n";
      exit(0);
      break;
    }

  bAdd->inter=new int* [bAdd->bSize];
  bAdd->cell=new int* [bAdd->bSize];
  bAdd->ifc=new int* [bAdd->bSize];
  bAdd->first=new int* [bAdd->bSize];
  bAdd->second=new int* [bAdd->bSize];
  for(i=0;i<bAdd->bSize;i++)
    {
      bAdd->inter[i]=new int [2];
      bAdd->cell[i]=new int [2];
      bAdd->ifc[i]=new int [2];
      bAdd->first[i]=new int [2];
      bAdd->second[i]=new int [2];
    }

  int inc=2*!side-1,ims;

  for(i=start;i<end;i++)
    {
      ims=i-start;
      bAdd->inter[ims][align]=bAdd->cell[ims][align]=bAdd->ifc[ims][align]=bAdd->first[ims][align]=bAdd->second[ims][align]=i;
      bAdd->cell[ims][!align]=1+side*(nxp2[!align]-1);
      bAdd->first[ims][!align]=bAdd->cell[ims][!align]+inc;
      bAdd->second[ims][!align]=bAdd->cell[ims][!align]+2*inc;
      bAdd->inter[ims][!align]=bAdd->cell[ims][!align]-inc;
      bAdd->ifc[ims][!align]=(side==0)?bAdd->first[ims][!align]:bAdd->cell[ims][!align];
    }
}

void flowDomain::resetBoundary()
{
  double **veloc[2]={u,v};
  int fir[2],sec[2];
  bZone* bAdd=allBoun;
  int i;

  while(bAdd)
    {
      for(i=0;i<bAdd->bSize;i++)
	switch(bAdd->bcType)
	  {
	  case(0):;

	  case(1):
	    veloc[!bAdd->aln][bAdd->ifc[i][0]][bAdd->ifc[i][1]]=bAdd->vel[i][!bAdd->aln];
	    veloc[bAdd->aln][bAdd->cell[i][0]][bAdd->cell[i][1]]=2.0*bAdd->vel[i][bAdd->aln]-veloc[bAdd->aln][bAdd->first[i][0]][bAdd->first[i][1]];
	    //TVD addition
	    if(bAdd->sid==0)
	      veloc[!bAdd->aln][bAdd->cell[i][0]][bAdd->cell[i][1]]=2.0*veloc[!bAdd->aln][bAdd->first[i][0]][bAdd->first[i][1]]-veloc[!bAdd->aln][bAdd->second[i][0]][bAdd->second[i][1]];
	    veloc[!bAdd->aln][bAdd->inter[i][0]][bAdd->inter[i][1]]=2.0*veloc[!bAdd->aln][bAdd->cell[i][0]][bAdd->cell[i][1]]-veloc[!bAdd->aln][bAdd->first[i][0]][bAdd->first[i][1]];
	    veloc[bAdd->aln][bAdd->inter[i][0]][bAdd->inter[i][1]]=2.0*veloc[bAdd->aln][bAdd->cell[i][0]][bAdd->cell[i][1]]-veloc[bAdd->aln][bAdd->first[i][0]][bAdd->first[i][1]];
	    break;

	  case(2):
	    fir[bAdd->aln]=sec[bAdd->aln]=bAdd->cell[i][bAdd->aln];
	    fir[!bAdd->aln]=bAdd->ifc[i][!bAdd->aln]+bAdd->sidFac;
	    sec[!bAdd->aln]=bAdd->ifc[i][!bAdd->aln]+2*bAdd->sidFac;
	    veloc[!bAdd->aln][bAdd->ifc[i][0]][bAdd->ifc[i][1]]=2.0*veloc[!bAdd->aln][fir[0]][fir[1]]-veloc[!bAdd->aln][sec[0]][sec[1]];
	    veloc[bAdd->aln][bAdd->cell[i][0]][bAdd->cell[i][1]]=2.0*veloc[bAdd->aln][bAdd->first[i][0]][bAdd->first[i][1]]-veloc[bAdd->aln][bAdd->second[i][0]][bAdd->second[i][1]];
	    p[bAdd->first[i][0]][bAdd->first[i][1]]=(p[bAdd->second[i][0]][bAdd->second[i][1]]+2.0*bAdd->pr[i])/3.0;
	    //TVD addition
	    if(bAdd->sid==0)
	      veloc[!bAdd->aln][bAdd->cell[i][0]][bAdd->cell[i][1]]=2.0*veloc[!bAdd->aln][bAdd->first[i][0]][bAdd->first[i][1]]-veloc[!bAdd->aln][bAdd->second[i][0]][bAdd->second[i][1]];
	    veloc[!bAdd->aln][bAdd->inter[i][0]][bAdd->inter[i][1]]=2.0*veloc[!bAdd->aln][bAdd->cell[i][0]][bAdd->cell[i][1]]-veloc[!bAdd->aln][bAdd->first[i][0]][bAdd->first[i][1]];
	    veloc[bAdd->aln][bAdd->inter[i][0]][bAdd->inter[i][1]]=2.0*veloc[bAdd->aln][bAdd->cell[i][0]][bAdd->cell[i][1]]-veloc[bAdd->aln][bAdd->first[i][0]][bAdd->first[i][1]];
	    break;
	  }
      bAdd=bAdd->next;
    }
}

void flowDomain::resetSolverCoef(int varIde)
{
  double **veloc[2]={u,v};
  int fir[2],sec[2];
  int is,js,i;
  double alFac[2]={1.0,cellAR};
  bZone *bAdd;
  bAdd=allBoun;

  switch(varIde)
    {
    case(0):;//xvel
    case(1)://yvel
      while(bAdd)
	{
	  switch(bAdd->bcType)
	    {
	    case(0):;//WALL
	    case(1)://VINLET
	      if(varIde==bAdd->aln)
		for(i=0;i<bAdd->bSize;i++)
		  {
		    mcof[varIde][bAdd->first[i][0]][bAdd->first[i][1]]-=sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]];
		    source[varIde][bAdd->first[i][0]][bAdd->first[i][1]]-=2.0*bAdd->vel[i][bAdd->aln]*sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]];
		    sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]]=0.0;
		  }
	      break;
	    case(2)://POTLET
	      if(varIde==bAdd->aln)
		for(i=0;i<bAdd->bSize;i++)
		  {
		    mcof[varIde][bAdd->first[i][0]][bAdd->first[i][1]]+=2*sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]];
		    sc[bAdd->aln][!bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]]-=sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]];
		    sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]]=0.0;
		  }
	      else
		for(i=0;i<bAdd->bSize;i++)
		  {
		    fir[bAdd->aln]=sec[bAdd->aln]=bAdd->cell[i][bAdd->aln];
		    fir[!bAdd->aln]=bAdd->ifc[i][!bAdd->aln]+bAdd->sidFac;
		    mcof[varIde][fir[0]][fir[1]]+=2*sc[bAdd->aln][bAdd->sid][fir[0]][fir[1]];
		    sc[bAdd->aln][!bAdd->sid][fir[0]][fir[1]]-=sc[bAdd->aln][bAdd->sid][fir[0]][fir[1]];
		    sc[bAdd->aln][bAdd->sid][fir[0]][fir[1]]=0.0;
		  }
	      break;
	    }
	  bAdd=bAdd->next;
	}
      break;
    case(2)://pressure
      while(bAdd)
	{
	  switch(bAdd->bcType)
	    {
	    case(0):;//xvel
	    case(1)://yvel
	      for(i=0;i<bAdd->bSize;i++)
		{
		  sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]]=0.0;
		  is=bAdd->first[i][0];
		  js=bAdd->first[i][1];
		  mcp[is][js]=-(scw[is][js]+sce[is][js]+scs[is][js]+scn[is][js])+pstab;
		}
	      break;
	    case(2)://POUTLET
	      for(i=0;i<bAdd->bSize;i++)
		{
		  mcp[bAdd->first[i][0]][bAdd->first[i][1]]=1.0;
		  soup[bAdd->first[i][0]][bAdd->first[i][1]]=0.0;
		  sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]]=0.0;
		  sc[bAdd->aln][!bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]]=-1.0/3;
		  sc[!bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]]=0.0;
		  sc[!bAdd->aln][!bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]]=0.0;
		}
	      break;
	    }
	  bAdd=bAdd->next;
	}
      break;
    default:
      std::cout<<"\nIncorrect solver coefficient option (0 for u, 1 for v, 2 for p). Terminating.\n\n";
      exit(0);
      break;
    }
}

double flowDomain::fluxLim(double r)
{
  if(r!=r)return(0.0);
  if(r<0.0)return(0.0);
  if(r<1.0)return(r);
  return(1.0);
}

void flowDomain::renewVariables()
{
  int i,j;
  for(i=3;i<nxp2[0];i++)
    for(j=2;j<nxp2[1];j++)
      un[i][j]+=relU*(pc[i-1][j]-pc[i][j])/mcu[i][j];

  for(i=2;i<nxp2[0];i++)
    for(j=3;j<nxp2[1];j++)
      vn[i][j]+=relV*cellAR*(pc[i][j-1]-pc[i][j])/mcv[i][j];

  for(i=2;i<nxp2[0];i++)
    for(j=2;j<nxp2[1];j++)
      p[i][j]+=relP*pc[i][j];
}

void flowDomain::swapVarMems()
{
  double **temp;
  temp=u;
  u=un;
  un=temp;

  temp=v;
  v=vn;
  vn=temp;
}

bool flowDomain::computeResiduals()
{
  int i,j;
  double tv;

  errorU=errorV=errorP=-1.0e40;

  for(i=3;i<nxp2[0];i++)
    for(j=2;j<nxp2[1];j++)
      {
	tv=fabs(un[i][j]-u[i][j]);
	if(tv>errorU)errorU=tv;
      }

  for(i=2;i<nxp2[0];i++)
    for(j=3;j<nxp2[1];j++)
      {
	tv=fabs(vn[i][j]-v[i][j]);
	if(tv>errorV)errorV=tv;
      }


  for(i=2;i<nxp2[0];i++)
    for(j=2;j<nxp2[1];j++)
      {
	tv=fabs(pc[i][j]);
	if(tv>errorP)errorP=tv;
      }

  if(errorU<sTolU && errorV<sTolV && errorP<sTolP)
    {
      std::cout<<"Solution converged.\n";
      return(true);
    }
  return(false);
}

void flowDomain::storeResiduals(int iter)
{
  std::ofstream F;

  F.open("res_u",std::ios::app);
  F<<iter<<" "<<log10(errorU)<<std::endl;
  F.close();

  F.open("res_v",std::ios::app);
  F<<iter<<" "<<log10(errorV)<<std::endl;
  F.close();

  F.open("res_p",std::ios::app);
  F<<iter<<" "<<log10(errorP)<<std::endl;
  F.close();
}

void flowDomain::displayResiduals()
{
  //std::cout<<sTolU<<" "<<sTolV<<" "<<sTolP<<std::endl;
  std::cout<<errorU<<" "<<errorV<<" "<<errorP<<std::endl;//exit(0);
}

void flowDomain::resetBoundaryT()
{
  double **Te=T;
  int i;
  bZone *bAdd;
  bAdd=allBoun;
  numDat=4;

  while(bAdd)
    {
      switch(bAdd->bcType)
	{
	case(0):;
	case(1):
	  for(i=0;i<bAdd->bSize;i++)
	    {
	      Te[bAdd->cell[i][0]][bAdd->cell[i][1]]=bAdd->mulFac0*Te[bAdd->first[i][0]][bAdd->first[i][1]]+bAdd->mulFac1*bAdd->Te[i];
	      Te[bAdd->inter[i][0]][bAdd->inter[i][1]]=2.0*Te[bAdd->cell[i][0]][bAdd->cell[i][1]]-Te[bAdd->first[i][0]][bAdd->first[i][1]];
	    }
	  break;
	case(2):
	  for(i=0;i<bAdd->bSize;i++)
	    {
	      Te[bAdd->cell[i][0]][bAdd->cell[i][1]]=2.0*Te[bAdd->first[i][0]][bAdd->first[i][1]]-Te[bAdd->second[i][0]][bAdd->second[i][1]];
	      Te[bAdd->inter[i][0]][bAdd->inter[i][1]]=2.0*Te[bAdd->cell[i][0]][bAdd->cell[i][1]]-Te[bAdd->first[i][0]][bAdd->first[i][1]];
	    }
	  break;
	}
      bAdd=bAdd->next;
    }
}

void flowDomain::resetSolverCoefT()
{
  int fir[2],sec[2];
  int is,js,i;
  bZone *bAdd;
  bAdd=allBoun;

  while(bAdd)
    {
      switch(bAdd->bcType)
	{
	case(0):
	  for(i=0;i<bAdd->bSize;i++)
	    {
	      is=bAdd->first[i][0];
	      js=bAdd->first[i][1];
	      mcT[is][js]+=sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]]*bAdd->mulFac0;
	      souT[is][js]-=bAdd->mulFac1*bAdd->Te[i]*sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]];
	      sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]]=0.0;
	    }
	  break;
	case(1):
	  for(i=0;i<bAdd->bSize;i++)
	    {
	      is=bAdd->first[i][0];
	      js=bAdd->first[i][1];
	      mcT[is][js]-=sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]];
	      souT[is][js]-=2.0*bAdd->Te[i]*sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]];
	      sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]]=0.0;
	    }
	  break;
	case(2):
	  for(i=0;i<bAdd->bSize;i++)
	    {
	      mcT[bAdd->first[i][0]][bAdd->first[i][1]]+=2.0*sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]];
	      sc[bAdd->aln][!bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]]-=sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]];
	      sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]]=0.0;
	    }
	  break;
	}
      bAdd=bAdd->next;
    }
}


void flowDomain::setMulFacs(double mFac0,double mFac1)
{
  allBoun->mulFac0=mFac0;
  allBoun->mulFac1=mFac1;
}

bool flowDomain::computeResidualsT()
{
  int i,j;
  double errT;
  errorT=-1.0e40;
  for(i=2;i<nxp2[0];i++)
    for(j=2;j<nxp2[1];j++)
      {
	errT=fabs(Tn[i][j]-T[i][j]);
	if(errorT<errT)errorT=errT;
      }
  if(errorT<1.0e-8)return(true);
  return(false);
}

void flowDomain::displayResidualsT()
{
  std::cout<<errorT<<std::endl;
}

void flowDomain::swapVarMemsT()
{
  double **temp;
  temp=T;
  T=Tn;
  Tn=temp;
}

void flowDomain::initialize(int dec,double **val)
{
  double **var;
  switch(dec)
    {
    case(0):
      var=u;
      break;
    case(1):
      var=v;
      break;
    case(2):
      var=p;
      break;
    case(3):
      var=T;
      break;
    default:
      std::cout<<"\nInvalid initialize option. Terminating from initialize().\n\n";
      exit(0);
      break;
    }
  int i,j;
  for(i=0;i<nxa[0];i++)
    for(j=0;j<nxa[1];j++)
      var[i][j]=val[i][j];
}

double flowDomain::getT(const int& i,const int& j)
{
  return(T[i][j]);
}

//class functions henceforth are for postprocessing and data storage

void flowDomain::initialize(std::string file)
{
  std::ifstream F(file.c_str(),std::ios::in|std::ios::binary);
  int i,j,nXa[2];

  F.read((char*)&nXa[0],sizeof(int));
  F.read((char*)&nXa[1],sizeof(int));
  F.read((char*)&numDat,sizeof(int));

  if(!(nXa[0]==nxa[0] && nXa[1]==nxa[1]))
    {
      std::cout<<"Discrepancy in input file. Dimension mismatch. Terminating.\n";
      exit(0);
    }

  switch(numDat)
    {
    case(3):
      std::cout<<"\nX-velocity, Y-velocity and Pressure read from unput file.\n\n";
      break;
    case(4):
      std::cout<<"\nX-velocity, Y-velocity, Pressure and Temperature read from unput file.\n\n";
      break;
    default:
      std::cout<<"\nIncorrect input file. Insufficient/incorrect number of data variables. Terminating.\n\n";
      exit(0);
      break;
    }

  for(i=0;i<nxa[0];i++)
    for(j=0;j<nxa[1];j++)
      F.read((char*)&u[i][j],sizeof(double));

  for(i=0;i<nxa[0];i++)
    for(j=0;j<nxa[1];j++)
      F.read((char*)&v[i][j],sizeof(double));

  for(i=0;i<nxa[0];i++)
    for(j=0;j<nxa[1];j++)
      F.read((char*)&p[i][j],sizeof(double));

  if(numDat==4)
    for(i=0;i<nxa[0];i++)
      for(j=0;j<nxa[1];j++)
	F.read((char*)&T[i][j],sizeof(double));

  F.close();
}

void flowDomain::writeData(int tstep)
{
  std::string filename=getFilename(tstep,domainName);
  std::ofstream F(filename.c_str(),std::ios::out|std::ios::binary);
  int i,j;

  F.write((char*)&nxa[0],sizeof(int));
  F.write((char*)&nxa[1],sizeof(int));
  F.write((char*)&numDat,sizeof(int));

  for(i=0;i<nxa[0];i++)
    for(j=0;j<nxa[1];j++)
      F.write((char*)&u[i][j],sizeof(double));

  for(i=0;i<nxa[0];i++)
    for(j=0;j<nxa[1];j++)
      F.write((char*)&v[i][j],sizeof(double));

  for(i=0;i<nxa[0];i++)
    for(j=0;j<nxa[1];j++)
      F.write((char*)&p[i][j],sizeof(double));

  if(numDat==4)
    for(i=0;i<nxa[0];i++)
      for(j=0;j<nxa[1];j++)
	F.write((char*)&T[i][j],sizeof(double));

  F.close();
}

std::string flowDomain::getFilename(int tstep,std::string filename,int m)
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

void flowDomain::dataStore(int dec)
{
  int i,j;
  std::ofstream F((domainName+varName[dec]).c_str());
  double **var;
  switch(dec)
    {
    case(0):
      var=u;
      break;
    case(1):
      var=v;
      break;
    case(2):
      var=p;
      break;
    case(3):
      var=T;
      break;
    default:
      std::cout<<"\nInvalid data storing option. Returning from dataStore().\n\n";
      break;
    }
  for(i=iStart[dec];i<iEnd[dec];i++)
    for(j=jStart[dec];j<jEnd[dec];j++)
      F<<(i-iCor[dec])*dx+xOri<<" "<<(j-jCor[dec])*dy+yOri<<" "<<var[i][j]<<std::endl;
  F.close();
}

void flowDomain::copyVars()
{
  int i,j;
  for(i=0;i<nxa[0];i++)
    for(j=0;j<nxa[1];j++)
      {
	uRes[i][j]=u[i][j];
	vRes[i][j]=v[i][j];
      }
}
