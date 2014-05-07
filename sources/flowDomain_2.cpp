#include"../headers/flowDomain_2.h"

void flowDomain_2::allocateMemoryFor(double ***var) //already correct
{
  int i,j;
  *var=new double* [nxa[0]];
  for(i=0;i<nxa[0];i++)
    (*var)[i]=new double [nxa[1]];
  for(i=0;i<nxa[0];i++)
    for(j=0;j<nxa[1];j++)
      (*var)[i][j]=0.0;
}

void flowDomain_2::allocateMemoryForAll() //modifying as required
{
  allocateMemoryFor(&u);
  allocateMemoryFor(&v);
  allocateMemoryFor(&p);
  allocateMemoryFor(&Tn);

  allocateMemoryFor(&un);
  allocateMemoryFor(&vn);
  allocateMemoryFor(&pc);
  allocateMemoryFor(&Tn);

  allocateMemoryFor(&fx);
  allocateMemoryFor(&fy);

  //allocateMemoryFor(&mcuv);
  allocateMemoryFor(&mcp);

  allocateMemoryFor(&mcuv);
  allocateMemoryFor(&scw);
  allocateMemoryFor(&sce);
  allocateMemoryFor(&scs);
  allocateMemoryFor(&scn);
  allocateMemoryFor(&souu);
  allocateMemoryFor(&souv);
  allocateMemoryFor(&soup);
  allocateMemoryFor(&extu);
  allocateMemoryFor(&extv);
  allocateMemoryFor(&mcT);
  allocateMemoryFor(&souT);

  allocateMemoryFor(&otu);
  allocateMemoryFor(&otv);
  allocateMemoryFor(&dAvgU);
  allocateMemoryFor(&dAvgV);

  allocateMemoryFor(&uFace);
  allocateMemoryFor(&vFace);

  sc[0][0]=scs;
  sc[0][1]=scn;
  sc[1][0]=scw;
  sc[1][1]=sce;

  source[0]=souu;
  source[1]=souv;
  source[2]=soup;
}

flowDomain_2::flowDomain_2(std::string nam) //ok until now
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

  xOri=yOri=0.0;
}

void flowDomain_2::setGridProp(int Nx,int Ny,double dxx,double dyy,double dtt) //ok until now
{
  dx[0]=dxx;
  dx[1]=dyy;
  dt=dtt;
  cellAR=dx[0]/dx[1];
  dxbdt=dx[0]/dt;

  nx[0]=Nx;
  nx[1]=Ny;
  nxa[0]=Nx+4;
  nxa[1]=Ny+4;
  nxp1[0]=Nx+1;
  nxp1[1]=Ny+1;
  nxp2[0]=Nx+2;
  nxp2[1]=Ny+2;
  nxp3[0]=Nx+3;
  nxp3[1]=Ny+3;

  allocateMemoryForAll();
}

void flowDomain_2::setFlowProp(double Rey,double Pra,int simtype) //corrected
{
  Re=Rey;

  dxReInv=1.0/(dx[0]*Re);
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

void flowDomain_2::setIterTol(double tolu,double tolv,double tolp) //corrected
{
  iTolU=tolu;
  iTolV=tolv;
  iTolP=tolp;
}

void flowDomain_2::setSoluTol(double tolu,double tolv,double tolp) //corrected
{
  sTolU=tolu;
  sTolV=tolv;
  sTolP=tolp;
}


void flowDomain_2::setRelPar(double relu,double relv,double relp) //corrected
{
  relU=relu;
  relV=relv;
  relP=relp;
}

void flowDomain_2::setPstab(double stab) //corrected
{
  pstab=stab;
}

void flowDomain_2::computeFlux() //corrected
{
  int i,j;

  for(i=2;i<nxp3[0];i++)
    for(j=2;j<nxp2[1];j++)
      fx[i][j]=0.5*(u[i-1][j]+u[i][j]);

  for(i=2;i<nxp2[0];i++)
    for(j=2;j<nxp3[1];j++)
      fy[i][j]=0.5*(v[i][j]+v[i][j-1]);
}

void flowDomain_2::computeCoefficients() //corrected
{
  int i,j,is,js;
  double Fw,Fe,Fs,Fn,WCTu,WCTv;
  int aw,ae,as,an;
  double veldif,rat,pfac,con0,con1,con2;
  double **var,**mc,**sou;
  int pdecI,pdecJ;

  mc=mcuv;
  con0=c0;
  con1=c1;
  con2=dxReInv;

  for(i=2;i<nxp2[0];i++)
    for(j=2;j<nxp2[1];j++)
      {
	Fw=fx[i][j];
	Fe=fx[i+1][j];
	Fs=fy[i][j];
	Fn=fy[i][j+1];

	aw=Fw>0.0;
	ae=Fe>0.0;
	as=Fs>0.0;
	an=Fn>0.0;

	mc[i][j]=Fe*ae-Fw*!aw+cellAR*(Fn*an-Fs*!as)+con0+dxbdt;
	scw[i][j]=-Fw*aw-con2;
	sce[i][j]=Fe*!ae-con2;
	scs[i][j]=-cellAR*Fs*as-con1;
	scn[i][j]=cellAR*Fn*!an-con1;

	WCTu=0.0;
	WCTv=0.0;
	//west face
	is=i-1+2*!aw;

	veldif=u[i][j]-u[i-1][j];
	rat=(u[is][j]-u[is-1][j])/veldif;
	WCTu-=0.5*Fw*fluxLim(rat)*(2.0*aw-1)*veldif;

	veldif=v[i][j]-v[i-1][j];
	rat=(v[is][j]-v[is-1][j])/veldif;
	WCTv-=0.5*Fw*fluxLim(rat)*(2.0*aw-1)*veldif;

	//east face
	is=i+2*!ae;

	veldif=u[i+1][j]-u[i][j];
	rat=(u[is][j]-u[is-1][j])/veldif;
	WCTu+=0.5*Fe*fluxLim(rat)*(2.0*ae-1)*veldif;

	veldif=v[i+1][j]-v[i][j];
	rat=(v[is][j]-v[is-1][j])/veldif;
	WCTv+=0.5*Fe*fluxLim(rat)*(2.0*ae-1)*veldif;

	//south face
	js=j-1+2*!as;

	veldif=u[i][j]-u[i][j-1];
	rat=(u[i][js]-u[i][js-1])/veldif;
	WCTu-=cellAR*0.5*Fs*fluxLim(rat)*(2.0*as-1)*veldif;

	veldif=v[i][j]-v[i][j-1];
	rat=(v[i][js]-v[i][js-1])/veldif;
	WCTv-=cellAR*0.5*Fs*fluxLim(rat)*(2.0*as-1)*veldif;

	//north face
	js=j+2*!an;

	veldif=u[i][j+1]-u[i][j];
	rat=(u[i][js]-u[i][js-1])/veldif;
	WCTu+=cellAR*0.5*Fn*fluxLim(rat)*(2.0*an-1)*veldif;

	veldif=v[i][j+1]-v[i][j];
	rat=(v[i][js]-v[i][js-1])/veldif;
	WCTv+=cellAR*0.5*Fn*fluxLim(rat)*(2.0*an-1)*veldif;

	souu[i][j]=0.5*(p[i-1][j]-p[i+1][j])-WCTu+dxbdt*u[i][j];
	souv[i][j]=0.5*cellAR*(p[i][j-1]-p[i][j+1])-WCTv+dxbdt*v[i][j];
      }
}

void flowDomain_2::computePCoefReq() //added
{
  int i,j;
  double **davg[2]={dAvgU,dAvgV};

  for(i=1;i<nxp3[0];i++)
    for(j=1;j<nxp3[1];j++)
      {
	dAvgU[i][j]=0.5/mcuv[i][j]+0.5/mcuv[i-1][j];
	dAvgV[i][j]=0.5/mcuv[i][j]+0.5/mcuv[i][j-1];

	otu[i][j]=u[i][j]-(p[i-1][j]-p[i+1][j])/(2.0*mcuv[i][j]);
	otv[i][j]=v[i][j]-(p[i][j-1]-p[i][j+1])*cellAR/(2.0*mcuv[i][j]);

	//if(i==nxp1[0]){dAvgU[i][j]=1.0/mcuv[i][j];otu[i+1][j]=otu[i][j];}
      }
  //return;
  bZone *bAdd=allBoun;
  int sec[2];
  while(bAdd)
    {
      //std::cout<<"btype "<<bAdd->bcType<<std::endl;
      switch(bAdd->bcType)
	{
	case(0):;//wall
	case(1)://vinelt
	  for(i=0;i<bAdd->bSize;i++)
	    {
	      if(bAdd->sid==0)
		davg[!bAdd->aln][bAdd->first[i][0]][bAdd->first[i][1]]=0.0;
	      else
		davg[!bAdd->aln][bAdd->cell[i][0]][bAdd->cell[i][1]]=0.0;

	      otu[bAdd->cell[i][0]][bAdd->cell[i][1]]=-otu[bAdd->first[i][0]][bAdd->first[i][1]]+2.0*bAdd->vel[i][0];
	      otv[bAdd->cell[i][0]][bAdd->cell[i][1]]=-otv[bAdd->first[i][0]][bAdd->first[i][1]]+2.0*bAdd->vel[i][1];
	    }
	  break;
	case(2)://pressure
	  for(i=0;i<bAdd->bSize;i++)
	    {
	      if(bAdd->sid==0)
		{
		  sec[bAdd->aln]=bAdd->second[i][bAdd->aln];
		  sec[!bAdd->aln]=bAdd->second[i][!bAdd->aln]+bAdd->sidFac;
		  davg[!bAdd->aln][bAdd->first[i][0]][bAdd->first[i][1]]=1.5*davg[!bAdd->aln][bAdd->second[i][0]][bAdd->second[i][1]]-0.5*davg[!bAdd->aln][sec[0]][sec[1]];
		}
	      else
		davg[!bAdd->aln][bAdd->cell[i][0]][bAdd->cell[i][1]]=1.5*davg[!bAdd->aln][bAdd->first[i][0]][bAdd->first[i][1]]-0.5*davg[!bAdd->aln][bAdd->second[i][0]][bAdd->second[i][1]];

	      otu[bAdd->cell[i][0]][bAdd->cell[i][1]]=2.0*otu[bAdd->first[i][0]][bAdd->first[i][1]]-otu[bAdd->second[i][0]][bAdd->second[i][1]];
	      otv[bAdd->cell[i][0]][bAdd->cell[i][1]]=2.0*otv[bAdd->first[i][0]][bAdd->first[i][1]]-otv[bAdd->second[i][0]][bAdd->second[i][1]];
	    }
	  break;
	case(3):
	  break;
	}
      bAdd=bAdd->next;
    }

  for(i=2;i<nxp3[0];i++)
    for(j=2;j<nxp3[1];j++)
      {
	uFace[i][j]=dAvgU[i][j]*(p[i-1][j]-p[i][j])+0.5*(otu[i-1][j]+otu[i][j]);
	vFace[i][j]=dAvgV[i][j]*(p[i][j-1]-p[i][j])*cellAR+0.5*(otv[i][j-1]+otv[i][j]);
	//uFace[i][j]=0.5*(u[i-1][j]+u[i][j]);
	//vFace[i][j]=0.5*(v[i][j-1]+v[i][j]);continue;
      }

  if(nxp2[0]==102)
    {
      i=101;
      std::ofstream F("fcu0",std::ios::out);
      for(j=2;j<nxp2[1];j++)
	F<<j<<" "<<otu[i][j]<<std::endl;
      //F<<j<<" "<<dAvgU[i][j]*(p[i-1][j]-p[i][j])+0.5*(otu[i-1][j]+otu[i][j])<<std::endl;
      F.close();
    }
  if(nxp2[0]==82)
    {
      i=1;
      std::ofstream F("fcu1",std::ios::out);
      for(j=22;j<nxp2[1];j++)
	F<<j-20<<" "<<otu[i][j]<<std::endl;
      //F<<j-20<<" "<<dAvgU[i][j]*(p[i-1][j]-p[i][j])+0.5*(otu[i-1][j]+otu[i][j])<<std::endl;
      F.close();
    }
}

void flowDomain_2::computePCoefficients() //corrected
{
  int i,j;
  for(i=2;i<nxp2[0];i++)
    for(j=2;j<nxp2[1];j++)
      {
	scw[i][j]=-dAvgU[i][j];
	sce[i][j]=-dAvgU[i+1][j];
	scs[i][j]=-cellAR*dAvgV[i][j];
	scn[i][j]=-cellAR*dAvgV[i][j+1];
	soup[i][j]=-(uFace[i+1][j]-uFace[i][j]+cellAR*(vFace[i][j+1]-vFace[i][j]));
	mcp[i][j]=-(scw[i][j]+sce[i][j]+scs[i][j]+scn[i][j])+pstab;
      }

  return;
  double mx=-1.0e40,mn=1.0e40;
  int mxi,mxj,mni,mnj;
  double **var=uFace;
  for(i=2;i<nxp2[0];i++)
    for(j=2;j<nxp2[1];j++)
      {
	if(mx<var[i][j]){mx=var[i][j];mxi=i;mxj=j;}
	if(mn>var[i][j]){mn=var[i][j];mni=i;mnj=j;}
      }
  std::cout<<"mx = "<<mx<<" mn = "<<mn<<std::endl;
  std::cout<<"mxi mxj "<<mxi<<" "<<mxj<<" mni mnj "<<mni<<" "<<mnj<<std::endl<<std::endl;
}

void flowDomain_2::getPCorrectionFromNeibs()
{
  bZone *bAdd=allBoun;
  int nind[2],aln,sid,offSet,i;
  flowDomain_2 *neib;
  while(bAdd)
    {
      if(bAdd->bcType==3)
	{
	  aln=bAdd->aln;
	  sid=bAdd->sid;
	  neib=bAdd->neibDom;
	  nind[!aln]=2+!sid*(neib->nxp2[!aln]-3);
	  offSet=bAdd->neibStart-bAdd->cell[0][aln];
	  for(i=0;i<bAdd->bSize;i++)
	    {
	      nind[aln]=bAdd->cell[i][aln]+offSet;
	      pc[bAdd->cell[i][0]][bAdd->cell[i][1]]=neib->pc[nind[0]][nind[1]];
	    }
	}
      bAdd=bAdd->next;
    }
}

void flowDomain_2::solveGS(int dec,double rel) //corrected
{
  int i,j;
  double newVal,absVal,maxDif,inc;
  double **mc,**var,**sou;
  double toler;
  switch(dec)
    {
    case(0):
      mc=mcuv;
      var=un;
      sou=souu;
      toler=iTolU;
      break;
    case(1):
      mc=mcuv;
      var=vn;
      sou=souv;
      toler=iTolV;
      break;
    case(2):
      mc=mcp;
      var=pc;
      sou=soup;
      toler=iTolP;
      break;
    }

  do
    {
      maxDif=-1.0e40;
      for(i=2;i<nxp2[0];i++)
	for(j=2;j<nxp2[1];j++)
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

void flowDomain_2::solvePSingle(double rel,double &maxDif)
{
  int i,j;
  double newVal,absVal,inc;
  double **mc,**var,**sou;
  double toler;

  mc=mcp;
  var=pc;
  sou=soup;
  toler=iTolP;

  for(i=2;i<nxp2[0];i++)
    for(j=2;j<nxp2[1];j++)
      {
	newVal=(sou[i][j]-(scw[i][j]*var[i-1][j]+sce[i][j]*var[i+1][j]+scs[i][j]*var[i][j-1]+scn[i][j]*var[i][j+1]))/mc[i][j];
	inc=rel*(newVal-var[i][j]);
	absVal=fabs(inc);
	var[i][j]+=inc;
	if(maxDif>absVal)maxDif=absVal;
      }
}

void flowDomain_2::defineBoundary(int align,int side,int start,int end,int bctype,double **velSource,double *pSource,double *TSource,double rk,flowDomain_2 *neib,int nStart) //corrected
{
  bZone *bAdd=new bZone;
  int i;
  bAdd->next=allBoun;
  allBoun=bAdd;
  bAdd->bSize=end-start;
  bAdd->bcType=bctype;
  bAdd->vel=velSource;
  bAdd->pr=pSource;
  bAdd->neibDom=neib;
  bAdd->neibStart=nStart;
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
    case(2):;//POUTLET
    case(3):
      break;
    default:
      std::cout<<"\nInvalid BC type. 0-wall, 1-velocity inlet, 2-pressure outlet. Terminating.\n\n";
      exit(0);
      break;
    }

  //std::cout<<"kkhello "<<bAdd->bSize<<"\n";
  //if(bAdd->bSize==0)return;

  bAdd->inter=new int* [bAdd->bSize];
  bAdd->cell=new int* [bAdd->bSize];
  bAdd->first=new int* [bAdd->bSize];
  bAdd->second=new int* [bAdd->bSize];
  for(i=0;i<bAdd->bSize;i++)
    {
      bAdd->inter[i]=new int [2];
      bAdd->cell[i]=new int [2];
      bAdd->first[i]=new int [2];
      bAdd->second[i]=new int [2];
    }

  int inc=2*!side-1,ims;

  for(i=start;i<end;i++)
    {
      ims=i-start;
      bAdd->inter[ims][align]=bAdd->cell[ims][align]=bAdd->first[ims][align]=bAdd->second[ims][align]=i;
      bAdd->cell[ims][!align]=1+side*(nxp2[!align]-1);
      bAdd->first[ims][!align]=bAdd->cell[ims][!align]+inc;
      bAdd->second[ims][!align]=bAdd->cell[ims][!align]+2*inc;
      bAdd->inter[ims][!align]=bAdd->cell[ims][!align]-inc;
    }
}

void flowDomain_2::resetBoundary() //corrected
{
  bZone* bAdd=allBoun;
  int i,aln,sid,nind[2],offSet,nind2[2];
  flowDomain_2 *neib;

  while(bAdd)
    {
      aln=bAdd->aln;
      sid=bAdd->sid;
      neib=bAdd->neibDom;
      if(neib)nind[!aln]=2+!sid*(neib->nxp2[!aln]-3);
      nind2[!aln]=nind[!aln]+2*sid-1;
      if(bAdd->bSize>0)offSet=bAdd->neibStart-bAdd->cell[0][aln];
      for(i=0;i<bAdd->bSize;i++)
	{
	  switch(bAdd->bcType)
	    {
	    case(0):;

	    case(1):
	      u[bAdd->cell[i][0]][bAdd->cell[i][1]]=2.0*bAdd->vel[i][0]-u[bAdd->first[i][0]][bAdd->first[i][1]];
	      v[bAdd->cell[i][0]][bAdd->cell[i][1]]=2.0*bAdd->vel[i][1]-v[bAdd->first[i][0]][bAdd->first[i][1]];
	      p[bAdd->cell[i][0]][bAdd->cell[i][1]]=2.0*p[bAdd->first[i][0]][bAdd->first[i][1]]-p[bAdd->second[i][0]][bAdd->second[i][1]];
	      break;

	    case(2):
	      u[bAdd->cell[i][0]][bAdd->cell[i][1]]=2.0*u[bAdd->first[i][0]][bAdd->first[i][1]]-u[bAdd->second[i][0]][bAdd->second[i][1]];
	      v[bAdd->cell[i][0]][bAdd->cell[i][1]]=2.0*v[bAdd->first[i][0]][bAdd->first[i][1]]-v[bAdd->second[i][0]][bAdd->second[i][1]];
	      p[bAdd->cell[i][0]][bAdd->cell[i][1]]=2.0*bAdd->pr[i]-p[bAdd->first[i][0]][bAdd->first[i][1]];
	      break;

	    case(3):
	      nind[aln]=nind2[aln]=bAdd->cell[i][aln]+offSet;
	      u[bAdd->cell[i][0]][bAdd->cell[i][1]]=neib->u[nind[0]][nind[1]];
	      v[bAdd->cell[i][0]][bAdd->cell[i][1]]=neib->v[nind[0]][nind[1]];
	      p[bAdd->cell[i][0]][bAdd->cell[i][1]]=neib->p[nind[0]][nind[1]];
	      u[bAdd->inter[i][0]][bAdd->inter[i][1]]=neib->u[nind2[0]][nind2[1]];
	      v[bAdd->inter[i][0]][bAdd->inter[i][1]]=neib->v[nind2[0]][nind2[1]];
	      p[bAdd->inter[i][0]][bAdd->inter[i][1]]=neib->p[nind2[0]][nind2[1]];
	      //std::cout<<"Copying to ["<<bAdd->cell[i][0]<<"]["<<bAdd->cell[i][1]<<"], copying from ["<<nind[0]<<"]["<<nind[1]<<"].\n";
	      //std::cout<<"Copying to ["<<bAdd->inter[i][0]<<"]["<<bAdd->inter[i][1]<<"], copying from ["<<nind2[0]<<"]["<<nind2[1]<<"].\n\n";
	      break;
	    }
	  //TVD addition
	  if(bAdd->bcType!=3)
	    {
	      u[bAdd->inter[i][0]][bAdd->inter[i][1]]=2.0*u[bAdd->cell[i][0]][bAdd->cell[i][1]]-u[bAdd->first[i][0]][bAdd->first[i][1]];
	      v[bAdd->inter[i][0]][bAdd->inter[i][1]]=2.0*v[bAdd->cell[i][0]][bAdd->cell[i][1]]-v[bAdd->first[i][0]][bAdd->first[i][1]];
	    }
	}

      bAdd=bAdd->next;
    }
}

void flowDomain_2::resetSolverCoef(int varIde) //woi - pressure part varIde = 2, at present guessing wont require varIde = 2
{
  double **veloc[2]={u,v};
  int fir[2],sec[2];
  int is,js,i;
  double alFac[2]={1.0,cellAR};
  bZone *bAdd;
  bAdd=allBoun;

  flowDomain_2 *neib;
  int nind[2],nind2[2],aln,sid,offSet,ci,cj;

  switch(varIde)
    {
    case(0):;//xvel
    case(1)://yvel
      while(bAdd)
	{
	  aln=bAdd->aln;
	  sid=bAdd->sid;
	  neib=bAdd->neibDom;
	  if(neib)nind[!aln]=2+!sid*(neib->nxp2[!aln]-3);
	  nind2[!aln]=nind[!aln]+2*sid-1;
	  if(bAdd->bSize>0)offSet=bAdd->neibStart-bAdd->cell[0][aln];
	  switch(bAdd->bcType)
	    {
	    case(0):;//WALL
	    case(1)://VINLET
	      for(i=0;i<bAdd->bSize;i++)
		{
		  mcuv[bAdd->first[i][0]][bAdd->first[i][1]]-=sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]];
		  source[0][bAdd->first[i][0]][bAdd->first[i][1]]-=2.0*bAdd->vel[i][0]*sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]];
		  source[1][bAdd->first[i][0]][bAdd->first[i][1]]-=2.0*bAdd->vel[i][1]*sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]];
		  sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]]=0.0;
		}
	      break;
	    case(2)://POUTLET
	      for(i=0;i<bAdd->bSize;i++)
		{
		  mcuv[bAdd->first[i][0]][bAdd->first[i][1]]+=2*sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]];
		  sc[bAdd->aln][!bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]]-=sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]];
		  sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]]=0.0;
		}
	      break;
	    case(3)://INTERFACE
	      break;//////////////////////
	      for(i=0;i<bAdd->bSize;i++)
		{
		  nind[aln]=nind2[aln]=bAdd->cell[i][aln]+offSet;
		  ci=bAdd->cell[i][0];
		  cj=bAdd->cell[i][1];
		  mcuv[ci][cj]=neib->mcuv[nind[0]][nind[1]];
		  //std::cout<<"Copying mcuv["<<ci<<"]["<<cj<<"] in "<<domainName<<" from mcuv["<<nind[0]<<"]["<<nind[1]<<"] in "<<neib->domainName<<std::endl;
		  //scw[ci][cj]=neib->scw[nind[0]][nind[1]];
		  //sce[ci][cj]=neib->sce[nind[0]][nind[1]];
		  //scs[ci][cj]=neib->scs[nind[0]][nind[1]];
		  //scn[ci][cj]=neib->scn[nind[0]][nind[1]];
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
	      // for(i=0;i<bAdd->bSize;i++)
	      // 	{
	      // 	  sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]]=0.0;
	      // 	  is=bAdd->first[i][0];
	      // 	  js=bAdd->first[i][1];
	      // 	  mcp[is][js]=-(scw[is][js]+sce[is][js]+scs[is][js]+scn[is][js])+pstab;
	      // 	}
	      break;
	    case(2)://POUTLET
	      for(i=0;i<bAdd->bSize;i++)
		{
		  mcp[bAdd->first[i][0]][bAdd->first[i][1]]-=sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]];
		  sc[bAdd->aln][bAdd->sid][bAdd->first[i][0]][bAdd->first[i][1]]=0.0;
		}
	      break;
	    case(3):
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

void flowDomain_2::manageInterface() //woi - pressure part varIde = 2, at present guessing wont require varIde = 2
{
  double **veloc[2]={u,v};
  int fir[2],sec[2];
  int is,js,i;
  double alFac[2]={1.0,cellAR};
  bZone *bAdd;
  bAdd=allBoun;

  flowDomain_2 *neib;
  int nind[2],nind2[2],aln,sid,offSet,ci,cj;
  while(bAdd)
    {
      if(bAdd->bcType==3)
	{
	  aln=bAdd->aln;
	  sid=bAdd->sid;
	  neib=bAdd->neibDom;
	  if(neib)nind[!aln]=2+!sid*(neib->nxp2[!aln]-3);
	  nind2[!aln]=nind[!aln]+2*sid-1;
	  if(bAdd->bSize>0)offSet=bAdd->neibStart-bAdd->cell[0][aln];

	  for(i=0;i<bAdd->bSize;i++)
	    {
	      nind[aln]=nind2[aln]=bAdd->cell[i][aln]+offSet;
	      ci=bAdd->cell[i][0];
	      cj=bAdd->cell[i][1];
	      mcuv[ci][cj]=neib->mcuv[nind[0]][nind[1]];
	      //std::cout<<"Copying mcuv["<<ci<<"]["<<cj<<"] in "<<domainName<<" from mcuv["<<nind[0]<<"]["<<nind[1]<<"] in "<<neib->domainName<<std::endl;
	    }
	}
      bAdd=bAdd->next;
    }
}

double flowDomain_2::fluxLim(double r) //already correct
{
  //return(0.0);/////
  if(r!=r)return(0.0);
  if(r<0.0)return(0.0);
  if(r<1.0)return(r);
  return(1.0);
}

void flowDomain_2::renewVariables() //corrected
{
  int i,j;

  bZone *bAdd=allBoun;
  while(bAdd)
    {
      switch(bAdd->bcType)
	{
	case(0):;
	case(1):
	  for(i=0;i<bAdd->bSize;i++)
	    pc[bAdd->cell[i][0]][bAdd->cell[i][1]]=2.0*pc[bAdd->first[i][0]][bAdd->first[i][1]]-pc[bAdd->second[i][0]][bAdd->second[i][1]];
	  break;
	case(2):
	  for(i=0;i<bAdd->bSize;i++)
	    pc[bAdd->cell[i][0]][bAdd->cell[i][1]]=-pc[bAdd->first[i][0]][bAdd->first[i][1]];
	  break;
	}
      bAdd=bAdd->next;
    }

  for(i=2;i<nxp2[0];i++)
    for(j=2;j<nxp2[1];j++)
      {
	un[i][j]+=relU*(pc[i-1][j]-pc[i+1][j])/(2.0*mcuv[i][j]);
	vn[i][j]+=relV*cellAR*(pc[i][j-1]-pc[i][j+1])/(2.0*mcuv[i][j]);
	p[i][j]+=relP*pc[i][j];
      }
}

void flowDomain_2::swapVarMems() //already correct
{
  double **temp;
  temp=u;
  u=un;
  un=temp;

  temp=v;
  v=vn;
  vn=temp;
}

bool flowDomain_2::computeResiduals() //corrected
{
  int i,j;
  double tv;

  errorU=errorV=errorP=-1.0e40;

  for(i=2;i<nxp2[0];i++)
    for(j=2;j<nxp2[1];j++)
      {
	tv=fabs(un[i][j]-u[i][j]);
	if(tv>errorU)errorU=tv;

	tv=fabs(vn[i][j]-v[i][j]);
	if(tv>errorV)errorV=tv;

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

void flowDomain_2::storeResiduals(int iter) //already correct
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

void flowDomain_2::displayResiduals() //already correct
{
  std::cout<<errorU<<" "<<errorV<<" "<<errorP<<std::endl;
}

void flowDomain_2::dataStore(int dec) //corrected
{
  int i,j;
  std::ofstream F;//((domainName+varName[dec]).c_str());
  double **var;
  switch(dec)
    {
    case(0):
      F.open((domainName+"u").c_str());
      var=u;
      break;
    case(1):
      F.open((domainName+"v").c_str());
      var=v;
      break;
    case(2):
      F.open((domainName+"p").c_str());
      var=p;
      break;
    case(3):
      F.open((domainName+"T").c_str());
      var=T;
      break;
    default:
      std::cout<<"\nInvalid data storing option. Returning from dataStore().\n\n";
      break;
    }
  for(i=2;i<nxp2[0];i++)
    for(j=2;j<nxp2[1];j++)
      F<<(i-1.5)*dx[0]+xOri<<" "<<(j-1.5)*dx[1]+yOri<<" "<<var[i][j]<<std::endl;
  F.close();
}

void flowDomain_2::setOrigin(double x0,double y0)
{
  xOri=x0;
  yOri=y0;
}

void flowDomain_2::readForPP(std::string file)
{
  std::ifstream F(file.c_str(),std::ios::in|std::ios::binary);
  int i,j;
 
  F.read((char*)&nxa[0],sizeof(int));
  F.read((char*)&nxa[1],sizeof(int));
  F.read((char*)&numDat,sizeof(int));
 
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
 
void flowDomain_2::writeData(int tstep)
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
 
std::string flowDomain_2::getFilename(int tstep,std::string filename,int m)
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

void flowDomain_2::setMulFacs(double mFac0,double mFac1)
{
  allBoun->mulFac0=mFac0;
  allBoun->mulFac1=mFac1;
}

double flowDomain_2::getT(const int& i,const int& j)
{
  return(T[i][j]);
}

void flowDomain_2::swapVarMemsT()
{
  double **temp;
  temp=T;
  T=Tn;
  Tn=temp;
}

void flowDomain_2::resetBoundaryT()
{
  //double **Te=T;
  int i;
  bZone *bAdd;
  bAdd=allBoun;
  numDat=4;
  int aln,sid,nind[2],offSet,nind2[2];
  flowDomain_2 *neib;
 
  while(bAdd)
    {
      aln=bAdd->aln;
      sid=bAdd->sid;
      neib=bAdd->neibDom;
      if(neib)nind[!aln]=2+!sid*(neib->nxp2[!aln]-3);
      nind2[!aln]=nind[!aln]+2*sid-1;
      if(bAdd->bSize>0)offSet=bAdd->neibStart-bAdd->cell[0][aln];
      switch(bAdd->bcType)
	{
	case(0):;
	case(1):
	  for(i=0;i<bAdd->bSize;i++)
	    {
	      T[bAdd->cell[i][0]][bAdd->cell[i][1]]=bAdd->mulFac0*T[bAdd->first[i][0]][bAdd->first[i][1]]+bAdd->mulFac1*bAdd->Te[i];
	      T[bAdd->inter[i][0]][bAdd->inter[i][1]]=2.0*T[bAdd->cell[i][0]][bAdd->cell[i][1]]-T[bAdd->first[i][0]][bAdd->first[i][1]];
	    }
	  break;
	case(2):
	  for(i=0;i<bAdd->bSize;i++)
	    {
	      T[bAdd->cell[i][0]][bAdd->cell[i][1]]=2.0*T[bAdd->first[i][0]][bAdd->first[i][1]]-T[bAdd->second[i][0]][bAdd->second[i][1]];
	      T[bAdd->inter[i][0]][bAdd->inter[i][1]]=2.0*T[bAdd->cell[i][0]][bAdd->cell[i][1]]-T[bAdd->first[i][0]][bAdd->first[i][1]];
	    }
	  break;

	case(3):
	  for(i=0;i<bAdd->bSize;i++)
	    {
	      nind[aln]=nind2[aln]=bAdd->cell[i][aln]+offSet;
	      T[bAdd->cell[i][0]][bAdd->cell[i][1]]=neib->T[nind[0]][nind[1]];
	      T[bAdd->inter[i][0]][bAdd->inter[i][1]]=neib->T[nind2[0]][nind2[1]];
	    }
	  break;
	}
      bAdd=bAdd->next;
    }
}
 
void flowDomain_2::resetSolverCoefT()
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

void flowDomain_2::computeCoefficientsT()
{
  int i,j,is,js;
  double Fw,Fe,Fs,Fn,WCTT;
  int aw,ae,as,an;
  double veldif,rat,pfac,con0,con1,con2;
  double **var,**mc,**sou;
  int pdecI,pdecJ;

  mc=mcT;
  con0=c2;
  con1=c3;
  con2=dxRePrInv;

  for(i=2;i<nxp2[0];i++)
    for(j=2;j<nxp2[1];j++)
      {
	Fw=fx[i][j];
	Fe=fx[i+1][j];
	Fs=fy[i][j];
	Fn=fy[i][j+1];

	aw=Fw>0.0;
	ae=Fe>0.0;
	as=Fs>0.0;
	an=Fn>0.0;

	mc[i][j]=Fe*ae-Fw*!aw+cellAR*(Fn*an-Fs*!as)+con0+dxbdt;
	scw[i][j]=-Fw*aw-con2;
	sce[i][j]=Fe*!ae-con2;
	scs[i][j]=-cellAR*Fs*as-con1;
	scn[i][j]=cellAR*Fn*!an-con1;

	WCTT=0.0;
	//west face
	is=i-1+2*!aw;

	veldif=T[i][j]-T[i-1][j];
	rat=(T[is][j]-T[is-1][j])/veldif;
	WCTT-=0.5*Fw*fluxLim(rat)*(2.0*aw-1)*veldif;

	//east face
	is=i+2*!ae;

	veldif=T[i+1][j]-T[i][j];
	rat=(T[is][j]-T[is-1][j])/veldif;
	WCTT+=0.5*Fe*fluxLim(rat)*(2.0*ae-1)*veldif;

	//south face
	js=j-1+2*!as;

	veldif=T[i][j]-T[i][j-1];
	rat=(T[i][js]-T[i][js-1])/veldif;
	WCTT-=cellAR*0.5*Fs*fluxLim(rat)*(2.0*as-1)*veldif;

	//north face
	js=j+2*!an;

	veldif=T[i][j+1]-T[i][j];
	rat=(T[i][js]-T[i][js-1])/veldif;
	WCTT+=cellAR*0.5*Fn*fluxLim(rat)*(2.0*an-1)*veldif;

	souT[i][j]=-WCTT+dxbdt*T[i][j];
      }
}
