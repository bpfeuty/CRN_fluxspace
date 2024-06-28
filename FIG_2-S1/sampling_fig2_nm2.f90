!! gfortran sampling_fig2.f90 -lseulex2 -lminpack -llapack


module GLOB
  implicit none
  integer,parameter::nv=8,nr=12,nm=2,npop=10000,nc=2,kn=1
  real(kind=kind(0.d0)),parameter::RT=2.5

  !! nv is the number of species
  !! nv=8 in default simulation nv=24 in figS1
  !! nr is the number of reactions
  !! nr=12 in default simulation 28 or 48 in figS1
  !! nm is the molecularity of reaction =1 or 2
  !! nc is the number of chemostat
  !! npop is the sample size
  !! npop=10^5 for figure 2B and FIGS1
  !! npop=10^4 for figure 2C
  !! For vizualization npop=1000 in figure 2C
  !! kn is the normalization used for kinetic rate
  !! kn=0 (no normalization) kn=1 (normalization)
  
  integer,parameter::nplot=500000;
  real(kind=kind(0.d0)),dimension(nv,nplot) :: xplot(nv,nplot),tplot(nplot);

  real(kind=kind(0.d0)),dimension(nv,nr)::stoch
  real(kind=kind(0.d0)),dimension(nr,nv)::stocht
  real(kind=kind(0.d0)),dimension(nv+nc,nr)::stoch2
  real(kind=kind(0.d0)),dimension(nr,nv+nc)::stoch2t 
  real(kind=kind(0.d0)),dimension(nv,nv)::jac
  integer,dimension(nv,nv)::conn
  integer,dimension(nr,4)::add 
  integer::nr_max,nmax
  real(kind=kind(0.d0))::kcat(nr),keq(nr),dg0(nr),yy(nc)
 
contains
  
  !! SUBROUTINE FOR THE ODE SYSTEM
   subroutine sys_rand(v,dvdt,t) 
    real(8),dimension(nv),intent(in)::v
    real(8),dimension(nv),intent(out)::dvdt
    real(8),intent(in)::t
    real(8)::reac(nr),knorm
    integer::i,j,ki
    dvdt=0.;
    do i=1,nc
       dvdt(i)=dvdt(i)+kcat(i)*(yy(i)-v(i));
    enddo
    do i=1,nv
       do j=nc+1,nr;
          if (kn==0) knorm=1D0
          if (kn==1) knorm=sqrt(keq(j))
          
          if (nm==1) then
             reac(j)=kcat(j)*knorm*(v(add(j,2))-v(add(j,1))/keq(j)); !1:prod 2:sub
          else
             reac(j)=kcat(j)*knorm*(v(add(j,3))*v(add(j,4))-v(add(j,1))*v(add(j,2))/keq(j)); 
          endif
          dvdt(i)=reac(j)*stoch(i,j)+dvdt(i);
       enddo
    enddo
  end subroutine sys_rand
  
  !! SUBROUTINE FOR THE SEULEX-BASED ODE INTEGRATION
  subroutine evolSEU(tbeg,tend,X,idid,tol) !! numerical integration
    real(kind=kind(0.d0)), intent(in) :: tol
    real(kind=kind(0.d0))             :: tbeg,tend
    real(kind=kind(0.d0)), intent(inout) :: X(nv)
    integer, intent(out) :: idid
    real(kind=kind(0.d0)),save::h
    integer,parameter :: lrcont=20,licont=10
    integer :: itol,ijac,ifcn,iout
    integer :: mljac=nv,mujac=nv,imas=0
    integer :: MUMAS=nv,MLMAS=nv
    integer, parameter :: lwork=2*(nv+10)*nv+61+20
    real(kind=kind(0.d0)) :: work(lwork),rpar
    integer, parameter :: liwork=2*nv+32+20
    integer :: iwork(liwork),ipar
    real(kind=kind(0.d0)) :: Rtol,Atol
    external fcn_HW,solout,mass,jac_HW !! code externe utilisé
    itol=0; rtol=tol;
    atol=1d-7;  !! précision abolue
    iout=1; 	!! appelle soulout
    ijac=0; 	!! difference finie pour le jacobien -> =1  utilise JAC
    ifcn=0;	!! F(X,Y) INDEPENDENT OF X (AUTONOMOUS) -> =1 nonautonomous
    work=0d0   !! work(2)=0.5D0 !! mettre en commentaire si adaptatif
    iwork=0
    iwork(1)=0 		!! pas de transfo hessenberg du jac -> =1 sinon
    iwork(2)=10**8 	!! nb de step max
    iwork(4)=3  	!! step size sequence
    call SEULEX(nv,FCN_HW,IFCN,tbeg,X,tend,H,&
         &RTOL,ATOL,ITOL,&
         &JAC_HW ,IJAC,MLJAC,MUJAC,&
         &MASS ,IMAS,MLMAS,MUMAS,&
         &SOLOUT,IOUT,&
         &WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
    nr_max=min(iwork(17),nplot) !! nombre de pas
  end subroutine evolSEU
end  module GLOB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_steady_state   
  use GLOB; implicit none;

contains

  ! compute the solutions of a NLODE system
  subroutine syssol(x,info) 
    use glob;  implicit none
    integer::n,lwa,ik
    parameter (n=nv, lwa=n*(3*n+13))
    real(kind=kind(0.d0))::wa(lwa),tol,fvec(n)
    real(kind=kind(0.d0)),intent(inout)::x(n)
    integer,intent(out)::info
    external FCN
    tol=1d-8;
    call hybrd1(FCN,n,x,fvec,tol,info,wa,lwa)
    do ik=1,nv
       if (x(ik)<0) info=5
    enddo
    !! INFO = 0  Improper input parameters.
    !! INFO = 1  Algorithm estimates that the relative error between X and the solution is at most TOL.
    !! INFO = 2  Number of calls to FCN has reached or exceeded 200 *(N+1).
    !! INFO = 3  TOL is too small.  No further improvement in the approximate solution X is possible.
    !! INFO = 4  Iteration is not making good progress.
  end subroutine syssol

  !! compute the steady state quantities from the concentration solution
  subroutine flux(xx,flu,dg,ent,ente,ibib) ! state -> flux properties
    !! ente is a three dimensional array containing Tsig,Tsig^I,Tsig^E
    integer,intent(in)::ibib
    real(kind=kind(0.d0)),intent(in)::xx(nv)
    real(kind=kind(0.d0)),intent(out)::flu(nr),ent(nr),ente(3),dg(nr)
    real(kind=kind(0.d0))::tmp,knorm
    integer::j,ki
    !! Tsig_r=-J_r DG_r
    !! DG=-RT*ln(P/S)=-RT*ln(J^+/J^-)

    ente=0D0;
    do j=1,nc !exchanged species
       flu(j)=kcat(j)*(yy(j)-xx(j));
       dg(j)=RT*log(xx(j)/yy(j));
       ent(j)=-flu(j)*dg(j);
    enddo

    do j=nc+1,nr !internal species
       if (kn==0) knorm=1D0
       if (kn==1)  knorm=sqrt(keq(j));
       
       if (nm==1) then
          flu(j)=kcat(j)*knorm*(xx(add(j,2))-xx(add(j,1))/keq(j))
          dg(j)=RT*log(xx(add(j,1))/keq(j)/xx(add(j,2))); !ln(S/P)
          ent(j)=-flu(j)*dg(j)
       else !amender
          flu(j)=kcat(j)*knorm*(xx(add(j,3))*xx(add(j,4))-xx(add(j,1))*xx(add(j,2))/keq(j));
          dg(j)=RT*log(xx(add(j,1))*xx(add(j,2))/keq(j)/xx(add(j,3))/xx(add(j,4)));
          ent(j)=-flu(j)*dg(j)
       endif
    enddo
    
    ente(1)=sum(ent(1:nr)); !total EPR
    ente(2)=sum(ent(1:nc)); !exchanged EPR
    ente(3)=ente(1)-ente(2) !internal EPR

    if (ibib==1) then  
       write(*,'(A,15f8.3)') ' J= ',flu(:)
       write(*,'(A,15f8.3)') ' x= ',xx
       write(*,'(A,15f8.3)') 'DG= ',dg
       write(*,'(A,15f8.3)') ' S= ',ent
    endif
    
  end subroutine flux
  
end module module_steady_state

!!!!!!!!!!!!!!! MAIN PROGRAM!!!!!!!!!!!!!!!!!

program principal 
  use glob; use module_steady_state; implicit none;

  real(kind=kind(0.d0)),dimension(nc)::muy0,chiyy,chimu,muy

  real(kind=kind(0.d0))::mu1,mu2,emux,smax
  real(kind=kind(0.d0)),dimension(nv)::mun,mu0,vo,xx,xx1,xx2,v1,v2
  real(kind=kind(0.d0)),dimension(nr)::ent,dg,dgs,ents
  
  integer::ipop,it,i1,i2,i3,idid,n0,ik,inf,nd(6)
  integer::nn(4),nnn(nr,4)
  real(kind=kind(0.d0))::r0,ti,tf,ente(4),entes(4),ff(nr)
  real(kind=kind(0.d0))::entex,flu(nr),dis(51)

  
  !! Compute the distribution of EPR (Ts) over npop (FIG.2B and FIG.S1)
  open(12,file='fig2S1_dis.dat');

  !! Compute thermodynamic variables versus Tsig over a subsample npop=10^4 (FIG.2C)
  open(33,file='fig2C-mub_nm2.dat'); !exchanged species only
  open(34,file='fig2C-mu_nm2.dat'); !all species
  open(38,file='fig2C-reac_nm2.dat'); !all species
  open(37,file='fig2C-reacb_nm2.dat'); !exchange species only
  open(39,file='fig2C-mean_nm2.dat'); 

  dis=0.;
  do ipop=1,npop !loop on the population sample
9    stoch=0;
     
     !! generate a stoiciometry matrix
     call stochio(nm,nnn,nd);

     !! exclude the case of no emergent cycle (EQUILIBRIUM WITH CHEMOSTAT j1=j2=0)
     if ((nd(1)-nd(3))==0) goto 9  
     
     !! exclude connection between exchange species to avoid oversampling of trivial topologies
     if (conn(1,2)>0.5) goto 9
     
     !! kinetic parameters from random uniform distribution within a range (default 10^[-2:1])
     call random_number(kcat); kcat=10**(3.*(kcat)-2.);  !!DEFAULT RANGE
     !call random_number(kcat); kcat=10**(3.*(kcat)-1);  !!RANGE USED IN FIG.S1D_bottom
     !call random_number(kcat); kcat=10**(3.*(kcat)-3);  !!RANGE USED IN FIG.S1D_bottom
     !call random_number(kcat); kcat=10**(5.*(kcat)-3.); !!RANGE USED IN FIG.S1D_top
     kcat(1:nc)=1D0;

     !! thermodynamics parameters mu° from random uniform distribution within a range (default [-6:6]=
     call random_number(mu0); mu0=12D0*(mu0-0.5);  !!DEFAULT RANGE
     !call random_number(mu0); mu0=24D0*(mu0-0.5); !!RANGE USED IN FIG.S1C_top
     !call random_number(mu0); mu0=0.D0*(mu0-0.5); !!!RANGE USED IN FIG.S1C_bottom
   
     !! chemostat parameters and computation of TS^max (eq.11)
     muy0=(/5.,-5./); yy=1D0; mu0(1:nc)=muy0(1:nc); muy=muy0+RT*log(yy);
     chiyy(1:nc)=kcat(1:nc)*yy(1:nc); chimu(1:nc)=kcat(1:nc)*exp(-muy0(1:nc)/RT); 
     emux=sum(chiyy)/sum(chimu);
     smax=sum(kcat(1:nc)*muy(1:nc)*exp(-muy0(1:nc)/RT)*(exp(muy(1:nc)/RT)-emux)); ! EQ.11

     !computation of DG0 and equilibrium constants
     do i1=nc+1,nr 
        if (nm==1) then !global variable
           dg0(i1)=mu0(nnn(i1,1))-mu0(nnn(i1,2));
        else
           dg0(i1)=mu0(nnn(i1,1))+mu0(nnn(i1,2))-mu0(nnn(i1,3))-mu0(nnn(i1,4))
        endif
        keq(i1)=exp(-dg0(i1)/RT);
     enddo
  
     !Find steady state concentration/potential solution
     xx=1.*exp(-mu0/RT); call syssol(xx,inf); xx2=xx;  call sys_rand(xx2,v2,0D0);
     !if syssol does not converge backup method with ODE integration
     if (inf.ne.1) then
        ti=0D0; tf=1D7; xx=exp(-mu0/RT); call evolSEU(ti,tf,xx,idid,1.d-6);
        xx1=xx; call sys_rand(xx1,v1,0D0);
     endif

     ! compute flux solutions
     call sys_rand(xx,ff,0D0);
     call flux(xx,flu,dg,ent,ente,0);
     entex=ente(1);
     mun=mu0+RT*log(xx);

     !distribution of steady state
     ik=int(50.*entex/smax)+1;
     if (ik>0.5) dis(ik)=dis(ik)+1;
     if (ik>50) dis(50)=dis(50)+1;

     !distribution of thermodynamic variables
     do i2=nc+1,nr
        write(38,*) entex,ent(i2),abs(dg(i2)),abs(flu(i2))
     enddo
     do i2=1,nc 
        write(37,*) entex,ent(i2)/smax,abs(dg(i2)),abs(flu(i2))
     enddo
     write(39,*) entex,ente(3),mun(1)-mun(2)
     do i3=1,nv !metabolites
        if (abs(mun(i3))<1D-5) goto 8
        write(34,*) entex,mun(i3),mu0(i3),xx(i3)
8       continue
     enddo
     do i3=1,nc !metabolites
        write(33,*) entex,mun(i3),xx(i3)*exp(mu0(i3)/RT)
     enddo
  enddo
  continue
  write(*,*) 'fin',smax
  do i1=1,51
     write(12,*) smax*real(i1-1)/5D1,5D1*dis(i1)/sum(dis(:))/smax
  enddo
  
contains

  !! Singular value decomposition function
  !! used to evaluate kernel and cokernel properties
  
  subroutine SVD(n1,n2,A,IT,nulldim) 
    implicit none
    integer,intent(in)::n1,n2,it
    real(kind=kind(0.d0)),intent(in)::A(n1,n2)
    real(kind=kind(0.d0))::S(n2)!,EV(nv,nr)
    integer,intent(out)::nulldim
    !integer:: LDA, LDU, LDVT
    integer,parameter:: LWMAX=1000
    integer::INFO,LWORK,i,i1
    real(kind=kind(0.d0))::U(n1,n1),VT(n2,n2),WORK(LWMAX)
    !LDA=M; LDU=M; LDVT=N
    LWORK=-1; VT=0.; S=0;
    call DGESVD('A','A',n1,n2,A,n1,S,U,n1,VT,n2,WORK,LWORK,INFO) 
    LWORK=MIN(LWMAX,INT(WORK(1)))
    call DGESVD('A','A',n1,n2,A,n1,S,U,n1,VT,n2,WORK,LWORK,INFO)
    !info>0 The SVD algorithm failed to converge
    nulldim=0; i1=0; 
    do i=1,n2 !n2 valeurs propres
       if (abs(S(i))<1d-2) then ! S(i) in singular value
          i1=i1+1; nulldim=nulldim+1;
       endif
    enddo
  end subroutine SVD

  subroutine jacobian(elas1,jac1) !valable uniquement pour uni=2
    integer::i1
    real(kind=kind(0.d0)),intent(out),dimension(nr,nv)::elas1
    real(kind=kind(0.d0)),intent(out),dimension(nv,nv)::jac1
    elas1=0D0; jac1=0D0;
    do i1=1,nc
       elas1(i1,i1)=-kcat(i1);
    enddo
    do i1=nc+1,nr
       elas1(i1,add(i1,1))=-kcat(i1);
       elas1(i1,add(i1,2))=kcat(i1)/keq(i1); !XXX-verif
    enddo
    jac1=matmul(stoch,elas1);
  end subroutine jacobian

  !! Function to generate random stoichiometric matrix
  subroutine stochio(uni,nnn,nd)
    integer,intent(in)::uni
    integer,intent(out)::nnn(nr,4),nd(6)
    integer::nl,nd2,nl2
    real(kind=kind(0.d0))::sto1(nv,nr)
    real(kind=kind(0.d0))::r0
    integer::nn(4),i1,i2,it,n0,ntmp(2),iz
    real(kind=kind(0.d0)),dimension(nr,nv)::wstocht
    real(kind=kind(0.d0)),dimension(nv+nc,nr)::wstoch2 !introduce chemostat 
    real(kind=kind(0.d0)),dimension(nr,nv+nc)::wstoch2t !transpose

    stoch=0; stoch2=0.; stocht=0.; stoch2t=0.;
    conn=0; nnn=0; 
    do i1=nc+1,nr
3      continue
       i2=1; nn=0;     
       do while (i2<(real(2*uni)+0.5))
1         continue;
          call random_number(r0)
          n0=1+int(r0*real(nv)); !random choice of an index
          if (i2>1.5) then 
             do iz=1,i2-1
                if (n0==nn(iz)) goto 1 !do not use the same metabolite
             enddo
          endif
          nn(i2)=n0; i2=i2+1;
       enddo
       
       nnn(i1,:)=nn;
       do iz=nc+1,i1-1
          if ((nn(1).eq.nnn(iz,1)).and.(nn(2).eq.nnn(iz,2))) goto 3 
          if ((nn(1).eq.nnn(iz,2)).and.(nn(2).eq.nnn(iz,1))) goto 3 
          if (uni==2) then
             if ((nn(3).eq.nnn(iz,3)).and.(nn(4).eq.nnn(iz,4))) goto 3 
             if ((nn(4).eq.nnn(iz,3)).and.(nn(3).eq.nnn(iz,4))) goto 3 
          endif
       enddo
       if (uni==1) then
          it=nn(1); stoch(it,i1)=1.; it=nn(2); stoch(it,i1)=-1.; add(i1,1:2)=nn(1:2);
          conn(nn(1),nn(2))=i1;   conn(nn(2),nn(1))=i1;
       else 
          it=nn(1); stoch(it,i1)=1.; it=nn(3); stoch(it,i1)=-1.;
          it=nn(2); stoch(it,i1)=1.; it=nn(4); stoch(it,i1)=-1.;
          add(i1,1:4)=nn(1:4); 
          conn(nn(1),nn(3))=i1;conn(nn(3),nn(1))=i1;conn(nn(1),nn(4))=i1;conn(nn(4),nn(1))=i1;
          conn(nn(2),nn(4))=i1;conn(nn(4),nn(2))=i1;conn(nn(2),nn(3))=i1;conn(nn(3),nn(2))=i1;
          conn(nn(2),nn(1))=i1;conn(nn(4),nn(3))=i1;
       endif
    enddo !end nr reaction
    do i1=1,nc !exchanged reactions
       add(i1,1)=i1; stoch(i1,i1)=1;
       stoch2(nv+i1,i1)=-1.;
    enddo ! stoch(nv,nr)

    !
    sto1=stoch; stoch2(1:nv,:)=stoch(1:nv,:); 
    stoch2t=transpose(stoch2); stocht(:,1:nv)=stoch2t(:,1:nv);
    wstoch2=stoch2; wstoch2t=stoch2t; wstocht=stocht;
    nl=0.; call SVD(nv+nc,nr,wstoch2,0,nl); nd(3)=nl; !cycles
    nl=0.; call SVD(nr,nv+nc,wstoch2t,0,nl); nd(4)=nl; !conservation
    nl=0;; call SVD(nv,nr,sto1,0,nl);  nd(1)=nl ! nd(1)-nd(3)=affinities
    nl=0.; call SVD(nr,nv,wstocht,0,nl); nd(2)=nl
    nd(5)=nv+2-nd(4)-(nr-nd(3)); nd(6)=nc-(nd(1)-nd(3)+nd(4)-nd(2))
    !'nd=',nd,   s-l     =  r-c   ;nc= [nd(1)-c]+ [l-nd(2)}
  end subroutine stochio

end program principal



!!!!!!!!!!!EXTERNAL SUBROUTINES FOR SEULEX !!!!!!!!!!!!!!!

subroutine solout (NRR,XOLD,X,Y,RC,LRC,IC,LIC,N,RPAR,IPAR,IRTRN)
  use glob
  implicit none
  integer :: N,NRR,LRC,LIC,IRTRN,ipar                               
  real(kind=kind(0.d0)) :: X,Y(N),RC(LRC),IC(LIC),rpar,XOLD
  real(4) :: t,val,diff
  real(kind=kind(0.d0)) :: F(N)

!!! temporal evolution of the ODE system
  if (NRR.LE.nplot) then
     tplot(NRR)=X;
     xplot(:,NRR)=Y(:);
  endif
  nmax=NRR
  !  if ((X>0).and.(Y(nv).LT.0.5)) then !! --- CRITERE D'ARRET de evol_SEU---!!YO!!
  !     Y=-1; IRTRN=-1; return;
  ! endif
  return
end subroutine solout

!!$ ---------------------------------------
!!$   Routines required for communication with numerical solver/integration routines
!!$----------------------------------------

SUBROUTINE FCN(n,x,fvec,iflag) ! used by sys_sol
  use glob; implicit none;
  integer::iflag,n
  real(kind=kind(0.d0))::t,x(n),fvec(n)
  call sys_rand(x,fvec,t) 
end SUBROUTINE FCN

SUBROUTINE FCN_HW(N,t,xHW,F,rpar,ipar) ! used by seulex
  use glob
  implicit none
  integer :: n,ipar
  real(kind=kind(0.d0)) :: t,xHW(n),f(n),rpar
  call sys_rand(XHW,f,t) !
end SUBROUTINE FCN_HW

subroutine jac_HW(N,X,Y,DFY,LDFY,RPAR,IPAR) ! used by seulex
  use glob
  implicit none
  integer :: N,LDFY,LRC,LIC,IRTRN,ipar                                       
  real(kind=kind(0.d0)) :: X,Y(N),DFY(LDFY,N),rpar
  !call jacobien(Y,DFY,X)
  return
end subroutine jac_HW

subroutine mass
  return
end subroutine mass

