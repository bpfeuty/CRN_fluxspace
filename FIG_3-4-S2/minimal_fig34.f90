!! gfortran sampling_fig34_nm2.f90 -lseulex2 -lminpack -llapack

module GLOB ! DEFINITION OF THE VARIABLE AND NETWORK
  implicit none
  
  integer,parameter::ns=5,nr=5; nm=2,nc=3  ! hyperparameters
  real(kind=kind(0.d0)),parameter::RT=2.5

  !! ns is the number of species (=8 in default simulations but varied in FIGS2)
  !! nr is the number of reactions (=12 in default simulations but varied in FIGS2)
  !! nm is the molecularity of reaction (=2 except panels B and some distributions).
  !! nc=3 is the number of chemostat 
  
  integer,parameter::nplot=500000; 
  real(kind=kind(0.d0)),dimension(ns,nplot) :: xplot(ns,nplot),tplot(nplot);    
  real(kind=kind(0.d0)),dimension(ns,nr)::stoch
  real(kind=kind(0.d0)),dimension(nr,ns)::stocht
  real(kind=kind(0.d0)),dimension(ns+nc,nr)::stoch2 
  real(kind=kind(0.d0)),dimension(nr,ns+nc)::stoch2t 
  real(kind=kind(0.d0)),dimension(ns,ns)::jac   
  integer,dimension(ns,ns)::conn
  integer,dimension(nr,4)::add 
  integer::nmax,nr_max
  real(kind=kind(0.d0))::yy(nc),kcat(nr),keq(nr)
  
contains

  subroutine sys_rand(v,dvdt,t) ! ODE system
    
    real(kind=kind(0.d0)),dimension(ns),intent(in)::v
    real(kind=kind(0.d0)),dimension(ns),intent(out)::dvdt
    real(kind=kind(0.d0)),intent(in)::t
    real(kind=kind(0.d0))::reac(nr),knorm
    integer::i,j,ki
    
    dvdt=0.;
    do i=1,nc
       dvdt(i)=dvdt(i)+kcat(i)*(yy(i)-v(i)); 
    enddo
    
    do i=1,ns
       do j=nc+1,nr; 
          knorm=sqrt(keq(j))
          if (nm==1) then
             reac(j)=kcat(j)*knorm*(v(add(j,2))-v(add(j,1))/keq(j)); !1:prod 2:sub!
          else
             reac(j)=kcat(j)*knorm*(v(add(j,3))*v(add(j,4))-v(add(j,1))*v(add(j,2))/keq(j));
          endif
          dvdt(i)=reac(j)*stoch(i,j)+dvdt(i);
       enddo
    enddo
  end subroutine sys_rand
  
  subroutine evolSEU(tbeg,tend,X,idid,tol) !! numerical integration
    real(kind=kind(0.d0)), intent(in) :: tol
    real(kind=kind(0.d0))             :: tbeg,tend
    real(kind=kind(0.d0)), intent(inout) :: X(ns)
    integer, intent(out) :: idid
    real(kind=kind(0.d0)),save::h
    integer,parameter :: lrcont=20,licont=10
    integer :: itol,ijac,ifcn,iout
    integer :: mljac=ns,mujac=ns,imas=0
    integer :: MUMAS=ns,MLMAS=ns
    integer, parameter :: lwork=2*(ns+10)*ns+61+20
    real(kind=kind(0.d0)) :: work(lwork),rpar
    integer, parameter :: liwork=2*ns+32+20
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
    call SEULEX(ns,FCN_HW,IFCN,tbeg,X,tend,H,&
         &RTOL,ATOL,ITOL,&
         &JAC_HW ,IJAC,MLJAC,MUJAC,&
         &MASS ,IMAS,MLMAS,MUMAS,&
         &SOLOUT,IOUT,&
         &WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
    nr_max=min(iwork(17),nplot) !! nombre de pas
  end subroutine evolSEU
end  module GLOB

!!!!!!!!!!!!!!!!!!!!!!!!!!! PROGRAMME PRINCIPAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program principal !!REF
  use glob; implicit none;

  !dn is the number of emergent cycles
  integer::ik,i0,ipop,it,i1,i2,i3,ifig   !loop integer
  integer::ib,idid                       ! warning integer
  integer::nsamp,nk,nec,ncc,nd(6),nn(4),nnn(nr,4) 

  !diverse system variable and parameters!
  !chemical potential mu, concentration xx, flux J, entropy production rate sig, gibbs free energy dg
  real(kind=kind(0.d0)),dimension(nc)::muy,chiyy,chimu,muy0
  real(kind=kind(0.d0)),dimension(ns)::xx,xxs,mu,mu0
  real(kind=kind(0.d0)),dimension(nr)::ent,dg,dg0  ! sortie de call flux et SVD
  real(kind=kind(0.d0))::ti,tf,ente(3),ff(nr),sigtot,flu(nr)
  real(kind=kind(0.d0))::emux,sigmax ! for computing  maximum entropy sig_max

  !specific outputs for determining optimal topologies
  real(kind=kind(0.d0))::maxj,disj1(50),disj2(50)    !store maximal values of specific fluxes and their distribution (FIG3)
  real(kind=kind(0.d0))::vmax,disj(12,12),disurf(100) !feasible space volume and its distribution (FIG4)

  open(13,file='fig3g.dat');
  open(14,file='fig4g.dat');  

  !i0 from 0 to 7 corresponds to a given hypeparameter associated to a specific figure (3 or 4) and panel
  !nk is the number of kinetic sampling for one given topology
  !nk=500 is used to test the sampled space of given topologies
  !nk=1 is used to sample together topologies and kinetics
  nsamp=1000; nk=2000; 

  
  ifig=3;  !for figure 3 nr=5 is selected (as the minimal number of reaction achieving maximum J3ex)
  !ifig=4; !for figure 4 nr=7 is selected (as the minimal number of reaction achieving maximum Vmax)
  !though not necessary range of chemical has been set to mu=0 for reproducibility of optimal state.
  
  !Chemostatting 
  muy0=0D0;
  if (ifig==3) then
     muy0=(/2.,-2.,3./)*RT; yy=(/1.,1.,0.5/); kcat(1:nc)=1d0; 
  endif
  if (ifig==4) then 
     muy0=(/1.,1.,-1./)*RT;  yy=(/1.,1.,1./); kcat(1:nc)=1d0; 
  endif
  muy=muy0+RT*log(yy);

  !Equations 10 and 11 to comptute mu^* and Ts^max
  chiyy(1:nc)=kcat(1:nc)*yy(1:nc); chimu(1:nc)=kcat(1:nc)*exp(-muy0(1:nc)/RT); 
  emux=sum(chiyy)/sum(chimu);
  sigmax=sum(kcat(1:nc)*muy(1:nc)*exp(-muy0(1:nc)/RT)*(exp(muy(1:nc)/RT)-emux)); 
  write(*,*) 'sig_max=',sigmax

  disj1=0D0;disj2=0D0; disurf=0.;
  do ipop=1,nsamp
     stoch=0.;
     !generate a stoichiometric matrix of multimolecularity
     !! (output are the mapping function nnn and the dimensionality of kernel and cokernel)
9    call stochio(nm,nnn,nd); 

     !! Selection of topologies with specific dimensionalities of the kernel
     !! (nbc, nec and ncc keeping in mind that nbc+nec=nc)
     ncc=nd(3); nec=nd(1)-ncc;  !ncc and nec are the number of closed and emergent reaction cycles
     !if (conn(1,2)>0.5) goto 9 !exclude connections between exchange species
     if (nec==0) goto 9

     !! sampling of thermodynamic parameters
     call random_number(mu0);
     !mu0=12D0*(mu0-0.5); !random chemical potential
     mu0=0D0*(mu0-0.5); !random chemical potential
     mu0(1:nc)=muy0(1:nc); 
     do i1=nc+1,nr ! from mu0 -> dg -> keq 
        if (nm==1) then 
           dg0(i1)=mu0(nnn(i1,1))-mu0(nnn(i1,2));
        else
           dg0(i1)=mu0(nnn(i1,1))+mu0(nnn(i1,2))-mu0(nnn(i1,3))-mu0(nnn(i1,4))
        endif
        keq(i1)=exp(-dg0(i1)/RT); 
     enddo
     write(*,*) ipop
     !kinetic sampling
     maxj=-1.;  disj=0D0
     do ik=1,nk 
        call random_number(kcat); kcat=10**(4.*kcat-2); 
        kcat(1:nc)=1d0;
        ! simulations
        ti=0D0; tf=1D4; xx=1.*exp(-mu0/RT);
        call evolSEU(ti,tf,xx,idid,1.d-6);
        call syssol(xx); !cross-check
        call sys_rand(xx,ff,0D0); mu=mu0+RT*log(xx); 
        ib=0;
        call flux(xx,flu,dg,ent,ente,ib);
        sigtot=ente(1);
        if (ib==2) goto 9
        write(10+ifig,*) sigtot/sigmax,flu(1:nc),ente(3)/ente(1),mu(4)-mu(5),mu(5)-mu(4)
        if (-flu(3)>maxj) maxj=-flu(3);
        disj(int((flu(1)+0.3)*1D1+1),int((flu(2)+0.3)*1D1)+1)=1.;
     enddo
     vmax=sum(disj); vmax=vmax/1D2/0.96;
     ! exit when the optimal state is sampled
     if (ifig==3) then
        if (maxj>0.37) goto 90
     endif
     if (ifig==4) then
        if (vmax>0.8) goto 90
     endif
     !distribution of flux space volume figure fig4E         
     write(*,*) ipop,vmax,maxj,nc-nec
     rewind(10+ifig)
  enddo !ipop loop
90 continue

  if (ifig==3) then
     write(*,*) '-J3max=',maxj
     write(*,*) 'nbc=',nc-nec
     !write(*,*) 'mu=',mu
     !write(*,*) 'J=',flu
  endif
  if (ifig==4) then
     write(*,*) 'Vmax=',vmax
  endif

contains

  subroutine stochio(nm,nnn,nd)
    integer,intent(in)::nm
    integer,intent(out)::nnn(nr,4),nd(6)
    integer::nl,nd2,nl2
    real(kind=kind(0.d0))::sto1(ns,nr)
    !real(kind=kind(0.d0)),intent(out)::ev(ns,nr)
    real(kind=kind(0.d0))::r0!
    integer::nn(4),i1,i2,it,n0,ntmp(2),iz
    real(kind=kind(0.d0)),dimension(nr,ns)::wstocht
    real(kind=kind(0.d0)),dimension(ns+nc,nr)::wstoch2 !introduce chemostat 
    real(kind=kind(0.d0)),dimension(nr,ns+nc)::wstoch2t !transpose
    stoch=0; stoch2=0.; stocht=0.; stoch2t=0.;
    conn=0; nnn=0; 
    do i1=nc+1,nr
3      continue
       i2=1; nn=0;     
       do while (i2<(real(nm*2)+0.5))
1         continue; call random_number(r0)
          n0=1+int(r0*real(ns)); !nb entier
          if (i2>1.5) then
             do iz=1,i2-1
                if (n0==nn(iz)) goto 1 !evite reutilise meme metabolite
             enddo
          endif
          nn(i2)=n0; i2=i2+1;
       enddo
       nnn(i1,:)=nn;
       do iz=nc+1,i1-1
          if ((nn(1).eq.nnn(iz,1)).and.(nn(2).eq.nnn(iz,2))) goto 3 !! 2 fois meme reac
          if ((nn(1).eq.nnn(iz,2)).and.(nn(2).eq.nnn(iz,1))) goto 3 !! 
          if (nm==2) then
             if ((nn(3).eq.nnn(iz,3)).and.(nn(4).eq.nnn(iz,4))) goto 3 !! bibi cas ?? xx
             if ((nn(4).eq.nnn(iz,3)).and.(nn(3).eq.nnn(iz,4))) goto 3 !! bibi cas ?? xx
          endif
       enddo
       if (nm==1) then
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
       stoch2(ns+i1,i1)=-1.;
    enddo ! stoch(ns,nr)
    sto1=stoch; stoch2(1:ns,:)=stoch(1:ns,:); 
    stoch2t=transpose(stoch2); stocht(:,1:ns)=stoch2t(:,1:ns);
    wstoch2=stoch2; wstoch2t=stoch2t; wstocht=stocht;
    nl=0.; call SVD(ns+nc,nr,wstoch2,0,nl); nd(3)=nl;  !cycles
    nl=0.; call SVD(nr,ns+nc,wstoch2t,0,nl); nd(4)=nl; !conservation
    nl=0;; call SVD(ns,nr,sto1,0,nl);  nd(1)=nl        ! nd(1)-nd(3)=affinities
    nl=0.; call SVD(nr,ns,wstocht,0,nl); nd(2)=nl      
    nd(5)=ns+2-nd(4)-(nr-nd(3)); !checks
    nd(6)=nc-(nd(1)-nd(3)+nd(4)-nd(2)); !checks
    !'nd=',nd,   s-l     =  r-c   ;nc= [nd(1)-c]+ [l-nd(2)}
  end subroutine stochio

  subroutine syssol(x) !compute jacobian analytically
    use glob;  implicit none
    integer::n,info,lwa
    parameter (n=ns, lwa=n*(3*n+13))
    real(kind=kind(0.d0))::wa(lwa),tol,fvec(n)
    real(kind=kind(0.d0))::x(n)
    external FCN
    tol=1d-8;
    call hybrd1(FCN,n,x,fvec,tol,info,wa,lwa)
    ! HYBRD - HYBRD1 HYBRJ HYBRJ1
    ! subroutine hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)
  end subroutine syssol

  subroutine flux(xx,flu,dg,ent,ente,ibib) ! state -> flux properties
    integer,intent(inout)::ibib
    real(kind=kind(0.d0)),intent(in)::xx(ns)
    real(kind=kind(0.d0)),intent(out)::flu(nr),ent(nr),ente(3),dg(nr)
    real(kind=kind(0.d0))::tmp,knorm,ktm(nr)
    integer::j,ki

    ! we use the following formula for computing DG (gibbs free energy) and sigma (entropy production rate)
    ! sigma_r = -J_r DG_r = J_r log(J_r^+/J_r^-)
    ! DG_r = -RT*log(J_r^+/J_^-)

    ente=0D0; 
    do j=1,nc
       flu(j)=kcat(j)*(yy(j)-xx(j));
       dg(j)=RT*log(xx(j)/yy(j));
       ent(j)=-flu(j)*dg(j);
    enddo
    do j=nc+1,nr
       knorm=sqrt(keq(j))
       if (nm==1) then
          flu(j)=kcat(j)*knorm*(xx(add(j,2))-xx(add(j,1))/keq(j))
          dg(j)=RT*log(xx(add(j,1))/keq(j)/xx(add(j,2)));
          ent(j)=-flu(j)*dg(j)
       else 
          flu(j)=kcat(j)*knorm*(xx(add(j,3))*xx(add(j,4))-xx(add(j,1))*xx(add(j,2))/keq(j));
          dg(j)=RT*log(xx(add(j,1))*xx(add(j,2))/keq(j)/xx(add(j,3))/xx(add(j,4)));
          ent(j)=-flu(j)*dg(j)
          if (ent(j)<0) ibib=2; !error detection
       endif
    enddo
    ente(1)=sum(ent(1:nr)); !s_tot
    ente(2)=sum(ent(1:nc)); !s_ext
    ente(3)=ente(1)-ente(2) !_sint
    if (ibib==1) then 
       write(*,'(A,15f8.3)') ' J= ',flu(:)
       write(*,'(A,15f8.3)') ' x= ',xx
       write(*,'(A,15f8.3)') ' DG=',dg
       write(*,'(A,15f8.3)') ' S= ',ent
    endif
  end subroutine flux

end program principal

!!!!!!EXTERNAL SUBROUTINES FOR SEULEX !!!!!

subroutine solout (NRR,XOLD,X,Y,RC,LRC,IC,LIC,N,RPAR,IPAR,IRTRN)
  use glob
  implicit none
  integer :: N,NRR,LRC,LIC,IRTRN,ipar                               
  real(kind=kind(0.d0)) :: X,Y(N),RC(LRC),IC(LIC),rpar,XOLD
  real(4) :: t,val,diff
  real(kind=kind(0.d0)) :: F(N)

  !!! stocke l'évolution temporelle dans un tableau
  if (NRR.LE.nplot) then
     tplot(NRR)=X;
     xplot(:,NRR)=Y(:);
  endif
  nmax=NRR
  !  if ((X>0).and.(Y(ns).LT.0.5)) then !! --- CRITERE D'ARRET de evol_SEU---!!YO!!
  !     Y=-1; IRTRN=-1; return;
  ! endif
  return
end subroutine solout

!!$ ---------------------------------------
!!$   Routines nécessaires à la communication avec 
!!$   l'intégrateur numérique
!!$----------------------------------------

SUBROUTINE FCN(n,x,fvec,iflag) ! utilisé par seulex 
  use glob; implicit none;
  integer::iflag,n
  real(kind=kind(0.d0))::t,x(n),fvec(n)
  call sys_rand(x,fvec,t) 
end SUBROUTINE FCN

SUBROUTINE FCN_HW(N,t,xHW,F,rpar,ipar) ! ??
  use glob
  implicit none
  integer :: n,ipar
  real(kind=kind(0.d0)) :: t,xHW(n),f(n),rpar
  call sys_rand(XHW,f,t) !
end SUBROUTINE FCN_HW

subroutine jac_HW(N,X,Y,DFY,LDFY,RPAR,IPAR) ! utilisé par seulex
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

subroutine SVD(n1,n2,A,IT,nulldim) !travaille sur la stoch
    implicit none
    integer,intent(in)::n1,n2,it
    real(kind=kind(0.d0)),intent(in)::A(n1,n2)
    real(kind=kind(0.d0))::S(n2)!,EV(ns,nr)
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
