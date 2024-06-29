! coarse grained model of carbohydrate metabolism
! fortran code that generates the figure 6cd
! compile with gfortran -llapack -lminpack -lseulex2  
! compilation time is about minutes

! important parameters
! nv=10 metabolite species (nc=7 exchanged and 3 internal)
! nr=12 metabolic reactions
! PT is the total concentration of ATP+ADP=10 mM
! num is the ,numericla scheme to compute solutions
! num=0 combine seulex and minpack for sol
! num=1 seulex
! num=2 minpack
! seulex method ensure that steady state solution are stable
! mu0 and DG0 values are used to compute keq in mass-action kinetics

! OPTIMIZATION
! perform optimization through random optimization with nsamp=500


module GLOB 
  implicit none
  integer,parameter::nv=10,nr=12,nc=7,num=0,nplot=500000;
  real(8)::RT=2.5D0,PT=0.01D0 !some constants

  !chemostat parameters
  real(kind=kind(0.d0)),dimension(nc)::muy,chiyy,chimu,muy0,yy

  !stoichiometric matrix
  real(kind=kind(0.d0)),dimension(nv,nr)::stochx(nv,nr),stochxt(nr,nv),stoch(nv+nc,nr),stocht(nr,nv+nc)
  real(kind=kind(0.d0))::mp(nv+nc,nr),mpt(nr,nv+nc),mpx(nv,nr),mptx(nr,nv); ! stoichiometric matrix

  !diverse global parameters/variables
  integer::nb(3),nbb(3),nr_max,nmax
  real(kind=kind(0.d0)),dimension(nv,nv)::jac
  real(kind=kind(0.d0)),dimension(nv,nplot) :: xplot(nv,nplot),tplot(nplot);
  real(kind=kind(0.d0)),dimension(nv)::mun,mu0 !chemical potentials
  real(kind=kind(0.d0)),dimension(nr)::kcat,keq,dg0,dg,ent,kcats 
  
contains

   subroutine sys_rand(v,dvdt,t) ! ODE system used by seulex and syssol
    real(kind=kind(0.d0)),dimension(nv-1),intent(in)::v
    real(kind=kind(0.d0)),dimension(nv-1),intent(out)::dvdt
    real(kind=kind(0.d0)),intent(in)::t
    real(kind=kind(0.d0))::reac(nr),tmp,tmp2
    integer::i,j
    dvdt=0.;
    do i=1,nc
       reac(i)=kcat(i)*(yy(i)-v(i)) 
       dvdt(i)=dvdt(i)+reac(i); ! exchange reactions
    enddo
    tmp=PT-v(9); tmp2=v(7);
    reac(8)=v(1)*(v(9)/PT)**nbb(1)-v(8)**2*v(4)**2*(tmp/PT*v(5))**nbb(1)*tmp2**8/keq(8); !gly
    reac(9)=v(8)*v(9)/PT-v(2)*tmp*tmp2/PT/keq(9); !pka
    reac(10)=v(8)*v(6)**3*tmp2**4*(v(9)/PT)**nbb(2)-v(4)**2*v(5)**3*(tmp*v(5)/PT)**nbb(2)/keq(10); !ares
    reac(11)=(tmp*v(5)/PT)**nbb(3)*v(8)-(v(9)/PT)**nbb(3)*v(3)/keq(11); !bs
    reac(12)=tmp/PT*v(5)-v(9)/PT/keq(12); !ah

    do i=1,nv-1
       do j=nc+1,nr; 
          dvdt(i)=kcat(j)*reac(j)*stochx(i,j)+dvdt(i);
       enddo
    enddo
  end subroutine sys_rand
  
  subroutine evolSEU(tbeg,tend,X,idid,tol) !! numerical integration
    real(kind=kind(0.d0)), intent(in) :: tol
    real(kind=kind(0.d0))             :: tbeg,tend
    real(kind=kind(0.d0)), intent(inout) :: X(nv-1)
    integer, intent(out) :: idid
    real(kind=kind(0.d0)),save::h
    integer,parameter :: lrcont=20,licont=10
    integer :: itol,ijac,ifcn,iout
    integer :: mljac=nv-1,mujac=nv-1,imas=0
    integer :: MUMAS=nv-1,MLMAS=nv-1
    integer, parameter :: lwork=2*(nv-1+10)*(nv-1)+61+20
    real(kind=kind(0.d0)) :: work(lwork),rpar
    integer, parameter :: liwork=2*(nv-1)+32+20
    integer :: iwork(liwork),ipar
    real(kind=kind(0.d0)) :: rtol,Atol
    external fcn_HW,solout,mass,jac_HW !! code externe utilisé
    itol=0; rtol=tol; idid=1; h=0D0;
    atol=1d-8;  !! précision abolue
    iout=1; 	!! appelle soulout
    ijac=0; 	!! difference finie pour le jacobien -> =1  utilise JAC
    ifcn=0;	!! F(X,Y) INDEPENDENT OF X (AUTONOMOUS) -> =1 nonautonomous
    work=0d0   !! work(2)=0.5D0 !! mettre en commentaire si adaptatif
    iwork=0
    iwork(1)=0 		!! pas de transfo hessenberg du jac -> =1 sinon
    iwork(2)=10**8 	!! nb de step max
    iwork(4)=3  	!! step size sequence
    call SEULEX(nv-1,FCN_HW,IFCN,tbeg,X,tend,H,&
         &RTOL,ATOL,ITOL,&
         &JAC_HW ,IJAC,MLJAC,MUJAC,&
         &MASS ,IMAS,MLMAS,MUMAS,&
         &SOLOUT,IOUT,&
         &WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
    nr_max=min(iwork(17),nplot) !! nombre de pas
  end subroutine evolSEU
end  module GLOB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_algebra   !! EXPM, EIGEN, JACSOL
  use GLOB; implicit none;
contains
  
  subroutine syssol(x) !compute jacobian analytically
    use glob;  implicit none
    integer::n,info,lwa
    parameter (n=nv-1, lwa=n*(3*n+13))
    real(kind=kind(0.d0))::wa(lwa),tol,fvec(n)
    real(kind=kind(0.d0))::x(n)
    external FCN
    tol=1d-15;
    call hybrd1(FCN,n,x,fvec,tol,info,wa,lwa)
    ! HYBRD - HYBRD1 HYBRJ HYBRJ1
    ! subroutine hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)
  end subroutine syssol

  subroutine flux(v,xfl,xdg,xsent) ! compute steady state quantities such as DG and Sent
    use glob; implicit none
    real(kind=kind(0.d0)),intent(in),dimension(nv-1)::v
    real(kind=kind(0.d0)),intent(out),dimension(nr)::xfl,xdg,xsent
    real(kind=kind(0.d0)),dimension(nr)::xflp,xflm
    real(kind=kind(0.d0))::tmp,tmp2
    integer::i
    xfl=0.; xflm=0.; xflp=0.;
    do i=1,nc
       xfl(i)=kcat(i)*(yy(i)-v(i))
       xdg(i)=-RT*log(yy(i)/v(i));
       xsent(i)=-xfl(i)*xdg(i)
       !if (kcat(i)<1D-6) xsent(i)=0D0
    enddo
    tmp=PT-v(9); tmp2=v(7)
    xflp(8)=v(1)*(v(9)/PT)**nbb(1) !glyc
    xflm(8)=tmp2**8*v(8)**2*v(4)**2*(v(5)*tmp/PT)**nbb(1)/keq(8)
    
    xflp(9)=v(8)*v(9)/PT !8:ferm
    xflm(9)=tmp2*v(2)*tmp/PT/keq(9);

    xflp(10)=tmp2**4*v(8)*v(6)**3*(v(9)/PT)**nbb(2) !9:resp
    xflm(10)=v(4)**2*v(5)**(3+nbb(2))*(tmp/PT)**nbb(2)/keq(10)
    
    xflp(11)=(tmp*v(5)/PT)**nbb(3)*v(8) !10:biomass
    xflm(11)=(v(9)/PT)**nbb(3)*v(3)/keq(11)

    xflp(12)=tmp*v(5)/PT
    xflm(12)=v(9)/keq(12)/PT; !ATPhydro consumes atp

    do i=nc+1,nr
       xfl(i)=kcat(i)*(xflp(i)-xflm(i));
       xdg(i)=-RT*log(xflp(i)/xflm(i));
       xsent(i)=-xfl(i)*xdg(i);
    enddo

  end  subroutine flux
end module module_algebra

!!!!!!!!!!!!!!!!!!!!!!!!!!! PROGRAMME PRINCIPAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program principal !!REF
  use glob; use module_algebra; implicit none;
  integer::it,i1,i2,i3,idid,dn,itx(4),isamp,nsamp
  real(kind=kind(0.d0))::thr,ti,tf,h_gly,h_ares,randd
  real(kind=kind(0.d0))::by,bys,eprs(2),sct,scts,etas
  real(kind=kind(0.d0)),dimension(nv-1)::xx,xeq,xxs
  real(kind=kind(0.d0)),dimension(nv)::mus
  integer,dimension(5,nv+nc)::bl
  real(kind=kind(0.d0)),dimension(nr)::flu,flus,ff,rrand,dgs
  integer,dimension(nr,5)::ec
  real(kind=kind(0.d0)),dimension(3)::entex 

  open(11,file='fig6d.dat');

  do i3=1,24 !stoch_gly,atp=0->23
     do i2=1,28 !stoch_ares,atp=0->56
        nsamp=500
        nb=(/10,10,10/);
        nb(1)=(i3-1)*1;
        nb(2)=(i2-1)*2;
        nbb=nb;
        !STOICHIOMETRIC MATRIX
        !species:gluc,ace,biom,co2,h2o,o2,acoa,adp,atp
        !reactions:gluc,ace,biom,co2,h2o,o2,glyc,ferm,resp,biom,atpc
        mp=0.; 
        do i1=1,6 
           mp(i1,i1)=1.; mp(nv+i1,i1)=-1.;
        enddo
        mp((/1,4,7,8/),8)=(/-1.,2.,8.,2./); mp((/5,9,10/),8)=nb(1)*(/1.,-1.,1./);
        mp((/2,7,8,9,10/),9)=(/1.,1.,-1.,-1.,1./);
        mp((/4,6,7,8/),10)=(/2.,-3.,-4.,-1./); mp((/9,10/),10)=nb(2)*(/-1.,1./); mp(5,10)=nb(2)+3; 
        mp((/3,8/),11)=(/1.,-1./);  mp((/5,9,10/),11)=nb(3)*(/-1.,1.,-1./);
        mp((/5,9,10/),12)=(/-1.,1.,-1./);  
        mpt=transpose(mp); mpx(1:nv,:)=mp(1:nv,:); mptx(:,1:nv)=mpt(:,1:nv); stoch=mp; stochx=mpx;

!!! algebraic analysis of stochiometric matrices
        if (0>1) then 
           call SVD(nv+nc,nr,mp,0,it); itx(1)=it; !total cycles
           call SVD(nr,nv+nc,mpt,0,it); itx(2)=it; 
           call SVD(nv,nr,mpx,0,it); itx(3)=it; !itx(1)-itx(3) gives affinities or emergent cycles
           call SVD(nr,nv,mptx,0,it); itx(4)=it;
           write(*,*) 'numbers of reaction cycles and conservation constraints'
           write(*,*) 'closed cycles=',itx(1)
           write(*,*) 'emergent cycle=',itx(3)-itx(1)
           write(*,*) 'broken conservation laws=',itx(2)-itx(4)
           write(*,*) 'unbroken conservation laws',itx(4)
           write(*,*) 'Dimension consistency (=0)',itx(3)-itx(1)+itx(2)-itx(4)-nc,nv+nc-itx(2)-[(nr)-itx(1)]
           ec(:,1)=(/1,-2,0,-2,2,0,-10,1,2,0,0,12/);    ! ferm-mode
           ec(:,2)=(/1,0,-2,-2,0,0,-8,1,0,0,2,-10/);   ! bio-mode 
           ec(:,3)=(/1,0,0,-6,-6,6,0,1,0,2,0,30/);    ! resp-mode
           !write(*,*) 'checks ec and bl'
           write(*,'(15f7.2)') matmul(stochx,ec(1:nr,1))
           write(*,'(15f7.2)') matmul(stochx,ec(1:nr,2))
           write(*,'(15f7.2)') matmul(stochx,ec(1:nr,3))
           !write(*,'(15f7.2)') matmul(stoch,ec(1:nr,4))
           bl=0.; bl(1,(/9,10/))=(/1,1/); ! phosphate cons
           bl(2,1:nv)=(/6,2,2,1,0,0,0,2,0,0/);  ! c conservation 
           bl(3,1:nv)=(/6,2,1,2,1,2,0,1,1,0/);  ! o conservation
           bl(4,1:nv)=(/12,3,2,0,2,0,1,2,2,0/);  ! o conservation
           bl(5,1:nv)=(/6,4,3,0,1,0,0,3,1,0/);  ! o conservation
           bl(:,nv+1:nv+nc)=bl(:,1:nc)
           ! check that conservation vectors spans left-null space of S
           write(*,'(15f7.2)') matmul(bl(1,1:nv+nc),stoch);
           write(*,'(15f7.2)') matmul(bl(2,1:nv+nc),stoch);
           write(*,'(15f7.2)') matmul(bl(3,1:nv+nc),stoch);
           write(*,'(15f7.2)') matmul(bl(4,1:nv+nc),stoch);
           write(*,'(15f7.2)') matmul(bl(5,1:nv+nc),stoch);
        endif

!!! THERMODYNAMIC PROPERTIES
        muy0(1:nc)=(/-392.,-238.,200.,-386.,-151.,16.,0./);
        mu0(1:nc)=muy0(1:nc);
        mu0(nc+1:nv)=(/-53.,-181.,0./) !formes fusionné -> 181 (30+151)  !8:ADP/9:ATP
        dg0(nc+1:nr)=matmul(mu0(1:nv),stochx(1:nv,nc+1:nr))
        !gly,   ferm,    resp,      biom   atp   gng 
        ! -176.  -4.00   -907.00    -57.00    -31.00     83.00 
        keq(nc+1:nr)=exp(-dg0(nc+1:nr)/2.5);
        !write(*,'(1A,18f11.2)') 'mu0=',mu0(1:nv)
        !write(*,'(1A,18f11.2)') 'dg0=',dg0(nc+1:nr)
        !write(*,'(1A,18e11.2)') 'keq=',keq(nc+1:nr) 
        thr=150.; 
        do i1=nc+1,nr !voir si j'enleve 
           if (keq(i1)<10**(-thr)) keq(i1)=10**(-thr)
           if (keq(i1)>10**(thr)) keq(i1)=10**thr
        enddo

!!! CHEMOSTATTING CONDTION
        kcat(1:nc)=(/1.,1.,1.,20.,20.,1.,20./); !j'ai changé kox=1
        !case of glucose nutrient source only
        yy=(/0.02,0.001,0.001,0.001,1.,1.,1./); 
        !case of acetate nutrient source
        muy=muy0+RT*log(yy); 

!!! SAMPLING PROCESS
        kcats(nc+1:nr)=(/1D6,1D0,1D0,1D6,1D0/); scts=-1D0; 
        do isamp=1,nsamp
1          xxs=0.01D0;
           call random_number(rrand);
           kcat(nc+1:nr)=kcats(nc+1:nr)*10**(0.2*(rrand(nc+1:nr)-0.5)); !optimization
           if (num==0) then
              ti=0D0; tf=1D10; xx=xxs; 
              call evolSEU(ti,tf,xx,idid,1d-10);
              mun(1:9)=mu0(1:9)+RT*log(xx); mun(10)=mu0(10)+RT*log(PT-xx(9));
              call sys_rand(xx,ff,0D0);  
              call syssol(xx); call sys_rand(xx,ff,0D0);

              if (sum(ff**2)>10**(-5.)) then
                 write(*,*) 'NO STEADY STATE!!!'
                 ff=0D0; xx=0D0; goto 1
              endif
              if (idid.ne.1) goto 1 !check
           endif
           call flux(xx,flu,dg,ent);
           if (isnan(ent(1))) goto 1 !check
           entex(1)=sum(ent);  entex(2)=sum(ent(1:nc)); entex(3)=sum(ent(nc+1:nr));
           h_ares=nb(2)*(-mun(9)+mun(10)+mun(5))/(mun(8)+3D0*mun(6)-3D0*mun(5)-2D0*mun(4))
           h_gly=nb(1)*(-mun(9)+mun(10)+mun(5))/(mun(1)-2D0*mun(8)-2D0*mun(4)) !glyco
           !write(10+i3,*) entex(1),flu,h_gly,h_ares,entex(1)/6D1
           sct=-flu(3);
           if (sct>scts) then
              scts=sct; 
              eprs(1)=entex(1); xxs=xx;flus=flu; kcats=kcat; etas=(h_ares+h_gly)/2.;
              eprs(2)=entex(3)
           endif
        enddo
        by=-flus(3)/flus(1)/3.;
        if (by<0) by=0
        write(11,*) nb(1:2),by,h_gly,h_ares,h_gly+h_ares,eprs
        write(*,*)  nb(1:2),by,h_gly,h_ares,h_gly+h_ares,eprs
     enddo
     write(11,*) ' '

  enddo

contains

  subroutine SVD(n1,n2,A,IT,nulldim) !travaille sur la stoch
    implicit none
    integer,intent(in)::n1,n2,it
    real(8),intent(in)::A(n1,n2)
    real(8)::S(n2)!,EV(nv,nr)
    integer,intent(out)::nulldim
    !integer:: LDA, LDU, LDVT
    integer,parameter:: LWMAX=1000
    integer::INFO,LWORK,i,i1
    real(8)::U(n1,n1),VT(n2,n2),WORK(LWMAX)
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

end program principal

!!!!!!!!!!!EXTERNAL SUBROUTINES FOR SEULEX !!!!!!!!!!!!!!!

subroutine solout (NRR,XOLD,X,Y,RC,LRC,IC,LIC,N,RPAR,IPAR,IRTRN)
  use glob
  implicit none
  integer :: N,NRR,LRC,LIC,IRTRN,ipar                               
  real(8) :: X,Y(N),RC(LRC),IC(LIC),rpar,XOLD
  real(4) :: t,val,diff
  real(8) :: F(N)

!!! stocke l'évolution temporelle dans un tableau
  if (NRR.LE.nplot) then
     tplot(NRR)=X;
     xplot(:,NRR)=Y(:);
  endif
  nmax=NRR
  if ((X>0).and.(NRR>3000)) then !! --- CRITERE D'AR1RET de evol_SEU---!!YO!!
     Y=-1; IRTRN=-1; return;
  endif
  return
end subroutine solout

!!$ ---------------------------------------
!!$   Routines nécessaires à la communication avec 
!!$   l'intégrateur numérique
!!$----------------------------------------

SUBROUTINE FCN(n,x,fvec,iflag) ! utilisé par seulex 
  use glob; implicit none;
  integer::iflag,n
  real(8)::t,x(n),fvec(n)
  call sys_rand(x,fvec,t) 
end SUBROUTINE FCN

SUBROUTINE FCN2(n,x,fvec,iflag) ! utilisé par seulex 
  use glob; implicit none;
  integer::iflag,n
  real(8)::t,x(n),fvec(n)
  call sys_rand(x,fvec,t) 
end SUBROUTINE FCN2

SUBROUTINE FCN_HW(N,t,xHW,F,rpar,ipar) ! ??
  use glob
  implicit none
  integer :: n,ipar
  real(8) :: t,xHW(n),f(n),rpar
  call sys_rand(XHW,f,t) !
end SUBROUTINE FCN_HW

subroutine jac_HW(N,X,Y,DFY,LDFY,RPAR,IPAR) ! utilisé par seulex
  use glob
  implicit none
  integer :: N,LDFY,LRC,LIC,IRTRN,ipar                                       
  real(8) :: X,Y(N),DFY(LDFY,N),rpar
  !call jacobien(Y,DFY,X)
  return
end subroutine jac_HW

subroutine mass
  return
end subroutine mass
