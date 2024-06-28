! check kernel and co-kernel for different stoichiometric matrixes 

program principal !!REF
  implicit none;
  integer,parameter::nv=10,nr=12,nc=7; !ih=feedback
  real(kind=kind(0.d0)),dimension(nv,nr)::stoch,mp6x
  real(kind=kind(0.d0)),dimension(nr,nv)::stocht,mp6tx
  real(kind=kind(0.d0)),dimension(nv+nc,nr)::stoch2,mp6
  real(kind=kind(0.d0)),dimension(nr,nv+nc)::stoch2t,mp6t
  integer::itx(4),i1,it,nb(4)
  integer,dimension(5,nv+nc)::bl
  integer,dimension(nr,4)::ec
 
  nb=(/10,10,10,1/); !glyc,tca,biom,gng   
  !algebraic analysis of the stoichiometric matrix 
  if (2>1) then !with atp-adp
     !species: 1gluc,2ace,3biom,4co2,5h2o,6o2,7acoa,8adp,9atp
     !reactions: 1gluc,2ace,3biom,4co2,5h2o,6o2,7glyc,8ferm,9resp,10biom,11atpc
     mp6=0.; 
     do i1=1,7
        mp6(i1,i1)=1.; mp6(nv+i1,i1)=-1.;
     enddo
     mp6((/1,4,7,8/),8)=(/-1.,2.,8.,2./); mp6((/5,9,10/),8)=nb(1)*(/1.,-1.,1./);
     mp6((/2,7,8,9,10/),9)=(/1.,1.,-1.,-1.,1./);
     mp6((/4,6,7,8/),10)=(/2.,-3.,-4.,-1./); mp6((/9,10/),10)=nb(2)*(/-1.,1./); mp6(5,10)=nb(2)+3; 
     mp6((/3,8/),11)=(/1.,-1./);  mp6((/5,9,10/),11)=nb(3)*(/-1.,1.,-1./);
     mp6((/5,9,10/),12)=(/-1.,1.,-1./); 
     mp6t=transpose(mp6); mp6x(1:nv,:)=mp6(1:nv,:); mp6tx(:,1:nv)=mp6t(:,1:nv); stoch2=mp6; stoch=mp6x;

     call SVD(nv+nc,nr,mp6,0,it); itx(1)=it; 
     call SVD(nr,nv+nc,mp6t,0,it); itx(2)=it; 
     call SVD(nv,nr,mp6x,0,it); itx(3)=it; 
     call SVD(nr,nv,mp6tx,0,it); itx(4)=it;
     write(*,*) 'cc=',itx(1)
     write(*,*) 'ec=',itx(3)-itx(1) 
     write(*,*) 'bl=',itx(2)-itx(4)
     write(*,*) 'ul=',itx(4)
     write(*,*) 'checks(=0)',itx(3)-itx(1)+itx(2)-itx(4)-nc,nv+nc-itx(2)-[(nr)-itx(1)]

     ec(:,1)=(/1,-2,0,-2,2,0,-10,1,2,0,0,12/);    ! ferm-mode 12
     ec(:,2)=(/1,0,-2,-2,0,0,-8,1,0,0,2,-10/);   ! bio-mode OK
     ec(:,3)=(/1,0,0,-6,-6,6,0,1,0,2,0,30/);    ! resp-mode
     !ec(:,4)=(/0,0,0,0,12,-6,1,0,-2,0,-11,-1/); ! gluconeo-mode no closed
     write(*,*) 'checks ec and bl'
     write(*,'(15f7.2)') matmul(stoch,ec(1:nr,1))
     write(*,'(15f7.2)') matmul(stoch,ec(1:nr,2))
     write(*,'(15f7.2)') matmul(stoch,ec(1:nr,3))
     bl=0.; bl(1,(/9,10/))=(/1.,1./); ! P conservation
     bl(2,1:nv)=(/6.,2.,2.,1.,0.,0.,0.,2.,0.,0./);  ! C conservation 
     bl(3,1:nv)=(/6.,2.,1.,2.,1.,2.,0.,1.,1.,0./);  ! O conservation
     bl(4,1:nv)=(/12.,3.,2.,0.,2.,0.,1.,2.,2.,0./); ! H conservation 
     bl(5,1:nv)=(/6.,4.,3.,0.,1.,0.,0.,3.,1.,0./); ! H conservation 
     
     bl(:,nv+1:nv+nc)=bl(:,1:nc)
     write(*,'(15f7.2)') matmul(bl(1,1:nv+nc),stoch2);
     write(*,'(15f7.2)') matmul(bl(2,1:nv+nc),stoch2);
     write(*,'(15f7.2)') matmul(bl(3,1:nv+nc),stoch2);
     write(*,'(15f7.2)') matmul(bl(4,1:nv+nc),stoch2);
     write(*,'(15f7.2)') matmul(bl(5,1:nv+nc),stoch2);
  endif
  ! s-l=r-c : n_cc-n_l = n_r-n_s  & nc=nec+nbl => n_cc=n_r-n_s+n_l
  ! 6=3+3 ; nl=4?  & n_cc=12-9+4 (unbroken law versus closed cycle)

contains
  
  subroutine SVD(n1,n2,A,IT,nulldim) 
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

subroutine mass
  return
end subroutine mass
