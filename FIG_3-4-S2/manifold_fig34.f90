! Compute the boundary manifold of the thermodynamically-feasible flux space of the generic CRN model
! gfortran manifold_fig34.f90 (no library is linked)

program principal
  implicit none
  real(kind=kind(0.d0)),parameter::RT=2.5d0
  real(kind=kind(0.d0)),dimension(3)::j,j_s,yy,muy0,muy,k,kt,l
  real(kind=kind(0.d0))::sig,sigmax,emux,sig_s,sigi,sigi_s,sigi_ss
  integer::i1,i2,i3,ifig,ns

  open(13,file='manifold_fig3.dat'); !plot j3 as function of sig (column#3 against column#5 )
  open(14,file='manifold_fig4.dat'); !plot j1 as function of j2 (column#1 against column#2 )

  ns=20000  ! grid parameter
  do ifig=3,4,1
     k=1D0;
     l=1D0; !broken mass conservation law
     if (ifig==3) then
        yy=(/1D0,1D0,0.5D0/);  muy0=RT*(/2D0,-2D0,3D0/); 
     endif
     if (ifig==4) then
        yy=(/1D0,1D0,1D0/); muy0=RT*(/1D0,1D0,-1D0/); 
     endif
     kt=k*exp(-muy0/RT); muy=muy0+RT*log(yy); !normalized kinetic constant
     emux=sum(k*yy)/sum(kt); sigmax=sum(kt*muy*(exp(muy/RT)-emux));
     write(*,*) "Figure",ifig
     write(*,*) "Tsig^max=",sigmax
     
     !iterative method
     do i2=1,ns
        sigi_ss=10.; sigi_s=10.;
        do i1=1,ns
           j(1)=4*real(i1-ns/2)/real(ns); !from -2 to 2
           j(2)=4.*real(i2-ns/2)/real(ns); !from -2 to 2
           j(3)=(-l(1)*j(1)-l(2)*j(2))/l(3);
           
           sigi_ss=sigi_s; sigi_s=sigi; sigi=0.;
           do i3=1,3 !internal entropy
              !sig_i=j_i*(RT*log(y_i-j_i/k_i)+muy0_i)
              sigi=sigi-j(i3)*(RT*log(yy(i3)-j(i3)/k(i3))+muy0(i3));
           enddo

           !! Eq.5 of the manuscript
           sig=sum(j*(muy0+RT*log(yy)))
           
           !detect intersection sigs=0
           if ((sigi_ss*sigi<0.).and.((abs(sigi_s)<abs(sigi)).and.(abs(sigi_s)<abs(sigi_ss)))) then
              write(10+ifig,*) j_s,sig_s,sig_s/sigmax
           endif
           j_s=j; sig_s=sig;
        enddo
     enddo
  enddo
  
end program principal
