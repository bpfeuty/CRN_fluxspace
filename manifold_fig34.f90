! compute the boundary manifold of the thermodynamically-feasible flux space
! computational time is about 10seconds in CPU.
! n_c=3=n_bc+n_ec

program principal
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp),dimension(3)::j,j_s,yy,muy0,muy,k,kt,l
  real(dp)::thr,sig,sigmax,RT,emux,sig_s,sigi,sigi_s,sigi_ss
  integer::i1,i2,i3,ifig,ns

  open(13,file='manifold_fig3.dat'); !plot j3 as function of sig (column#3 against column#4 )
  open(15,file='test.dat'); !plot j3 as function of sig (column#3 against column#4 )
  open(14,file='manifold_fig4.dat'); !plot j1 as function of j2 (column#1 against column#2 )
  RT=2.5D0; ! gas constant times temperature
  ns=20000  ! grid parameter
  thr=2D-1
  do ifig=3,4,1
     k=1D0; l=1D0;
     if (ifig==3) then
        yy=(/1D0,1D0,0.5D0/);  muy0=RT*(/2D0,-2D0,3D0/); 
     endif
     if (ifig==4) then
        yy=(/1D0,1D0,1D0/); muy0=RT*(/1D0,1D0,-1D0/); 
     endif
     kt=k*exp(-muy0/RT); muy=muy0+RT*log(yy); !normalized kinetic constant
     
     !emux=sum(kt*exp(muy/RT))/sum(kt);
     emux=sum(k*yy)/sum(kt);         
     sigmax=sum(kt*muy*(exp(muy/RT)-emux));
     write(*,*) emux,sigmax
     !iterative methods
     do i2=1,ns
        sigi_ss=10.; sigi_s=10.;
        do i1=1,ns
           j(1)=4*real(i1-ns/2)/real(ns); !from -2 to 2
           j(2)=4.*real(i2-ns/2)/real(ns); !from -2 to 2
           j(3)=(-l(1)*j(1)-l(2)*j(2))/l(3)
           
           sigi_ss=sigi_s; sigi_s=sigi; sigi=0.;
           do i3=1,3 !internal entropy
              !sig_i=j_i*(RT*log(y_i-j_i/k_i)+muy0_i)
              sigi=sigi-j(i3)*(RT*log(yy(i3)-j(i3)/k(i3))+muy0(i3));
           enddo
           !write(*,*) i2,i1,sigi

           !! Eq.5 of the manuscript
           sig=sum(j*(muy0+RT*log(yy)))
           if ((i1>1.5).and.(i2>1.5)) then
              if ((abs(sigi_s)<thr).and.((abs(sigi_s)<abs(sigi)).and.(abs(sigi_s)<abs(sigi_ss)))) then
                 write(10+ifig,*) j_s,sig_s,sig_s/sigmax
                 if (ifig==3) write(15,*) sig_s/sigmax,j_s(3)
              endif
           endif
           j_s=j; sig_s=sig;
        enddo
     enddo
  enddo
  
end program principal
