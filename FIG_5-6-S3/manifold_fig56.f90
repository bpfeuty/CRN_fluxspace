! Compute the boundary manifold of the thermodynamically-feasible flux space of the coarse grained model
! gfortran manifold_fig56.f90; (no library need to be linked)
! computational time is about 100 seconds on a "Intel Core i7-8665U"
! "feasible_fig5.dat" is a large file (72GO).
program principal
  real(16),parameter::RT=2.5D0
  integer,parameter::nc=7,n3=200,n2=40,n1=1100,n6=60 !grid parameter

  real(16),dimension(nc)::l1,l2,l3,l4 ! conservation vectors
  real(16),dimension(nc)::j,yy,mu0y,k

  real(16),dimension(3,n1,7)::jmaxglu
  real(16),dimension(3,n2,7)::jmaxace
  real(16),dimension(3,n6,7)::jmaxox

  real(16)::sigi,sav,tmp,thr,epr,zz1,zz2
  integer::i1,i2,i3,i6,i0,ii

  open(50,file='boundary_fig5.dat'); !compute the boundary manifold 
  open(51,file='feasible_fig5.dat'); !compute the feasible subspace
  open(52,file='jaexglu_fig5.dat'); ! compute the boundary projection: Jaex=f(Jgex)
  open(53,file='jbexglu_fig5.dat'); ! compute the boundary projection: Jbex=f(Jgex)
  open(54,file='jcexglu_fig5.dat'); ! compute the boundary projection: Jcex=f(Jgex)
  open(55,file='jaexace_fig5.dat'); ! compute the boundary projection: Jaex=f(Jaex)
  open(56,file='jbexace_fig5.dat'); ! compute the boundary projection: Jbex=f(Jaex)
  open(57,file='jcexace_fig5.dat'); ! compute the boundary projection: Jcex=f(Jaex)
  open(58,file='jaexox_fig5.dat'); ! compute the boundary projection:  Jaex=f(Joex)
  open(59,file='jbexox_fig5.dat'); ! compute the boundary projection:  Jbex=f(Joex)
  open(60,file='jcexox_fig5.dat'); ! compute the boundary projection:  Jcex=f(Joex)

  !chemostatting property 
  yy=(/2D-2,1D-3,1D-3,1D-3,1D0,1D0,1D0/);
  mu0y=(/-392.,-238.,200.,-386.,-150.,16.,0./);
  k=(/1D0,1D0,1D0,2D1,2D1,1D0,2D1/);

  !conservation vectors l^y
  l1=(/6,2,2,1,0,0,0/);  !carbon conservation
  l2=(/6,2,1,2,1,2,0/);  !oxygen conservation
  l3=(/12,3,2,0,2,0,1/); !hydrogen conservation
  l4=(/6,4,3,0,1,0,0/);  !topology-based conservation

  jmaxglu=0D0; jmaxox=0D0; jmaxace=0D0;  
  do i3=1,n3
     j(3)=k(3)*yy(3)-0.061D0*real(i3-1)/(real(n3-1))
     do i2=1,n2
        j(2)=k(2)*yy(2)-0.061D0*real(i2-1)/(real(n2-1)) 
        do i1=1,n1
           zz1=(real(i1)-1002)/2D0; zz2=10D0**zz1;          
           if (i1<998.5) then
              j(1)=k(1)*yy(1)*(1.-zz2);  !yy(1)-zz2;
           else
              j(1)=k(1)*yy(1)*(1.-real(i1-998)/5D1) !yy(1)-real(i1-999)/2D2
           endif
           j(5)=-(l4(1)*j(1)+l4(3)*j(3)+l4(2)*j(2))/l4(5); ! additional broken law
           j(4)=-(l1(1)*j(1)+l1(3)*j(3)+l1(2)*j(2))/l1(4); ! Carbon conservation law
           j(6)=-(l2(1)*j(1)+l2(3)*j(3)+l2(4)*j(4)+l2(5)*j(5)+l2(2)*j(2))/l2(6); !Ox conservation law
           j(7)=-(l3(1)*j(1)+l3(2)*j(2)+l3(3)*j(3)+l3(5)*j(5))/l3(7);

           i6=int(j(6)/0.16D0*5D1)+10
           sigi=0.; !important
           do i0=1,6
              if (i0==1) then
                 if (i1<998.5) then
                    sigi=sigi+j(i0)*(mu0y(i0)+RT*(log(yy(i0))+zz1*log(1D1)));
                 else
                    sigi=sigi+j(i0)*(mu0y(i0)+RT*log(yy(i0)-j(i0)/k(i0)));
                 endif
              else
                 sigi=sigi+j(i0)*(mu0y(i0)+RT*log(yy(i0)-j(i0)/k(i0)));
              endif
           enddo
           epr=sum(j*(mu0y+RT*log(yy)));
           thr=0.3;
           if (sigi>0) then
              write(51,'(10f11.4)') j,epr
              do ii=1,3
                 if (j(ii+1)<jmaxglu(ii,i1,ii+1)) then
                    jmaxglu(ii,i1,1:6)=j(1:6);
                    jmaxglu(ii,i1,7)=epr;
                 endif

                 if (j(ii+1)<jmaxace(ii,i2,ii+1)) then
                    jmaxace(ii,i2,1:6)=j(1:6);
                    jmaxace(ii,i2,7)=epr;
                 endif
                 if ((i6<n6+0.5).and.(i6>0.5)) then
                    if (j(ii+1)<jmaxox(ii,i6,ii+1)) then
                       jmaxox(ii,i6,1:6)=j(1:6);
                       jmaxox(ii,i6,7)=epr;
                    endif
                 endif
              enddo
           endif
           if (abs(sigi)<thr) then
              write(50,*) j(1:6),epr !glucimp,aceimp
           endif
        enddo
     enddo
  enddo

  do i1=1,n1
     do ii=1,3
        write(51+ii,'(17f13.5)') jmaxglu(ii,i1,1:7)
     enddo
  enddo
  do i1=1,n2
     do ii=1,3
        write(54+ii,'(17f13.5)') jmaxace(ii,i1,1:7)
     enddo
  enddo
  do i1=1,n6
     do ii=1,3
        write(57+ii,'(17f13.5)') jmaxox(ii,i1,1:7)
     enddo
  enddo
end program principal
