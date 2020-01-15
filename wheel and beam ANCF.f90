!ANCF simply supported beam 2 dimension model
!Sophia university
!Michael Lim

!*********Calculation tool module************
module tool
implicit none
!dx（Extract X component vector）
real(8),parameter:: dx(1:2)=(/1.0d0,0.0d0/)
!dy（Extract y component vector）
real(8),parameter:: dy(1:2)=(/0.0d0,1.0d0/)
!χ（注：Invers）
real(8),parameter:: kai(2,2)=reshape((/0.0d0,1.0d0,-1.0d0,0.0d0/),(/2,2/))
!trans_χ（注：Invers）
real(8),parameter:: tkai(2,2)=reshape((/0.0d0,-1.0d0,1.0d0,0.0d0/),(/2,2/))
!Identity matrix
real(8),parameter:: ide(2,2)=reshape((/1.0d0,0.0d0,0.0d0,1.0d0/),(/2,2/))

contains
!matrix transformation counter clockwise
function cm(th) result(cm_out)
real(8),intent(in):: th
real(8) cm_out(2,2)
cm_out(1:2,1:2)=0.0d0
cm_out(1,1)= dcos(th)
cm_out(1,2)=-dsin(th)
cm_out(2,1)= dsin(th)
cm_out(2,2)= dcos(th)
return
end function cm

!matrix transformation clockwise（転置）
function tcm(th) result(tcm_out)
real(8),intent(in):: th
real(8) tcm_out(2,2)
tcm_out(1:2,1:2)=0.0d0
tcm_out(1,1)= dcos(th)
tcm_out(1,2)= dsin(th)
tcm_out(2,1)=-dsin(th)
tcm_out(2,2)= dcos(th)
return
end function tcm

!vector・matrix・vector func
function vmv(v1,m1,v2) result(v1m1v2)
real(8),intent(in):: v1(:),m1(:,:),v2(:)
real(8) v1m1(size(m1,2)),v1m1v2
if(size(v1) /= size(m1,1).or.size(m1,2) /= size(v2))then
  write(*,*) 'dimension error vmv'
  stop
endif
v1m1(1:size(m1,2))=0.0d0
v1m1v2=0.0d0
v1m1=matmul(v1,m1)
v1m1v2=dot_product(v1m1,v2)
return
end function vmv

!vector.matrix.matrix.vector func
function vmmv(v1,m1,m2,v2) result(v1m1m2v2)
real(8),intent(in):: v1(:),m1(:,:),m2(:,:),v2(:)
real(8) v1m1(size(m1,2)),v1m1m2(size(m2,2)),v1m1m2v2
if(size(v1) /= size(m1,1).or.size(m1,2) /= size(m2,1).or.size(m2,2) /= size(v2))then
  write(*,*) 'dimension error vmmv'
  stop
endif
v1m1(1:size(m1,2))=0.0d0
v1m1m2(1:size(m2,2))=0.0d0
v1m1m2v2=0.0d0
v1m1=matmul(v1,m1)
v1m1m2=matmul(v1m1,m2)
v1m1m2v2=dot_product(v1m1m2,v2)
return
end function vmmv

!vector.matrix.matrix.matrix.vector fun
function vmmmv(v1,m1,m2,m3,v2) result(v1m1m2m3v2)
real(8),intent(in):: v1(:),m1(:,:),m2(:,:),m3(:,:),v2(:)
real(8) v1m1(size(m1,2)),v1m1m2(size(m2,2)),v1m1m2m3(size(m3,2)),v1m1m2m3v2
if(size(v1) /= size(m1,1).or.size(m1,2) /= size(m2,1)&
  &.or.size(m2,2) /= size(m3,1).or.size(m3,2) /= size(v2))then
  write(*,*) 'dimension error vmmmv'
  stop
endif
v1m1(1:size(m1,2))=0.0d0
v1m1m2(1:size(m2,2))=0.0d0
v1m1m2m3(1:size(m3,2))=0.0d0
v1m1m2m3v2=0.0d0
v1m1=matmul(v1,m1)
v1m1m2=matmul(v1m1,m2)
v1m1m2m3=matmul(v1m1m2,m3)
v1m1m2m3v2=dot_product(v1m1m2m3,v2)
return
end function vmmmv

!vector.matrix.matrix
function vmm(v1,m1,m2) result(v1m1m2)
real(8),intent(in):: v1(:),m1(:,:),m2(:,:)
real(8) v1m1(size(m1,2)),v1m1m2(size(m2,2))
if(size(v1) /= size(m1,1).or.size(m1,2) /= size(m2,1))then
  write(*,*) 'dimension error vmm'
  stop
endif
v1m1(1:size(m1,2))=0.0d0
v1m1m2(1:size(m2,2))=0.0d0
v1m1=matmul(v1,m1)
v1m1m2=matmul(v1m1,m2)
return
end function vmm

!vector.matrix.matrix.matrix
function vmmm(v1,m1,m2,m3) result(v1m1m2m3)
real(8),intent(in):: v1(:),m1(:,:),m2(:,:),m3(:,:)
real(8) v1m1(size(m1,2)),v1m1m2(size(m2,2)),v1m1m2m3(size(m3,2))
if(size(v1) /= size(m1,1).or.size(m1,2) /= size(m2,1).or.size(m2,2) /= size(m3,1))then
  write(*,*) 'dimension error vmmm'
  stop
endif
v1m1(1:size(m1,2))=0.0d0
v1m1m2(1:size(m2,2))=0.0d0
v1m1m2m3(1:size(m3,2))=0.0d0
v1m1=matmul(v1,m1)
v1m1m2=matmul(v1m1,m2)
v1m1m2m3=matmul(v1m1m2,m3)
return
end function vmmm

!matrix.matrix.vetor
function mmv(m1,m2,v1) result(m1m2v1)
real(8),intent(in):: m1(:,:),m2(:,:),v1(:)
real(8) m1m2(size(m1,1),size(m2,2)),m1m2v1(size(m1,1))
if(size(m1,2) /= size(m2,1).or.size(m2,2) /= size(v1))then
  write(*,*) 'dimension error mmv'
  stop
endif
m1m2=matmul(m1,m2)
m1m2v1=matmul(m1m2,v1)
return
end function mmv

!corrdinate distance between point（v1,v2 inertia vector）
function dist(v1,v2) result(dist_out)
real(8),intent(in):: v1(2),v2(2)
real(8) dist_out
dist_out=dsqrt(dot_product(v1-v2,v1-v2))
end function 

!v direction unit vector 
function evec(v) result(evec_out)
real(8),intent(in):: v(2)
real(8) evec_out(2)
evec_out(1:2)=0.0d0
evec_out=v(1:2)/dsqrt(dot_product(v,v))
end function
				
!１dimension tensor multiplication
function tens(v1,v2) result(tens_out)
real(8),intent(in):: v1(2),v2(2)
real(8) tens_out(2,2)
tens_out(1:2,1:2)=0.0d0
tens_out(1,1:2)=v1(1)*v2(1:2)
tens_out(2,1:2)=v1(2)*v2(1:2)
end function

!i()) x coordinate contact point
function nr1(i) result(nr1_out)
integer,intent(in):: i
integer nr1_out
nr1_out=4*i-3
end function nr1

!(i) y coordinate contact point
function nr2(i) result(nr2_out)
integer,intent(in):: i
integer nr2_out
nr2_out=4*i-2
end function nr2

!(i)dr/dx(1) coordinate contact point
function nr3(i) result(nr3_out)
integer,intent(in):: i
integer nr3_out
nr3_out=4*i-1
end function nr3

!i/dx(2) coordinate contact point
function nr4(i) result(nr4_out)
integer,intent(in):: i
integer nr4_out
nr4_out=4*i
end function nr4


end module tool

!****************Parameter setting module*********************
module wheel_data
implicit none
real(8),parameter:: dt=1.0d-6    !1.0d-6 !2.5d-7        !-----Time Interval[s]
real(8),parameter:: ga=9.80665d0                        !-----Gravity acceleration
real(8),parameter:: pi=3.1415926535897932d0             !-----Pi
real(8),parameter:: eps1=1.0d-6 !1.0d-7(best) !5.0d-8   !-----Tangent constraint error parameter
real(8),parameter:: eps2=1.0d-6 !1.0d-7(best) !5.0d-8   !-----Position constrant error parameter
real(8),parameter:: epsI=5.0d-11                        !-----Initial value condition
real(8),parameter:: v0=4.0d0                            !-----initial speed[km/h]
!wheel parameter
integer,parameter:: nwh=1						        !-----amount of Wheels
real(8),parameter:: ma=1.3461d0!1.3461d0 !5.45d0        !-----wheel mass
real(8),parameter:: r0=0.043d0   !0.43d0	            !-----wheel diameter
real(8),parameter:: ja=ma*(r0**2)/2.0d0      !46.37d0   !-----Wheel moment of inertia




!wheel and rail parameter(rigidity,spring,elasticity)
real(8),parameter:: rora=7850.0d0    !2700.0d0          !-----Rail density
real(8),parameter:: era=2.07d11           !2.07d11      !-----Rail young modulus
real(8),parameter:: pora=0.3d0      !0.35d0             !-----Rail poasson ratio



real(8),parameter:: kc=2.27875d10 !2.6023d10  !2.2895d10  !2.3233d10  !-----contact rigidity
real(8),parameter:: cc=0!1.08d4   !2guzai*rute.mkc'       !-----contact damper guzai 0.344


!軌道パラメータ
real(8),parameter:: ks=0.9d7  !0.9d7  !0.6d7  !1.0d7          !-----Sleeper spring constant
real(8),parameter:: cs=1.0d2!1.0d2                            !-----Sleeper damping constant
real(8),parameter:: kb=1.3d7  !1.3d7  !2.6d7                  !---Foundation constant
real(8),parameter:: cb=1.0d3!1.0d3                            !---Foundation spring　　　　　
real(8),parameter:: ms=0.51d0                                  !---Sleeper mass

real(8),parameter:: pitch=0.08d0             !0.581d0          !-----Sleeper interval

integer,parameter:: nsl=4 !10 !8  !2                   !-----Amount of sleeper

real(8),parameter:: hra=17.0d-3 !0.2312d0               !-----Rail height
real(8),parameter:: wra=6.5d-3  !3.0d-2				    !-----Rail width
real(8),parameter:: ira=(wra*hra**3)/12.0d0	            !-----Rail moment of inertia
real(8),parameter:: ara=hra*wra				            !-----Rail cross section
integer,parameter:: ne=2 	  !7                        !-----Rail element（inside sleeper interval）
real(8),parameter:: dl=pitch/real(ne)                   !-----1 element length
integer,parameter:: nn=ne*(nsl-1)+1		                !-----Number of rail nodes
integer,parameter:: nrc=(nn-1)/2+1		                !-----Rail center nodes
integer,parameter:: nra=4*nn	                        !-----Rail matrix size
integer,parameter:: rm= ne*4                            !------ sleeper location on the rail nodes

integer,parameter:: jq=2*nsl
integer,parameter:: tm=nra+jq

integer,parameter:: nq=nra+jq+3*nwh  	                !-----Generalized coordinates
integer,parameter:: ns=2					            !-----Not Generalized coordinates(wheel angle/ril position)
integer,parameter:: nz=12*nn-10                         !-----mass matrix (For sparse matrix)
integer,parameter:: n=nq		    		            !-----Matrix size

integer,parameter:: sl_sta=2                             !wheel starting point（～th sleeper）
integer,parameter:: e_sta2=ne*(sl_sta-1)+1              !wheel starting points nodes
real(8),parameter:: e_staF=pitch*(sl_sta-1)             !wheel starting points position

integer,parameter:: sl_end=3                            !wheel end point（～th sleeper）
integer,parameter:: e_end=ne*(sl_end-1)+1               !wheel end point
real(8),parameter:: e_endF=pitch*(sl_end-1)             !wheel end point position


!Lap parameter
integer,parameter:: start=1                          !start at 1
integer,parameter:: round=1      !2000               !Total lap
integer,parameter:: interval=1                       !output frequency

end module wheel_data

!**************Main program********************
program main
    use wheel_data
    use tool
    implicit none
    character filename*128
    integer i,j,k,cou,el_cou,z,ie1,st,mode,lap,shape

    integer t1,t2,t_rate2,t_max2,diff2 !計算時間測定
    real*8 t,q(2*n),qdd(n),el_q(n),el_s(4),sW1,sR1
    real*8 mass(tm,tm),stiff(nra,nra),Kt(nra,nra)!,damp(nra,nra)
    real*8 const(2),lra,rr1(2),rw1(2),tr1(2),nr_1(2)
    real*8 ef1(2),dn1,jt1(2),nef1,nra2
   !
    real*8 ee1(8),d_x

    integer irow(nz),jcol(nz)
    real*8 as(nz)

    call system_clock(t1)

!Rail Length
lra=dl*dble(nn-1)

!Initial condition initialization
q(1:2*n)=0.0d0
el_q(1:n)=0.0d0
el_s(1:4)=0.0d0

!Rail contact point initial condition(Displacement)

do i=1,nn
    q(nr1(i))=dble(i-1)*dl
	q(nr2(i))=0.0d0
	q(nr3(i))=1.0d0
	q(nr4(i))=0.0d0
enddo

!Wheel initial condition(Position)
    q(nra+jq+1)=e_staF    
    q(nra+jq+2)=r0 !r0      
    q(nra+jq+3)=0.0d0



!Not generalized coordinate
    sW1=0.0d0             !wheel angle
	sR1=0.0d0             !wheel length
    
    
!Rail matrix・Rigid matrix calculation
 
    call mass_matrix(mass,irow,jcol,as,0)
    call stiff_matrix(q,stiff,Kt,0)

    
    !Numerical calculatin
      cou=100
      d_x=dl/10.0d0
      el_cou=100  !100
      t=0.0d0
      ie1=e_sta2           !----Initial element beam
      
      st=1
    
      do !(start at static state)
        if(el_cou.eq.100)then   
 !	    
          el_cou=0
        endif
          el_cou=el_cou+1 
       
 
        call search_cp(q,sW1,sR1,ie1,start)
                 
        call runge(t,q,qdd,ie1,st,sW1,sR1,rr1,rw1,tr1,ef1,const,start,nef1,jt1,dn1) !Start not lap!
 
           t=t+dt
        if(t>0.05d0)exit !0.3  
      enddo

     !Save initial condition
       el_q(1:n)=q(1:n)
       el_s(1)=sW1
       el_s(2)=sR1

 
      write(*,*) 'Balance end'
     
     
     !make it dynamic
       st=0
     !Time initialization
       t=0.0d0
     !shape initialization
       shape=0
 
     !Speed initialization
       q(n+1:2*n)=0.0d0

 !Intitial speed
       q(n+nra+jq+1)=v0/3.6d0  !*(q(nr3(e_sta))*dcos(atan(dsy))-q(nr4(e_sta))*dsin(atan(dsy)))
       q(n+nra+jq+2)=0.0d0     !v0/3.6d0*(q(nr4(e_sta))*dsin(atan(dsy))+q(nr3(e_sta))*dsin(atan(dsy)))
       q(n+nra+jq+3)=-q(n+nra+jq+1)/r0

       !********Output ************

       write(*,*) 'round and interval',round,interval

       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
       do lap=start,round

   
        !****Built Dat file for every lap****
         if(lap==1.OR.mod(lap,interval)==0)then    
          
         
          write (filename, '( i4.4,"Wheel Displacement.dat")') lap
          open (50, file=filename, status='replace')
          
          write (filename, '( i4.4,"Wheel Velocity.dat")') lap
          open (60, file=filename, status='replace')
          
    
          write (filename, '( i4.4,"component.dat")') lap
          open (90, file=filename, status='replace')
          
          write(filename, '( i4.4,"Rail displacement at center.dat")') lap
          open (100, file=filename, status='replace')
     
       
      
         endif    
        
        if(lap.ne.1)then
          if(mod(lap,interval)==0)then 
        
   
          endif
         endif
       !***********************************
     
       if(lap.ne.1)then  !
        
         !initial condition
         q(1:2*n)=0.0d0
     
         !Rail contact point initial condition(Displacement)
     
        do i=1,nn
         q(nr1(i))=dble(i-1)*dl
         q(nr2(i))=0.0d0
         q(nr3(i))=1.0d0
         q(nr4(i))=0.0d0
        enddo
     
       !Whell initial condition(Position)
         q(nra+jq+1)=e_staF    
         q(nra+jq+2)=r0 !r0      
         q(nra+jq+3)=0.0d0
         

  !Not generalized coordinate initial condition
    sW1=0.0d0             !wheel angle
   
    sR1=0  !wheel length
    
    !Intial condition(Need more research)
    q(nra+jq+2)=q(nra+jq+2)



     
    !****Create initial status output file****	 
  	
   do !(Static state)
       
	   if(lap==1.OR.mod(lap,interval)==0)then

		 if(el_cou.eq.100)then

           el_cou=0
		 		 
		 endif
	       el_cou=el_cou+1 
	   
	   else
	   endif
	 
		 call search_cp(q,sW1,sR1,ie1,lap)
   		
         call runge(t,q,qdd,ie1,st,sW1,sR1,rr1,rw1,tr1,ef1,const,lap,nef1,jt1,dn1)

		  t=t+dt
	   if(t>0.05d0)exit
  enddo

  
  !Dynamic
  st=0
    
  !Time initialization
      t=0.0d0
  !shape initialization
    shape=0

  !Speed initialization
   q(n+1:2*n)=0.0d0

    q(n+nra+jq+1)=v0/3.6d0 
    q(n+nra+jq+2)=0.0d0    
    q(n+nra+jq+3)=-q(n+nra+jq+1)/r0
 
    endif

!---Output file
    do 
   
        call search_cp(q,sW1,sR1,ie1,lap)
        
        call  runge(t,q,qdd,ie1,st,sW1,sR1,rr1,rw1,tr1,ef1,const,lap,nef1,jt1,dn1)


    if(lap==1.OR.mod(lap,interval)==0)then

	  if(cou.eq.100)then

		!write(10,"(100e30.15)") t,q(nr1(1)),q(nr1(nn)),q(nr1(nrc))
		!write(20,"(100e30.15)") t,q(nr2(1)),q(nr2(nn)),q(nr2(nrc)),-qdd(nr2(nrc))       
		!write(30,"(100e30.15)") t,q(nr3(1)),q(nr3(nrc)),q(nr3(nn))
		!write(40,"(100e30.15)") t,q(nr4(1)),q(nr4(nrc)),q(nr4(nn))

		write(50,"(100e30.15)") t,q(nra+jq+1:nra+jq+2)
		write(60,"(100e30.15)") t,q(n+nra+jq+1:n+nra+jq+2)

		write(90,"(100e30.15)") t,ef1,jt1,nef1,dn1,rr1,rw1,rw1-rr1!,rr1(1),rr1(2),
        write(100,"(100e30.15)") t,q(nr2(4))


		cou=0
	  endif
	  cou=cou+1
	endif
	
	t=t+dt
	
	if(rr1(1)>=q(nr1(e_end)))then
	
	exit 
    endif
	
  enddo
    write(*,*) 'This model ran',lap,'laps'

enddo


close(50)
close(60)

close(90)


write(*,*) 'finish'
stop
end


!========Below this part  are subroutine ==========

!#######################################
! Contact point search
!#######################################
subroutine search_cp(q,sW1,sR1,ie1,lap)
  use wheel_data
  use tool
  implicit none
  integer,intent(in)::lap
  integer,intent(inout)::ie1
  real*8,intent(in)::q(2*n)
  real*8,intent(inout)::sW1,sR1

  integer i
  real*8 roa1(2),ee1(8),dee1(8),zeta1,xs1,sf1(2,8),dsf1(2,8),d2sf1(2,8)
  real*8 twlsW1(2),trlsR1(2)
  real*8 sWst1,sRst1
  real*8 rsin,rcos  
  real*8 ds1(2)
  real*8 rr1(2),tr1(2),nr_1(2)                       !Rail contact point　position,tangent,normal
  real*8 rw1(2),tw1(2),nw1(2)                        !Wheel contact point　position,tangent,normal
  real*8 E1(2)                                       !E(1) contact condition
  real*8 rwr1(2)                                     !Rail and wheel (rw-rr)
  real*8 Es1(2,2),Esinv1(2,2)
  
 
  !Wheel center of gravity displacement
    roa1(1:2)=q(nra+jq+1:nra+jq+2)

 
  !wheel contact element node coordinates
    ee1(1:8)=q((ie1-1)*4+1:(ie1-1)*4+8) 
    dee1(1:8)=q((ie1-1)*4+1+n:(ie1-1)*4+8+n)
 

   !sR1 for switching element
    if(sR1>dl)then
   ie1=ie1+1
   sR1=sR1-dl
    ee1(1:8)=q((ie1-1)*4+1:(ie1-1)*4+8) 
    dee1(1:8)=q((ie1-1)*4+1+n:(ie1-1)*4+8+n)
    else
     if(sR1<0.0d0)then
     ie1=ie1-1
     sR1=sR1+dl
     ee1(1:8)=q((ie1-1)*4+1:(ie1-1)*4+8) 
     dee1(1:8)=q((ie1-1)*4+1+n:(ie1-1)*4+8+n)
     endif
    endif
 
   
 !-----Contact point search---------------
   !Save initial state
   sWst1=sW1

   sRst1=sR1

 
  do i=1,1000
     zeta1=sR1/dl

 
    if(zeta1 < 0.0d0.or.zeta1 > 1.0d0)exit
    

 
   !ANCF Shape matrix
   sf1(1:2,1:8)=0.0d0
   sf1(1,1)=1.0d0-3.0d0*(zeta1**2)+2.0d0*(zeta1**3)
   sf1(1,3)=dl*(zeta1-2.0d0*(zeta1**2)+(zeta1**3))
   sf1(1,5)=3.0d0*(zeta1**2)-2.0d0*(zeta1**3)
   sf1(1,7)=dl*(-(zeta1**2)+(zeta1**3))
   sf1(2,2)=sf1(1,1)
   sf1(2,4)=sf1(1,3)
   sf1(2,6)=sf1(1,5)
   sf1(2,8)=sf1(1,7)
 
   
 
   !dsf/ds
   dsf1(1:2,1:8)=0.0d0
   dsf1(1,1)=(-6.0d0*zeta1+6.0d0*(zeta1**2))/dl
   dsf1(1,3)=1.0d0-4.0d0*zeta1+3.0d0*(zeta1**2)
   dsf1(1,5)=(6.0d0*zeta1-6.0d0*(zeta1**2))/dl
   dsf1(1,7)=-2.0d0*zeta1+3.0d0*(zeta1**2)
   dsf1(2,2)=(-6.0d0*zeta1+6.0d0*(zeta1**2))/dl
   dsf1(2,4)=1.0d0-4.0d0*zeta1+3.0d0*(zeta1**2)
   dsf1(2,6)=(6.0d0*zeta1-6.0d0*(zeta1**2))/dl
   dsf1(2,8)=-2.0d0*zeta1+3.0d0*(zeta1**2)
 
   
 
   !d2sf/ds2
   d2sf1(1:2,1:8)=0.0d0
   d2sf1(1,1)=(-6.0d0+12.0d0*zeta1)/(dl**2)
   d2sf1(1,3)=(-4.0d0+6.0d0*zeta1)/dl
   d2sf1(1,5)=(6.0d0-12.0d0*zeta1)/(dl**2)
   d2sf1(1,7)=(-2.0d0+6.0d0*zeta1)/dl
   d2sf1(2,2)=(-6.0d0+12.0d0*zeta1)/(dl**2)
   d2sf1(2,4)=(-4.0d0+6.0d0*zeta1)/dl
   d2sf1(2,6)=(6.0d0-12.0d0*zeta1)/(dl**2)									
   d2sf1(2,8)=(-2.0d0+6.0d0*zeta1)/dl
 

   
 
!===wheel contact point position===!
   rw1(1)=roa1(1)-r0*dsin(sW1)  !x position vector
   rw1(2)=roa1(2)-r0*dcos(sW1)  !y position vector
   tw1(1)=-r0*dcos(sW1)         !Tangent x component vector
   tw1(2)=r0*dsin(sW1)          !Tangent y component vector
   nw1(1)=-r0*dsin(sW1)         !Normal x component vector
   nw1(2)=-r0*cos(sW1)          !Normal y component vector

 

!===Rail contact point position===!
   rr1=matmul(sf1,ee1)
   tr1=matmul(dsf1,ee1)    ! rail tangent vector      
   nr_1=matmul(kai,tr1)    ! rail normal vector              



!===Vector difference between rail and wheel===!
   rwr1=rw1-rr1



!===Partial derivative of tw with respect to sw===!
   twlsW1(1)=r0*dsin(sW1)
   twlsW1(2)=r0*dcos(sW1)


!===Partial derivative of tr with respect to sr===!
   trlsR1=matmul(d2sf1,ee1)



!=====  wheel contact conditionE(s) ======!
   E1(1)=dot_product(tw1,nr_1)
   E1(2)=dot_product(tr1,rwr1)
  


 !==== determine the error ====!
   if (abs(E1(1))<=eps1.and.abs(E1(2))<=eps2)then  !eps1 and eps2 are the allowed error value, 
   exit 

   endif


 !=====  Contact condition (partial derivative with respec ti sr/sw) Es(2,2) =====!
   Es1(1,1)=dot_product(twlsW1,nr_1)
   Es1(1,2)=dot_product(tr1,tw1)
   Es1(2,1)=dot_product(tw1,matmul(kai,trlsR1))
   Es1(2,2)=dot_product(trlsR1,rwr1)-dot_product(tr1,tr1)

 call minver(Es1,2)             !----Inverse matrix

 ds1(1:2)=matmul(Es1,E1)  

   sW1=sW1-ds1(1)
   sR1=sR1-ds1(2)
   
   !wheel switchig element
   if(sR1>dl)then
	ie1=ie1+1
	sR1=sR1-dl
    ee1(1:8)=q((ie1-1)*4+1:(ie1-1)*4+8)
    dee1(1:8)=q((ie1-1)*4+1+n:(ie1-1)*4+8+n)
   else
    if(sR1<0.0d0)then
	ie1=ie1-1
    sR1=sR1+dl
    ee1(1:8)=q((ie1-1)*4+1:(ie1-1)*4+8)
    dee1(1:8)=q((ie1-1)*4+1+n:(ie1-1)*4+8+n)
    else
	endif
  endif


 enddo

!  Determine the error!!
 if (abs(E1(1))>eps1.or.abs(E1(2))>eps2)then
  sW1=sWst1
  sR1=sRst1
  write(*,*) 'tansaku_error1',E1(1),E1(2),lap
  !stop
 endif



 return
 end


 !************************************************************
 !************************************************************

!***********Runge kutta 4th order**************
subroutine runge(t,q,qdd,ie1,st,sW1,sR1,rr1,rw1,tr1,ef1,const,lap,nef1,jt1,dn1)

use wheel_data
implicit none
integer,intent(in)::ie1,st,lap
real*8,intent(in):: sW1,sR1
real*8,intent(inout):: q(2*n)
real*8 rr1(2),tr1(2),const(4),slip1,pmax1,rw1(2)
real*8 ef1(2),nef1,jt1(2),dn1

real*8,intent(out):: qdd(n)
real*8,dimension(2*n):: b
real*8,dimension(2*n,4):: k
real*8 t

call fafv(q,qdd,t,ie1,st,sW1,sR1,rr1,rw1,tr1,ef1,const,lap,nef1,jt1,dn1)         

k(1:n,1)=q(n+1:2*n)
k(n+1:2*n,1)=qdd(1:n)

b(1:2*n)=q(1:2*n)+0.5d0*dt*k(1:2*n,1)

call fafv(b,qdd,t,ie1,st,sW1,sR1,rr1,rw1,tr1,ef1,const,lap,nef1,jt1,dn1)

k(1:n,2)=b(n+1:2*n)
k(n+1:2*n,2)=qdd(1:n)

b(1:2*n)=q(1:2*n)+0.5d0*dt*k(1:2*n,2)

call fafv(b,qdd,t,ie1,st,sW1,sR1,rr1,rw1,tr1,ef1,const,lap,nef1,jt1,dn1)

k(1:n,3)=b(n+1:2*n)
k(n+1:2*n,3)=qdd(1:n)

b(1:2*n)=q(1:2*n)+dt*k(1:2*n,3)

call fafv(b,qdd,t,ie1,st,sW1,sR1,rr1,rw1,tr1,ef1,const,lap,nef1,jt1,dn1)    

k(1:n,4)=b(n+1:2*n)
k(n+1:2*n,4)=qdd(1:n)

q(1:2*n)=q(1:2*n)+dt*((k(1:2*n,1)+2.0d0*(k(1:2*n,2)+k(1:2*n,3))+k(1:2*n,4))/6.0d0)
qdd(1:n)=(k(n+1:2*n,1)+2.0d0*(k(n+1:2*n,2)+k(n+1:2*n,3))+k(n+1:2*n,4))/6.0d0 
return
end

!*************MAss matrix with constraint condition***************
subroutine mass_matrix(mass,irow,jcol,as,mode)
use wheel_data
implicit none
integer i,j,k,mode,d
real*8,save:: mass0(tm,tm)      
real*8,save:: mass_v(tm,tm)     
real*8,intent(out):: mass(tm,tm)
real*8 me(8,8)
real*8 a,b,c
integer,save:: irow0(nz),jcol0(nz)
real*8,save:: as0(nz)
integer,intent(out):: irow(nz),jcol(nz)
real*8,intent(out):: as(nz)
real*8,save:: massinv(tm,tm)


select case (mode)
case (0)

a=ara*rora*dl
b=a*dl
c=b*dl

me(1:8,1:8)=0.0d0
me(1,1)= a*0.371428571428571d0
me(1,3)= b*0.052380952380952d0
me(1,5)= a*0.128571428571428d0
me(1,7)=-b*0.030952380952380d0
me(2,2)= me(1,1)
me(2,4)= me(1,3)
me(2,6)= me(1,5)
me(2,8)= me(1,7)

me(3,3)= c*0.009523809523809d0
me(3,5)= -me(1,7)
me(3,7)=-c*0.007142857142857d0
me(4,4)= me(3,3)
me(4,6)= me(3,5)
me(4,8)= me(3,7)

me(5,5)= me(1,1)
me(5,7)=-me(1,3)
me(6,6)= me(1,1)
me(6,8)= me(5,7)

me(7,7)= me(3,3)
me(8,8)= me(7,7)

do i=3,8
	do j=1,i-2
		me(i,j)=me(j,i)
	enddo
enddo


mass0(1:tm,1:tm)=0.0d0
do k=1,nn-1
	do i=1,8
		do j=1,8
			mass0((k-1)*4+i,(k-1)*4+j)=mass0((k-1)*4+i,(k-1)*4+j)+me(i,j)
		enddo
	enddo
enddo


do k=1,nsl
  !constraint condition for simply supported beam
        mass0((2*k-1)+nra,((k-1)*rm+1))= 1d0   !kousoku x axis   row
        mass0(2*k+nra,((k-1)*rm+2))= 1d0       !kousoku y axis   row
        mass0((k-1)*rm+1,(2*k-1)+nra)=1d0     !kousoku x axis    column
        mass0((k-1)*rm+2,2*k+nra)=1d0          !kousoku y axis   column
 
enddo



!Contraint condition
mass_v=mass0

mass_v(1,1:nra)=0.0d0
mass_v(1:nra,1)=0.0d0
mass_v(1,1)=1.0d0
mass_v(nra-3,1:nra)=0.0d0
mass_v(1:nra,nra-3)=0.0d0
mass_v(nra-3,nra-3)=1.0d0


!Sparse mode
call sparse_mode(mass_v,nra,irow0,jcol0,as0) 

call minver(mass_v,tm)             !----Inverse matrix

case (1)
 mass=mass_v
 irow=irow0
 jcol=jcol0
 as=as0

case (2)
 mass=mass_v
 
case (3) !mass call!!!!
 mass=mass0

end select

return
end

!! ============simple stiffness matrix of ANC ================================
subroutine stiff_matrix(q,stiff,Kt,mode)
use wheel_data
implicit none
integer i,j,k,mode
!***** mode0 variable definition *****!
real*8,save:: Kle(8,8),Kte(8,8),Kt0(nra,nra),Ktb(15,nra) !15=上バンド幅7+下バンド幅7+1
real*8 a,b
!***** mode1 variable definition*****!
real*8,intent(in):: q(2*n)
real*8,intent(inout):: stiff(nra,nra)
real*8 ipsiron,ee(8)
!***** mode2 variable definitin *****!
real*8,intent(out):: Kt(nra,nra)

select case (mode)

!******************** mode(0) ********************!
case (0)   !element matrix definition

Kle(1:8,1:8)=0.0d0
Kte(1:8,1:8)=0.0d0

a=era*ara/dl
b=era*ira/(dl**3)

Kle(1,1)=1.0d0
Kle(1,5)=-1.0d0
Kle(2,2)=1.0d0
Kle(2,6)=-1.0d0
Kle(5,1)=-1.0d0
Kle(5,5)=1.0d0
Kle(6,2)=-1.0d0
Kle(6,6)=1.0d0
Kle=Kle*a

Kte(1,1:8)=(/12.0d0 , 0.0d0 ,6.0d0*dl, 0.0d0, -12.0d0, 0.0d0, 6.0d0*dl, 0.0d0 /)
Kte(2,1:8)=(/0.0d0, 12.0d0, 0.0d0, 6.0d0*dl, 0.0d0, -12.0d0, 0.0d0, 6.0d0*dl /)
Kte(3,1:8)=(/6.0d0*dl, 0.0d0, 4.0d0*dl**2, 0.0d0, -6.0d0*dl, 0.0d0, 2.0d0*dl**2, 0.0d0/)
Kte(4,1:8)=(/0.0d0, 6.0d0*dl, 0.0d0, 4.0d0*dl**2, 0.0d0, -6.0d0*dl, 0.0d0, 2.0d0*dl**2/)
Kte(5,1:8)=(/-12.0d0, 0.0d0, -6.0d0*dl, 0.0d0, 12.0d0, 0.0d0, -6.0d0*dl, 0.0d0 /)
Kte(6,1:8)=(/0.0d0, -12.0d0, 0.0d0, -6.0d0*dl, 0.0d0, 12.0d0, 0.0d0, -6.0d0*dl /)
Kte(7,1:8)=(/6.0d0*dl, 0.0d0, 2.0d0*dl**2, 0.0d0, -6.0d0*dl, 0.0d0, 4*dl**2, 0.0d0/)
Kte(8,1:8)=(/0.0d0, 6.0d0*dl, 0.0d0, 2.0d0*dl**2, 0.0d0, -6.0d0*dl, 0.0d0, 4*dl**2/)
Kte=Kte*b

Kt0(1:nra,1:nra)=0.0d0

do k=1,nn-1
	do i=1,8
		do j=1,8
			Kt0((k-1)*4+i,(k-1)*4+j)=Kt0((k-1)*4+i,(k-1)*4+j)+Kte(i,j)
		enddo
	enddo
enddo


!******************** end mode(0) ********************!

!******************** mode(1) ********************!
case(1)   !rigidity matrix definition

Stiff(1:nra,1:nra)=0.0d0

do i=1,nn-1
 ee(1:8)=q(4*(i-1)+1:4*(i-1)+8)
 ipsiron=1.0d0-dl/sqrt((ee(5)-ee(1))**2+(ee(6)-ee(2))**2)

 !ipsiron=sqrt( (ee(5)-ee(1))**2 + (ee(6)-ee(2))**2 )/dl - 1.0d0 
 Stiff(4*(i-1)+1:4*(i-1)+8,4*(i-1)+1:4*(i-1)+8)=Stiff(4*(i-1)+1:4*(i-1)+8,4*(i-1)+1:4*(i-1)+8)+ipsiron*Kle(1:8,1:8) 
enddo

 Stiff(1:nra,1:nra)=Stiff(1:nra,1:nra)+Kt0(1:nra,1:nra)

!******************** end mode(1) ********************!

!******************** mode(2) ********************!
case(2)   !call elastic rigidity matrix

Kt(1:nra,1:nra)=Kt0(1:nra,1:nra)

!******************** end mode(2) ********************!

end select

return
end

!***********************************************************

!*************(Equation of motion)**************
subroutine fafv(q,qdd,t,ie1,st,sW1,sR1,rr1,rw1,tr1,ef1,const,lap,nef1,jt1,dn1)
                 
use wheel_data
use tool

implicit none
integer i,j
integer,parameter:: ITWKSP=0
integer,intent(in)::ie1,st,lap
real*8 t

real*8,intent(in):: q(2*n),sW1,sR1,rr1(2),tr1(2)
real*8,intent(out):: qdd(n),const(2)
real*8 mass(tm,tm),stiff(nra,nra),Kt(nra,nra)
real*8 c1,c2,qk(nra),qc(nra)
real*8 exf(nq),nef1,jt1,dn1(2),rw1(2)

real*8 zeta1,xs
real*8 ef1(2),cpf1(2),Ref1(8)

real*8 mw1(3,3)!,mw2(3,3)
integer irow(nz),jcol(nz),IPARAM(6)
real*8 as(nz),RPARAM(5)


 call  contact_creep_force(q,sW1,sR1,ie1,st,ef1,rr1,rw1,Ref1,tr1,lap,nef1,jt1,dn1)

!external force
exf(1:nq)=0.0d0


!---***Rail external force definition***---
  !Rail gravity
c1=0.5d0*(-rora*ara*dl*ga)
c2=dl*0.125d0*(-rora*ara*dl*ga)

do i=1,nn-1
   exf(4*i-2)=exf(4*i-2)+c1
 exf(4*i)  =exf(4*i)  +c2
 exf(4*i+2)=exf(4*i+2)+c1
 exf(4*i+4)=exf(4*i+4)-c2
enddo

!Rail elasticity
qk(1:nra)=0.0d0
qc(1:nra)=0.0d0

call stiff_matrix(q,stiff,Kt,1)
qk=qk+matmul(stiff,q(1:nra))         !Elasticity


exf(1:nra)=exf(1:nra)-qk(1:nra) !-qc(1:nra)

!normal force and tangent force from wheel
exf(nr1(ie1):nr4(ie1+1))=exf(nr1(ie1):nr4(ie1+1))+Ref1    !normal force




!External force from sleeper
!do i=1,nsl 
 !exf(nr2(ne*(i-1)+1))=exf(nr2(ne*(i-1)+1))&
  ! &-ks*(q(nr2(ne*(i-1)+1))-q(nra+jq+3*nwh+i))&
 !&-cs*(q(n+nr2(ne*(i-1)+1))-q(n+nra+jq+3*nwh+i))!-----sleeper rigidity and damping，
!enddo

!Finite element boundary condition
exf(nr1(1))=0.0d0                      
exf(nr1(nn))=0.0d0  

!---***End rail external force definition***---!



!---***Wheel mass matrix definition***---!
 mw1(1:3,1:3)=0.0d0
 
 mw1(1,1)=ma
 mw1(2,2)=ma
 mw1(3,3)=ja

!***Wheel mass inverse matrix definition***!
call minver(mw1,3)             !----inverse matrix

!---***Wheel external force Definition***---!


!***wheel external force**
  exf(nra+jq+1:nra+jq+2)=exf(nra+jq+1:nra+jq+2)+ef1!+cpf1!+qwf1 !normal + (tangent) no tangent in this case
  exf(nra+jq+2)=exf(nra+jq+2)-ma*ga                     !gravity of the wheel
  exf(nra+jq+3)= 0                    !assuming no creep force

!***end wheel external force definition***!


!--*******second derivative calculation*******--
  qdd(1:nq)=0.0d0
 
  
 !***Rail contacts point***!
  call mass_matrix(mass,irow,jcol,as,2)  !1:sparse 2:inverse

  qdd(1:tm)=matmul(mass,exf(1:tm))
     
  const(1)=q(nr1(1))
  const(2)=q(nr1(nn))
 
 !***Wheel***!
  qdd(nra+jq+1:nra+jq+3)=matmul(mw1,exf(nra+jq+1:nra+jq+3))  !wheel q


 !***sleeper***!
  !do i=1,nsl
    !qdd(nra+3*nwh+i)=(ks*(q(nr2(ne*(i-1)+1))-q(nra+3*nwh+i))&
                   ! &+cs*(q(n+nr2(ne*(i-1)+1))-q(n+nra+3*nwh+i))&
                    !&-kb*q(nra+3*nwh+i)-cb*q(n+nra+3*nwh+i))/ms-ga
 !  enddo  

return
end

!*************Creep force***************
subroutine contact_creep_force(q,sW1,sR1,ie1,st,ef1,rr1,rw1,Ref1,tr1,lap,nef1,jt1,dn1)
  
use wheel_data
use tool
implicit none
integer i,j,k
integer,intent(in):: ie1,st,lap
real*8,intent(in)::q(2*n),sW1,sR1

real*8,intent(out)::ef1(2),rr1(2),tr1(2),Ref1(8)
real*8 roa1(2),tha1,voa1(2),oma1,ee1(8),dee1(8)
real*8 zeta1,xs1
real*8 sf1(2,8),dsf1(2,8)

real*8 nef1
real*8 nr_1(2)                    
real*8 rw1(2),tw1(2),nw1(2)            !Contact point on wheel Position, tangent, normal
real*8 it1(2),jt1(2),alph1(2,2),ltrl1,djt1(2)

real*8 vw1(2),vr1(2)
real*8 dn1,ddn1
integer,save::flag



!wheel center of gravity displacement / speed
roa1(1:2)=q(nra+jq+1:nra+jq+2)
tha1=q(nra+jq+3)
voa1(1:2)=q(n+nra+jq+1:n+nra+jq+2)
oma1=q(n+nra+jq+3)



!Contact elemnt node coordinates
ee1(1:8)=q((ie1-1)*4+1:(ie1-1)*4+8)
dee1(1:8)=q((ie1-1)*4+1+n:(ie1-1)*4+8+n)



!ζcalculation
zeta1=sR1/dl
!zeta2=sR2/dl

!ANC matrix shape
sf1(1:2,1:8)=0.0d0
sf1(1,1)=1.0d0-3.0d0*(zeta1**2)+2.0d0*(zeta1**3)
sf1(1,3)=dl*(zeta1-2.0d0*(zeta1**2)+(zeta1**3))
sf1(1,5)=3.0d0*(zeta1**2)-2.0d0*(zeta1**3)
sf1(1,7)=dl*(-(zeta1**2)+(zeta1**3))
sf1(2,2)=sf1(1,1)
sf1(2,4)=sf1(1,3)
sf1(2,6)=sf1(1,5)
sf1(2,8)=sf1(1,7)



!dsf/ds
dsf1(1:2,1:8)=0.0d0
dsf1(1,1)=(-6.0d0*zeta1+6.0d0*(zeta1**2))/dl
dsf1(1,3)=1.0d0-4.0d0*zeta1+3.0d0*(zeta1**2)
dsf1(1,5)=(6.0d0*zeta1-6.0d0*(zeta1**2))/dl
dsf1(1,7)=-2.0d0*zeta1+3.0d0*(zeta1**2)
dsf1(2,2)=(-6.0d0*zeta1+6.0d0*(zeta1**2))/dl
dsf1(2,4)=1.0d0-4.0d0*zeta1+3.0d0*(zeta1**2)
dsf1(2,6)=(6.0d0*zeta1-6.0d0*(zeta1**2))/dl
dsf1(2,8)=-2.0d0*zeta1+3.0d0*(zeta1**2)





!===wheel contact point position===!
rw1(1)=roa1(1)-r0*dsin(sW1)
rw1(2)=roa1(2)-r0*dcos(sW1)
tw1(1)=-r0*dcos(sW1)
tw1(2)=r0*dsin(sW1)
nw1(1)=-r0*dsin(sW1)
nw1(2)=-r0*cos(sW1)



!===wheel rail contact point position===!
rr1=matmul(sf1,ee1)
tr1=matmul(dsf1,ee1)             
nr_1=matmul(kai,tr1)   

ltrl1=dsqrt(tr1(1)**2+tr1(2)**2)             



!===wheel contact point speed===!
!vw(1)=voa(1);vw(2)=voa(2)
vw1(1)=voa1(1)+r0*oma1*dcos(sW1)
vw1(2)=voa1(2)-r0*oma1*dsin(sW1)     



!===wheel rail contact point speed===!
vr1=matmul(sf1,dee1)



!前輪接線・法線方向単位ベクトル
it1=tr1/ltrl1
jt1=matmul(kai,it1)



!前輪法線方向単位ベクトルの時間微分(固定点)
alph1=(ide-tens(tr1,tr1)/dot_product(tr1,tr1))/ltrl1 
djt1=matmul(kai,mmv(alph1,dsf1,dee1))



!===Elastic force calculation===!
!dn(Contact amount)
dn1=dot_product(jt1,rw1-rr1)


!ddn(derivative of dn1)
ddn1=dot_product(jt1,vw1-vr1)+dot_product(djt1,rw1-rr1)


!wheel normal contact force
if(dn1<0.0d0)then
nef1=-kc*dn1*dsqrt(dabs(dn1))!-cc*ddn1*dabs(dn1)
else
nef1=0.0d0
endif



!Contact force in absolute coordinate system
ef1=jt1*nef1


!for Rail!
Ref1=matmul(transpose(sf1),-ef1)    !Normal force






return
end

!*****Matrix inverse subroutine*****
subroutine minver(a,n)
  implicit real*8 (a-z)
  integer i,iw,j,k,l,lr,m,err,n
  integer work(n)
  real*8 a(n,n)
     l=n
     m=n
     eps=1.0d-12  
       
       if (m.lt.2 .or. m.gt.n .or. eps.le.0.) then
         err=999
       else
         err=0
         det=1.0d0
         do 10 i=1,m
           work(i)=i
  10     continue
         do 70 k=1,m
           wmax=0.0d0
           do 20 i=k,m
             w=dabs(a(i,k))
             if (w.gt.wmax) then
               wmax=w
               lr=i
             endif
  20       continue
           pivot=a(lr,k)
           api=dabs(pivot)
           if (api.le.eps) then
             err=1
             return
           endif
           det=det*pivot
           if (lr.ne.k) then
             det=-det
             iw=work(k)
             work(k)=work(lr)
             work(lr)=iw
             do 30 j=1,m
               w=a(k,j)
               a(k,j)=a(lr,j)
               a(lr,j)=w
  30         continue
           endif
           do 40 i=1,m
             a(k,i)=a(k,i)/pivot
  40       continue
           do 60 i=1,m
             if (i.ne.k) then
         w=a(i,k)
               if (w.ne.0.) then
                 do 50 j=1,m
                   if (j.ne.k) a(i,j)=a(i,j)-w*a(k,j)
  50             continue
                 a(i,k)=-w/pivot
               endif
             endif
  60       continue
           a(k,k)=1.0d0/pivot
  70     continue
         do 100 i=1,m
  80       k=work(i)
           if (k.ne.i) then
             iw=work(k)
             work(k)=work(i)
             work(i)=iw
             do 90 j=1,m
               w=a(j,i)
               a(j,i)=a(j,k)
               a(j,k)=w
  90         continue
             goto 80
           endif
  100    continue
       endif
       return
       end
 
 
 
 !don't TOuch!!!!
 !=========IMSL用 Sparse mode generation module========!
 !
 !INPUT a:matrix
 !      n:matrix size
 !OUTPUT  
 !      irow:非ゼロ要素　row
 !      jcol:非ゼロ要素　column
 !        as:非ゼロ要素リスト
 !        nz:非ゼロ要素数
 
 subroutine sparse_mode(a,n,irow,jcol,as)
! Do not Touch
 use wheel_data,only:nz
 implicit none
 integer i,j,k,n
 real*8 a(n,n)
 integer irow(nz),jcol(nz)
 real*8 as(nz)
 
 !**** initial value ****
  irow(1:nz)=0
  jcol(1:nz)=0
    as(1:nz)=0.0d0
 !****************
 
 !irow,jcol,as produce
 
 do i=1,n
   do j=1,i
      if(a(i,j)/=0.0d0)then
         do k=1,nz
         if(as(k)==0.0d0)then
             irow(k)=i
         jcol(k)=j
         as(k)=a(i,j)
             exit
         endif
       enddo
      endif 
   enddo
 enddo
 
 return
 end