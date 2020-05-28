program project
integer i,j,pgopen
!integer env
character*32 str_lon1,str_lat1,str_lon2,str_lat2
integer lat_d,lon_d
real lat_m,lon_m,depth,d,p,di
real lon1,lat1,lon2,lat2,x1,y1,delta
real y(263584),x(263584),z(263584),cosr(263584),plotx(263584),distance(263584)
real x2(263584),y2(263584),delta2(263584),tw_x(1483),tw_y(1483)
real tx(3145728),ty(3145728),tz(3145728),tcosr(3145728),tplotx(3145728),tdistance(3145728)
real tx2(3145728),ty2(3145728),tdelta2(3145728)

istat=pgopen('profile.ps/vcps')
 if(istat<=0) stop 'ERR opening for PS file!'
 open(1,file='15.txt',status='old')
 open(2,file='output_test.txt',status='unknown')
 open(3,file='Taiwan.txt',status='old')
 Open(10,file='new.txt',status='old')

 10 read(1,'(18x,i2,f5.2,i3,f5.2,f6.2)',err=99)lat_d,lat_m,lon_d,lon_m,depth
	write(2,'(f9.5,1x,f9.5,1x,f6.2)')lat_d+lat_m/60,lon_d+lon_m/60,depth
 goto 10
 99 close(1)
 close(2)
 
 Print*,'start point:"lon1" "lat1",end point:"lon2" "lat2"'
 Read(*,*) lon1,lat1,lon2,lat2
 Print*,'set "the distance to the profile" "depth"'
 Read(*,*) d,di
 !Print*,'fixed the env,input "1", or not "0"'
 !Read(*,*) env

 call delaz(lat1,lon1,lat2,lon2,x1,y1,delta)
 
  write(str_lon1,'(f6.2)')lon1
  write(str_lat1,'(f5.2)')lat1
  write(str_lon2,'(f6.2)')lon2
  write(str_lat2,'(f5.2)')lat2
  
	do j=1,1483
	read(3,'(f7.3,1x,f6.3)') tw_x(j),tw_y(j)
	end do
	close(3)

 do i=1,263584
 open(2,file='output_test.txt',status='old')
 read(2,'(f9.5,1x,f9.5,1x,f6.2)')y(i),x(i),z(i)
 call delaz(lat1,lon1,y(i),x(i),x2(i),y2(i),delta2(i))
 cosr(i)=((x2(i)*x1)+(y2(i)*y1))/(delta2(i)*delta)
 plotx(i)=delta2(i)*cosr(i)
 distance(i)=sqrt((delta2(i)*delta2(i))-(plotx(i)*plotx(i)))
 end do
 
 !print*,plotx
 
	call pgsubp(1,3)
 
 	!TW
	call pgsci(1)
	call pgslw(1)
	call pgscf(3)
	call pgenv(119.,123.,21.,26.,1,0)
	call pglab('longitude','latitude','Taiwan')
	open(3,file='Taiwan.txt',status='old')

	do j=1,1483
	read(3,'(f7.3,1x,f6.3)')tw_x(j),tw_y(j)
	end do
	
	call pgsci(2)
	call pgtext(lon1,lat1-.2,str_lon1)
	call pgtext(lon1,lat1-.4,str_lat1)
	call pgtext(lon2,lat2-.2,str_lon2)
	call pgtext(lon2,lat2-.4,str_lat2)
	call pgsci(1)
	!call PGAXIS('N1',121.5,21.5,122.,21.5,0.0,50.0,50.0,1.0,1.0,1.0,0.0,0.0,0.0)


	call pgsci(1)
	call pgslw(2)
	call pgsls(1)
	call pgline(1483,tw_x,tw_y)

	call pgsci(3)
	call pgslw(3)
	call pgsls(3)
	call pgmove(lon1,lat1)
	call pgdraw(lon2,lat2)
	
	
	!topo
	 do i=1,3145728
 read(10,'(f9.5,1x,f8.5,1x,f10.4)')tx(i),ty(i),tz(i)
 call delaz(lat1,lon1,ty(i),tx(i),tx2(i),ty2(i),tdelta2(i))
 tcosr(i)=((tx2(i)*x1)+(ty2(i)*y1))/(tdelta2(i)*delta)
 tplotx(i)=tdelta2(i)*tcosr(i)
 tdistance(i)=sqrt((tdelta2(i)*tdelta2(i))-(tplotx(i)*tplotx(i)))
 end do
 
 call pgsci(1)
 call pgslw(1)
 call pgscf(3)
 call pgsch(1.0)
 call pgenv(0.,delta,minval(tz),maxval(tz),0,0)
 call pglab('distance(km)','hight(km)','topo profile')
 
 do i=1,3145728
 if(tdistance(i).le.1.)then
 call pgsch(1.0)
 call pgpt1(tplotx(i),tz(i),4)
 call pgsci(4)
 end if
 end do
	call pgsci(1)
	call pgslw(3)
	call pgsls(2)
	call pgmove(0.0,0.0)
	call pgdraw(delta,0.0)


 !profile
 call pgsci(1)
 call pgslw(1)
 call pgscf(3)
 call pgsch(1.0)
 call pgenv(0.,delta,di,0.,0,0)
 call pglab('distance(km)','depth(km)','profile')
 
 
 p=delta/20
 call pgtext(delta/2-3.5*p,di-10.,str_lon1)
 call pgtext(delta/2-1.5*p,di-10.,str_lat1)
 call pgtext(delta/2,di-10.,'to')
 call pgtext(delta/2+1.*p,di-10.,str_lon2)
 call pgtext(delta/2+3.*p,di-10.,str_lat2)
 
 
 do i=1,263584
 if(distance(i).le.d)then
 call pgpt1(plotx(i),z(i),4)
 call pgsci(2)
 end if
 end do
 
 call pgend
 

end program

subroutine delaz(elat,elon,slat,slon,dx,dy,delta)
      avlat=0.5*(elat+slat)
      a=1.840708+avlat*(.0015269+avlat*(-.00034+avlat*(1.02337e-6)))
      b=1.843404+avlat*(-6.93799e-5+avlat*(8.79993e-6+avlat*(-6.47527e-8)))
      dlat=slat-elat
      dlon=slon-elon
      dx=a*dlon*60.
      dy=b*dlat*60.
      delta=sqrt(dx*dx+dy*dy)
end	subroutine delaz