

!==========================================================================
      subroutine getcon(ifile,lenfile)

      use msflib
      character*(*) ifile
      integer*2    n,n1

      narg = nargs()

      if (narg.lt.2) then
        write(*,'(//,a,/)') ' Usage: exefile inputfile'
        stop
      end if
      n = 1
      call getarg(n,ifile,n1)
      lenfile = n1
      return
      end



!==========================================================================
      subroutine cmpdate(iyr,imo,idy, ju,jyr,jmo,jdy)

      if(imo.eq.jmo .and. idy.eq.jdy .and. iyr.eq.jyr) return

      write(*,'(a)') ' Date problem in input data file'
      write(*,'(i7,2i2,a   )') iyr,imo,idy,'  Simulation Date'
      write(*,'(i7,2i2,a,i3)') jyr,jmo,jdy,'  Date in file unit',ju
      stop
      end



!==========================================================================
      function     mndays(iyr,imn)

      ! returns number of days for month (iyr,imn)

      dimension    ndays(12,0:1)
      data ndays/31,28,31,30,31,30,31,31,30,31,30,31,
     &           31,29,31,30,31,30,31,31,30,31,30,31/

      leapyr = 0
      if(mod(iyr,4).eq.0) leapyr=1
      if(mod(iyr,100).eq.0) leapyr=0    ! Excel has error for 1900
      if(mod(iyr,400).eq.0) leapyr=1
      mndays = ndays(imn,leapyr)
      return
      end


!==========================================================================
      function     mndayex(iyr,imn)

      ! this subroutine is the version for excel date

      dimension    ndays(12,0:1)
      data ndays/31,28,31,30,31,30,31,31,30,31,30,31,
     &           31,29,31,30,31,30,31,31,30,31,30,31/

      leapyr = 0
      if(mod(iyr,4).eq.0) leapyr=1
c      if(mod(iyr,100).eq.0) leapyr=0
c      if(mod(iyr,400).eq.0) leapyr=1
      mndayex = ndays(imn,leapyr)
      return
      end


!==========================================================================
      subroutine adddays(iyr,imo,idy, nad, jyr,jmo,jdy)

c     this subroutine advances nad days from (iyr,imo,idy) and
c     returns (jyr,jmo,jdy)

      jyr = iyr
      jmo = imo
      jdy = idy + nad

 10   mday = mndays(jyr,jmo)
      if(jdy.le.mday) goto 20

      jmo = jmo+1
      jdy = jdy-mday
      if(jmo.gt.12) then
        jyr = jyr + 1
        jmo = 1
      end if
      goto 10

 20   write(*,'(1x,3i4,i7,2x,3i4)') iyr,imo,idy,nad,jyr,jmo,jdy
      return
      end



!==========================================================================
      subroutine skipline(iu,ch1)

      character*1 ch1,c1

      do while (.true.)                             ! skip comment lines
        read(iu,'(a1)',end=20) c1
        if(ch1.ne.c1) goto 10
      end do
 10   backspace(iu)
      return
 20   write(*,'(a,i5)') ' End of file in unit = ',iu
      stop
      end



!==========================================================================
c     This subroutine removes any character with,
c     ASCII code less than or equal to 32 (blank is ACHAR(32)), and
c     ASCII code greater or equal to 127, and
c     converts lower case to upper case,
c     and insert trailing blanks.
c     to 8 character variable
c     Jeffrey H. Hwang, Ph.D.    5/2/01

      subroutine ctrim8r(cb)

      character  ca*8, cb*8, ch1*1

      j = 0
      ca = cb
      do i=1,8
        read(ca(i:i),'(a1)') ch1
        k = iachar(ch1)
        if (k.gt.32  .and. k.lt.127) then
          if(k.ge.97 .and. k.le.122) k = k - 32   ! convert to upper case
          j=j+1
          write(cb(j:j),'(a1)') achar(k)
        end if
      end do
      if (j.lt.8) then
        j=j+1
        do i=j,8
          write(cb(i:i),'(a1)') achar(32)
        end do
      end if
      return
      end



!==========================================================================
      subroutine openfile(iunit,filename,nskip)
      character filename*50, ch1*1

      write(*,'(i4,1x,a)') iunit,filename
      open(iunit,file=filename,status='OLD')
      if(nskip.gt.0) then
        do i=1,nskip
          read(iunit,'(a1)') ch1
        end do
      end if
      return
      end



!==========================================================================
! open input data file
      subroutine openifile(iu1,iu2,ifile)
      character ifile*50, bl50*50, ch1*1

      bl50 = '                                                  '

      read(iu1,'(a50,i2)') ifile,nskip
      if(ifile.eq.bl50) return

      write(*,'(/,a,a,i3)') ' Opening ',ifile,iu2
      open(iu2,file=ifile,status='OLD')
      if(nskip.gt.0) then
        do i=1,nskip
          read(iu2,'(a1)') ch1
        end do
      end if
      return
      end

!=======================================================================

      subroutine openhdat(iu1,iu2,nsy,ismonth,ifile,fmt)

      character*50 ifile
      character*40 fmt
      character*1 ch1
      
      ch1 = '*'
      call OpnIF2(iu1,iu2,ifile,ch1)
      read(iu2,'(a40)') fmt

      idy = 1
      call LocTS1(iu2,nsy,ismonth,idy)
      return
      end

      
!==========================================================================
      subroutine OpnIF1(iu,ifile,cskip)
      character ifile*50, cskip*1

      write(*,'(i4,1x,a)') iu,ifile
      open(iu,file=ifile,status='OLD')
      if(cskip.eq.' ') return
      call skipline(iu,cskip)
      return
      end
      

!==========================================================================
      subroutine OpnIF2(iu1,iu2,ifile,cskip)
      character ifile*50, bl50*50, cskip*1

      bl50 = '                                                  '

      read(iu1,'(a50)') ifile
      if(ifile.eq.bl50) return
      
      open(iu2,file=ifile,status='OLD')
      write(*,'(i4,1x,a,a)') iu2,ifile,' is open'
      if(cskip.eq.' ') return
      call skipline(iu2,cskip)
      return
      end
      

      
!==========================================================================
      subroutine OpnOF2(iu1,iu2,ofile,cha,chb)
      character*(*)  cha, chb      
      character ofile*50, bl50*50

      bl50 = '                                                  '

      read(iu1,'(a50)') ofile
      if(ofile.eq.bl50) return
      
      open(iu2,file=ofile,status='UNKNOWN')
      write(iu2,'(a,a)') '* ',cha
      write(iu2,'(a,a)') '* Response file = ',chb
      write(*,'(i4,1x,a,a)') iu2,Ofile,' is open'
      return
      end
c      if(cha.ne.bl50) write(iu2,'(a)') cha
c      if(chb.ne.bl50) write(iu2,'(a,a)') '* Response file = ',chb


      
!==========================================================================
! locate time series data file
      subroutine LocTS1(iu1,iyr,imn,idy)
c     character rfmt*40

      ! assume first three positions of rfmt is for year,month,day
 10   read(iu1,'(i4,2i2)',err=10,end=20) jyr,jmn,jdy
      if(jyr.ne.iyr .or. jmn.ne.imn .or. jdy.ne.idy) goto 10
      backspace(iu1)
      write(*,'(a,i3,a)') ' File with iunit=',iu1,
     &                    ' is located at start of simulation.'
      return
        
 20   write(*,'(a,i3,a)') ' File with iunit=',iu1,
     &                    ' has problem in data.  Stop execution.'
      stop
      end
      

!==========================================================================
! locate monthly time series data file
      subroutine LocTS2(iu1,rfmt,iyr,imn)
      character rfmt*40

      ! assume first three positions of rfmt is for year,month,day
 10   read(iu1,rfmt,err=20,end=20) jyr,jmn
      if(jyr.ne.iyr .or. jmn.ne.imn) goto 10
      backspace(iu1)
      write(*,'(a,i3,a)') ' File with iunit=',iu1,
     &                    ' is located at start of simulation.'
      return
        
 20   write(*,'(a,i3,a)') ' File with iunit=',iu1,
     &                    ' has problem in data.  Stop execution.'
      stop
      end
      
!==========================================================================
! locate monthly time series data file
      subroutine LocTS3(iu1,iyr,imn)

      ! each line starts with yyyymmdd
 10   read(iu1,'(i4,2i2)',err=20,end=20) jyr,jmn
      if(jyr.ne.iyr .or. jmn.ne.imn) goto 10
      backspace(iu1)
      write(*,'(a,i3,a)') ' File with iunit=',iu1,
     &                    ' is located at start of simulation.'
      return
        
 20   write(*,'(a,i3,a)') ' File with iunit=',iu1,
     &                    ' has problem in data.  Stop execution.'
      stop
      end

!==========================================================================

      subroutine stagec(slop,xn,b,z,dpb,dpz,q,d,uperc)
c
c     q=k*s**.5
c
c     solve (q/s**.5) - K = 0.  with Newton's method
c     F =  (q * xn) / (1.486 * sqrt(slop)) - xk = 0.
c
c     k=1.486 * area * (hydualic radius) **.6667
c
c    initial  guess at D
c

      if(b.le.0.) b = 1.  !!!!!!!!!!!!!!!!!!!!!!!
      v=2.
      area=q/v
      d = (sqrt(b*b + 4.*z*area) - b)/(2.*z)
      dsave = d
 1    rm=d*(b+d*z)
      rn=b+2*d*sqrt(z*z+1)
      xk1 = (b+2*d) * (rm/rn)**.667
      xk2 = .667 * rn**.667 * rm**(-0.333) * (b+2.*d*z)
      xk3 = .667 * rn**(-0.333) * rm**.667 * (2.*sqrt(z*z+1))

      F1 = -(xk1 + rm * (xk2 - xk3) / rn**1.333)
      F = (q*xn)/(1.486*sqrt(slop)) - rm * (rm/rn)**.667
      delta = - F/F1
      if(delta.gt.0.) delta=delta/2.
      d = d + delta
c      write(*,'(a,3f10.5)') ' StageC',q,d,delta !!!!!!!!!!!delete
c      if (abs(delta).le.0.00001) go to 10
      if (abs(delta).le.0.001) go to 10
      if(d.gt.0.) goto 1

      write(*,'(//,a)') ' Problem in StageC'                           ! 5/4/05
      write(*,'(f10.0,7f10.4)') q,dsave,d,slop,xn,b,z

      d = dsave/2.
      dsave = d
      goto 1
      write(*,'(a)') ' Increase bottom width and/or reduce side slope' ! 5/4/05
      uperc = -99999.                                                  ! 5/4/05
      return                                                           ! 5/4/05

 10   continue
c
c     compute hydraulic elements
c
      area=d*(b+d*z)
      wp=b+2*d*sqrt(z*z+1)
      sp=wp-b
      r=area/wp
      uperc=dpb*b+dpz*sp
      return
      end


!==========================================================================

      subroutine stagep()
c     subroutine stagep(slop,xn,diam,q,d)
      return
      end



cccccccc!==========================================================================
ccccccccc this function interpolates from multiple tables
cccccccc
cccccccc      function rintp(it,jt,xv)
cccccccc      ! xv    - x value to interpolate
cccccccc      ! it    - type of tables, 0 for flow, 1 for TDS, 2 for TIN
cccccccc      ! jt    - table number
cccccccc
cccccccc
cccccccc      include 'ro_dim_3.max'
cccccccc
cccccccc      include 'tbprpt.var'
cccccccc
cccccccc
cccccccc      do k=2,ntab(it,jt)
cccccccc        if(xv.le.table(it,jt,k,1)) go to 10
cccccccc      end do
cccccccc      k = ntab(it,jt)
cccccccc 10   continue
cccccccc
cccccccc      yv = table(it,jt,k-1,2) + (table(it,jt,k,2)-table(it,jt,k-1,2))
cccccccc     &                        / (table(it,jt,k,1)-table(it,jt,k-1,1))
cccccccc     &                        * (xv              -table(it,jt,k-1,1))
cccccccc
ccccccccc     write(*,'(3i5,3(/,10x,2f10.3))') it,jt,k,                   !deubg-delete
ccccccccc    &                   table(it,jt,k-1,1),table(it,jt,k-1,2),   !deubg-delete
ccccccccc    &                   table(it,jt,k  ,1),table(it,jt,k  ,2),   !deubg-delete
ccccccccc    &                   xv,yv                                    !deubg-delete
cccccccc
cccccccc      rintp = yv
cccccccc      return
cccccccc      end
cccccccc


!==========================================================================
c this function interpolates from multiple tables

      function rintpm(iftn,ktab,xv)
      ! xv    - x value to interpolate
      ! iftn  - type of tables, 0 for flow, 1 for TDS, 2 for TIN
      ! ktab  - table number

      include 'ro_dim_3.max'

      include 'tbprpt.var'

      if(iftn.eq.0) then
        do k=2,ntabf(ktab)
          if(xv.le.fwtab(ktab,k,1)) go to 10
        end do
        k = ntabf(ktab)
 10     continue
        yv = fwtab(ktab,k-1,2) + (fwtab(ktab,k,2)-fwtab(ktab,k-1,2))
     &                         / (fwtab(ktab,k,1)-fwtab(ktab,k-1,1))
     &                         * (xv             -fwtab(ktab,k-1,1))
      else
        do k=2,ntabq(ktab,iftn)      
          if(xv.le.fqtab(ktab,iftn,k,1)) goto 20
        end do
        k = ntabq(ktab,iftn)
 20     continue
        yv = fqtab(ktab,iftn,k-1,2) 
     &           + (fqtab(ktab,iftn,k,2)-fqtab(ktab,iftn,k-1,2))
     &           / (fqtab(ktab,iftn,k,1)-fqtab(ktab,iftn,k-1,1))
     &           * (xv                  -fqtab(ktab,iftn,k-1,1))
      end if

      rintpm = yv
      return
      end



!==========================================================================
      subroutine xlsdate(iopt,iyr,imn,idy,nexday,julday)
      ! this subroutine mimics excel date informaiton storage
      ! nexday counts from 1/1/1900
      ! julday  is julian day, counting from Jan 1 each year
      ! if iopt = 0, calculate nexday from iyr,imn,idy
      !           1, calculate iyr,imn,idy from nexday
      ! note that Excel treats 1900 as a leap year, which is an error
      ! this program mimics excel and treats 1900 as a leap year

      dimension    ndays(12,0:1),nsad(0:3)
      data ndays/31,28,31,30,31,30,31,31,30,31,30,31,
     &           31,29,31,30,31,30,31,31,30,31,30,31/
      data nsad  /0,366,731,1096/


      if(iopt.eq.0) then
        iyr1 = int((iyr-1900)/4)
        iyr2 = iyr-1900-iyr1*4

        leapyr = 0
        if(mod(iyr,4).eq.0) leapyr=1
c       if(mod(iyr,100).eq.0) leapyr=0   for excel
c       if(mod(iyr,400).eq.0) leapyr=1

        nexday = iyr1*1461 + nsad(iyr2)

        julday = 0
        if(imn.gt.1) then
          do jmn=1,imn-1
            julday = julday + ndays(jmn,leapyr)
          end do
        end if
        julday = julday + idy
        nexday = nexday + julday
        return
      else
        iyr1 = int(nexday/1461)
        idy = nexday - iyr1*1461
        do jyr=1,3
          iyr2=jyr-1
          if(idy.lt.nsad(jyr)) goto 10
        end do
        iyr2=3
 10     idy = idy - nsad(iyr2)
        julday = idy
        iyr = 1900+iyr1*4+iyr2

        leapyr = 0
        if(mod(iyr,4).eq.0) leapyr=1
        do imn=1,12
          if(idy.lt.ndays(imn,leapyr)) goto 20
          idy = idy - ndays(imn,leapyr)
        end do

 20     if(idy.eq.0) then
          imn = imn-1
          if(imn.gt.0) then
            idy = ndays(imn,leapyr)
          else
            iyr = iyr - 1
            imn = 12
            idy = 31
          end if
        end if
        return
      end if
      end




!==========================================================================

      subroutine vrinit()

      include 'RO_DIM_3.max'
c      include 'lkprpt.var'
      include 'lkprpt2.var'
      include 'rvprpt.var'

c     initialize variables
      lkname = ' NO LINK'            ! (MXLK)

      usn = 0
      dsn = 0                        ! (MXLK)
      convtype = 0                   ! (MXLK)
      incoml = 0                     ! (MXLK,MXIL)
      rxn = 0.                       ! (MXLK)
      cqlout = 0.                    ! (MXLK,4,0:MXPS)
      jqlmax = 0                     ! (MXLK)
      qlmax = 0.                     ! (MXLK)
      linkrch = 0                    ! (MXLK)
      lkrch = 0                      ! (MXLK)
      rvrch = 0                      ! (MXTB,3)
      ptrch = 0                      ! (MXLK)
      ptrrv = 0                      ! (MXLK)
      ptrdv = 0                      ! (MXLK)

      ndname = ' NO NODE'            ! (MXLK)
      nha = 0                        ! (MXLK)
      nodha = 0                      ! (MXLK,MXHN)
      rmult = 0.                     ! (MXLK,MXHN)
      nincl = 0                      ! (MXLK)
      qout = 0.                      ! (MXLK,0:MXPS)
      cqout = 0.                     ! (MXLK,4,0:MXPS)

      rvelev = 0.                    ! (MXRV,MXTB)
      rvarea = 0.                    ! (MXRV,MXTB)
      rvstor = 0.                    ! (MXRV,MXTB)
      rvout = 0.                     ! (MXRV,MXTB,4)
      rvpr = 0.                      ! (MXRV,MXTB)
      alimit = 0.                    ! (MXRV)
      slimit = 0.                    ! (MXRV,3)
      evwts = 0.                     ! (MXRV,MXEV)
      qin = 0.                       ! (MXLK)
      percl = 0.                     ! (MXRV,0:MXPS)
      storl = 0.                     ! (MXRV,0:MXPS)
      evapl = 0.                     ! (MXRV)
      out1l = 0.                     ! (MXRV)
      out2l = 0.                     ! (MXRV)
      cstorl = 0.                    ! (MXRV,4,0:MXPS)
      npts = 0                       ! (MXRV)
      lkatrv = 0.                    ! (MXRV)
      rdlink = '        '       !0.  ! (MXRV,2)
      linkrv = 0                     ! (MXRV,2)

      return
      end
