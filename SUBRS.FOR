CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                          C
C  Rainfall-Runoff-Router-Rootzone Model                                   C
C  Copy Right (C) 1996 - 2017 Wildermuth Environmental, Inc.               C
C                                                                          C
C  This program is free software: you can redistribute it and/or modify    C
C  it under the terms of the GNU General Public License as published by    C
C  the Free Software Foundation, either version 3 of the License, or       C
C  (at your option) any later version.                                     C
C                                                                          C
C  This program is distributed in the hope that it will be useful,         C
C  but WITHOUT ANY WARRANTY; without even the implied warranty of          C
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           C
C  GNU General Public License for more details.                            C
C                                                                          C
C  You should have received a copy of the GNU General Public License       C
C  along with this program.  If not, see <http://www.gnu.org/licenses/>.   C
C                                                                          C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!=======================================================================
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



!=======================================================================
      subroutine cmpdate(iyr,imo,idy, ju,jyr,jmo,jdy)

      if(imo.eq.jmo .and. idy.eq.jdy .and. iyr.eq.jyr) return

      write(*,'(a)') ' Date problem in input data file'
      write(*,'(i7,2i2,a   )') iyr,imo,idy,'  Simulation Date'
      write(*,'(i7,2i2,a,i3)') jyr,jmo,jdy,'  Date in file unit',ju
      stop
      end



!=======================================================================
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


!=======================================================================
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


!=======================================================================
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



!=======================================================================
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



!=======================================================================

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


!=======================================================================

      subroutine openhdat(iu1,iu2,nsy,ismonth,ifile,fmt)

      character*50 ifile
      character*40 fmt

      call openifile(iu1,iu2,ifile)
      read(iu2,'(a40)') fmt

      do while (.true.)
        read(iu2,fmt,err=1,end=2) m1,m2,m3
        if(m1.eq.nsy .and. m2.eq.ismonth .and. m3.eq.1) then   ! CY
          backspace(iu2)
          goto 3
        end if
 1      continue ! if daily data is in multiple line, skip error reading date
      end do
 2    write(*,'(//,a,i5/)') ' Check data for input file',iu2
      stop

 3    write(6,'(a,a)') ifile,'is ready'
      return
      end


!=======================================================================
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



!=======================================================================
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
      
      subroutine openofile(iu1,iu2,ofile,cha,chb)
      character*(*)  cha, chb
      character  ofile*50, bl50*50

      bl50 = '                                                  '
      read(iu1,'(a50)') ofile
      if(ofile.eq.bl50) then
        write(*,'(a,i2)') ' Blank file name for unit ',iu2
        stop
      end if
      write(*,'(/,a,a)') ' Opening ',ofile
      open(iu2,file=ofile,status='UNKNOWN')
      if(cha.ne.bl50) write(iu2,'(a)') cha
      if(chb.ne.bl50) write(iu2,'(a,a)') 'Response file = ',chb
      return
      end

!=======================================================================

      subroutine openmfile(iu1,iu2,nha,ofile,cha,chb)

      include 'rf_dim_3.max'  

      common /havar/haname(MXHA),haarea(MXHA)
      character*8  haname
      character*(*)  cha, chb
      character*50   ofile

      call openofile(iu1,iu2,ofile,cha,chb)
      write(iu2,'(a8,25a8,/,(8x,25a8))') 'HANAME',
     &                                  (haname(iha),iha=1,nha)
      write(iu2,'(a8,25f8.0,/,(8x,25f8.0))') 'HAAREA',
     &                                   (haarea(iha),iha=1,nha)
      write(6,'(a,a)') ofile,'is ready'
      return
      end
   
!=======================================================================

      ! subroutine to get timeseries file 
      subroutine tsfpos(iu,isyr,ismo,isdy,fname)
      character*50 fname 
      
      do while (.true.)
        read(iu,'(i4,2i2)') iyr,imo,idy
        if(iyr.eq.isyr .and. imo.eq.ismo .and. idy.eq.isdy) then
          backspace(iu)
          write(*,'(a,a)') ' Ready file --- ',fname
          goto 1
        end if
      end do
      write(*,'(a,a)') ' Failed to ready file --- ',fname
      stop
  1   return
      end

!=======================================================================
      
      subroutine openqfile(iu1,iu2,ofile,yesno)
c     to op optional file
c     useful to put at the end of main file for optional output
      character*50 ofile,bl50
      logical*1    yesno

      bl50 = '                                                  '
      read(iu1,'(a50)',end=10,err=10) ofile
      if(ofile.eq.bl50) goto 10
      open(iu2,file=ofile,status='UNKNOWN')
      write(*,'(a,a)') ' Opening ',ofile
      yesno = .TRUE.
      return
 10   yesno = .FALSE.
      return
      end      
      
!=======================================================================

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

      write(*,'(//,a)') ' Problem in StageC'  
      write(*,'(f10.0,7f10.4)') q,dsave,d,slop,xn,b,z

      d = dsave/2.
      dsave = d
      goto 1
      write(*,'(a)') ' Increase bottom width and/or reduce side slope'
      uperc = -99999.
      return 

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


!=======================================================================

      subroutine stagep()
c     subroutine stagep(slop,xn,diam,q,d)
      return
      end



!=======================================================================
c this function interpolates from multiple tables

      function rintp(it,jt,xv)
      ! xv    - x value to interpolate
      ! it    - type of tables, 0 for flow, 1 for TDS, 2 for TIN
      ! jt    - tables number


      include 'ro_dim_4.max'

      include 'tbprpt.var'


      do k=2,ntab(it,jt)
        if(xv.le.table(it,jt,k,1)) go to 10
      end do
      k = ntab(it,jt)
 10   continue

      yv = table(it,jt,k-1,2) + (table(it,jt,k,2)-table(it,jt,k-1,2))
     &                        / (table(it,jt,k,1)-table(it,jt,k-1,1))
     &                        * (xv              -table(it,jt,k-1,1))

      rintp = yv
      return
      end





!=======================================================================
      subroutine xlsdate(iopt,iyr,imn,idy,nexday,julday)
      ! this subroutine mimics excel date informaiton storage
      ! nexday counts from 1/1/1900
      ! julday  is julian day, counting from Jan 1 each year
      ! if iopt = 0, calculate nexday from iyr,imn,idy
      !           1, calculate iyr,imn,idy from isdate
      ! note that Excel treats 1900 as leap year, which is an error
      ! this program mimics excel and treats 1900 as leap year

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
c20     return
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




!=======================================================================

      subroutine vrinit()

      include 'RO_DIM_4.max'
      include 'lkprpt.var'
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
      rdlink = 0.                    ! (MXRV,2)
      linkrv = 0                     ! (MXRV,2)

      return
      end
