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

c This routine reads and sets up output files and variables
c
      subroutine rdoutp()

      include 'RO_DIM_4.MAX'

      include 'GenParm.VAR'
      include 'OVprpt.VAR'
      include 'lkprpt.VAR'

      character*8 ch8(100)
      character*50 ffile

c     read location ids for daily flow

      read(1,'(i8)') nplk
      write(*,'(a,i5)') ' NPLK =',nplk

      if(nplk.gt.0) then
        write(6,'(/,a)') 'Daily flow at following locations, NPLK'
        write(6,'(a)')   'LinkName  Number'
        do i=1,nplk
          read(1,'(a8)') ch8(i)
          call ctrim8r(ch8(i))      ! trim the name and convert to upper case
          do j=1,nlinks
            if(ch8(i).eq.lkname(j)) then
              iplk(i,1) = j
              iplk(i,2) = 1
              goto 1462
            end if
          end do
          do j=1,nnodes
            if(ch8(i).eq.ndname(j)) then
              iplk(i,1) = j
              iplk(i,2) = 2
              goto 1462
            end if
          end do

          write(*,'(a,a8,a)')' Link or node ',ch8(i),' is not available'
          ierror = .true.
 1462     continue
          write(6,'(a8,2i4)') ch8(i),iplk(i,1),iplk(i,2)
        end do

c       read file name to print flow at iplk(i, ) and open file
c       print at each link using
c
        read(1,'(a40)') ffile
        open(61,file=ffile,status='unknown')

        write(61,'(i8)') nplk
        write(61,'(8x,30a8)') (ch8(i), i=1,nplk)
        write(61,'(8x,30i8)') (iplk(i,1),i=1,nplk)
        write(61,'(8x,30i8)') (iplk(i,2),i=1,nplk)
        if(qopt) then
          read(1,'(a40)') ffile   ! TDS output file
          open(62,file=ffile,status='unknown')
          write(62,'(i8)') nplk
          write(62,'(8x,30a8)') (ch8(i), i=1,nplk)
          write(62,'(8x,30i8)') (iplk(i,1),i=1,nplk)
          write(62,'(8x,30i8)') (iplk(i,2),i=1,nplk)
          read(1,'(a40)') ffile   ! TIN output file
          open(63,file=ffile,status='unknown')
          write(63,'(i8)') nplk
          write(63,'(8x,30a8)') (ch8(i), i=1,nplk)
          write(63,'(8x,30i8)') (iplk(i,1),i=1,nplk)
          write(63,'(8x,30i8)') (iplk(i,2),i=1,nplk)
        end if
      end if !nplk

c     read location ids for monthly flow

      read(1,'(i8)') npnd
      write(*,'(a,i5)') ' NPND =',npnd

      if (npnd.gt.MXPL) then
        write(*,'(a)') ' npnd > MXPL, modify dimension file'
        stop
      else if(npnd.gt.0) then
        write(6,'(/,a)') 'Monthly flow at following locations, NPND'
        write(6,'(a)') 'printed for total and points loads, NPND'
        write(6,'(a)')   'NodeName  Number'
        do i=1,npnd
          read(1,'(a8)') ch8(i)
          call ctrim8r(ch8(i))      ! trim the name and convert to upper case
          do in=1,nnodes
            if(ch8(i).eq.ndname(in)) then
              ipnd(i) = in
              goto 1463
            end if
          end do
          write(*,'(a,a8,a)')' Node ',ch8(i),' is not available'
          ierror = .true.
 1463     continue
          write(6,'(a8,2i4)') ch8(i),ipnd(i)
        end do

c       read file name to print flow at ipnd(i) and open file
c       print at each link using
c
        read(1,'(a40)') ffile
        open(64,file=ffile,status='unknown')
        write(64,'(a)') ' Total monthly flow passing the nodes'
        write(64,'(a)') ' Units (ac-ft)'
        write(64,'(i6)') npnd
        write(64,'(6x,30a8)') (ch8(i), i=1,npnd)
        write(64,'(1x)')
        if(qopt) then
          read(1,'(a40)') ffile   ! TDS output file
          open(65,file=ffile,status='unknown')
          read(1,'(a40)') ffile   ! TIN output file
          open(66,file=ffile,status='unknown')
        end if
      end if !npnd

c     read reach names
      read(1,'(i8)') nprch
      write(*,'(a,i5)') ' NPRCH =',nprch

      if (nprch.gt.MXRH) then
        write(*,'(a)') ' nprch > MXRH'
        stop
      else if(nprch.gt.0) then
        do i=1,nprch
          read(1,'(a8)') rchname(i)
          call ctrim8r(rchname(i))
        end do
        read(1,'(a40)') ffile
        open(67,file=ffile,status='unknown')
c modified to be same format as file 81
        write(67,'(a)') ' Total monthly percolation'
        write(67,'(a)') ' Units|(ac-ft)'
        write(67,'(i5)') nprch
        write(67,'(6x,(50a8))') (rchname(i), i=1,nprch)
        write(67,'(1x)')
c end

        if(qopt) then
          read(1,'(a40)') ffile   ! TDS output file
          open(68,file=ffile,status='unknown')
          read(1,'(a40)') ffile   ! TIN output file
          open(69,file=ffile,status='unknown')
        end if

        write(6,'(/,a)') 'At following reaches, monthly percolation'
        write(6,'(a)') ' will be printed, NPRCH'
        write(6,'(a)')   'LinkName  Number'
        do i=1,nprch
          write(6,'(a8)') rchname(i)
        end do
      end if

c     read recharge basins ids for monthly percolation and evaporation
      read(1,'(i8)') nqls
      write(*,'(a,i5)') ' NQLS =',nqls

      if(nqls.gt.0) then
        write(6,'(/,a,a)') 'At following links, monthly percolation',
     &                     ' will be printed, NQLS'
        write(6,'(a)')   'LinkName  Number'
        do i=1,nqls
          read(1,'(a8)') ch8(i)
          call ctrim8r(ch8(i))      ! trim the name and convert to upper case
          do il=1,nlinks
            if(ch8(i).eq.lkname(il)) then
              ipql(i) = il
              goto 81
            end if
          end do
          write(*,'(a,a8,a)')' Link ',ch8(i),' is not available'
          ierror = .true.
          xqevp(i) = 0.
 81       xqloss(i) = 0.
          write(6,'(a8,2i4)') ch8(i),ipql(i)
        end do
        read(1,'(a40)') ffile
        open(81,file=ffile,status='unknown')
        write(81,'(a)') ' Total monthly percolation'
        write(81,'(a)') ' Units|(ac-ft)'
        write(81,'(i5)') nqls
        write(81,'(6x,(50a8))') (ch8(i), i=1,nqls)
        write(81,'(6x,(50i8))') (ipql(i),i=1,nqls)
        read(1,'(a40)') ffile
        open(82,file=ffile,status='unknown')
        write(82,'(a)') ' Total monthly evaporation'
        write(82,'(a)') ' Units|(ac-ft)'
        write(82,'(i5)') nqls
        write(82,'(6x,(50a8))') (ch8(i), i=1,nqls)
        write(82,'(6x,(50i8))') (ipql(i),i=1,nqls)
      end if

      return
      end
