      BLOCK DATA ARDATA
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
cc      SAVE /ARPRNT/
      SAVE   
      DATA ARPAR1/1.19, 99 * 0.0/
      DATA IAPAR2/3, 49 * 0/
      DATA ARINT1/100 * 0.0/
      DATA IAINT2/50 * 0/
      END
c=======================================================================
c.....Routine borrowed from ZPC.
c.....double precision  is modified to real*4.
cbz1/29/99
c      subroutine index1(n, m, arrin, indx)
