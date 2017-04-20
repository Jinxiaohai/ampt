      double precision function ran1(idum)
*     return a uniform random deviate between 0.0 and 1.0. set idum to 
*     any negative value to initialize or reinitialize the sequence.
      implicit double precision (a-h, o-z)
      dimension r(97)
      common /rndm1/ number
cc      SAVE /rndm1/
      parameter (m1 = 259200, ia1 = 7141, ic1 = 54773, rm1 = 1d0 / m1)
      parameter (m2 = 134456, ia2 = 8121, ic2 = 28411, rm2 = 1d0 / m2)
      parameter (m3 = 243000, ia3 = 4561, ic3 = 51349)
clin-6/23/00 save ix1-3:
clin-10/30/02 r unsaved, causing wrong values for ran1 when compiled with f77:
cc      SAVE ix1,ix2,ix3,r
      SAVE   
      data iff/0/
      if (idum .lt. 0 .or. iff .eq. 0) then
         iff = 1
         ix1 = mod(ic1 - idum, m1)
         ix1 = mod(ia1 * ix1 + ic1, m1)
         ix2 = mod(ix1, m2)
         ix1 = mod(ia1 * ix1 + ic1, m1)
         ix3 = mod(ix1, m3)
         do 11 j = 1, 97
            ix1 = mod(ia1 * ix1 + ic1, m1)
            ix2 = mod(ia2 * ix2 + ic2, m2)
            r(j) = (dble(ix1) + dble(ix2) * rm2) * rm1
 11         continue
         idum = 1
      end if
      ix1 = mod(ia1 * ix1 + ic1, m1)
      ix2 = mod(ia2 * ix2 + ic2, m2)
      ix3 = mod(ia3 * ix3 + ic3, m3)
clin-7/01/02       j = 1 + (97 * i x 3) / m3
      j=1+(97*ix3)/m3
clin-4/2008:
c      if (j .gt. 97 .or. j .lt. 1) pause
      if (j .gt. 97 .or. j .lt. 1) print *, 'In zpc ran1, j<1 or j>97',j
      ran1 = r(j)
      r(j) = (dble(ix1) + dble(ix2) * rm2) * rm1
clin-6/23/00 check random number generator:
      number = number + 1
c      if(number.le.100000) write(99,*) 'number, ran1=', number,ran1
      return
      end
