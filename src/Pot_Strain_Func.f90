
SUBROUTINE Pot_Strain_Func(ntide, theta, phi)
   use param_Y
   use prec
   implicit none
   integer :: ntide, i
   real(kind=dp) :: theta, phi, P20, P22, P20_t, P22_t, P20_tt, P22_tt
   real(kind=dp) :: a, b, fact_e, Rp, phi2
   complex(kind=dp) :: Pe, Pe_t, Pe_p, Pe_tp, Pe_tt, Pe_pp, rc, OO1, ciw
   complex(kind=dp), dimension(ntide) :: ur,uth,uph,Pot,sigrr,sigtt,sigpp,sigtp, sigrt, sigrp 
   complex(kind=dp), dimension(ntide) :: epsrr,epstt,epspp,epstp, epsrp, epsrt
   complex(kind=dp), parameter :: c0=(0._dp,0._dp), c2=(2._dp,0._dp)

OO1        = CMPLX(O*(O+1._dp),KIND=dp)
ciw        = (0._dp,1._dp)*CMPLX(w,KIND=dp)

!
!--- LEGENDRE POLYNOMIAL FUNCTIONS AND DERIVATIVE

phi2=phi*2._dp
Rp=r(nc)
P20   = .5_dp*(3._dp*COS(theta)**2-1._dp)
P22   = 3._dp*SIN(theta)**2
fact_e = CMPLX((Rp*w)**2*ecc,KIND=dp)

!Fourier coeff. for:
!-->potential
a = -1.5_dp*(P20 - .5_dp*P22*COS(phi2))
b = P22*SIN(phi2)
Pe = fact_e*CMPLX(a,-b,KIND=dp)

!-->colatitude-derivative of the potential:
a = 2.25_dp*SIN(2._dp*theta)*(1._dp+COS(phi2))
b = 3._dp*SIN(2._dp*theta)*SIN(phi2)
Pe_t = fact_e*CMPLX(a,-b,KIND=dp)

!-->colatitude-double derivative of the potential:
a = 4.5_dp*COS(2._dp*theta)*(1._dp+COS(phi2)) 
b = 6._dp*COS(2._dp*theta)*SIN(phi2)
Pe_tt = fact_e*CMPLX(a,-b,KIND=dp)

!-->longitude-derivative of the potential:
a = -1.5_dp*P22*SIN(phi2) 
b =   2._dp*P22*COS(phi2)
Pe_p = fact_e*CMPLX(a,-b,KIND=dp)

!-->longitude-double derivative of the potential:
a = -3._dp*P22*COS(phi2)
b = -4._dp*P22*SIN(phi2)
Pe_pp = fact_e*CMPLX(a,-b,KIND=dp)

!-->colatitude-longitude-cross derivative of the potential:
a = -4.5_dp*SIN(2._dp*theta)*SIN(phi2)
b = 6._dp*SIN(2._dp*theta)*COS(phi2) 
Pe_tp = fact_e*CMPLX(a,-b,KIND=dp)


!
! Displacements, potential, stress and strain
!

do i=2,ntide

rc = CMPLX(r(i),KIND=dp)

ur(i)=Yi(1,i)*Pe
uth(i)=Yi(3,i)*Pe_t
uph(i)=Yi(3,i)*Pe_p/CMPLX(sin(theta),KIND=kind(dp))

Pot(i)=Yi(5,i)*Pe

sigrr(i)=Yi(2,i)*Pe

sigtt(i) =(lamc(i)*dy1dr(i) + (lamc(i)+c2*muc(i))/rc*(c2*Yi(1,i)-OO1*Yi(3,i)) &
                                         - c2*muc(i)/rc*Yi(1,i)) *Pe &
                  - c2*muc(i)/rc*Yi(3,i) * (Pe_t/CMPLX(TAN(theta),KIND=dp) &
                         +Pe_pp/CMPLX(SIN(theta)**2,KIND=dp) )

sigpp(i) =( lamc(i)*dy1dr(i) + (lamc(i)+c2*muc(i))/rc*(c2*Yi(1,i)-OO1*Yi(3,i)) - c2*muc(i)/rc*Yi(1,i) )*Pe &
                  - c2*muc(i)/rc*Yi(3,i) * Pe_tt

sigrt(i)  = Yi(4,i) * Pe_t

sigtp(i) = c2*muc(i)/rc*Yi(3,i) * ( Pe_tp/CMPLX(SIN(theta),KIND=dp) -Pe_p/CMPLX(TAN(theta)*SIN(theta),KIND=dp) )

sigrp(i)  = Yi(4,i)/CMPLX(SIN(theta),KIND=dp) * Pe_p

epsrr(i) = ciw * dy1dr(i) * Pe
epstt(i) = ciw * (Yi(3,i)*Pe_tt + Yi(1,i)*Pe)/rc
epspp(i) = ciw * (Yi(3,i)*Pe_pp/CMPLX(SIN(theta)**2,KIND=dp) &
                 +Yi(3,i)*Pe_t/CMPLX(TAN(theta),KIND=dp) + Yi(1,i)*Pe)/rc
epsrt(i)  = ciw * Yi(4,i) * Pe_t/muc(i)
epstp(i) = ciw * c2/rc*Yi(3,i) * (Pe_tp/CMPLX(SIN(theta),KIND=dp)  &
                     -Pe_p/CMPLX(TAN(theta)*SIN(theta),KIND=dp) )
epsrp(i)  = ciw * Yi(4,i)/CMPLX(SIN(theta),KIND=dp) * Pe_p/muc(i)

strain_energy(i)=CMPLX(.5_dp,KIND=dp) * ( CONJG(sigrr(i))*epsrr(i)+CONJG(sigtt(i))*epstt(i)+CONJG(sigpp(i))*epspp(i) &
                                         +CONJG(sigrt(i))*epsrt(i)+CONJG(sigrp(i))*epsrp(i)+CONJG(sigtp(i))*epstp(i))


htide=21._dp/10._dp*(w**5)*(Rp**4)*(ecc**2)/(rc**2)*REAL(HH(2,i))*IMAG(muc(i))


enddo





END SUBROUTINE Pot_Strain_Func 


SUBROUTINE Pot_Strain_Func_out(ntide, theta, phi)
   use param_Y
   use prec
   implicit none
   integer :: ntide, i
   real(kind=dp) :: theta, phi, P20, P22, P20_t, P22_t, P20_tt, P22_tt
   real(kind=dp) :: a, b, fact_e, Rp, phi2
   complex(kind=dp) :: Pe, Pe_t, Pe_p, Pe_tp, Pe_tt, Pe_pp, rc, OO1, ciw
   complex(kind=dp), dimension(ntide) :: ur,uth,uph,Pot,sigrr,sigtt,sigpp,sigtp, sigrt, sigrp 
   complex(kind=dp), dimension(ntide) :: epsrr,epstt,epspp,epstp, epsrp, epsrt
   complex(kind=dp), parameter :: c0=(0._dp,0._dp), c2=(2._dp,0._dp)

OO1        = CMPLX(O*(O+1._dp),KIND=dp)
ciw        = (0._dp,1._dp)*CMPLX(w,KIND=dp)

!
!--- LEGENDRE POLYNOMIAL FUNCTIONS AND DERIVATIVE

phi2=phi*2._dp
Rp=r(nc)
P20   = .5_dp*(3._dp*COS(theta)**2-1._dp)
P22   = 3._dp*SIN(theta)**2
fact_e = CMPLX((Rp*w)**2*ecc,KIND=dp)

!Fourier coeff. for:
!-->potential
a = -1.5_dp*(P20 - .5_dp*P22*COS(phi2))
b = P22*SIN(phi2)
Pe = fact_e*CMPLX(a,-b,KIND=dp)

!-->colatitude-derivative of the potential:
a = 2.25_dp*SIN(2._dp*theta)*(1._dp+COS(phi2))
b = 3._dp*SIN(2._dp*theta)*SIN(phi2)
Pe_t = fact_e*CMPLX(a,-b,KIND=dp)

!-->colatitude-double derivative of the potential:
a = 4.5_dp*COS(2._dp*theta)*(1._dp+COS(phi2)) 
b = 6._dp*COS(2._dp*theta)*SIN(phi2)
Pe_tt = fact_e*CMPLX(a,-b,KIND=dp)

!-->longitude-derivative of the potential:
a = -1.5_dp*P22*SIN(phi2) 
b =   2._dp*P22*COS(phi2)
Pe_p = fact_e*CMPLX(a,-b,KIND=dp)

!-->longitude-double derivative of the potential:
a = -3._dp*P22*COS(phi2)
b = -4._dp*P22*SIN(phi2)
Pe_pp = fact_e*CMPLX(a,-b,KIND=dp)

!-->colatitude-longitude-cross derivative of the potential:
a = -4.5_dp*SIN(2._dp*theta)*SIN(phi2)
b = 6._dp*SIN(2._dp*theta)*COS(phi2) 
Pe_tp = fact_e*CMPLX(a,-b,KIND=dp)


!
! Displacements, potential, stress and strain
!

do i=2,ntide

rc = CMPLX(r(i),KIND=dp)

ur(i)=Yi(1,i)*Pe
uth(i)=Yi(3,i)*Pe_t
uph(i)=Yi(3,i)*Pe_p/CMPLX(sin(theta),KIND=kind(dp))

Pot(i)=Yi(5,i)*Pe

sigrr(i)=Yi(2,i)*Pe

sigtt(i) =(lamc(i)*dy1dr(i) + (lamc(i)+c2*muc(i))/rc*(c2*Yi(1,i)-OO1*Yi(3,i)) &
                                         - c2*muc(i)/rc*Yi(1,i)) *Pe &
                  - c2*muc(i)/rc*Yi(3,i) * (Pe_t/CMPLX(TAN(theta),KIND=dp) &
                         +Pe_pp/CMPLX(SIN(theta)**2,KIND=dp) )

sigpp(i) =( lamc(i)*dy1dr(i) + (lamc(i)+c2*muc(i))/rc*(c2*Yi(1,i)-OO1*Yi(3,i)) - c2*muc(i)/rc*Yi(1,i) )*Pe &
                  - c2*muc(i)/rc*Yi(3,i) * Pe_tt

sigrt(i)  = Yi(4,i) * Pe_t

sigtp(i) = c2*muc(i)/rc*Yi(3,i) * ( Pe_tp/CMPLX(SIN(theta),KIND=dp) -Pe_p/CMPLX(TAN(theta)*SIN(theta),KIND=dp) )

sigrp(i)  = Yi(4,i)/CMPLX(SIN(theta),KIND=dp) * Pe_p

epsrr(i) = ciw * dy1dr(i) * Pe
epstt(i) = ciw * (Yi(3,i)*Pe_tt + Yi(1,i)*Pe)/rc
epspp(i) = ciw * (Yi(3,i)*Pe_pp/CMPLX(SIN(theta)**2,KIND=dp) &
                 +Yi(3,i)*Pe_t/CMPLX(TAN(theta),KIND=dp) + Yi(1,i)*Pe)/rc
epsrt(i)  = ciw * Yi(4,i) * Pe_t/muc(i)
epstp(i) = ciw * c2/rc*Yi(3,i) * (Pe_tp/CMPLX(SIN(theta),KIND=dp)  &
                     -Pe_p/CMPLX(TAN(theta)*SIN(theta),KIND=dp) )
epsrp(i)  = ciw * Yi(4,i)/CMPLX(SIN(theta),KIND=dp) * Pe_p/muc(i)

strain_energy(i)=CMPLX(.5_dp,KIND=dp) * ( CONJG(sigrr(i))*epsrr(i)+CONJG(sigtt(i))*epstt(i)+CONJG(sigpp(i))*epspp(i) &
                                         +CONJG(sigrt(i))*epsrt(i)+CONJG(sigrp(i))*epsrp(i)+CONJG(sigtp(i))*epstp(i))


htide=21._dp/10._dp*(w**5)*(Rp**4)*(ecc**2)/(rc**2)*REAL(HH(2,i))*IMAG(muc(i))

enddo



write(206,*) 90.d0-theta/pi*180.D0, phi*180.d0/pi, real(sigtt(ntide)), real(sigpp(ntide)), real(sigtp(ntide))
write(206,*) 90.d0-theta/pi*180.D0, phi*180.d0/pi, dimag(sigtt(ntide)), dimag(sigpp(ntide)), dimag(sigtp(ntide))

write(207,*) 90.d0-theta/pi*180.D0, phi*180.d0/pi, real(ur(ntide)), real(uth(ntide)), real(uph(ntide))
write(207,*) 90.d0-theta/pi*180.D0, phi*180.d0/pi, dimag(ur(ntide)), dimag(uth(ntide)), dimag(uph(ntide))



END SUBROUTINE Pot_Strain_Func_out 




