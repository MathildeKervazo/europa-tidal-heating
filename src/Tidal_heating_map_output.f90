SUBROUTINE Tidal_heating_map_output

    use Param_Y
    use param_rheo
    use prec

implicit none

integer:: i, i_output
real(kind=dp) :: R_output

open(600,file='Htide_Output/Tot/Htide_R_colat_per.txt')

lont=0/180._dp*pi
do i=1,nmant
colatt=0.5_dp/180._dp*pi
do while (colatt.lt.(179.5_dp/180._dp*pi))
if (vs(i).gt. 0._dp) then
call Pot_Strain_Func(i, colatt, lont)
write(600,*) R(i), colatt/pi*180._dp, lont/pi*180._dp, real(strain_energy(i))+hrad
endif
colatt=colatt+5._dp/180._dp*pi
enddo
colatt=179.5_dp/180._dp*pi
if (vs(i).gt. 0._dp) then
call Pot_Strain_Func(i, colatt, lont)
write(600,*) R(i), colatt/pi*180._dp, lont/pi*180._dp, real(strain_energy(i))+hrad
endif
enddo

open(700,file='Htide_Output/Tot/Htide_colat_lon_per.txt')
R_output=1302.d3
do i=1,nmant
if (R(i).le.R_output) i_output=i
enddo
if ((R_output-R(i_output)).gt.(R(i_output+1)-R_output)) i_output=i_output+1
write (*,*) 'i_output=', i_output, 'R=', R(i_output)/1.D3
if (vs(i_output).eq.0._dp) then
write(*,*) 'R_output in fluid layer, no dissipation computed there'
endif
i=i_output
colatt=0.5_dp/180._dp*pi
do while (colatt.lt.(179.5_dp/180._dp*pi))
lont=0._dp
do while (lont.lt.(2._dp*pi))
call Pot_Strain_Func(i, colatt, lont)
write(700,*) R(i), colatt/pi*180._dp, lont/pi*180._dp, real(strain_energy(i))+hrad
lont=lont+5._dp/180._dp*pi
enddo
write(*,*) i, 'colat=', colatt/pi*180._dp, 'R_output=', R(i)
colatt=colatt+5._dp/180._dp*pi
enddo

open(800,file='Htide_Output/Tot/Htide_colatt_lont_radius.dat')
do i=1,nmant
colatt=0.5_dp/180._dp*pi
do while (colatt.lt.(179.5_dp/180._dp*pi))
lont=0._dp
do while (lont.lt.(2._dp*pi))
if (vs(i).gt. 0._dp) then
call Pot_Strain_Func(i, colatt, lont)
write(800,*) R(i), colatt/pi*180._dp, lont/pi*180._dp, real(strain_energy(i))+hrad
endif
lont=lont+5._dp/180._dp*pi
enddo
colatt=colatt+5._dp/180._dp*pi
enddo
enddo

!!!!!!!!!!!!!!!!!!!Total!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(840,file='Htide_Output/Tot/Htide_colatt_lont_radius_tot.dat')
do i=ncore+1, nmant
colatt=0.5_dp/180._dp*pi
do while (colatt.lt.(179.5_dp/180._dp*pi))
lont=0._dp
do while (lont.lt.(2._dp*pi))
call Pot_Strain_Func(i, colatt, lont)
write(840,*) R(i), colatt/pi*180._dp, lont/pi*180._dp, real(strain_energy(i))+hrad
lont=lont+5._dp/180._dp*pi
enddo
colatt=colatt+5._dp/180._dp*pi
enddo
enddo


END SUBROUTINE Tidal_heating_map_output
