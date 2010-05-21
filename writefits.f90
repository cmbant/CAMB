!subroutine to export Cls in FITS format for HEALPix 1.2
!Antony Lewis July 2003


 subroutine WriteFitsCls(Clsfile, lmx)
  use CAMB
  use head_fits, ONLY : add_card 
  use fitstools, only : write_asctab
  implicit none
  character(LEN=*), INTENT(IN) ::  Clsfile
  integer, INTENT(IN) :: lmx
  CHARACTER(LEN=80), DIMENSION(1:120) :: header
  INTEGER nlheader,i, j
  real, allocatable, dimension (:,:) :: clout,allcl
  real(dl) :: fac,  PowerVals(20)
  character(Len=40) :: unitstr
  character(Len=8) :: PowerKeys(20)
  logical COBEnorm


  if (CP%InitPower%nn>1) write(*,*) &
       'Warning: FITS file contains result for first power spectrum only'

  allocate(clout(2:lmx,1:4))
   
  call CAMB_GetCls(clout, lmx, 1, .false.)
  !HealPix 1.2 uses E-B conventions

  if (CP%OutputNormalization == outCOBE) then
     fac=2*pi*CP%tcmb**2
  else
     if (CP%OutputNormalization >=2) then
      fac=1
     else
      fac=OutputDenominator*CP%tcmb**2
     end if
  end if

!FITS file has Cls without l(l+1)/twopi factors
  do i=2,lmx
     clout(i,:)=clout(i,:)/i/dble(i+1)*fac
  end do

  allocate(allcl(0:lmx,1:4))
  allcl(2:lmx,1:4) = clout
  allcl(0:1,1:4) = 0
  deallocate(clout)

  header = ''

  if (CP%OutputNormalization == outCOBE) then
     unitstr='Kelvin-squared'
  else
     unitstr='unknown'
  end if

  call add_card(header,'COMMENT','-----------------------------------------------')
  call add_card(header,'COMMENT','     CMB power spectrum C(l) keywords          ')
  call add_card(header,'COMMENT','-----------------------------------------------')
  call add_card(header,'EXTNAME','''COMPUTED POWER SPECTRUM''')
  call add_card(header,'COMMENT',' POWER SPECTRUM : C(l) ')
  call add_card(header)
  call add_card(header,'CREATOR','CAMB',        'Software creating the FITS file')
  call add_card(header,'VERSION',version,     'Version of the simulation software')
  call add_card(header,'POLAR',.true.,'Polarisation included (True/False)')
  call add_card(header,'POLNORM','CMBFAST','Uses E-B conventions')
  call add_card(header)
  call add_card(header)
  call add_card(header,'TTYPE1', 'TEMPERATURE','Temperature C(l)')
  call add_card(header,'TUNIT1', unitstr,'unit')
  call add_card(header)

     call add_card(header,'TTYPE2', 'E-mode C_l','ELECTRIC polarisation C(l)')
     call add_card(header,'TUNIT2', unitstr,'unit')
     call add_card(header)

     call add_card(header,'TTYPE3', 'B-mode C_l','MAGNETIC polarisation C(l)')
     call add_card(header,'TUNIT3', unitstr,'unit')
     call add_card(header)

     call add_card(header,'TTYPE4', 'E-T cross corr.','Gradient-Temperature cross terms')
     call add_card(header,'TUNIT4', unitstr,'unit')
     call add_card(header)

 call add_card(header,'COMMENT','-----------------------------------------------')
 call add_card(header,'COMMENT','     Cosmological parameters')
 call add_card(header,'COMMENT','-----------------------------------------------')
 call add_card(header,'OMEGAB',CP%omegab, 'Omega in baryons')
 call add_card(header,'OMEGAC',CP%omegac, 'Omega in CDM')
 call add_card(header,'OMEGAV',CP%omegav, 'Omega in cosmological constant')
 call add_card(header,'OMEGAN',CP%omegan, 'Omega in neutrinos')
 call add_card(header,'HUBBLE', CP%h0, 'Hublle constant in km/s/Mpc')
 call add_card(header,'NNUNR',CP%Num_Nu_massive, 'number of massive neutrinos')
 call add_card(header,'NNUR',CP%Num_Nu_massless, 'number of massless neutrinos')
 call add_card(header,'TCMB',CP%tcmb, 'CMB temperature in Kelvin')
 call add_card(header,'HELFRACT',CP%yhe, 'Helium fraction')
 call add_card(header,'OPTDLSS',CP%Reion%optical_depth, 'reionisation optical depth')
 call add_card(header,'IONFRACT',CP%Reion%fraction, 'ionisation fraction')
 call add_card(header,'ZREION',CP%reion%redshift, 'reionisation redshift')
 call add_card(header,'COMMENT','-----------------------------------------------')
 call add_card(header,'COMMENT','     Other parameters')
 call add_card(header,'COMMENT','-----------------------------------------------')
 call add_card(header,'SCALARS',CP%WantScalars, 'includes scalar modes')
 call add_card(header,'TENSORS',CP%WantTensors, 'includes tensor modes')
 call add_card(header,'INITFLAG',CP%Scalar_initial_condition, 'initial condition flag') 
 COBEnorm = CP%outputNormalization==outCOBE
 call add_card(header,'COBENORM',COBEnorm, 'COBE normalized') 
 call add_card(header,'KETA_MAX',CP%Max_eta_k, 'Max wavenumber') 
 call add_card(header,'PRECIS',AccuracyBoost, 'Relative computation accuracy') 
 call add_card(header,'EQS_FILE',Eqns_name, 'Gauge-dependent and background equations') 
 call add_card(header,'POW_FILE',Power_Name, 'Initial power spectrum file') 
 i = Power_Descript(1,CP%WantScalars,CP%WantTensors,PowerKeys,PowerVals)
 do j=1,i
 call add_card(header,PowerKeys(j),PowerVals(j), 'Initial power spectrum details') 
 end do
  
  nlheader = SIZE(header)
  call write_asctab (allcl, lmx, 4, header, nlheader, Clsfile)
  deallocate(allcl)

end subroutine WriteFitsCls
