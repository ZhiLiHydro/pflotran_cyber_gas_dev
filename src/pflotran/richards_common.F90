module Richards_Common_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Richards_Aux_module
  use Global_Aux_module
  
  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal
  
  implicit none
  
  private 

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-8
  PetscReal, parameter :: floweps   = 1.D-24
  PetscReal, parameter :: perturbation_tolerance = 1.d-6
  
  public RichardsAccumulation, &
         RichardsAccumDerivative, &
         RichardsFlux, &
         RichardsFluxDerivative, &
         RichardsBCFlux, &
         RichardsBCFluxDerivative
  
contains

! ************************************************************************** !

subroutine RichardsAccumDerivative(rich_auxvar,global_auxvar, &
                                   material_auxvar, &
                                   option, &
                                   characteristic_curves, &
                                   J)
  ! 
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 

  use Option_module
  use Characteristic_Curves_module
  use Material_Aux_class, only : material_auxvar_type, &
                                 soil_compressibility_index, &
                                 MaterialAuxVarInit, &
                                 MaterialAuxVarCopy, &
                                 MaterialAuxVarStrip, &
                                 MaterialCompressSoil
  
  implicit none

  type(richards_auxvar_type) :: rich_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  type(characteristic_curves_type) :: characteristic_curves
  PetscReal :: J(option%nflowdof,option%nflowdof)
     
  PetscInt :: ispec 
  PetscReal :: vol_over_dt
  PetscReal :: tempreal
  PetscReal :: derivative

  PetscInt :: iphase, ideriv
  type(richards_auxvar_type) :: rich_auxvar_pert
  type(global_auxvar_type) :: global_auxvar_pert
  ! leave as type
  type(material_auxvar_type) :: material_auxvar_pert
  PetscReal :: x(1), x_pert(1), pert, res(1), res_pert(1), J_pert(1,1)

  vol_over_dt = material_auxvar%volume/option%flow_dt
  
  ! accumulation term units = dkmol/dp
  J(1,1) = (material_auxvar%dporosity_dp*global_auxvar%sat(1)* &
            global_auxvar%den(1) + &
            (global_auxvar%sat(1)*rich_auxvar%dden_dp + &
             rich_auxvar%dsat_dp*global_auxvar%den(1)) * &
            material_auxvar%porosity) * &
            vol_over_dt

  if (option%flow%numerical_derivatives) then
    call GlobalAuxVarInit(global_auxvar_pert,option)  
    call MaterialAuxVarInit(material_auxvar_pert,option)  
    call RichardsAuxVarCopy(rich_auxvar,rich_auxvar_pert,option)
    call GlobalAuxVarCopy(global_auxvar,global_auxvar_pert,option)
    call MaterialAuxVarCopy(material_auxvar,material_auxvar_pert,option)
    x(1) = global_auxvar%pres(1)
    call RichardsAccumulation(rich_auxvar,global_auxvar,material_auxvar, &
                              option,res)
    ideriv = 1
    pert = max(dabs(x(ideriv)*perturbation_tolerance),0.1d0)
    x_pert = x
    if (x_pert(ideriv) < option%reference_pressure) pert = -1.d0*pert
    x_pert(ideriv) = x_pert(ideriv) + pert
    
    call RichardsAuxVarCompute(x_pert(1),rich_auxvar_pert,global_auxvar_pert, &
                               material_auxvar_pert, &
                               characteristic_curves, &
                               -999, &
                               PETSC_TRUE,option)
    call RichardsAccumulation(rich_auxvar_pert,global_auxvar_pert, &
                              material_auxvar_pert, &
                              option,res_pert)
    J_pert(1,1) = (res_pert(1)-res(1))/pert
    J = J_pert
    call GlobalAuxVarStrip(global_auxvar_pert)  
    call MaterialAuxVarStrip(material_auxvar_pert)  
  endif
   
end subroutine RichardsAccumDerivative

! ************************************************************************** !

subroutine RichardsAccumulation(rich_auxvar,global_auxvar, &
                                material_auxvar, &
                                option,Res)
  ! 
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 

  use Option_module
  use Material_Aux_class, only : material_auxvar_type, &
                                 soil_compressibility_index, &
                                 MaterialCompressSoil
  
  implicit none

  type(richards_auxvar_type) :: rich_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: Res(1:option%nflowdof)
  
  PetscReal :: vol_over_dt
       
  vol_over_dt = material_auxvar%volume/option%flow_dt

  ! accumulation term units = kmol/s
  Res(1) = global_auxvar%sat(1) * global_auxvar%den(1) * &
           material_auxvar%porosity * vol_over_dt
  
end subroutine RichardsAccumulation

! ************************************************************************** !

subroutine RichardsFluxDerivative(rich_auxvar_up,global_auxvar_up, &
                                  material_auxvar_up, & 
                                  rich_auxvar_dn,global_auxvar_dn, &
                                  material_auxvar_dn, &
                                  area, dist, &
                                  aniso_data, &
                                  option, &
                                  characteristic_curves_up, &
                                  characteristic_curves_dn, &
                                  Jup,Jdn)
  ! 
  ! Computes the derivatives of the internal flux terms
  ! for the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 
  use Option_module 
  use Characteristic_Curves_module
  use Material_Aux_class
  use Connection_module
  use Anisotropy_Data_module
  
  implicit none
  
  type(richards_auxvar_type) :: rich_auxvar_up, rich_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(aniso_data_type), pointer :: aniso_data
  type(option_type) :: option
  PetscReal :: v_darcy, v_darcy_aniso, area, dist(-1:3)
  type(characteristic_curves_type) :: characteristic_curves_up
  type(characteristic_curves_type) :: characteristic_curves_dn
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
     
  PetscReal :: q
  PetscReal :: ukvr,Dq
  PetscReal :: dd_up, dd_dn, perm_up, perm_dn
  PetscReal :: dist_gravity
  PetscReal :: upweight,density_ave,cond,gravity,dphi
  
  PetscReal :: dden_ave_dp_up, dden_ave_dp_dn
  PetscReal :: dgravity_dden_up, dgravity_dden_dn
  PetscReal :: dphi_dp_up, dphi_dp_dn
  PetscReal :: dukvr_dp_up, dukvr_dp_dn
!  PetscReal :: dukvr_x_dp_up, dukvr_x_dp_dn
!  PetscReal :: dukvr_y_dp_up, dukvr_y_dp_dn
!  PetscReal :: dukvr_z_dp_up, dukvr_z_dp_dn
  PetscReal :: dq_dp_up, dq_dp_dn
  PetscReal :: dv_darcy_dp_up, dv_darcy_dp_dn

  PetscReal :: dq_dp_up_adjXm, dq_dp_dn_adjXm, dq_dp_up_adjXp, dq_dp_dn_adjXp
  PetscReal :: dq_dp_up_adjYm, dq_dp_dn_adjYm, dq_dp_up_adjYp, dq_dp_dn_adjYp
  PetscReal :: dq_dp_up_adjZm, dq_dp_dn_adjZm, dq_dp_up_adjZp, dq_dp_dn_adjZp

  PetscReal, pointer :: dq_dp_up_adj(:)
  PetscReal, pointer :: dq_dp_dn_adj(:)
  PetscInt :: iadj

  PetscInt :: iphase, ideriv
  type(richards_auxvar_type) :: rich_auxvar_pert_up, rich_auxvar_pert_dn
  type(global_auxvar_type) :: global_auxvar_pert_up, global_auxvar_pert_dn
  ! leave as type
  type(material_auxvar_type) :: material_auxvar_pert_up, material_auxvar_pert_dn
  PetscReal :: x_up(1), x_dn(1), x_pert_up(1), x_pert_dn(1), pert_up, pert_dn, &
            res(1), res_pert_up(1), res_pert_dn(1), J_pert_up(1,1), J_pert_dn(1,1)
  
  v_darcy = 0.d0  
  v_darcy_aniso = 0.d0
  ukvr = 0.d0
  
  Jup = 0.d0
  Jdn = 0.d0 
  
  dden_ave_dp_up = 0.d0
  dden_ave_dp_dn = 0.d0
  dgravity_dden_up = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_up = 0.d0
  dphi_dp_dn = 0.d0
  dukvr_dp_up = 0.d0
  dukvr_dp_dn = 0.d0
  dq_dp_up = 0.d0
  dq_dp_dn = 0.d0
  
  call ConnectionCalculateDistances(dist,option%gravity,dd_up,dd_dn, &
                                    dist_gravity,upweight)
  if(.not. associated(aniso_data)) then
    call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
    call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  else
    ! Flux component due head gradients directly across the interface can be
    ! calculated using the isotropic approach provided that permeabilty is set
    ! consistent with anisotropic calculation (avoid collapsing to scalar).
    ! Additional anisotropic terms due to tangential head gradients are added
    ! by later call to RichardsFluxAniso()
    if(associated(aniso_data%struct)) then
      ! sjb - can do same calls as isotropic case here, should just extract
      !       k_xx, k_yy or k_zz depending on orientation
      ! sjb - would be better to be specific here rather than rely on
      !       TensorToScalar call to do the right thing
      call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
      call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
    else
      ! u-direction is chosen to be normal to the interface, so
      perm_up = aniso_data%unstruct%up%perm_uu
      perm_dn = aniso_data%unstruct%dn%perm_uu
    endif
  endif  

  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
  
! Flow term
  if (rich_auxvar_up%kvr > eps .or. &
      rich_auxvar_dn%kvr > eps) then
   
    if (global_auxvar_up%sat(1) <eps) then 
      upweight=0.d0
    else if (global_auxvar_dn%sat(1) <eps) then 
      upweight=1.d0
    endif
    density_ave = upweight*global_auxvar_up%den(1)+ &
                  (1.D0-upweight)*global_auxvar_dn%den(1)
    dden_ave_dp_up = upweight*rich_auxvar_up%dden_dp
    dden_ave_dp_dn = (1.D0-upweight)*rich_auxvar_dn%dden_dp

    gravity = (upweight*global_auxvar_up%den(1) + &
              (1.D0-upweight)*global_auxvar_dn%den(1)) &
              * FMWH2O * dist_gravity
    dgravity_dden_up = upweight*FMWH2O*dist_gravity
    dgravity_dden_dn = (1.d0-upweight)*FMWH2O*dist_gravity

    dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1)  + gravity
    dphi_dp_up = 1.d0 + dgravity_dden_up*rich_auxvar_up%dden_dp
    dphi_dp_dn = -1.d0 + dgravity_dden_dn*rich_auxvar_dn%dden_dp

    if (dphi>=0.D0) then
      ukvr = rich_auxvar_up%kvr
      dukvr_dp_up = rich_auxvar_up%dkvr_dp
    else
      ukvr = rich_auxvar_dn%kvr
      dukvr_dp_dn = rich_auxvar_dn%dkvr_dp
    endif      
   
    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi

      dv_darcy_dp_up = Dq*(dukvr_dp_up*dphi+ukvr*dphi_dp_up)
      dv_darcy_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)


      ! calculate anisotropic contributions
      if(associated(aniso_data)) then
        ! Call RichardsFluxAniso with calc_derivs=PETSC_TRUE
        Call RichardsFluxAniso(aniso_data,material_auxvar_up,material_auxvar_dn,&
                                global_auxvar_up,global_auxvar_dn,&
                                option,dist,&
                                ukvr,dukvr_dp_up,dukvr_dp_dn,&
                                v_darcy, dv_darcy_dp_up, dv_darcy_dp_dn,&
                                PETSC_TRUE,v_darcy_aniso)
        v_darcy = v_darcy + v_darcy_aniso
        ! dv_darcy_aniso_dp_... terms are stored in aniso_data - used below
      endif

      q = v_darcy * area

      dq_dp_up = dv_darcy_dp_up*area
      dq_dp_dn = dv_darcy_dp_dn*area
      
      if(associated(aniso_data)) then
        ! (note ukvr and dukvr_dp_up/dn already taken account of in RichardsFluxAniso())
        if(associated(aniso_data%struct)) then
          ! structured grid case
          ! add aniso terms to dq_dp_up/dn
          dq_dp_up = dq_dp_up + aniso_data%struct%dv_darcy_aniso_dp_up * area
          dq_dp_dn = dq_dp_dn + aniso_data%struct%dv_darcy_aniso_dp_dn * area
          ! additional aniso terms for dv_darcy_aniso_dp_up/dn_adjX/Y/Z/pm
          dq_dp_up_adjXm = aniso_data%struct%dv_darcy_aniso_dp_up_adjXm * area
          dq_dp_dn_adjXm = aniso_data%struct%dv_darcy_aniso_dp_dn_adjXm * area
          dq_dp_up_adjXp = aniso_data%struct%dv_darcy_aniso_dp_up_adjXp * area
          dq_dp_dn_adjXp = aniso_data%struct%dv_darcy_aniso_dp_dn_adjXp * area
          dq_dp_up_adjYm = aniso_data%struct%dv_darcy_aniso_dp_up_adjYm * area
          dq_dp_dn_adjYm = aniso_data%struct%dv_darcy_aniso_dp_dn_adjYm * area
          dq_dp_up_adjYp = aniso_data%struct%dv_darcy_aniso_dp_up_adjYp * area
          dq_dp_dn_adjYp = aniso_data%struct%dv_darcy_aniso_dp_dn_adjYp * area
          dq_dp_up_adjZm = aniso_data%struct%dv_darcy_aniso_dp_up_adjZm * area
          dq_dp_dn_adjZm = aniso_data%struct%dv_darcy_aniso_dp_dn_adjZm * area
          dq_dp_up_adjZp = aniso_data%struct%dv_darcy_aniso_dp_up_adjZp * area
          dq_dp_dn_adjZp = aniso_data%struct%dv_darcy_aniso_dp_dn_adjZp * area
        endif
        if(associated(aniso_data%unstruct)) then
          ! unstructured grid case
          ! add aniso terms to dq_dp_up/dn

          allocate(dq_dp_up_adj(aniso_data%unstruct%up%num_adj))
          allocate(dq_dp_dn_adj(aniso_data%unstruct%dn%num_adj))

          dq_dp_up = dq_dp_up + aniso_data%unstruct%dv_darcy_aniso_dp_up * area
          dq_dp_dn = dq_dp_dn + aniso_data%unstruct%dv_darcy_aniso_dp_dn * area
          ! additional aniso terms for dv_darcy_aniso_dp_up/dn_adj
          do iadj=1,aniso_data%unstruct%up%num_adj
            dq_dp_up_adj(iadj) = aniso_data%unstruct%dv_darcy_aniso_dp_up_adj(iadj) * area
          enddo
          do iadj=1,aniso_data%unstruct%dn%num_adj
            dq_dp_dn_adj(iadj) = aniso_data%unstruct%dv_darcy_aniso_dp_dn_adj(iadj) * area
          enddo
        endif
      endif
      
      ! flux = q*density_ave, hence density_ave terms below
      Jup(1,1) = (dq_dp_up*density_ave+q*dden_ave_dp_up)
      Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)

      if(associated(aniso_data)) then
        ! Note dden_ave_dp_up_adj etc. assumed = 0 so no dden_... terms below
        ! (The quantities below are used in RichardsJacobianInternalConn when constructing J)
        if(associated(aniso_data%struct)) then
          aniso_data%struct%Jup_adjXm = dq_dp_up_adjXm * density_ave
          aniso_data%struct%Jdn_adjXm = dq_dp_dn_adjXm * density_ave
          aniso_data%struct%Jup_adjXp = dq_dp_up_adjXp * density_ave
          aniso_data%struct%Jdn_adjXp = dq_dp_dn_adjXp * density_ave
          aniso_data%struct%Jup_adjYm = dq_dp_up_adjYm * density_ave
          aniso_data%struct%Jdn_adjYm = dq_dp_dn_adjYm * density_ave
          aniso_data%struct%Jup_adjYp = dq_dp_up_adjYp * density_ave
          aniso_data%struct%Jdn_adjYp = dq_dp_dn_adjYp * density_ave
          aniso_data%struct%Jup_adjZm = dq_dp_up_adjZm * density_ave
          aniso_data%struct%Jdn_adjZm = dq_dp_dn_adjZm * density_ave
          aniso_data%struct%Jup_adjZp = dq_dp_up_adjZp * density_ave
          aniso_data%struct%Jdn_adjZp = dq_dp_dn_adjZp * density_ave
        endif
        if(associated(aniso_data%unstruct)) then
          do iadj=1,aniso_data%unstruct%up%num_adj
            aniso_data%unstruct%Jup_adj(iadj) = dq_dp_up_adj(iadj) * density_ave
          enddo
          do iadj=1,aniso_data%unstruct%dn%num_adj
            aniso_data%unstruct%Jdn_adj(iadj) = dq_dp_dn_adj(iadj) * density_ave
          enddo
        endif
      endif
    endif
  endif

 ! note: Res is the flux contribution, for node up J = J + Jup
 !                                              dn J = J - Jdn  

  if (option%flow%numerical_derivatives) then
    call GlobalAuxVarInit(global_auxvar_pert_up,option)
    call GlobalAuxVarInit(global_auxvar_pert_dn,option)  
    call MaterialAuxVarInit(material_auxvar_pert_up,option)
    call MaterialAuxVarInit(material_auxvar_pert_dn,option)  
    call RichardsAuxVarCopy(rich_auxvar_up,rich_auxvar_pert_up,option)
    call RichardsAuxVarCopy(rich_auxvar_dn,rich_auxvar_pert_dn,option)
    call GlobalAuxVarCopy(global_auxvar_up,global_auxvar_pert_up,option)
    call GlobalAuxVarCopy(global_auxvar_dn,global_auxvar_pert_dn,option)
    call MaterialAuxVarCopy(material_auxvar_up,material_auxvar_pert_up,option)
    call MaterialAuxVarCopy(material_auxvar_dn,material_auxvar_pert_dn,option)
    x_up(1) = global_auxvar_up%pres(1)
    x_dn(1) = global_auxvar_dn%pres(1)
    call RichardsFlux(rich_auxvar_up,global_auxvar_up,material_auxvar_up, &
                      rich_auxvar_dn,global_auxvar_dn,material_auxvar_dn, &
                      area, dist, &
                      aniso_data, &
                      option,v_darcy,res)
    ideriv = 1
!    pert_up = x_up(ideriv)*perturbation_tolerance
    pert_up = max(dabs(x_up(ideriv)*perturbation_tolerance),0.1d0)
    if (x_up(ideriv) < option%reference_pressure) pert_up = -1.d0*pert_up
!    pert_dn = x_dn(ideriv)*perturbation_tolerance
    pert_dn = max(dabs(x_dn(ideriv)*perturbation_tolerance),0.1d0)
    if (x_dn(ideriv) < option%reference_pressure) pert_dn = -1.d0*pert_dn
    x_pert_up = x_up
    x_pert_dn = x_dn
    x_pert_up(ideriv) = x_pert_up(ideriv) + pert_up
    x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
    call RichardsAuxVarCompute(x_pert_up(1),rich_auxvar_pert_up, &
                               global_auxvar_pert_up, &
                               material_auxvar_pert_up, &
                               characteristic_curves_up, &
                               -999, &
                               PETSC_TRUE,option)
    call RichardsAuxVarCompute(x_pert_dn(1),rich_auxvar_pert_dn, &
                               global_auxvar_pert_dn, &
                               material_auxvar_pert_dn, &
                               characteristic_curves_dn, &
                               -999, &
                               PETSC_TRUE,option)
    call RichardsFlux(rich_auxvar_pert_up,global_auxvar_pert_up, &
                      material_auxvar_pert_up, &
                      rich_auxvar_dn,global_auxvar_dn, &
                      material_auxvar_dn, &
                      area, dist, &
                      aniso_data, &
                      option,v_darcy,res_pert_up)
    call RichardsFlux(rich_auxvar_up,global_auxvar_up, &
                      material_auxvar_up, &
                      rich_auxvar_pert_dn,global_auxvar_pert_dn, &
                      material_auxvar_pert_dn, &
                      area, dist, &
                      aniso_data, &
                      option,v_darcy,res_pert_dn)
    J_pert_up(1,ideriv) = (res_pert_up(1)-res(1))/pert_up
    J_pert_dn(1,ideriv) = (res_pert_dn(1)-res(1))/pert_dn
    Jup = J_pert_up
    Jdn = J_pert_dn
    call GlobalAuxVarStrip(global_auxvar_pert_up)
    call GlobalAuxVarStrip(global_auxvar_pert_dn)    
    call MaterialAuxVarStrip(material_auxvar_pert_up)
    call MaterialAuxVarStrip(material_auxvar_pert_dn)    
  endif

end subroutine RichardsFluxDerivative

! ************************************************************************** !

subroutine RichardsFluxAniso(aniso_data, &
                            material_auxvar_up,material_auxvar_dn, &
                            global_auxvar_up,global_auxvar_dn, &
                            option, &
                            dist, & 
                            ukvr, dukvr_dp_up, dukvr_dp_dn, & ! kvr = krel/visc (derivs only needed if calculating derivatives, o/w set zero)
                            v_darcy_iso, dv_darcy_iso_dp_up, dv_darcy_iso_dp_dn, & ! needed in unstruct case for dphi/du (derivs only needed if calculating derivatives, o/w set zero)
                            calc_derivs, &
                            v_darcy_aniso)
!
! Computes the contributions to internal flux terms and Jacobian arising
! from pressure gradients perpendicular to cell-centre connections due
! to anisotropic permeability terms (structured grids only)
!
! Author: Steve Benbow
! Date: 01/06/20
!
  use Option_module
  use Material_Aux_class
  use Anisotropy_Data_module
  use Variables_module
  use Connection_module


  ! inputs
  type(aniso_data_type), pointer :: aniso_data
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(option_type) :: option
  PetscReal :: dist(-1:3)
  PetscReal :: ukvr

  PetscReal :: dukvr_dp_up, dukvr_dp_dn
  PetscBool :: calc_derivs ! terms are stored in aniso_data%dv_darcy_aniso_...

  PetscReal :: v_darcy_iso ! v_darcy = Dq * ukvr * dphi, as calculated by RichardsFlux()
  PetscReal :: dv_darcy_iso_dp_up,dv_darcy_iso_dp_dn ! derivs of above

  ! outputs
  PetscReal :: v_darcy_aniso

  ! local
  PetscReal :: mu_up, mu_dn, rho_up, rho_dn, dummy
  PetscReal :: v_darcy_aniso_up, v_darcy_aniso_dn
  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: gravity

  ! local - structured grid case
  type(aniso_data_structured_type), pointer :: aniso_ds
  type(aniso_deriv_data_structured_type), pointer :: aniso_dds
  PetscReal :: k_up_xx, k_up_yy, k_up_zz, k_up_xy, k_up_yz, k_up_xz
  PetscReal :: k_dn_xx, k_dn_yy, k_dn_zz, k_dn_xy, k_dn_yz, k_dn_xz
  PetscReal :: dphi_x_up, dphi_x_dn, dphi_y_up, dphi_y_dn, dphi_z_up, dphi_z_dn

  ! local - unstructured grid case
  type(aniso_data_unstructured_type), pointer :: aniso_du
  type(aniso_deriv_data_unstructured_type), pointer :: aniso_ddu
  type(aniso_data_unstructured_cell_type), pointer :: aniso_updn
  PetscReal :: pres_updn, rho_updn
  PetscReal, target :: dphi_v_up, dphi_v_dn, dphi_w_up, dphi_w_dn
  PetscReal, pointer :: dphi, ddphi_dp, ddphi_dp_adj
  PetscInt :: iupdn, iadj
  PetscReal :: LSrhs(2), dot, res, dres_dp_up, dres_dp_dn, dres_dp_adj
  PetscReal :: dphi_du, d2phi_du_dp_up, d2phi_du_dp_dn

  nullify(aniso_ds)
  nullify(aniso_dds)

  nullify(aniso_du)
  nullify(aniso_ddu)

  nullify(dphi)
  nullify(ddphi_dp)
  nullify(ddphi_dp_adj)

  ! Initialise output
  v_darcy_aniso = 0.d0

  ! Ideally would have these for adjacent cells too
  rho_up = global_auxvar_up%den(1) * FMWH2O
  rho_dn = global_auxvar_dn%den(1) * FMWH2O


  if(associated(aniso_data%struct)) then
    ! structured grid case
    aniso_ds => aniso_data%struct

    ! Permeability tensor to up/dn cells
    k_up_xx = MaterialAuxVarGetValue(material_auxvar_up,PERMEABILITY_X)
    k_up_yy = MaterialAuxVarGetValue(material_auxvar_up,PERMEABILITY_Y)
    k_up_zz = MaterialAuxVarGetValue(material_auxvar_up,PERMEABILITY_Z)
    k_up_xy = MaterialAuxVarGetValue(material_auxvar_up,PERMEABILITY_XY)
    k_up_yz = MaterialAuxVarGetValue(material_auxvar_up,PERMEABILITY_YZ)
    k_up_xz = MaterialAuxVarGetValue(material_auxvar_up,PERMEABILITY_XZ)

    k_dn_xx = MaterialAuxVarGetValue(material_auxvar_dn,PERMEABILITY_X)
    k_dn_yy = MaterialAuxVarGetValue(material_auxvar_dn,PERMEABILITY_Y)
    k_dn_zz = MaterialAuxVarGetValue(material_auxvar_dn,PERMEABILITY_Z)
    k_dn_xy = MaterialAuxVarGetValue(material_auxvar_dn,PERMEABILITY_XY)
    k_dn_yz = MaterialAuxVarGetValue(material_auxvar_dn,PERMEABILITY_YZ)
    k_dn_xz = MaterialAuxVarGetValue(material_auxvar_dn,PERMEABILITY_XZ)

    ! allow off-diagonal terms to be uninitialised, in which case set to zero
    if ( k_up_xy == UNINITIALIZED_DOUBLE ) then
      k_up_xy = 0.d0
    endif
    if ( k_up_yz == UNINITIALIZED_DOUBLE ) then
      k_up_yz = 0.d0
    endif
    if ( k_up_xz == UNINITIALIZED_DOUBLE ) then
      k_up_xz = 0.d0
    endif
    if ( k_dn_xy == UNINITIALIZED_DOUBLE ) then
      k_dn_xy = 0.d0
    endif
    if ( k_dn_yz == UNINITIALIZED_DOUBLE ) then
      k_dn_yz = 0.d0
    endif
    if ( k_dn_xz == UNINITIALIZED_DOUBLE ) then
      k_dn_xz = 0.d0
    endif

    ! components of grad P on connection, in the plane of the interface
    dphi_x_up = 0.d0
    dphi_x_dn = 0.d0
    dphi_y_up = 0.d0
    dphi_y_dn = 0.d0
    dphi_z_up = 0.d0
    dphi_z_dn = 0.d0

    if (calc_derivs .EQV. PETSC_TRUE) then
      allocate(aniso_dds)

      aniso_dds%ddphi_x_dp_up_adjXm = 0.d0
      aniso_dds%ddphi_x_dp_dn_adjXm = 0.d0
      aniso_dds%ddphi_x_dp_up_adjXp = 0.d0
      aniso_dds%ddphi_x_dp_dn_adjXp = 0.d0
      aniso_dds%ddphi_x_dp_up       = 0.d0
      aniso_dds%ddphi_x_dp_dn       = 0.d0

      aniso_dds%ddphi_y_dp_up_adjYm = 0.d0
      aniso_dds%ddphi_y_dp_dn_adjYm = 0.d0
      aniso_dds%ddphi_y_dp_up_adjYp = 0.d0
      aniso_dds%ddphi_y_dp_dn_adjYp = 0.d0
      aniso_dds%ddphi_y_dp_up       = 0.d0
      aniso_dds%ddphi_y_dp_dn       = 0.d0

      aniso_dds%ddphi_z_dp_up_adjZm = 0.d0
      aniso_dds%ddphi_z_dp_dn_adjZm = 0.d0
      aniso_dds%ddphi_z_dp_up_adjZp = 0.d0
      aniso_dds%ddphi_z_dp_dn_adjZp = 0.d0
      aniso_dds%ddphi_z_dp_up       = 0.d0
      aniso_dds%ddphi_z_dp_dn       = 0.d0
    endif

    ! P gradient in x-dir - upstream
    if( aniso_ds%id_up_adjXp .NE. -1 .AND. aniso_ds%id_up_adjXm .NE. -1 ) then ! has contribution in x-dir with adj p/m cells
      dphi_x_up = ((aniso_ds%pres_up_adjXp - rho_up * aniso_ds%loc_gravity_up_adjXp) &
                      - (aniso_ds%pres_up_adjXm - rho_up * aniso_ds%loc_gravity_up_adjXm)) &
                      / (aniso_ds%dist_up_adjXp + aniso_ds%dist_up_adjXm)
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%ddphi_x_dp_up_adjXm = -1.d0 / (aniso_ds%dist_up_adjXp + aniso_ds%dist_up_adjXm)
        aniso_dds%ddphi_x_dp_up_adjXp =  1.d0 / (aniso_ds%dist_up_adjXp + aniso_ds%dist_up_adjXm)
        aniso_dds%ddphi_x_dp_up       =  0.d0
      endif
    else if ( aniso_ds%id_up_adjXp .NE. -1 ) then ! has contribution in x-dir with adj p cell
      dphi_x_up = ((aniso_ds%pres_up_adjXp - rho_up * aniso_ds%loc_gravity_up_adjXp) &
                    - (global_auxvar_up%pres(1) - rho_up * aniso_ds%loc_gravity_up)) &
                    / aniso_ds%dist_up_adjXp
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%ddphi_x_dp_up_adjXm =  0.d0
        aniso_dds%ddphi_x_dp_up_adjXp =  1.d0 / aniso_ds%dist_up_adjXp
        aniso_dds%ddphi_x_dp_up       = -1.d0 / aniso_ds%dist_up_adjXp
      endif
    else if ( aniso_ds%id_up_adjXm .NE. -1 ) then ! has contribution in x-dir with adj m cell
      dphi_x_up = ((global_auxvar_up%pres(1) - rho_up * aniso_ds%loc_gravity_up) &
                    - (aniso_ds%pres_up_adjXm - rho_up * aniso_ds%loc_gravity_up_adjXm)) &
                    / aniso_ds%dist_up_adjXm
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%ddphi_x_dp_up_adjXm = -1.d0 / aniso_ds%dist_up_adjXm
        aniso_dds%ddphi_x_dp_up_adjXp =  0.d0
        aniso_dds%ddphi_x_dp_up       =  1.d0 / aniso_ds%dist_up_adjXm
      endif
    endif

    ! P gradient in x-dir - dnstream
    if( aniso_ds%id_dn_adjXp .NE. -1 .AND. aniso_ds%id_dn_adjXm .NE. -1 ) then ! has contribution in x-dir with adj p/m cells
      dphi_x_dn = ((aniso_ds%pres_dn_adjXp - rho_dn * aniso_ds%loc_gravity_dn_adjXp) &
                    - (aniso_ds%pres_dn_adjXm - rho_dn * aniso_ds%loc_gravity_dn_adjXm)) &
                    / (aniso_ds%dist_dn_adjXp + aniso_ds%dist_dn_adjXm)
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%ddphi_x_dp_dn_adjXm = -1.d0 / (aniso_ds%dist_dn_adjXp + aniso_ds%dist_dn_adjXm)
        aniso_dds%ddphi_x_dp_dn_adjXp =  1.d0 / (aniso_ds%dist_dn_adjXp + aniso_ds%dist_dn_adjXm)
        aniso_dds%ddphi_x_dp_dn       =  0.d0
      endif
    else if ( aniso_ds%id_dn_adjXp .NE. -1 ) then ! has contribution in x-dir with adj p cell
      dphi_x_dn = ((aniso_ds%pres_dn_adjXp - rho_dn * aniso_ds%loc_gravity_dn_adjXp) &
                  - (global_auxvar_dn%pres(1) - rho_dn * aniso_ds%loc_gravity_dn)) &
                  / aniso_ds%dist_dn_adjXp
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%ddphi_x_dp_dn_adjXm =  0.d0
        aniso_dds%ddphi_x_dp_dn_adjXp =  1.d0 / aniso_ds%dist_dn_adjXp
        aniso_dds%ddphi_x_dp_dn       = -1.d0 / aniso_ds%dist_dn_adjXp
      endif
    else if ( aniso_ds%id_dn_adjXm .NE. -1 ) then ! has contribution in x-dir with adj m cell
      dphi_x_dn = ((global_auxvar_dn%pres(1) - rho_dn * aniso_ds%loc_gravity_dn) &
                  - (aniso_ds%pres_dn_adjXm - rho_dn * aniso_ds%loc_gravity_dn_adjXm)) &
                  / aniso_ds%dist_dn_adjXm
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%ddphi_x_dp_dn_adjXm = -1.d0 / aniso_ds%dist_dn_adjXm
        aniso_dds%ddphi_x_dp_dn_adjXp =  0.d0
        aniso_dds%ddphi_x_dp_dn       =  1.d0 / aniso_ds%dist_dn_adjXm
      endif
    endif

    ! P gradient in y-dir - upstream
    if( aniso_ds%id_up_adjYp .NE. -1 .AND. aniso_ds%id_up_adjYm .NE. -1 ) then ! has contribution in y-dir with adj p/m cells
      dphi_y_up = ((aniso_ds%pres_up_adjYp - rho_up * aniso_ds%loc_gravity_up_adjYp) &
                    - (aniso_ds%pres_up_adjYm - rho_up * aniso_ds%loc_gravity_up_adjYm)) &
                    / (aniso_ds%dist_up_adjYp + aniso_ds%dist_up_adjYm)
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%ddphi_y_dp_up_adjYm = -1.d0 / (aniso_ds%dist_up_adjYp + aniso_ds%dist_up_adjYm)
        aniso_dds%ddphi_y_dp_up_adjYp =  1.d0 / (aniso_ds%dist_up_adjYp + aniso_ds%dist_up_adjYm)
        aniso_dds%ddphi_y_dp_up       =  0.d0
      endif
    else if ( aniso_ds%id_up_adjYp .NE. -1 ) then ! has contribution in y-dir with adj p cell
      dphi_y_up = ((aniso_ds%pres_up_adjYp - rho_up * aniso_ds%loc_gravity_up_adjYp) &
                  - (global_auxvar_up%pres(1) - rho_up * aniso_ds%loc_gravity_up)) &
                  / aniso_ds%dist_up_adjYp
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%ddphi_y_dp_up_adjYm =  0.d0
        aniso_dds%ddphi_y_dp_up_adjYp =  1.d0 / aniso_ds%dist_up_adjYp
        aniso_dds%ddphi_y_dp_up       = -1.d0 / aniso_ds%dist_up_adjYp
      endif
    else if ( aniso_ds%id_up_adjYm .NE. -1 ) then ! has contribution in y-dir with adj m cell
      dphi_y_up = ((global_auxvar_up%pres(1) - rho_up * aniso_ds%loc_gravity_up) &
                  - (aniso_ds%pres_up_adjYm - rho_up * aniso_ds%loc_gravity_up_adjYm)) &
                  / aniso_ds%dist_up_adjYm
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%ddphi_y_dp_up_adjYm = -1.d0 / aniso_ds%dist_up_adjYm
        aniso_dds%ddphi_y_dp_up_adjYp =  0.d0
        aniso_dds%ddphi_y_dp_up       =  1.d0 / aniso_ds%dist_up_adjYm
      endif
    endif

    ! P gradient in y-dir - dnstream
    if( aniso_ds%id_dn_adjYp .NE. -1 .AND. aniso_ds%id_dn_adjYm .NE. -1 ) then ! has contribution in y-dir with adj p/m cells
      dphi_y_dn = ((aniso_ds%pres_dn_adjYp - rho_dn * aniso_ds%loc_gravity_dn_adjYp) &
                    - (aniso_ds%pres_dn_adjYm - rho_dn * aniso_ds%loc_gravity_dn_adjYm)) &
                    / (aniso_ds%dist_dn_adjYp + aniso_ds%dist_dn_adjYm)
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%ddphi_y_dp_dn_adjYm = -1.d0 / (aniso_ds%dist_dn_adjYp + aniso_ds%dist_dn_adjYm)
        aniso_dds%ddphi_y_dp_dn_adjYp =  1.d0 / (aniso_ds%dist_dn_adjYp + aniso_ds%dist_dn_adjYm)
        aniso_dds%ddphi_y_dp_dn       =  0.d0
      endif
    else if ( aniso_ds%id_dn_adjYp .NE. -1 ) then ! has contribution in y-dir with adj p cell
      dphi_y_dn = ((aniso_ds%pres_dn_adjYp - rho_dn * aniso_ds%loc_gravity_dn_adjYp) &
                - (global_auxvar_dn%pres(1) - rho_dn * aniso_ds%loc_gravity_dn)) &
                / aniso_ds%dist_dn_adjYp
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%ddphi_y_dp_dn_adjYm =  0.d0
        aniso_dds%ddphi_y_dp_dn_adjYp =  1.d0 / aniso_ds%dist_dn_adjYp
        aniso_dds%ddphi_y_dp_dn       = -1.d0 / aniso_ds%dist_dn_adjYp
      endif
    else if ( aniso_ds%id_dn_adjYm .NE. -1 ) then ! has contribution in y-dir with adj m cell
      dphi_y_dn = ((global_auxvar_dn%pres(1) - rho_dn * aniso_ds%loc_gravity_dn) &
                - (aniso_ds%pres_dn_adjYm - rho_dn * aniso_ds%loc_gravity_dn_adjYm)) &
                / aniso_ds%dist_dn_adjYm
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%ddphi_y_dp_dn_adjYm = -1.d0 / aniso_ds%dist_dn_adjYm
        aniso_dds%ddphi_y_dp_dn_adjYp =  0.d0
        aniso_dds%ddphi_y_dp_dn       =  1.d0 / aniso_ds%dist_dn_adjYm
      endif
    endif

    ! P gradient in z-dir - upstream
    if( aniso_ds%id_up_adjZp .NE. -1 .AND. aniso_ds%id_up_adjZm .NE. -1 ) then ! has contribution in z-dir with adj p/m cells
      dphi_z_up = ((aniso_ds%pres_up_adjZp - rho_up * aniso_ds%loc_gravity_up_adjZp) &
                  - (aniso_ds%pres_up_adjZm - rho_up * aniso_ds%loc_gravity_up_adjZm)) &
                  / (aniso_ds%dist_up_adjZp + aniso_ds%dist_up_adjZm)
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%ddphi_z_dp_up_adjZm = -1.d0 / (aniso_ds%dist_up_adjZp + aniso_ds%dist_up_adjZm)
        aniso_dds%ddphi_z_dp_up_adjZp =  1.d0 / (aniso_ds%dist_up_adjZp + aniso_ds%dist_up_adjZm)
        aniso_dds%ddphi_z_dp_up       =  0.d0
      endif
    else if ( aniso_ds%id_up_adjZp .NE. -1 ) then ! has contribution in z-dir with adj p cell
      dphi_z_up = ((aniso_ds%pres_up_adjZp - rho_up * aniso_ds%loc_gravity_up_adjZp) &
                - (global_auxvar_up%pres(1) - rho_up * aniso_ds%loc_gravity_up)) &
                / aniso_ds%dist_up_adjZp
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%ddphi_z_dp_up_adjZm =  0.d0
        aniso_dds%ddphi_z_dp_up_adjZp =  1.d0 / aniso_ds%dist_up_adjZp
        aniso_dds%ddphi_z_dp_up       = -1.d0 / aniso_ds%dist_up_adjZp
      endif
    else if ( aniso_ds%id_up_adjZm .NE. -1 ) then ! has contribution in z-dir with adj m cell
      dphi_z_up = ((global_auxvar_up%pres(1) - rho_up * aniso_ds%loc_gravity_up) &
                - (aniso_ds%pres_up_adjZm - rho_up * aniso_ds%loc_gravity_up_adjZm)) &
                / aniso_ds%dist_up_adjZm
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%ddphi_z_dp_up_adjZm = -1.d0 / aniso_ds%dist_up_adjZm
        aniso_dds%ddphi_z_dp_up_adjZp =  0.d0
        aniso_dds%ddphi_z_dp_up       =  1.d0 / aniso_ds%dist_up_adjZm
      endif
    endif

    ! P gradient in z-dir - dnstream
    if( aniso_ds%id_dn_adjZp .NE. -1 .AND. aniso_ds%id_dn_adjZm .NE. -1 ) then ! has contribution in z-dir with adj p/m cells
      dphi_z_dn = ((aniso_ds%pres_dn_adjZp - rho_dn * aniso_ds%loc_gravity_dn_adjZp) &
                  - (aniso_ds%pres_dn_adjZm - rho_dn * aniso_ds%loc_gravity_dn_adjZm)) &
                  / (aniso_ds%dist_dn_adjZp + aniso_ds%dist_dn_adjZm)
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%ddphi_z_dp_dn_adjZm = -1.d0 / (aniso_ds%dist_dn_adjZp + aniso_ds%dist_dn_adjZm)
        aniso_dds%ddphi_z_dp_dn_adjZp =  1.d0 / (aniso_ds%dist_dn_adjZp + aniso_ds%dist_dn_adjZm)
        aniso_dds%ddphi_z_dp_dn       =  0.d0
      endif
    else if ( aniso_ds%id_dn_adjZp .NE. -1 ) then ! has contribution in z-dir with adj p cell
      dphi_z_dn = ((aniso_ds%pres_dn_adjZp - rho_dn * aniso_ds%loc_gravity_dn_adjZp) &
              - (global_auxvar_dn%pres(1) - rho_dn * aniso_ds%loc_gravity_dn)) &
              / aniso_ds%dist_dn_adjZp
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%ddphi_z_dp_dn_adjZm =  0.d0
        aniso_dds%ddphi_z_dp_dn_adjZp =  1.d0 / aniso_ds%dist_dn_adjZp
        aniso_dds%ddphi_z_dp_dn       = -1.d0 / aniso_ds%dist_dn_adjZp
      endif
    else if ( aniso_ds%id_dn_adjZm .NE. -1 ) then ! has contribution in z-dir with adj m cell
      dphi_z_dn = ((global_auxvar_dn%pres(1) - rho_dn * aniso_ds%loc_gravity_dn) &
              - (aniso_ds%pres_dn_adjZm - rho_dn * aniso_ds%loc_gravity_dn_adjZm)) &
              / aniso_ds%dist_dn_adjZm
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%ddphi_z_dp_dn_adjZm = -1.d0 / aniso_ds%dist_dn_adjZm
        aniso_dds%ddphi_z_dp_dn_adjZp =  0.d0
        aniso_dds%ddphi_z_dp_dn       =  1.d0 / aniso_ds%dist_dn_adjZm
      endif
    endif

    ! dbg
    ! write(*,"('dphi_x/y/z_up/dn =',e10.3,e10.3,e10.3,e10.3,e10.3,e10.3)") &
    !   dphi_x_up, dphi_x_dn, dphi_y_up, dphi_y_dn, dphi_z_up, dphi_z_dn

    ! Calculate anisotropic darcy contribution for given orientation
    if ( aniso_ds%orientation .eq. 'X' ) then
      v_darcy_aniso_up = -1.d0 * (k_up_xy * dphi_y_up + k_up_xz * dphi_z_up) * ukvr
      v_darcy_aniso_dn = -1.d0 * (k_dn_xy * dphi_y_dn + k_dn_xz * dphi_z_dn) * ukvr
      mu_up =       dist(-1) *dist(0)/k_up_xx  ! weight - up (dist part = dist to interface)
      mu_dn = (1.d0-dist(-1))*dist(0)/k_dn_xx  ! weight - dn (dist part = dist to interface)
      ! average up/dn components to obtain a consistent interface midpoint pressure
      v_darcy_aniso = (mu_up*v_darcy_aniso_up + mu_dn*v_darcy_aniso_dn) / (mu_up + mu_dn)
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%dv_darcy_aniso_up_dp_up_adjXm = 0.d0
        aniso_dds%dv_darcy_aniso_up_dp_up_adjXp = 0.d0

        aniso_dds%dv_darcy_aniso_up_dp_up_adjYm = -1.0d0 * k_up_xy * aniso_dds%ddphi_y_dp_up_adjYm * ukvr
        aniso_dds%dv_darcy_aniso_up_dp_up_adjYp = -1.0d0 * k_up_xy * aniso_dds%ddphi_y_dp_up_adjYp * ukvr

        aniso_dds%dv_darcy_aniso_up_dp_up_adjZm = -1.0d0 * k_up_xz * aniso_dds%ddphi_z_dp_up_adjZm * ukvr
        aniso_dds%dv_darcy_aniso_up_dp_up_adjZp = -1.0d0 * k_up_xz * aniso_dds%ddphi_z_dp_up_adjZp * ukvr

        aniso_dds%dv_darcy_aniso_up_dp_up = &
              -1.d0 * ((k_up_xy * aniso_dds%ddphi_y_dp_up + k_up_xz * aniso_dds%ddphi_z_dp_up) * ukvr &
              + (k_up_xy * dphi_y_up +  k_up_xz * dphi_z_up) * dukvr_dp_up)

        aniso_dds%dv_darcy_aniso_dn_dp_dn_adjXm = 0.d0
        aniso_dds%dv_darcy_aniso_dn_dp_dn_adjXp = 0.d0

        aniso_dds%dv_darcy_aniso_dn_dp_dn_adjYm = -1.0d0 * k_dn_xy * aniso_dds%ddphi_y_dp_dn_adjYm * ukvr
        aniso_dds%dv_darcy_aniso_dn_dp_dn_adjYp = -1.0d0 * k_dn_xy * aniso_dds%ddphi_y_dp_dn_adjYp * ukvr

        aniso_dds%dv_darcy_aniso_dn_dp_dn_adjZm = -1.0d0 * k_dn_xz * aniso_dds%ddphi_z_dp_dn_adjZm * ukvr
        aniso_dds%dv_darcy_aniso_dn_dp_dn_adjZp = -1.0d0 * k_dn_xz * aniso_dds%ddphi_z_dp_dn_adjZp * ukvr

        aniso_dds%dv_darcy_aniso_dn_dp_dn = &
              -1.d0 * ((k_dn_xy * aniso_dds%ddphi_y_dp_dn + k_dn_xz * aniso_dds%ddphi_z_dp_dn) * ukvr &
              + (k_dn_xy * dphi_y_dn +  k_dn_xz * dphi_z_dn) * dukvr_dp_dn)
      endif
    endif
    if ( aniso_ds%orientation .eq. 'Y' ) then
      v_darcy_aniso_up = -1.d0 * (k_up_xy * dphi_x_up + k_up_yz * dphi_z_up) * ukvr
      v_darcy_aniso_dn = -1.d0 * (k_dn_xy * dphi_x_dn + k_dn_yz * dphi_z_dn) * ukvr
      mu_up =       dist(-1) *dist(0)/k_up_yy  ! weight - up (dist part = dist to interface)
      mu_dn = (1.d0-dist(-1))*dist(0)/k_dn_yy  ! weight - dn (dist part = dist to interface)
      ! average up/dn components to obtain a consistent interface midpoint pressure
      v_darcy_aniso = (mu_up*v_darcy_aniso_up + mu_dn*v_darcy_aniso_dn) / (mu_up + mu_dn)
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%dv_darcy_aniso_up_dp_up_adjXm = -1.0d0 * k_up_xy * aniso_dds%ddphi_x_dp_up_adjXm * ukvr
        aniso_dds%dv_darcy_aniso_up_dp_up_adjXp = -1.0d0 * k_up_xy * aniso_dds%ddphi_x_dp_up_adjXp * ukvr

        aniso_dds%dv_darcy_aniso_up_dp_up_adjYm = 0.d0
        aniso_dds%dv_darcy_aniso_up_dp_up_adjYp = 0.d0

        aniso_dds%dv_darcy_aniso_up_dp_up_adjZm = -1.0d0 * k_up_yz * aniso_dds%ddphi_z_dp_up_adjZm * ukvr
        aniso_dds%dv_darcy_aniso_up_dp_up_adjZp = -1.0d0 * k_up_yz * aniso_dds%ddphi_z_dp_up_adjZp * ukvr

        aniso_dds%dv_darcy_aniso_up_dp_up = &
              -1.d0 * ((k_up_xy * aniso_dds%ddphi_x_dp_up + k_up_yz * aniso_dds%ddphi_z_dp_up) * ukvr &
              + (k_up_xy * dphi_x_up +  k_up_yz * dphi_z_up) * dukvr_dp_up)

        aniso_dds%dv_darcy_aniso_dn_dp_dn_adjXm = -1.0d0 * k_dn_xy * aniso_dds%ddphi_x_dp_dn_adjXm * ukvr
        aniso_dds%dv_darcy_aniso_dn_dp_dn_adjXp = -1.0d0 * k_dn_xy * aniso_dds%ddphi_x_dp_dn_adjXp * ukvr

        aniso_dds%dv_darcy_aniso_dn_dp_dn_adjYm = 0.d0
        aniso_dds%dv_darcy_aniso_dn_dp_dn_adjYp = 0.d0

        aniso_dds%dv_darcy_aniso_dn_dp_dn_adjZm = -1.0d0 * k_dn_yz * aniso_dds%ddphi_z_dp_dn_adjZm * ukvr
        aniso_dds%dv_darcy_aniso_dn_dp_dn_adjZp = -1.0d0 * k_dn_yz * aniso_dds%ddphi_z_dp_dn_adjZp * ukvr

        aniso_dds%dv_darcy_aniso_dn_dp_dn = &
              -1.d0 * ((k_dn_xy * aniso_dds%ddphi_x_dp_dn + k_dn_yz * aniso_dds%ddphi_z_dp_dn) * ukvr &
              + (k_dn_xy * dphi_x_dn +  k_dn_yz * dphi_z_dn) * dukvr_dp_dn)
      endif
    endif
    if ( aniso_ds%orientation .eq. 'Z' ) then
      v_darcy_aniso_up = -1.d0 * (k_up_xz * dphi_x_up + k_up_yz * dphi_y_up) * ukvr
      v_darcy_aniso_dn = -1.d0 * (k_dn_xz * dphi_x_dn + k_dn_yz * dphi_y_dn) * ukvr
      mu_up =       dist(-1) *dist(0)/k_up_zz  ! weight - up (dist part = dist to interface)
      mu_dn = (1.d0-dist(-1))*dist(0)/k_dn_zz  ! weight - dn (dist part = dist to interface)
      ! average up/dn components to obtain a consistent interface midpoint pressure
      v_darcy_aniso = (mu_up*v_darcy_aniso_up + mu_dn*v_darcy_aniso_dn) / (mu_up + mu_dn)
      if (calc_derivs .EQV. PETSC_TRUE) then
        aniso_dds%dv_darcy_aniso_up_dp_up_adjXm = -1.0d0 * k_up_xz * aniso_dds%ddphi_x_dp_up_adjXm * ukvr
        aniso_dds%dv_darcy_aniso_up_dp_up_adjXp = -1.0d0 * k_up_xz * aniso_dds%ddphi_x_dp_up_adjXp * ukvr

        aniso_dds%dv_darcy_aniso_up_dp_up_adjYm = -1.0d0 * k_up_yz * aniso_dds%ddphi_y_dp_up_adjYm * ukvr
        aniso_dds%dv_darcy_aniso_up_dp_up_adjYp = -1.0d0 * k_up_yz * aniso_dds%ddphi_y_dp_up_adjYp * ukvr

        aniso_dds%dv_darcy_aniso_up_dp_up_adjZm = 0.d0
        aniso_dds%dv_darcy_aniso_up_dp_up_adjZp = 0.d0

        aniso_dds%dv_darcy_aniso_up_dp_up = &
              -1.d0 * ((k_up_xz * aniso_dds%ddphi_x_dp_up + k_up_yz * aniso_dds%ddphi_y_dp_up) * ukvr &
              + (k_up_xz * dphi_x_up +  k_up_yz * dphi_y_up) * dukvr_dp_up)

        aniso_dds%dv_darcy_aniso_dn_dp_dn_adjXm = -1.0d0 * k_dn_xz * aniso_dds%ddphi_x_dp_dn_adjXm * ukvr
        aniso_dds%dv_darcy_aniso_dn_dp_dn_adjXp = -1.0d0 * k_dn_xz * aniso_dds%ddphi_x_dp_dn_adjXp * ukvr

        aniso_dds%dv_darcy_aniso_dn_dp_dn_adjYm = -1.0d0 * k_dn_yz * aniso_dds%ddphi_y_dp_dn_adjYm * ukvr
        aniso_dds%dv_darcy_aniso_dn_dp_dn_adjYp = -1.0d0 * k_dn_yz * aniso_dds%ddphi_y_dp_dn_adjYp * ukvr

        aniso_dds%dv_darcy_aniso_dn_dp_dn_adjZm = 0.d0
        aniso_dds%dv_darcy_aniso_dn_dp_dn_adjZp = 0.d0

        aniso_dds%dv_darcy_aniso_dn_dp_dn = &
              -1.d0 * ((k_dn_xz * aniso_dds%ddphi_x_dp_dn + k_dn_yz * aniso_dds%ddphi_y_dp_dn) * ukvr &
              + (k_dn_xz * dphi_x_dn +  k_dn_yz * dphi_y_dn) * dukvr_dp_dn)
      endif
    endif
    
    ! dbg
    ! call PrintAnisoData(aniso_data)
    ! write (*,"(' v_darcy_aniso=',e10.3)") v_darcy_aniso

    if (calc_derivs .EQV. PETSC_TRUE) then
      aniso_ds%dv_darcy_aniso_dp_up_adjXm = (mu_up * aniso_dds%dv_darcy_aniso_up_dp_up_adjXm) &
                                              / (mu_up + mu_dn)
      aniso_ds%dv_darcy_aniso_dp_dn_adjXm = (mu_dn * aniso_dds%dv_darcy_aniso_dn_dp_dn_adjXm) &
                                              / (mu_up + mu_dn)
      aniso_ds%dv_darcy_aniso_dp_up_adjXp = (mu_up * aniso_dds%dv_darcy_aniso_up_dp_up_adjXp) &
                                              / (mu_up + mu_dn)
      aniso_ds%dv_darcy_aniso_dp_dn_adjXp = (mu_dn * aniso_dds%dv_darcy_aniso_dn_dp_dn_adjXp) &
                                              / (mu_up + mu_dn)

      aniso_ds%dv_darcy_aniso_dp_up_adjYm = (mu_up * aniso_dds%dv_darcy_aniso_up_dp_up_adjYm) &
                                              / (mu_up + mu_dn)
      aniso_ds%dv_darcy_aniso_dp_dn_adjYm = (mu_dn * aniso_dds%dv_darcy_aniso_dn_dp_dn_adjYm) &
                                              / (mu_up + mu_dn)
      aniso_ds%dv_darcy_aniso_dp_up_adjYp = (mu_up * aniso_dds%dv_darcy_aniso_up_dp_up_adjYp) &
                                              / (mu_up + mu_dn)
      aniso_ds%dv_darcy_aniso_dp_dn_adjYp = (mu_dn * aniso_dds%dv_darcy_aniso_dn_dp_dn_adjYp) &
                                              / (mu_up + mu_dn)

      aniso_ds%dv_darcy_aniso_dp_up_adjZm = (mu_up * aniso_dds%dv_darcy_aniso_up_dp_up_adjZm) &
                                              / (mu_up + mu_dn)
      aniso_ds%dv_darcy_aniso_dp_dn_adjZm = (mu_dn * aniso_dds%dv_darcy_aniso_dn_dp_dn_adjZm) &
                                              / (mu_up + mu_dn)
      aniso_ds%dv_darcy_aniso_dp_up_adjZp = (mu_up * aniso_dds%dv_darcy_aniso_up_dp_up_adjZp) &
                                              / (mu_up + mu_dn)
      aniso_ds%dv_darcy_aniso_dp_dn_adjZp = (mu_dn * aniso_dds%dv_darcy_aniso_dn_dp_dn_adjZp) &
                                              / (mu_up + mu_dn)

      aniso_ds%dv_darcy_aniso_dp_up = (mu_up * aniso_dds%dv_darcy_aniso_up_dp_up) / (mu_up + mu_dn)
      aniso_ds%dv_darcy_aniso_dp_dn = (mu_dn * aniso_dds%dv_darcy_aniso_dn_dp_dn) / (mu_up + mu_dn)

      deallocate(aniso_dds)
    endif
  else
    ! unstructured grid case
    aniso_du => aniso_data%unstruct

    ! !dbg 
    ! call PrintAnisoData(aniso_data)

    ! v,w are the coordinate vectors in the plane of the interface over which flux is being calculated
    dphi_v_up = 0
    dphi_w_up = 0
    dphi_v_dn = 0
    dphi_w_dn = 0

    if (calc_derivs) then
      allocate(aniso_ddu)
      allocate(aniso_ddu%ddphi_v_up_dp_up_adj(aniso_du%up%num_adj)) 
      allocate(aniso_ddu%ddphi_w_up_dp_up_adj(aniso_du%up%num_adj)) 
      allocate(aniso_ddu%ddphi_v_dn_dp_dn_adj(aniso_du%dn%num_adj))
      allocate(aniso_ddu%ddphi_w_dn_dp_dn_adj(aniso_du%dn%num_adj))
      aniso_ddu%ddphi_v_up_dp_up     = 0 ! scalar
      aniso_ddu%ddphi_v_up_dp_dn     = 0 ! scalar
      aniso_ddu%ddphi_v_up_dp_up_adj = 0 ! array
      aniso_ddu%ddphi_w_up_dp_up     = 0 ! scalar
      aniso_ddu%ddphi_w_up_dp_dn     = 0 ! scalar
      aniso_ddu%ddphi_w_up_dp_up_adj = 0 ! array
      aniso_ddu%ddphi_v_dn_dp_dn     = 0 ! scalar
      aniso_ddu%ddphi_v_dn_dp_up     = 0 ! scalar
      aniso_ddu%ddphi_v_dn_dp_dn_adj = 0 ! array
      aniso_ddu%ddphi_w_dn_dp_dn     = 0 ! scalar
      aniso_ddu%ddphi_w_dn_dp_up     = 0 ! scalar
      aniso_ddu%ddphi_w_dn_dp_dn_adj = 0 ! array
    endif


    do iupdn=1,2

      if (iupdn==1) then
        ! upstream
        aniso_updn => aniso_du%up
        pres_updn = global_auxvar_up%pres(1)
        rho_updn = rho_up
      else
        ! downstream
        aniso_updn => aniso_du%dn
        pres_updn = global_auxvar_dn%pres(1)
        rho_updn = rho_dn
      endif

      ! Calculate ddphi_du from v_darcy (= Dq * ukvr * dphi) supplied by RichardsFlux(),
      ! where Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
      ! (n.b ddphi_du is basically -(P_up-P_dn)/dist) + gravity terms)
      dphi_du = -1.d0* v_darcy_iso / (ukvr*aniso_updn%perm_uu)
      d2phi_du_dp_up = -1.d0* dv_darcy_iso_dp_up / (ukvr*aniso_updn%perm_uu) &
                       + v_darcy_iso * dukvr_dp_up / (ukvr*ukvr*aniso_updn%perm_uu)
      d2phi_du_dp_dn = -1.d0* dv_darcy_iso_dp_dn / (ukvr*aniso_updn%perm_uu) &
                       + v_darcy_iso * dukvr_dp_dn / (ukvr*ukvr*aniso_updn%perm_uu)

      do iadj = 1,aniso_updn%num_adj
        ! Could form inv(M)*sum(rhs_terms) here, but it's easier to accumulate the derivative
        ! if written as sum(inv(M)*rhs_term)
        dot = aniso_updn%dists_vw(1,iadj)*aniso_updn%dists_vw(1,iadj) &
              + aniso_updn%dists_vw(2,iadj)*aniso_updn%dists_vw(2,iadj)
        res =  (aniso_updn%pres_adj(iadj) - rho_updn*aniso_updn%adu_cell_adj(iadj)%ptr%loc_gravity ) &
                 - (pres_updn - rho_updn*aniso_updn%loc_gravity) &
                 - dphi_du*aniso_updn%dists_u(iadj)

        LSrhs(1) = (1.d0/dot) * aniso_updn%dists_vw(1,iadj)
        LSrhs(2) = (1.d0/dot) * aniso_updn%dists_vw(2,iadj)

        if (iupdn==1) then
          ! upstream
          ! apply inverse of LS matrix-residual
          dphi_v_up = dphi_v_up + aniso_updn%LSMatInv(1,1)*LSrhs(1)*res + aniso_updn%LSMatInv(1,2)*LSrhs(2)*res
          dphi_w_up = dphi_w_up + aniso_updn%LSMatInv(2,1)*LSrhs(1)*res + aniso_updn%LSMatInv(2,2)*LSrhs(2)*res
          if (calc_derivs) then

            dres_dp_up  = -1.d0 - d2phi_du_dp_up*aniso_updn%dists_u(iadj)
            dres_dp_dn  =       - d2phi_du_dp_dn*aniso_updn%dists_u(iadj)
            dres_dp_adj =  1.d0
    
            aniso_ddu%ddphi_v_up_dp_up = aniso_ddu%ddphi_v_up_dp_up                         &
                                         + aniso_updn%LSMatInv(1,1) * LSrhs(1) * dres_dp_up &
                                         + aniso_updn%LSMatInv(1,2) * LSrhs(2) * dres_dp_up 
            aniso_ddu%ddphi_v_up_dp_dn = aniso_ddu%ddphi_v_up_dp_dn                         &
                                         + aniso_updn%LSMatInv(1,1) * LSrhs(1) * dres_dp_dn &
                                         + aniso_updn%LSMatInv(1,2) * LSrhs(2) * dres_dp_dn 

            aniso_ddu%ddphi_w_up_dp_up = aniso_ddu%ddphi_w_up_dp_up                         &
                                         + aniso_updn%LSMatInv(2,1) * LSrhs(1) * dres_dp_up &
                                         + aniso_updn%LSMatInv(2,2) * LSrhs(2) * dres_dp_up
            aniso_ddu%ddphi_w_up_dp_dn = aniso_ddu%ddphi_w_up_dp_dn                         &
                                         + aniso_updn%LSMatInv(2,1) * LSrhs(1) * dres_dp_dn &
                                         + aniso_updn%LSMatInv(2,2) * LSrhs(2) * dres_dp_dn

            aniso_ddu%ddphi_v_up_dp_up_adj(iadj) = aniso_updn%LSMatInv(1,1) * LSrhs(1) * dres_dp_adj &
                                                 + aniso_updn%LSMatInv(1,2) * LSrhs(2) * dres_dp_adj
            aniso_ddu%ddphi_w_up_dp_up_adj(iadj) = aniso_updn%LSMatInv(2,1) * LSrhs(1) * dres_dp_adj &
                                                 + aniso_updn%LSMatInv(2,2) * LSrhs(2) * dres_dp_adj
          endif
        else
          ! downstream
          ! apply inverse of LS matrix-residual
          dphi_v_dn = dphi_v_dn + aniso_updn%LSMatInv(1,1)*LSrhs(1)*res + aniso_updn%LSMatInv(1,2)*LSrhs(2)*res
          dphi_w_dn = dphi_w_dn + aniso_updn%LSMatInv(2,1)*LSrhs(1)*res + aniso_updn%LSMatInv(2,2)*LSrhs(2)*res
          if (calc_derivs) then

            dres_dp_up  =       - d2phi_du_dp_up*aniso_updn%dists_u(iadj)
            dres_dp_dn  = -1.d0 - d2phi_du_dp_dn*aniso_updn%dists_u(iadj)
            dres_dp_adj =  1.d0
    
            aniso_ddu%ddphi_v_dn_dp_dn = aniso_ddu%ddphi_v_dn_dp_dn                         &
                                         + aniso_updn%LSMatInv(1,1) * LSrhs(1) * dres_dp_dn &
                                         + aniso_updn%LSMatInv(1,2) * LSrhs(2) * dres_dp_dn
            aniso_ddu%ddphi_v_dn_dp_up = aniso_ddu%ddphi_v_dn_dp_up                         &
                                         + aniso_updn%LSMatInv(1,1) * LSrhs(1) * dres_dp_up &
                                         + aniso_updn%LSMatInv(1,2) * LSrhs(2) * dres_dp_up
                                         
            aniso_ddu%ddphi_w_dn_dp_dn = aniso_ddu%ddphi_w_dn_dp_dn                         &
                                         + aniso_updn%LSMatInv(2,1) * LSrhs(1) * dres_dp_dn &
                                         + aniso_updn%LSMatInv(2,2) * LSrhs(2) * dres_dp_dn
            aniso_ddu%ddphi_w_dn_dp_up = aniso_ddu%ddphi_w_dn_dp_up                         &
                                         + aniso_updn%LSMatInv(2,1) * LSrhs(1) * dres_dp_up &
                                         + aniso_updn%LSMatInv(2,2) * LSrhs(2) * dres_dp_up
                                        
            aniso_ddu%ddphi_v_dn_dp_dn_adj(iadj) = aniso_updn%LSMatInv(1,1) * LSrhs(1) * dres_dp_adj &
                                                 + aniso_updn%LSMatInv(1,2) * LSrhs(2) * dres_dp_adj
            aniso_ddu%ddphi_w_dn_dp_dn_adj(iadj) = aniso_updn%LSMatInv(2,1) * LSrhs(1) * dres_dp_adj &
                                                 + aniso_updn%LSMatInv(2,2) * LSrhs(2) * dres_dp_adj
          endif
        endif
      enddo

    enddo

    ! Calculate anisotropic darcy contribution in up/dn cells
    v_darcy_aniso_up = -1.d0 * (aniso_du%up%perm_uv * dphi_v_up + aniso_du%up%perm_uw * dphi_w_up) * ukvr
    v_darcy_aniso_dn = -1.d0 * (aniso_du%dn%perm_uv * dphi_v_dn + aniso_du%dn%perm_uw * dphi_w_dn) * ukvr
    mu_up =       dist(-1) *dist(0)/aniso_du%up%perm_uu  ! weight - up (dist to interface/Kuu)
    mu_dn = (1.d0-dist(-1))*dist(0)/aniso_du%dn%perm_uu  ! weight - dn (dist to interface/Kuu)
    ! average up/dn components to obtain a consistent interface midpoint pressure
    v_darcy_aniso = (mu_up*v_darcy_aniso_up + mu_dn*v_darcy_aniso_dn) / (mu_up + mu_dn)

    if (calc_derivs) then
      aniso_du%dv_darcy_aniso_dp_up = (mu_up/(mu_up+mu_dn)) * &
                                      (-1.0d0 * (  aniso_du%up%perm_uv * aniso_ddu%ddphi_v_up_dp_up &
                                                 + aniso_du%up%perm_uw * aniso_ddu%ddphi_w_up_dp_up) * ukvr &
                                       -1.0d0 * (aniso_du%up%perm_uv * dphi_v_up + aniso_du%up%perm_uw * dphi_w_up) * dukvr_dp_up) &
                                       + (mu_dn/(mu_up+mu_dn)) * &
                                       (-1.0d0 * (  aniso_du%dn%perm_uv * aniso_ddu%ddphi_v_dn_dp_up &
                                                  + aniso_du%dn%perm_uw * aniso_ddu%ddphi_w_dn_dp_up) * ukvr &
                                       -1.0d0 * (aniso_du%dn%perm_uv * dphi_v_dn + aniso_du%dn%perm_uw * dphi_w_dn) * dukvr_dp_up)
      aniso_du%dv_darcy_aniso_dp_dn = (mu_dn/(mu_up+mu_dn)) * &
                                      (-1.0d0 * (  aniso_du%dn%perm_uv * aniso_ddu%ddphi_v_dn_dp_dn &
                                                 + aniso_du%dn%perm_uw * aniso_ddu%ddphi_w_dn_dp_dn) * ukvr &
                                       -1.0d0 * (aniso_du%dn%perm_uv * dphi_v_dn + aniso_du%dn%perm_uw * dphi_w_dn) * dukvr_dp_dn) &
                                       + (mu_up/(mu_up+mu_dn)) * &
                                       (-1.0d0 * (  aniso_du%up%perm_uv * aniso_ddu%ddphi_v_up_dp_dn &
                                                  + aniso_du%up%perm_uw * aniso_ddu%ddphi_w_up_dp_dn) * ukvr &
                                        -1.0d0 * (aniso_du%up%perm_uv * dphi_v_up + aniso_du%up%perm_uw * dphi_w_up) * dukvr_dp_dn) 
      do iadj=1,aniso_du%up%num_adj
        aniso_du%dv_darcy_aniso_dp_up_adj(iadj) = (mu_up/(mu_up+mu_dn)) * &
                                               -1.0d0 * (  aniso_du%up%perm_uv * aniso_ddu%ddphi_v_up_dp_up_adj(iadj) &
                                                         + aniso_du%up%perm_uw * aniso_ddu%ddphi_w_up_dp_up_adj(iadj)) * ukvr
      enddo
      do iadj=1,aniso_du%dn%num_adj
        aniso_du%dv_darcy_aniso_dp_dn_adj(iadj) = (mu_dn/(mu_up+mu_dn)) * &
                                               -1.0d0 * (  aniso_du%dn%perm_uv * aniso_ddu%ddphi_v_dn_dp_dn_adj(iadj) &
                                                         + aniso_du%dn%perm_uw * aniso_ddu%ddphi_w_dn_dp_dn_adj(iadj)) * ukvr
      enddo
      ! dbg
      ! call PrintAnisoData(aniso_data)
      ! write (*,"(' v_darcy_aniso=',e10.3)") v_darcy_aniso
    endif

    if (calc_derivs) then
      deallocate(aniso_ddu)
    endif
  endif
end subroutine RichardsFluxAniso

! ************************************************************************** !

subroutine RichardsFlux(rich_auxvar_up,global_auxvar_up, &
                        material_auxvar_up, &
                        rich_auxvar_dn,global_auxvar_dn, &
                        material_auxvar_dn, &
                        area, dist, &
                        aniso_data, &
                        option,v_darcy,Res)
  ! 
  ! Computes the internal flux terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 
  use Option_module
  use Material_Aux_class
  use Connection_module
  use Anisotropy_Data_module
  
  implicit none
  
  type(richards_auxvar_type) :: rich_auxvar_up, rich_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(aniso_data_type), pointer :: aniso_data
  type(option_type) :: option
  PetscReal :: v_darcy, area, dist(-1:3)
  PetscReal :: Res(1:option%nflowdof) 
     
  PetscInt :: ispec
  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dd_up, dd_dn, perm_up, perm_dn
  PetscReal :: fluxm, q
  PetscReal :: ukvr,Dq
  PetscReal :: upweight, density_ave, cond, gravity, dphi
  
  PetscReal :: v_darcy_aniso

  fluxm = 0.d0
  v_darcy = 0.D0  
  ukvr = 0.d0
  
  call ConnectionCalculateDistances(dist,option%gravity,dd_up,dd_dn, &
                                    dist_gravity,upweight)
  if(.not. associated(aniso_data)) then
    call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
    call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  else
    ! Flux component due head gradients directly across the interface can be
    ! calculated using the isotropic approach provided that permeabilty is set
    ! consistent with anisotropic calculation (avoid collapsing to scalar).
    ! Additional anisotropic terms due to tangential head gradients are added
    ! by later call to RichardsFluxAniso()
    if(associated(aniso_data%struct)) then
      ! sjb - can do same calls as isotropic case here, should just extract
      !       k_xx, k_yy or k_zz depending on orientation
      ! sjb - would be better to be specific here rather than rely on
      !       TensorToScalar call to do the right thing
      call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
      call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
    else
      ! u-direction is chosen to be parallel to centre->centre direction, so
      perm_up = aniso_data%unstruct%up%perm_uu
      perm_dn = aniso_data%unstruct%dn%perm_uu
    endif
  endif

  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
  
! Flow term
  if (rich_auxvar_up%kvr > eps .or. &
      rich_auxvar_dn%kvr > eps) then
    if (global_auxvar_up%sat(1) <eps) then 
      upweight=0.d0
    else if (global_auxvar_dn%sat(1) <eps) then 
      upweight=1.d0
    endif
    

    density_ave = upweight*global_auxvar_up%den(1)+ &
                  (1.D0-upweight)*global_auxvar_dn%den(1)

    gravity = (upweight*global_auxvar_up%den(1) + &
              (1.D0-upweight)*global_auxvar_dn%den(1)) &
              * FMWH2O * dist_gravity

    dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1)  + gravity

    if (dphi>=0.D0) then
      ukvr = rich_auxvar_up%kvr
    else
      ukvr = rich_auxvar_dn%kvr
    endif      

    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi
   
      ! add anisotropic contributions
      if(associated(aniso_data)) then
        ! Call PrintAnisoData(aniso_data)
        Call RichardsFluxAniso(aniso_data,material_auxvar_up,material_auxvar_dn,&
                                  global_auxvar_up,global_auxvar_dn,&
                                  option,dist,&
                                  ukvr,0.d0,0.d0,&
                                  v_darcy,0.d0,0.d0,&
                                  PETSC_FALSE,v_darcy_aniso)
        v_darcy = v_darcy + v_darcy_aniso
      endif
    endif

    q = v_darcy * area

    fluxm = q*density_ave
  endif 

  Res(1) = fluxm
 ! note: Res is the flux contribution, for node 1 R = R + Res_FL
 !                                              2 R = R - Res_FL  

end subroutine RichardsFlux

! ************************************************************************** !

subroutine RichardsBCFluxDerivative(ibndtype,auxvars, &
                                    rich_auxvar_up,global_auxvar_up, &
                                    rich_auxvar_dn,global_auxvar_dn, &
                                    material_auxvar_dn, &
                                    area,dist,option, &
                                    characteristic_curves_dn, &
                                    Jdn)
  ! 
  ! Computes the derivatives of the boundary flux
  ! terms for the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 
  use Option_module
  use Characteristic_Curves_module
  use Material_Aux_class
  use EOS_Water_module
  use Utility_module
 
  implicit none
  
  PetscInt :: ibndtype(:)
  type(richards_auxvar_type) :: rich_auxvar_up, rich_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(option_type) :: option
  PetscReal :: auxvars(:) ! from aux_real_var array in boundary condition
  PetscReal :: area
  ! dist(-1) = fraction_upwind
  ! dist(0) = magnitude
  ! dist(1:3) = unit vector
  ! dist(0)*dist(1:3) = vector
  PetscReal :: dist(-1:3)
  type(characteristic_curves_type) :: characteristic_curves_dn
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  
  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: perm_dn
  
  PetscReal :: v_darcy
  PetscReal :: q,density_ave
  PetscReal :: ukvr,diffdp,Dq
  PetscReal :: upweight,cond,gravity,dphi

  PetscReal :: dden_ave_dp_dn
  PetscReal :: dgravity_dden_dn
  PetscReal :: dphi_dp_dn
  PetscReal :: dukvr_dp_dn
  PetscReal :: dq_dp_dn
  PetscInt :: pressure_bc_type
  PetscReal :: dphi_x_dp_dn,dphi_y_dp_dn,dphi_z_dp_dn
  PetscReal :: dHdn_x_dp_dn, dHdn_y_dp_dn, dHdn_z_dp_dn

  PetscInt :: iphase, ideriv
  type(richards_auxvar_type) :: rich_auxvar_pert_dn, rich_auxvar_pert_up
  type(global_auxvar_type) :: global_auxvar_pert_dn, global_auxvar_pert_up
  class(material_auxvar_type), allocatable :: material_auxvar_pert_dn, &
                                              material_auxvar_pert_up
  PetscReal :: perturbation
  PetscReal :: x_dn(1), x_up(1), x_pert_dn(1), x_pert_up(1), pert_dn, res(1), &
            res_pert_dn(1), J_pert_dn(1,1)
  PetscReal :: rho
  PetscReal :: dum1
  PetscReal :: dq_lin,dP_lin
  PetscReal :: q_approx, dq_approx
  PetscErrorCode :: ierr

  v_darcy = 0.d0
  ukvr = 0.d0
  density_ave = 0.d0
  q = 0.d0

  Jdn = 0.d0 
  
  dden_ave_dp_dn = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_dn = 0.d0
  dukvr_dp_dn = 0.d0
  dq_dp_dn = 0.d0

  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  
  ! Flow
  pressure_bc_type = ibndtype(RICHARDS_PRESSURE_DOF)
  select case(pressure_bc_type)
    ! figure out the direction of flow
    case(DIRICHLET_BC,DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
         HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC, &
         HET_SURF_HYDROSTATIC_SEEPAGE_BC, &
         HET_DIRICHLET_BC,HET_HYDROSTATIC_SEEPAGE_BC, &
         HET_HYDROSTATIC_CONDUCTANCE_BC)

      ! dist(0) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3) = vector(3) - unit vector
      dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))

      if (pressure_bc_type == HYDROSTATIC_CONDUCTANCE_BC .or. &
          pressure_bc_type == DIRICHLET_CONDUCTANCE_BC .or. &
          pressure_bc_type == HET_HYDROSTATIC_CONDUCTANCE_BC) then
        Dq = auxvars(RICHARDS_CONDUCTANCE_DOF)
      else
        Dq = perm_dn / dist(0)
      endif
      ! Flow term
      if (rich_auxvar_up%kvr > eps .or. &
          rich_auxvar_dn%kvr > eps) then
        upweight=1.D0
        if (global_auxvar_up%sat(1) < eps) then 
          upweight=0.d0
        else if (global_auxvar_dn%sat(1) < eps) then 
          upweight=1.d0
        endif
        
        density_ave = upweight*global_auxvar_up%den(1) + &
                      (1.D0-upweight)*global_auxvar_dn%den(1)
        dden_ave_dp_dn = (1.D0-upweight)*rich_auxvar_dn%dden_dp

        gravity = (upweight*global_auxvar_up%den(1) + &
                  (1.D0-upweight)*global_auxvar_dn%den(1)) &
                  * FMWH2O * dist_gravity
        dgravity_dden_dn = (1.d0-upweight)*FMWH2O*dist_gravity

        dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity
        dphi_dp_dn = -1.d0 + dgravity_dden_dn*rich_auxvar_dn%dden_dp

        select case(pressure_bc_type)
          case(HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC, &
               DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
               HET_SURF_HYDROSTATIC_SEEPAGE_BC, &
               HET_HYDROSTATIC_SEEPAGE_BC,HET_HYDROSTATIC_CONDUCTANCE_BC)
                ! flow in
            if (dphi > 0.d0 .and. &
                ! boundary cell is <= pref
                global_auxvar_up%pres(1)-option%reference_pressure < eps) then
              dphi = 0.d0
              dphi_dp_dn = 0.d0
            endif
        end select
        
        if (dphi>=0.D0) then
          ukvr = rich_auxvar_up%kvr
        else
          ukvr = rich_auxvar_dn%kvr
          dukvr_dp_dn = rich_auxvar_dn%dkvr_dp
        endif      
     
        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
          q = v_darcy * area
          dq_dp_dn = Dq*(dukvr_dp_dn*dphi + ukvr*dphi_dp_dn)*area

          ! If running with surface-flow model, ensure (darcy_velocity*dt) does
          ! not exceed depth of standing water.
          if (option%surf_flow_on) then
          if (Equal(rich_auxvar_dn%vars_for_sflow(11),0.d0)) then
            if (pressure_bc_type == HET_SURF_HYDROSTATIC_SEEPAGE_BC .and. &
                option%surf_flow_on) then
              call EOSWaterdensity(option%reference_temperature, &
                                   option%reference_pressure,rho,dum1,ierr)

              if (global_auxvar_dn%pres(1) <= &
                  rich_auxvar_dn%vars_for_sflow(1)) then

                ! Linear approximation
                call Interpolate(rich_auxvar_dn%vars_for_sflow(8), &
                                 rich_auxvar_dn%vars_for_sflow(7), &
                                 global_auxvar_dn%pres(1), &
                                 rich_auxvar_dn%vars_for_sflow(10), &
                                 rich_auxvar_dn%vars_for_sflow(9), &
                                 q_approx)
                v_darcy = q_approx/area
                q       = q_approx

                dP_lin = rich_auxvar_dn%vars_for_sflow(8) - &
                         rich_auxvar_dn%vars_for_sflow(7)
                dq_lin = rich_auxvar_dn%vars_for_sflow(10) - &
                         rich_auxvar_dn%vars_for_sflow(9)
                dq_dp_dn = dq_lin/dP_lin

              else
                if (global_auxvar_dn%pres(1) <= &
                    rich_auxvar_dn%vars_for_sflow(2)) then
                  ! Cubic approximation
                  call CubicPolynomialEvaluate( &
                        rich_auxvar_dn%vars_for_sflow(3:6), &
                        global_auxvar_dn%pres(1) - option%reference_pressure, &
                        q_approx, dq_approx)
                  v_darcy = q_approx/area
                  q = q_approx
                  dq_dp_dn = dq_approx
                endif
              endif
            endif
          endif
          endif

        endif

      endif

    case(NEUMANN_BC)
      if (dabs(auxvars(RICHARDS_PRESSURE_DOF)) > floweps) then
        v_darcy = auxvars(RICHARDS_PRESSURE_DOF)
        if (v_darcy > 0.d0) then 
          density_ave = global_auxvar_up%den(1)
        else 
          density_ave = global_auxvar_dn%den(1)
          dden_ave_dp_dn = rich_auxvar_dn%dden_dp
        endif 
        q = v_darcy * area
      endif

    case(UNIT_GRADIENT_BC)

      if(rich_auxvar_dn%kvr > eps) then
        dphi = dot_product(option%gravity,dist(1:3))* &
                           global_auxvar_dn%den(1)*FMWH2O      
        density_ave = global_auxvar_dn%den(1)

        dphi_dp_dn = dot_product(option%gravity,dist(1:3))* &
                                 rich_auxvar_dn%dden_dp*FMWH2O
        dden_ave_dp_dn = rich_auxvar_dn%dden_dp

        ! since boundary auxvar is meaningless (no pressure specified there), only use cell
        ukvr = rich_auxvar_dn%kvr
        dukvr_dp_dn = rich_auxvar_dn%dkvr_dp
     
        if (ukvr*perm_dn>floweps) then
          v_darcy = perm_dn * ukvr * dphi
          q = v_darcy*area
          dq_dp_dn = perm_dn*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
        endif 
     
      endif

  end select

  Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)

  if (option%flow%numerical_derivatives) then
    call GlobalAuxVarInit(global_auxvar_pert_up,option)
    call GlobalAuxVarInit(global_auxvar_pert_dn,option)  
    allocate(material_auxvar_pert_up,material_auxvar_pert_dn)
    call MaterialAuxVarInit(material_auxvar_pert_up,option)  
    call MaterialAuxVarInit(material_auxvar_pert_dn,option)  
    call RichardsAuxVarCopy(rich_auxvar_up,rich_auxvar_pert_up,option)
    call RichardsAuxVarCopy(rich_auxvar_dn,rich_auxvar_pert_dn,option)
    call GlobalAuxVarCopy(global_auxvar_up,global_auxvar_pert_up,option)
    call GlobalAuxVarCopy(global_auxvar_dn,global_auxvar_pert_dn,option)
    call MaterialAuxVarCopy(material_auxvar_dn,material_auxvar_pert_up, &
                            option)
    call MaterialAuxVarCopy(material_auxvar_dn,material_auxvar_pert_dn, &
                            option)
    x_up(1) = global_auxvar_up%pres(1)
    x_dn(1) = global_auxvar_dn%pres(1)
    ideriv = 1
    if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
      x_up(ideriv) = x_dn(ideriv)
    endif
    call RichardsBCFlux(ibndtype,auxvars, &
                        rich_auxvar_up,global_auxvar_up, &
                        rich_auxvar_dn,global_auxvar_dn, &
                        material_auxvar_dn, &
                        area,dist,option,v_darcy,res)
    if (pressure_bc_type == ZERO_GRADIENT_BC) then
      x_pert_up = x_up
    endif
    ideriv = 1
!    pert_dn = x_dn(ideriv)*perturbation_tolerance    
    pert_dn = max(dabs(x_dn(ideriv)*perturbation_tolerance),0.1d0)
    if (x_dn(ideriv) < option%reference_pressure) pert_dn = -1.d0*pert_dn
    x_pert_dn = x_dn
    x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
    x_pert_up = x_up
    if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
      x_pert_up(ideriv) = x_pert_dn(ideriv)
    endif   
    call RichardsAuxVarCompute(x_pert_dn(1),rich_auxvar_pert_dn, &
                               global_auxvar_pert_dn, &
                               material_auxvar_pert_dn, &
                               characteristic_curves_dn, &
                               -999, &
                               PETSC_TRUE,option)
    call RichardsAuxVarCompute(x_pert_up(1),rich_auxvar_pert_up, &
                               global_auxvar_pert_up, &
                               material_auxvar_pert_up, &
                               characteristic_curves_dn, &
                               -999, &
                               PETSC_TRUE,option)
    call RichardsBCFlux(ibndtype,auxvars, &
                        rich_auxvar_pert_up,global_auxvar_pert_up, &
                        rich_auxvar_pert_dn,global_auxvar_pert_dn, &
                        material_auxvar_pert_dn, &
                        area,dist,option,v_darcy,res_pert_dn)
    J_pert_dn(1,ideriv) = (res_pert_dn(1)-res(1))/pert_dn
    Jdn = J_pert_dn
    call GlobalAuxVarStrip(global_auxvar_pert_up)
    call GlobalAuxVarStrip(global_auxvar_pert_dn)   
    call MaterialAuxVarStrip(material_auxvar_pert_up)
    call MaterialAuxVarStrip(material_auxvar_pert_dn)
    deallocate(material_auxvar_pert_up,material_auxvar_pert_dn)
  endif

end subroutine RichardsBCFluxDerivative

! ************************************************************************** !

subroutine RichardsBCFlux(ibndtype,auxvars, &
                          rich_auxvar_up, global_auxvar_up, &
                          rich_auxvar_dn, global_auxvar_dn, &
                          material_auxvar_dn, &
                          area, dist, option,v_darcy,Res)
  ! 
  ! Computes the  boundary flux terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 
  use Option_module
  use Material_Aux_class
  use EOS_Water_module
  use Utility_module
 
  implicit none
  
  PetscInt :: ibndtype(:)
  type(richards_auxvar_type) :: rich_auxvar_up, rich_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(option_type) :: option
  PetscReal :: auxvars(:) ! from aux_real_var array
  PetscReal :: v_darcy, area
  ! dist(-1) = fraction_upwind - not applicable here
  ! dist(0) = magnitude
  ! dist(1:3) = unit vector
  ! dist(0)*dist(1:3) = vector
  PetscReal :: dist(-1:3)
  PetscReal :: Res(1:option%nflowdof) 
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscInt :: ispec
  PetscReal :: perm_dn
  PetscReal :: fluxm,q,density_ave
  PetscReal :: ukvr,diffdp,Dq
  PetscReal :: upweight,cond,gravity,dphi
  PetscInt :: pressure_bc_type
  PetscReal :: dphi_x,dphi_y,dphi_z
  PetscReal :: rho
  PetscReal :: dum1
  PetscReal :: q_approx, dq_approx
  PetscErrorCode :: ierr
  
  fluxm = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0
  ukvr = 0.d0

  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)  

  ! Flow  
  pressure_bc_type = ibndtype(RICHARDS_PRESSURE_DOF)
  select case(pressure_bc_type)
    ! figure out the direction of flow
    case(DIRICHLET_BC,DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
         HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC, &
         HET_SURF_HYDROSTATIC_SEEPAGE_BC, &
         HET_DIRICHLET_BC,HET_HYDROSTATIC_SEEPAGE_BC, &
         HET_HYDROSTATIC_CONDUCTANCE_BC)

      ! dist(0) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3) = vector(3) - unit vector
      dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))

      if (pressure_bc_type == HYDROSTATIC_CONDUCTANCE_BC .or. &
          pressure_bc_type == DIRICHLET_CONDUCTANCE_BC .or. &
          pressure_bc_type == HET_HYDROSTATIC_CONDUCTANCE_BC) then
        Dq = auxvars(RICHARDS_CONDUCTANCE_DOF)
      else
        Dq = perm_dn / dist(0)
      endif
      ! Flow term
      if (rich_auxvar_up%kvr > eps .or. &
          rich_auxvar_dn%kvr > eps) then
        upweight=1.D0
        if (global_auxvar_up%sat(1) < eps) then 
          upweight=0.d0
        else if (global_auxvar_dn%sat(1) < eps) then 
          upweight=1.d0
        endif
        density_ave = upweight*global_auxvar_up%den(1)+(1.D0-upweight)*global_auxvar_dn%den(1)
   
        gravity = (upweight*global_auxvar_up%den(1) + &
                  (1.D0-upweight)*global_auxvar_dn%den(1)) &
                  * FMWH2O * dist_gravity
       
        dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity

        select case(pressure_bc_type)
          case(HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC, &
               DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
               HET_SURF_HYDROSTATIC_SEEPAGE_BC, &
               HET_HYDROSTATIC_SEEPAGE_BC,HET_HYDROSTATIC_CONDUCTANCE_BC)
                ! flow in
            if (dphi > 0.d0 .and. &
                ! boundary cell is <= pref
                global_auxvar_up%pres(1)-option%reference_pressure < eps) then
              dphi = 0.d0
            endif
        end select
   
       if (dphi>=0.D0) then
         ukvr = rich_auxvar_up%kvr
       else
         ukvr = rich_auxvar_dn%kvr
       endif
        
       if (ukvr*Dq>floweps) then
        v_darcy = Dq * ukvr * dphi

        ! If running with surface-flow model, ensure (darcy_velocity*dt) does
        ! not exceed depth of standing water.
        if (pressure_bc_type == HET_SURF_HYDROSTATIC_SEEPAGE_BC .and. &
            option%surf_flow_on) then
          call EOSWaterdensity(option%reference_temperature, &
                               option%reference_pressure,rho,dum1,ierr)

          if (Equal(rich_auxvar_dn%vars_for_sflow(11),0.d0)) then
            if (global_auxvar_dn%pres(1) <= rich_auxvar_dn%vars_for_sflow(1)) then

              if (Equal(rich_auxvar_dn%vars_for_sflow(7),-99999.d0)) then
                call PrintErrMsg(option,'Coeffs for linear approx for darcy flux not set')
              endif

              ! Linear approximation
              call Interpolate(rich_auxvar_dn%vars_for_sflow(8), &
                               rich_auxvar_dn%vars_for_sflow(7), &
                               global_auxvar_dn%pres(1), &
                               rich_auxvar_dn%vars_for_sflow(2), &
                               rich_auxvar_dn%vars_for_sflow(1), &
                               q_approx)
              v_darcy = q_approx/area

            else if (global_auxvar_dn%pres(1) <= rich_auxvar_dn%vars_for_sflow(2)) then

              if (Equal(rich_auxvar_dn%vars_for_sflow(3),-99999.d0)) then
                call PrintErrMsg(option,'Coeffs for cubic approx for darcy flux not set')
              endif

              ! Cubic approximation
              call CubicPolynomialEvaluate(rich_auxvar_dn%vars_for_sflow(3:6), &
                                           global_auxvar_dn%pres(1) - option%reference_pressure, &
                                           q_approx, dq_approx)
              v_darcy = q_approx/area
            endif
          endif

        endif
       endif
      endif 

    case(NEUMANN_BC)
      if (dabs(auxvars(RICHARDS_PRESSURE_DOF)) > floweps) then
        v_darcy = auxvars(RICHARDS_PRESSURE_DOF)

        if (v_darcy > 0.d0) then 
          density_ave = global_auxvar_up%den(1)
        else 
          density_ave = global_auxvar_dn%den(1)
        endif 
      endif

    case(UNIT_GRADIENT_BC)

      dphi = dot_product(option%gravity,dist(1:3))*global_auxvar_dn%den(1)*FMWH2O
      density_ave = global_auxvar_dn%den(1)

      ! since boundary auxvar is meaningless (no pressure specified there), only use cell
      ukvr = rich_auxvar_dn%kvr
     
      if (ukvr*perm_dn>floweps) then
        v_darcy = perm_dn * ukvr * dphi
      endif 

  end select

  q = v_darcy * area

  fluxm = q*density_ave

  Res(1)=fluxm

end subroutine RichardsBCFlux

end module Richards_Common_module
