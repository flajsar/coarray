﻿module sim_step
    use basic_fnc
    use constants
    use interactions
    use motion
    use termo_barostat
    implicit none
    contains
        subroutine equilibration_step(x,v,a,particles,l,co_potencial_el,co_potencial_lj,step,pressure,temp,co_press_tensor,co_virial_tensor,coms,co_virial,LJ_spheres,LJ_spheres_force)
            real(8),intent(inout)::x(num_images(),mp*N/num_images(),3)[*],v(num_images(),3*n/num_images(),3)[*],a(num_images(),3*n/num_images(),3)[*],particles(mp*n,2),l[*],LJ_spheres(N_LJ_spheres,3),LJ_spheres_force(N_LJ_spheres,3)
            real(8),intent(inout)::co_potencial_el[*],co_potencial_lj[*],pressure[*],temp[*],co_press_tensor(3,3)[*],co_virial_tensor(3,3)[*],coms(num_images(),mp*N/num_images(),3)[*],co_virial[*]
            integer,intent(in)::step
            
            co_virial_tensor=0
            SYNC ALL
            !uporabimo rattle algoritem
            call rattle_polozaji(x,a,v,particles,l)
            if (mp==4) then
                call postavitev_m(x,l) !posodobi hitrosti na halfstep, zračuna g-je iterativno, zracuna nove popravljene položaje
            end if
            call synh_arrays(x)
            call synh_arrays(v)
            a=0
            LJ_spheres_force=0
            co_potencial_el=0
            co_potencial_lj=0
            call calc_coms(x,particles,l,coms)
            call ewald_sum(x,particles,a,l,co_virial,co_potencial_el,co_virial_tensor,coms)
            call lj_del_pot(x,particles,a,l,co_potencial_lj,co_virial,co_virial_tensor,coms)
            call LJ_sphere_correction(x,particles,a,l,co_potencial_lj,co_virial,co_virial_tensor,coms,LJ_spheres,LJ_spheres_force)
            !call sym_tensor(co_virial_tensor)
            sync all
            call CO_SUM(co_potencial_el)
            call CO_SUM(co_potencial_lj)
            call CO_SUM(co_virial_tensor)
            call CO_SUM(a)
            sync all
                !nove sile -> nove hitrosti
            call rattle_hit(x,a,v,particles,l)
            call synh_arrays(v)
            call calc_temp(v,particles,temp)
            call calc_pressure(co_virial_tensor,v,l,pressure,particles,x,a,co_press_tensor,temp)
            !if (mod(step,samplingst*10)==0) then
            call remove_comm(x,v,particles) !odstranimo center-of-mass motion
            !end if
            Sync all
            call berendsen_termostat(v,temp0,temp,particles,dt)!tau=dt
            call synh_arrays(v)
            !    !zapišemo rezultate
            if (this_image()==1) then
                !if (mod(step,samplingst*10)==0) then
                    print*,step,temp*temp_ref, pressure*press_ref/10**5
                !end if
            end if
        end subroutine
        
        subroutine nvt_nh_step(x,v,a,particles,l,co_virial,co_potencial_el,co_potencial_lj,step,pressure,temp,th_eta,co_press_tensor,coms,co_virial_tensor,LJ_spheres,LJ_spheres_force)
            real(8),intent(inout)::x(num_images(),mp*N/num_images(),3)[*],v(num_images(),3*n/num_images(),3)[*],a(num_images(),3*n/num_images(),3)[*],particles(mp*n,2),l[*],LJ_spheres(N_LJ_spheres,3),LJ_spheres_force(N_LJ_spheres,3)
            real(8),intent(inout)::co_potencial_el[*],co_potencial_lj[*],pressure[*],temp[*],co_press_tensor(3,3)[*],co_virial_tensor(3,3)[*],coms(num_images(),mp*N/num_images(),3)[*],co_virial[*],th_eta[*]
            integer,intent(in)::step
            
            co_virial_tensor=0
            SYNC ALL
            !uporabimo rattle algoritem
            call calc_th_eta(th_eta,v,particles)
            call nh_v_adjust(v,th_eta)
            call calc_th_eta(th_eta,v,particles)
            !uporabimo rattle algoritem
            
            call rattle_polozaji(x,a,v,particles,l)
            if (mp==4) then
                call postavitev_m(x,l)!posodobi hitrosti na halfstep, zraèuna g-je iterativno, zraèuna nove popravljene položaje
            end if
            call synh_arrays(x)
            call synh_arrays(v)
                !zracunamo sile na novo (z novimi položaji)
            a=0
            LJ_spheres_force=0
            co_potencial_el=0
            co_potencial_lj=0
            call calc_coms(x,particles,l,coms)
            call ewald_sum(x,particles,a,l,co_virial,co_potencial_el,co_virial_tensor,coms)
            call lj_del_pot(x,particles,a,l,co_potencial_lj,co_virial,co_virial_tensor,coms)
            call LJ_sphere_correction(x,particles,a,l,co_potencial_lj,co_virial,co_virial_tensor,coms,LJ_spheres,LJ_spheres_force)

                !nove sile -> nove hitrosti
            sync all
            call CO_SUM(co_potencial_el)
            call CO_SUM(co_potencial_lj)
            call CO_SUM(co_virial_tensor)
            call CO_SUM(a)
            sync all
                !nove sile -> nove hitrosti
            call rattle_hit(x,a,v,particles,l)
            call synh_arrays(v)
            !    !zraèunamo tlak
            call calc_temp(v,particles,temp)
            call calc_pressure(co_virial_tensor,v,l,pressure,particles,x,a,co_press_tensor,temp)
            
        
            call calc_th_eta(th_eta,v,particles)
            call  nh_v_adjust(v,th_eta)
            call calc_th_eta(th_eta,v,particles)
                !zapišemo rezultate
            if (this_image()==1) then
            !if (mod(step,samplingst*10)==0) then
                print*,step,temp*temp_ref, pressure*press_ref/10**5
            !end if
            end if
        end subroutine
!        
        subroutine npt_mtk_step(x,v,a,particles,l,co_virial,co_potencial_el,co_potencial_lj,step,pressure,temp,co_press_tensor,co_virial_tensor,th_eta,mtk_nib,coms,LJ_spheres,LJ_spheres_force)
            real(8),intent(inout)::x(num_images(),mp*N/num_images(),3)[*],v(num_images(),3*n/num_images(),3)[*],a(num_images(),3*n/num_images(),3)[*],particles(mp*n,2),l[*],LJ_spheres(N_LJ_spheres,3),LJ_spheres_force(N_LJ_spheres,3)
            real(8),intent(inout)::co_potencial_el[*],co_potencial_lj[*],pressure[*],temp[*],co_press_tensor(3,3)[*],co_virial_tensor(3,3)[*],coms(num_images(),mp*N/num_images(),3)[*],co_virial[*]
            real(8),intent(inout)::th_eta[*],mtk_nib[*]
            integer,intent(in)::step
            sync all
            !MTK termo/barostat
            call calc_temp(v,particles,temp)
            call calc_pressure(co_virial_tensor,v,l,pressure,particles,x,a,co_press_tensor,temp)
            call mtk_calc_th_eta(th_eta,particles,mtk_nib,temp)
            call calc_mtk_nib(th_eta,l,mtk_nib,pressure,temp)
            call mtk_v_adjust(v,th_eta,mtk_nib)
            call calc_temp(v,particles,temp)
            call calc_pressure(co_virial_tensor,v,l,pressure,particles,x,a,co_press_tensor,temp)
            call calc_mtk_nib(th_eta,l,mtk_nib,pressure,temp)
            call mtk_calc_th_eta(th_eta,particles,mtk_nib,temp)
            call mtk_box(l,mtk_nib)
            
            co_virial_tensor=0
            sync all
            !uporabimo rattle algoritem
            if (mp==4) then
                call shift_coms(x,mtk_nib,particles,l,coms)
            end if
            !
            call rattle_polozaji(x,a,v,particles,l)
            if (mp==4) then
                call postavitev_m(x,l)!posodobi hitrosti na halfstep, zraèuna g-je iterativno, zraèuna nove popravljene položaje
            end if
            call synh_arrays(x)
            call synh_arrays(v)
                !zracunamo sile na novo (z novimi položaji)
            a=0
            LJ_spheres_force=0
            co_potencial_el=0
            co_potencial_lj=0
            call calc_coms(x,particles,l,coms)
            call ewald_sum(x,particles,a,l,co_virial,co_potencial_el,co_virial_tensor,coms)
            call lj_del_pot(x,particles,a,l,co_potencial_lj,co_virial,co_virial_tensor,coms)
            call LJ_sphere_correction(x,particles,a,l,co_potencial_lj,co_virial,co_virial_tensor,coms,LJ_spheres,LJ_spheres_force)

                !nove sile -> nove hitrosti
            sync all
            call CO_SUM(co_potencial_el)
            call CO_SUM(co_potencial_lj)
            call CO_SUM(co_virial_tensor)
            call CO_SUM(a)
            sync all
                !nove sile -> nove hitrosti
            call rattle_hit(x,a,v,particles,l)
            call synh_arrays(v)
            call calc_temp(v,particles,temp)
            call calc_pressure(co_virial_tensor,v,l,pressure,particles,x,a,co_press_tensor,temp)
            !
            call mtk_calc_th_eta(th_eta,particles,mtk_nib,temp)
            call calc_mtk_nib(th_eta,l,mtk_nib,pressure,temp)
            call mtk_v_adjust(v,th_eta,mtk_nib)
            call calc_temp(v,particles,temp)
            call calc_pressure(co_virial_tensor,v,l,pressure,particles,x,a,co_press_tensor,temp)
            call calc_mtk_nib(th_eta,l,mtk_nib,pressure,temp)
            call mtk_calc_th_eta(th_eta,particles,mtk_nib,temp)
        
        
            !zapišemo rezultate
            if (this_image()==1) then
            !if (mod(step,samplingst*10)==0) then
                print*,step,temp*temp_ref, pressure*press_ref/10**5
            !end if
            end if
        end subroutine
!    
    
    
    
    
end module
    
