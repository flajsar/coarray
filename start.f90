module start
    use constants
    use basic_fnc
    use motion
    implicit none
    
    contains
    
    subroutine ref_units(l)
        real(8),intent(inout)::l
        real(8)::lj_parameter(2)=0
    
        ener_ref=1.0/avog_num
        time_ref=sqrt(atomic_mass*angstrom**2/ener_ref)
        temp_ref=ener_ref/k_boltz
        press_ref=ener_ref/(angstrom**3.)
        velocity_ref=sqrt(ener_ref/atomic_mass)
        force_ref=ener_ref/angstrom
        eps0=epsilon0 * angstrom * ener_ref /unit_charge**2
    
        tau=time_ref/TB_strength
        mtk_w=(6*N-3)*(tau/2)**2*10
        th_q=(6*N-3)*tau**2
        G=0.35
        dt=time_step/time_ref
        kb=k_boltz * temp_ref/(ener_ref)
        temp0=temperature0/temp_ref
        press0=pressure0* 10**5 / press_ref
        rdf_max=rdf_parameters(2)
        del_r=rdf_parameters(1)
        n_points=floor(rdf_max/del_r)+1
        count_Rdf=0
        
        switch_d=r_upper**2-r_lower**2
        switch_a=-10/switch_d**3
        switch_b=15/switch_d**4
        switch_c=-6/switch_d**5
    
        pi=atan(1.)*4
        models(1,:)=(/3.0,15.99903,-0.82  ,1.00811,0.41   ,0.0,0.0     ,1.0   ,1.632980862,109.47,0.0   ,649.8,3.1656/) !SPC
        models(2,:)=(/3.0,15.99903,-0.8476,1.00811,0.4238 ,0.0,0.0     ,1.0   ,1.632980862,109.47,0.0   ,649.8,3.1656/) !SPC/E
        models(3,:)=(/4.0,15.99903,0.0    ,1.00811,0.52   ,0.0,-1.04   ,0.9572,1.513900655,104.52,0.15  ,648.5,3.15365/)  !TIP4P
        models(4,:)=(/4.0,15.99903,0.0    ,1.00811,0.52422,0.0,-1.04844,0.9572,1.513900655,104.52,0.125 ,680.946,3.16435/) !TIP4P-Ew
        models(5,:)=(/4.0,15.99903,0.0    ,1.00811,0.5564 ,0.0,-1.1128 ,0.9572,1.513900655,104.52,0.1546,774.3,3.1589/) !TIP4P/2005
        lj_parameter(2)=2.782
        lj_parameter(1)=3.2135*unit_charge/(10**3)*avog_num
        LJ_sphere_parameters(1,:)=lj_parameter

        lj_parameter(2)=3.403
        lj_parameter(1)=119.7*k_boltz*avog_num
        LJ_sphere_parameters(2,:)=lj_parameter
        lj_parameter(2)=3.64
        lj_parameter(1)=168.7*k_boltz*avog_num       
        LJ_sphere_parameters(3,:)=lj_parameter
        lj_parameter(2)=4.07
        lj_parameter(1)=223.6*k_boltz*avog_num           
        LJ_sphere_parameters(4,:)=lj_parameter
        mass_o=models(model,2)
        charge_o=models(model,3)
        mass_h=models(model,4)
        charge_h=models(model,5)
        mass_m=models(model,6)
        charge_m=models(model,7)
        r_oh=models(model,8)
        r_hh=models(model,9)
        angle0=models(model,10)*pi/180
        r_om=models(model,11)
        a1=models(model,12)
        a2=models(model,13)
        
        LJ_atom_sigma=LJ_sphere_parameters(1,1)
        LJ_atom_eps=LJ_sphere_parameters(1,2)
        
 
            
    
    end subroutine
    
    
    subroutine calc_self_int(particles)
        real(8),intent(in)::particles(mp*N,2)
        integer::i
        self_int=0
        do i=1,mp*N
            self_int=self_int+particles(i,2)**2
        end do
        self_int=self_int*g/(4*eps0*pi**(3.0/2.0))
    end subroutine
    
    subroutine calc_lj_integrals(l)
        real(8),intent(in)::l
        lj_integral1=(1/(315*r_lower**3*r_upper**9))*(4*a1*a2**6*(a2**6*(switch_a*(-35*r_lower**9 + 135*r_lower**7*r_upper**2 - 189*r_lower**5*r_upper**4 + 105*r_lower**3*r_upper**6 - 16*r_upper**9) + r_lower**2*(r_lower - r_upper)**5*(-5*switch_c*r_lower*(r_lower - r_upper)*(7*r_lower**4 + 42*r_lower**3*r_upper + 102*r_lower**2*r_upper**2 + 122*r_lower*r_upper**3 + 63*r_upper**4) + switch_b*(35*r_lower**4 + 175*r_lower**3*r_upper + 345*r_lower**2*r_upper**2 + 325*r_lower*r_upper**3 + 128*r_upper**4))) + 3*r_lower**3*(r_lower - r_upper)**4*r_upper**6*(35*switch_a*(r_lower**2 + 4*r_lower*r_upper + r_upper**2) - (r_lower - r_upper)*(7*switch_b*(5*r_lower**3 + 25*r_lower**2*r_upper + 15*r_lower*r_upper**2 + 3*r_upper**3) + 5*switch_c*(-7*r_lower**5 - 35*r_lower**4*r_upper + 24*r_lower**2*r_upper**3 + 15*r_lower*r_upper**4 + 3*r_upper**5)))))
        lj_integral2=(4.0*a1*a2**6*(a2**6 - 3*r_upper**6))/(9.0*r_upper**9)
    end subroutine
    
    subroutine calc_v_mat(v_mat)
        implicit none
        real(8),intent(out)::v_mat(3,3)
        
        v_mat(1,1)=(1/mass_o + 1/mass_h)* r_oh**2
        v_mat(1,2)=(1/mass_o)* r_oh**2 * cos(angle0)
        v_mat(1,3)=-(1/mass_h)* r_oh*r_hh * sin(-angle0/2)
        
        v_mat(2,1)=(1/mass_o)* r_oh**2 * cos(angle0)
        v_mat(2,2)=(1/mass_o + 1/mass_h)* r_oh**2
        v_mat(2,3)=(1/mass_h)* r_oh*r_hh * sin(angle0/2)
        
        v_mat(3,1)=-(1/mass_h)*r_oh*r_hh * sin(-angle0/2)
        v_mat(3,2)=(1/mass_h)* r_oh*r_hh * sin(angle0/2)
        v_mat(3,3)=(1/mass_h + 1/mass_h)* r_hh**2
        
        call invert_33_matrix(v_mat)
        
    end subroutine
    
    subroutine set_particles_values(particles)
        real(8),intent(inout)::particles(mp*N,2)
        integer::i
        
        do i=1,mp*N,mp
            particles(i,1)=mass_o
            particles(i,2)=charge_o
            particles(i+1,1)=mass_h
            particles(i+1,2)=charge_h
            particles(i+2,1)=mass_h
            particles(i+2,2)=charge_h
            if (mp==4) then
                particles(i+3,1)=mass_m
                particles(i+3,2)=charge_m
            end if
        end do
    end subroutine
    
    subroutine stconx_from_file(x,v,a,l,particles,file_unitXCON)
        real(8),intent(in)::l,particles(mp*N,2)
        real(8),intent(out)::x(num_images(),mp*N/num_images(),3)[*],v(num_images(),3*n/num_images(),3)[*],a(num_images(),3*n/num_images(),3)[*]
        integer,intent(in)::file_unitXCON
        real(8)::tempread(3),xx(3),tempx(mp*N,3),tempv(3*N,3)
        integer::i,mol,imag, npi,npiv
        
        npi=mp*N/num_images()
        npiv=3*N/num_images()
        do i=1,mp*N
            read(file_unitXCON,*)tempread
            tempread=(tempread)/1.26+l/2
            tempx(i,:)=periodic_bound( tempread ,l)
        end do
        do i=1,3*N,3
            call random_number(xx)
            tempv(i,1) = dsqrt(temp0/(particles(i,1)+particles(i+1,1)*2))*sqrt(-2*log(xx(1)))*cos(2.*pi*xx(2))
            tempv(i,2) = dsqrt(temp0/(particles(i,1)+particles(i+1,1)*2))*sqrt(-2*log(xx(1)))*sin(2.*pi*xx(2))
            call random_number(xx)
            tempv(i,3) = dsqrt(temp0/(particles(i,1)+particles(i+1,1)*2))*sqrt(-2*log(xx(1)))*cos(2.*pi*xx(2))
            tempv(i+1,:)= tempv(i,:)
            tempv(i+2,:)= tempv(i,:)
        enddo
         
        do i=1,num_images()
            x(i,:,:)=tempx(1+npi*(i-1):i*npi,:)
            v(i,:,:)=tempv(1+npiv*(i-1):i*npiv,:)
            a(i,:,:)=0
        end do
        
    end subroutine
    
    subroutine random_starting_conditions(x,v,particles,a,l)
        real(8),intent(in)::particles(mp*N,2),l[*]
        real(8),intent(out)::x(num_images(),mp*N/num_images(),3)[*],v(num_images(),3*n/num_images(),3)[*],a(num_images(),3*n/num_images(),3)[*]
        real(8)::xx(3),deltax,deltay,r,mol(3,3),gama,vec(3)!,tempx(mp*N,3),tempv(3*N,3)
        logical::covering
        integer::i,j,npi,nvi,imag,imag2
        
        npi=mp*N/num_images()
        nvi=3*N/num_images()
        call random_seed()
        call random_number(xx)
        x(1,1,:)=(xx(:))*l
        do imag=1,num_images()
            do i=mp+1-mp*(1-kronecker(1,imag)),npi,mp !postavitev kisikov
                covering=.true.
                do while (covering)
                    call random_number(xx)
                    xx=(xx)*l
                    covering=.false.
                    do imag2=1,num_images()-imag
                        do j=1,i-mp,mp
                            r=dsqrt(image(x(imag2,j,1),xx(1),l)**2+image(x(imag2,j,2),xx(2),l)**2+image(x(imag2,j,3),xx(3),l)**2)
                            if (r<3) covering=.true.
                        enddo
                    end do
                enddo
                x(imag,i,:)=xx(:)
            
            enddo
        end do 
        !postavitev vodikov
        do imag=1,num_images()
            do i=1,npi,mp 
                call random_number(xx)
                xx(:)=floor(xx(:)*3)
                deltax=r_hh/2
                deltay=r_oh*cos(angle0/2)
                x(imag,i+1,1)=x(imag,i,1)+deltax
                x(imag,i+2,1)=x(imag,i,1)-deltax
                x(imag,i+1,2)=x(imag,i,2)+deltay
                x(imag,i+2,2)=x(imag,i,2)+deltay
                x(imag,i+1,3)=x(imag,i,3)
                x(imag,i+2,3)=x(imag,i,3)
                x(imag,i,:)=periodic_bound(x(imag,i,:),l)
                x(imag,i+1,:)=periodic_bound(x(imag,i+1,:),l)
                x(imag,i+2,:)=periodic_bound(x(imag,i+2,:),l)
            end do
        end do
        if (mp==4) then
            do imag=1,num_images()
                do i=1,npi,mp
                    mol(1,:)=x(imag,i,:)
                    mol(2,:)=x(imag,i+1,:)
                    mol(3,:)=x(imag,i+2,:)
                    mol(2,:)=mol(2,:)+l*int((mol(1,:)-mol(2,:))*2/l)
                    mol(3,:)=mol(3,:)+l*int((mol(1,:)-mol(3,:))*2/l)
                
                    vec=0.5*(mol(2,:)+mol(3,:))-mol(1,:)
                    gama=r_om/(dsqrt(len_sq(vec)))
                    x(imag,i+3,:)=mol(1,:)+gama*vec
                    x(imag,i+3,:)=periodic_bound(x(imag,i+3,:),l)

                end do
            end do
        end if
        do imag=1,num_images()
            do i=1,nvi,3
                call random_number(xx)
                v(imag,i,1) = dsqrt(temp0/(mass_o+mass_h*2))*sqrt(-2*log(xx(1)))*cos(2.*pi*xx(2))
                v(imag,i,2) = dsqrt(temp0/(mass_o+mass_h*2))*sqrt(-2*log(xx(1)))*sin(2.*pi*xx(2))
                call random_number(xx)
                v(imag,i,3) = dsqrt(temp0/(mass_o+mass_h*2))*sqrt(-2*log(xx(1)))*cos(2.*pi*xx(2))
                v(imag,i+1,:)= v(imag,i,:)
                v(imag,i+2,:)= v(imag,i,:)
            end do
        end do
        a=0
 
            
        
    end subroutine
end module