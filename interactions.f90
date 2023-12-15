
    
    module interactions
    use basic_fnc
    use constants
    implicit none
        !vsebuje subroutines za izračun potenciala, sil in virial
        !lj_del_pot -> zračuna lj del potenciala,sile,virial med kisiki 
        ! ewald_sum-> zračuna elektrostatkse sile kot ewaldovo vsoto
    
    contains
    
    subroutine lj_del_pot(x,particles,a,l,co_potencial_lj,co_virial,co_virial_tensor,coms)
    
        implicit none
        real(8),intent(in)::x(num_images(),mp*N/num_images(),3),particles(mp*N,2),coms(num_images(),mp*N/num_images(),3)
        real(8),intent(in)::l
        real(8),intent(inout)::a(num_images(),3*n/num_images(),3),co_potencial_lj,co_virial,co_virial_tensor(3,3)
        !temp spremenljivke:
        real(8),allocatable::f(:,:,:,:,:),vf(:,:,:)
        real(8)::x1,x2,r,r2,temppot,v1(3),v2(3),mol(mp,3),switch,switch_der,lj_long_range_corr,fcom
        real(8)::vt11,vt12,vt13,vt22,vt23,vt33,potencial_lj_t,virialt,vx12(3),coms_diff(3)
        integer::i,j,j1,cim2,cim,start,ending,npi,image
        allocate(f(num_images(),mp*N/num_images(),num_images(),mp*N/num_images(),3),vf(num_images(),mp*n/num_images(),3))
        f=0
        vf=0
        npi=mp*N/num_images()
        image=this_image()
        !zračunamo LJ potencial in sledeče sile med kisiki
        do i=1,npi,mp
            do cim2=image,image+num_images()/2
                cim=mod(cim2,num_images())+kronecker(cim2,num_images())*num_images()
                start=1+kronecker(cim2,image)*(mp*(i/mp)+mp)+  kronecker(cim2,image+num_images()/2)*((i-1)*ceiling(real(floor(real(cim)/real(image)))/(num_images()+1))+(i-1+mp)*ceiling(real(floor(real(image)/real(cim)))/(num_images()+1)))
                ending=npi+kronecker(cim2,image+num_images()/2)*(-npi+start+npi/2-1)
                do j1=start,ending,mp
                    j=mod(j1,npi)+kronecker(j1,npi)*npi
                    vx12=diff12(x(cim,j,:),x(image,i,:),l)
                    r2=len_sq(vx12)
                    r=dsqrt(r2)
                    call calc_switch(r,switch,switch_der)
                    switch=1
                    switch_der=0
                    x1=(a2/r)**3
                    x1=x1**2
                    temppot=4*a1*x1*(x1-1)
                    x2=-24*a1*x1*(2*x1-1)*switch+temppot*switch_der

                    co_potencial_lj=co_potencial_lj+temppot*switch
                    fcom=-x2/r2
                    f(image,i,cim,j,:)=fcom*vx12
                    f(cim,j,image,i,:)=-f(image,i,cim,j,:)
                    vf(image,i,:)=vf(image,i,:)+fcom*vx12
                    vf(cim,j,:)=vf(cim,j,:)-fcom*vx12
                    coms_diff=diff12(coms(cim,j,:),coms(image,i,:),l)
                    
                    co_virial_tensor(1,1)=co_virial_tensor(1,1)-f(image,i,cim,j,1)*coms_diff(1)
                    co_virial_tensor(1,2)=co_virial_tensor(1,2)-f(image,i,cim,j,1)*coms_diff(2)
                    co_virial_tensor(2,1)=co_virial_tensor(2,1)-f(image,i,cim,j,2)*coms_diff(1)
                    co_virial_tensor(1,3)=co_virial_tensor(1,3)-f(image,i,cim,j,1)*coms_diff(3)
                    co_virial_tensor(3,1)=co_virial_tensor(3,1)-f(image,i,cim,j,3)*coms_diff(1)
                    co_virial_tensor(2,2)=co_virial_tensor(2,2)-f(image,i,cim,j,2)*coms_diff(2)
                    co_virial_tensor(2,3)=co_virial_tensor(2,3)-f(image,i,cim,j,2)*coms_diff(3)
                    co_virial_tensor(3,2)=co_virial_tensor(3,2)-f(image,i,cim,j,3)*coms_diff(2)
                    co_virial_tensor(3,3)=co_virial_tensor(3,3)-f(image,i,cim,j,3)*coms_diff(3)
                end do
            end do
        end do
        if (mp==4) then
            call project_realf(vf,a,x,l,coms) !preračuna sile iz 'dummy' atoma na realne
        else
            a(:,:,:)=a(:,:,:)+vf(:,:,:)
        end if
        
    end subroutine
    !
    !subroutine calc_lj_corr(l,lj_long_range_corr)
    !    real(8),intent(in)::l
    !    real(8),intent(out)::lj_long_range_corr
    !    integer::i
    !    lj_long_range_corr=0
    !    lj_long_range_corr=(lj_integral1+lj_integral2)*2*pi*N**2/l**3
    !end subroutine
    !    
        
        
        
    subroutine ewald_sum(x,particles,a,l,co_virial,co_potencial_el,co_virial_tensor,coms)
    
        implicit none
        real(8),intent(in)::x(num_images(),mp*N/num_images(),3),particles(mp*N,2),coms(num_images(),mp*N/num_images(),3)
        real(8),intent(in)::l
        real(8),intent(inout)::a(num_images(),3*n/num_images(),3),co_potencial_el,co_virial,co_virial_tensor(3,3)
        !temp spremenljivke:
        real(8),allocatable::f(:,:,:,:,:),f_rec(:,:,:),vf(:,:,:)
        real(8)::x1,x2,r,r2,sqpi,fakt1,fakt2,temppot,potencial_el_t,vt11,vt12,vt13,vt22,vt23,vt33
        integer::i,j,j2,k1,k2,k3,image,j1
        real(8)::vec_k(3),rec_v(3,3),d_i(4,3),force(3),mol(4,3),vs_re(mp*N),vs_im(mp*N),fcom,vx12(3),coms_diff(3),switch,switch_der
        complex(8)::vs1,vs2,im,y1,t,vs1t
        integer::npi,cim,cim2,start,ending,k_index
        allocate(f(num_images(),mp*N/num_images(),num_images(),mp*N/num_images(),3),f_rec(num_images(),mp*n/num_images(),3),vf(num_images(),mp*n/num_images(),3))
        f=0
        vf=0
        image=this_image()
        npi=mp*N/num_images()
        
        sqpi=dsqrt(pi)
        co_potencial_el=0
        co_virial_tensor=0
        !zračunamo short-range del elektrostatskega potenciala za vse atome
        do i=1,npi
            do cim2=image,image+num_images()/2
                cim=mod(cim2,num_images())+kronecker(cim2,num_images())*num_images()
                
                !start=1+kronecker(cim2,image)*(i)+kronecker(cim2,image+num_images()/2)*ceiling(real(floor(real(image)/real(cim)))/(num_images()+1))*npi/2
                start=1+kronecker(cim2,image)*(i)+kronecker(cim2,image+num_images()/2)*((i-1)*ceiling(real(floor(real(cim)/real(image)))/(num_images()+1))+(i)*ceiling(real(floor(real(image)/real(cim)))/(num_images()+1)))
                
                !ending=npi-kronecker(cim2,image+num_images()/2)*ceiling(real(floor(real(cim)/real(image)))/(num_images()+1))*npi/2
                ending=npi+kronecker(cim2,image+num_images()/2)*(-npi+start+npi/2-1)
                do j1=start,ending
                    j=mod(j1,npi)+kronecker(j1,npi)*npi
                    vx12=(diff12(x(cim,j,:),x(image,i,:),l))
                    r2=len_sq(vx12) !kvadrat razdalje med atomoma
                    r=dsqrt(r2)
                    coms_diff=diff12(coms(cim,j,:),coms(image,i,:),l)
                    call calc_switch(dsqrt(len_sq(coms_diff)),switch,switch_der)
                    switch=1
                    switch_der=0
                    fakt2=kronecker(cim,image)*ceiling(real(floor(real(mp-4*kronecker(mod(i,mp)+1,1)-mod(i,mp))/real(abs(i-j)))/real((npi+1))))!fakt1 in fakt2 sta uporabljena zato da določita če sta atoma del iste molekule
                    fakt1=((fakt2+1.0)/2.0)**(-1)-1 !če sta del iste molekule: fakt1=0, fakt2=1 in obratno
                    x1=fakt1*erfc(G*r)-fakt2*erf(G*r)
                    x2=particles(i+(image-1)*npi,2)*particles(j+(cim-1)*npi,2) !naboj*naboj
                    temppot=x2*x1/(r*4.0*pi*eps0)*switch
                    co_potencial_el=temppot+co_potencial_el
                    fcom=switch*((temppot+(4.0*pi*eps0)**(-1)*  2 * G/sqpi * x2 * exp(-G**2 * r2))/r2)+temppot*switch_der
                    
                    
                    co_virial_tensor(1,1)=co_virial_tensor(1,1)-fcom*coms_diff(1)*vx12(1)*fakt1
                    co_virial_tensor(3,3)=co_virial_tensor(3,3)-fcom*coms_diff(3)*vx12(3)*fakt1
                    co_virial_tensor(2,2)=co_virial_tensor(2,2)-fcom*coms_diff(2)*vx12(2)*fakt1
                    co_virial_tensor(1,2)=co_virial_tensor(1,2)-fcom*(coms_diff(1)*vx12(2)+coms_diff(2)*vx12(1))*fakt1/2
                    co_virial_tensor(1,3)=co_virial_tensor(1,3)-fcom*(coms_diff(1)*vx12(3)+coms_diff(3)*vx12(1))*fakt1/2
                    co_virial_tensor(2,3)=co_virial_tensor(2,3)-fcom*(coms_diff(3)*vx12(2)+coms_diff(2)*vx12(3))*fakt1/2
                    
                    f(image,i,cim,j,:)=vx12*fcom
                    vf(image,i,:)=vf(image,i,:)+vx12*fcom
                    f(cim,j,image,i,:)=-f(image,i,cim,j,:)
                    vf(cim,j,:)=vf(cim,j,:)-vx12*fcom
                end do
            end do
        end do
        co_virial_tensor(2,1)=co_virial_tensor(1,2)
        co_virial_tensor(3,1)=co_virial_tensor(1,3)
        co_virial_tensor(3,2)=co_virial_tensor(2,3)
        
        !!recipročni del potenciala
        
        f_rec=0
        fakt2=1/(2*eps0*l**3)
        temppot=0
        vs1=(0.0,0.0)
        vs2=vs1
        im=(0.0,1.0)
        rec_v(1,:)=(/1.d0,0.d0,0.d0/)
        rec_v(2,:)=(/0.d0,1.d0,0.d0/)
        rec_v(3,:)=(/0.d0,0.d0,1.d0/)
        rec_v(1,:)=rec_v(1,:)*2*pi/l
        rec_v(2,:)=rec_v(2,:)*2*pi/l
        rec_v(3,:)=rec_v(3,:)*2*pi/l
        
        call get_k_bound(start,ending)
        do k_index=start,ending
            call get_ks(k1,k2,k3,k_index)
            if ((k1 .eq. 0) .and. (k2 .eq. 0) .and. (k3 .eq. 0)) then
            else
                vec_k=k1*rec_v(1,:)+k2*rec_v(2,:)+k3*rec_v(3,:)
                vs1=(0.0,0.0)
                vs2=vs1
                fakt1=exp(-len_sq(vec_k)/(4*g**2))/len_sq(vec_k)
                vs1t=0
                do i=1,npi
                    do j=1,num_images()
                        vs1t=particles(i+(j-1)*npi,2)*exp(im*dot_product(vec_k,x(j,i,:)))+vs1t
                    end do
                end do
                vs2=DCONJG(vs1t)
                t=vs1t*vs2
                temppot=temppot+fakt1*t%re
                do i=1,npi
                    do j=1,num_images()
                        y1=cdexp(im*dot_product(vec_k,x(j,i,:))) * vs2
                        force=vec_k(:)*(particles(i+(j-1)*npi,2)/(eps0*l**3) * fakt1 * y1%IM)
                        vf(j,i,:)=vf(j,i,:)+force
                        f_rec(j,i,:)=f_rec(j,i,:)+force
                    end do
                end do
                do j=1,3
                    do i=1,3
                        co_virial_tensor(j,i)=co_virial_tensor(j,i)-fakt1*fakt2*(kronecker(i,j)-vec_k(i)*vec_k(j)*(2/len_sq(vec_k) + 1/(2*g**2)) )*t%re
                    end do
                end do
                        
                end if
        end do
        temppot=temppot*fakt2
        co_potencial_el=co_potencial_el+temppot
        do image=1,num_images()
            do i=1,npi
                do j=1,3
                    do j2=1,3
                        vx12=diff12(coms(image,i,:),x(image,i,:),l)
                        co_virial_tensor(j,j2)=co_virial_tensor(j,j2)+f_rec(image,i,j)*vx12(j2)
                    end do
                end do
            end do
        end do
        image=this_image()
        
        !prištejemo self interaction potencial
        co_potencial_el=co_potencial_el-self_int/num_images()
        if (mp==4) then
            call project_realf(vf,a,x,l,coms) !preračuna sile iz 'dummy' atoma na realne
        else
            a(:,:,:)=a(:,:,:)+vf(:,:,:)
        end if
        
        end subroutine
        
        subroutine project_realf(vf,a,x,l,coms)
            implicit none
            real(8),intent(in)::vf(num_images(),mp*n/num_images(),3),x(num_images(),mp*N/num_images(),3),coms(num_images(),mp*N/num_images(),3),l
            real(8),intent(inout)::a(num_images(),3*n/num_images(),3)
            real(8)::gama,mol(4,3),f1(3),vec(3),vec1(3),fm(3)
            integer::i,imag,npi
            npi=mp*N/num_images()
            do imag=1,num_images()
                do i=1,npi,mp
                    a(imag,i-i/4,:)=a(imag,i-i/4,:)+vf(imag,i,:)
                    a(imag,i+1-(i+1)/4,:)=a(imag,i+1-(i+1)/4,:)+vf(imag,i+1,:)
                    a(imag,i+2-(i+2)/4,:)=a(imag,i+2-(i+2)/4,:)+vf(imag,i+2,:)
                    mol(1,:)=x(imag,i,:)
                    mol(2,:)=x(imag,i+1,:)
                    mol(3,:)=x(imag,i+2,:)
                    mol(4,:)=x(imag,i+3,:)
                    mol(2,:)=mol(2,:)+l*int((mol(1,:)-mol(2,:))*2/l)
                    mol(3,:)=mol(3,:)+l*int((mol(1,:)-mol(3,:))*2/l)
                    mol(4,:)=mol(4,:)+l*int((mol(1,:)-mol(4,:))*2/l)
                    
                    vec=0.5*(mol(2,:)+mol(3,:))-mol(1,:)
                    gama=r_om/(dsqrt(len_sq(vec)))
                    vec1=mol(4,:)-mol(1,:)
                    f1=vec1*(dot_product(vf(imag,i+3,:),vec1)/dot_product(vec1,vec1))
                    fm=vf(imag,i+3,:)
                    a(imag,i-i/4,:)=a(imag,i-i/4,:)+fm-gama*(fm-f1)
                    a(imag,i+1-(i+1)/4,:)=a(imag,i+1-(i+1)/4,:)+gama*(fm-f1)*0.5
                    a(imag,i+2-(i+2)/4,:)=a(imag,i+2-(i+2)/4,:)+gama*(fm-f1)*0.5
                end do
            end do
        end subroutine



        subroutine LJ_sphere_correction(x,particles,a,l,co_potencial_lj,co_virial,co_virial_tensor,coms,LJ_spheres,LJ_spheres_force)
            implicit none
            real(8),intent(in)::x(num_images(),mp*N/num_images(),3),particles(mp*N,2),coms(num_images(),mp*N/num_images(),3),LJ_spheres(N_LJ_spheres,3)
            real(8),intent(in)::l
            real(8),intent(inout)::a(num_images(),3*n/num_images(),3),co_potencial_lj,co_virial,co_virial_tensor(3,3),LJ_spheres_force(N_LJ_spheres,3)
            !temp spremenljivke:
            real(8),allocatable::vf(:,:)
            real(8)::x1,x2,r,r2,temppot,v1(3),v2(3),mol(mp,3),switch,switch_der,lj_long_range_corr,fcom,force(3)
            real(8)::vt11,vt12,vt13,vt22,vt23,vt33,potencial_lj_t,virialt,vx12(3),coms_diff(3),a1s,a2s
            integer::i,j,j1,cim2,cim,start,ending,npi,image
            allocate(vf(n/num_images(),3))
            vf=0
            npi=mp*N/num_images()
            image=this_image()
            !zračunamo LJ potencial in sledeče sile med kisiki
            a1s=sqrt(a1*LJ_sphere_parameters(lj_atom,1))
            a2s=(a2+LJ_sphere_parameters(lj_atom,2))/2            
            do j=1,N_LJ_spheres
                do i=1,npi,mp
                    vx12=diff12(LJ_spheres(j,:),x(image,i,:),l)
                    r2=len_sq(vx12)
                    r=dsqrt(r2)
                    x1=(a2s/r)**3
                    x1=x1**2
                    temppot=4*a1s*x1*(x1-1)
                    x2=-24*a1s*x1*(2*x1-1)
                    co_potencial_lj=co_potencial_lj+temppot
                    fcom=-x2/r2
                    force=fcom*vx12
                    vf(int(i/mp)+1,:)=vf(int(i/mp)+1,:)+force
                    LJ_spheres_force(j,:)=-force+LJ_spheres_force(j,:)
                    coms_diff=diff12(LJ_spheres(j,:),coms(image,i,:),l)                        
                    co_virial_tensor(1,1)=co_virial_tensor(1,1)-force(1)*coms_diff(1)
                    co_virial_tensor(1,2)=co_virial_tensor(1,2)-force(1)*coms_diff(2)
                    co_virial_tensor(2,1)=co_virial_tensor(2,1)-force(2)*coms_diff(1)
                    co_virial_tensor(1,3)=co_virial_tensor(1,3)-force(1)*coms_diff(3)
                    co_virial_tensor(3,1)=co_virial_tensor(3,1)-force(3)*coms_diff(1)
                    co_virial_tensor(2,2)=co_virial_tensor(2,2)-force(2)*coms_diff(2)
                    co_virial_tensor(2,3)=co_virial_tensor(2,3)-force(2)*coms_diff(3)
                    co_virial_tensor(3,2)=co_virial_tensor(3,2)-force(3)*coms_diff(2)
                    co_virial_tensor(3,3)=co_virial_tensor(3,3)-force(3)*coms_diff(3)
                end do
            end do
            if (mp==4) then
                do i=1,3*n/num_images(),3
                    a(image,i,:)=vf(int(i/3)+1,:)+a(image,i,:)
                end do
            else
                a(image,:,:)=a(image,:,:)+vf(:,:)
            end if
        end subroutine
end module
    