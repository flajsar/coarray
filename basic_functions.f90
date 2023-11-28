module basic_fnc
    use constants
    implicit none
    contains
    !image- izracuna razliko med dvema koordinatama po mininum image konveciji
    !len_sq- izracuna kvadrat dolžine vektorja
    !diff12- izracuna razliko dveh vektorjev z upoštevanjem periodnih mej
    !periodic_bound- preveri če je atom prek meje škatle in to popravi da ga premakne na drugo stran
    !calc_pressure- z virial in hitrostmi zračuna tlak v sistemu po enačbi P=1/(3*V) * (2*K-Virial)
    !calc_rdf_pict- zračuna za posamezno sliko dodatek k  rdf 
    !kronecker(n1,n2) vrne 0 če sta n1 in n2 enaka, drugače 0
    
    
    !function image(x1,x2,l)
    !    implicit none
    !    real(8) :: x1,x2,l,x12,image
    !    x12=x2-x1
    !    x12=x12-l*int(2*x12/l) !če je x12 večji od l/2, se upošteva bližja slika
    !    image=x12
    !end function image
    !
    !
    !function len_sq(vektor)
    !    implicit none
    !    real(8):: vektor(3),len_sq
    !    len_sq=vektor(1)**2+vektor(2)**2+vektor(3)**2
    !end function len_sq
    !
    !
    !function diff12(v1,v2,l)
    !    implicit none
    !    real(8):: v1(3),v2(3),diff12(3),l
    !    diff12(1)=image(v1(1),v2(1),l)
    !    diff12(2)=image(v1(2),v2(2),l)
    !    diff12(3)=image(v1(3),v2(3),l)
    !end function diff12
    !
    !
    !function periodic_bound(v,l)
    !    implicit none
    !    real(8)::periodic_bound(3),v(3),l
    !    v(1)=v(1)-l*int(v(1)/l)-floor(v(1)/(1000*l))*l !del z int preveri če je večje od l in l odšteje
    !    v(2)=v(2)-l*int(v(2)/l)-floor(v(2)/(1000*l))*l !del z floor preveri če je negativno in l prišteje
    !    v(3)=v(3)-l*int(v(3)/l)-floor(v(3)/(1000*l))*l
    !    periodic_bound(:)=v(:)
    !end function
    !
    !function kronecker(i,j)
    !    integer::kronecker,i,j
    !    
    !    kronecker=int((float((i+j)-abs(i-j)))/(float((i+j)+abs(i-j))))
    !end function
    !
    !
    !subroutine calc_switch(r,switch,switch_der)
    !    real(8),intent(in)::r
    !    real(8),intent(out)::switch,switch_der
    !    real(8)::z,fs
    !    
    !    z=ceiling(real(floor(r/r_lower))/real(10d10))*(r**2-r_lower**2)
    !    fs=switch_a*z**2+switch_b*z**3+switch_c*z**4
    !    switch=(1+fs*z)*(1-ceiling(real(floor(r/r_upper))/real(10d10)))
    !    switch_der=2*r*(3*fs+switch_b*z**3+2*switch_c*z**4)*(1-ceiling(real(floor(r/r_upper))/real(10d10)))
    !end subroutine
    !
    !subroutine fill_virial_tensor(virial_tensor,vt11,vt12,vt13,vt22,vt23,vt33)
    !    real(8),intent(in)::vt11,vt12,vt13,vt22,vt23,vt33
    !    real(8),intent(inout)::virial_tensor(3,3)
    !        virial_tensor(1,1)=vt11+virial_tensor(1,1)
    !        virial_tensor(1,2)=vt12+virial_tensor(1,2)
    !        virial_tensor(2,1)=vt12+virial_tensor(2,1)
    !        virial_tensor(3,1)=vt13+virial_tensor(3,1)
    !        virial_tensor(1,3)=vt13+virial_tensor(1,3)
    !        virial_tensor(2,2)=vt22+virial_tensor(2,2)
    !        virial_tensor(3,3)=vt33+virial_tensor(3,3)
    !        virial_tensor(2,3)=vt23+virial_tensor(2,3)
    !        virial_tensor(3,2)=vt23+virial_tensor(3,2)
    !end subroutine
    !
    !subroutine sym_tensor(co_virial_tensor)
    !    real(8),intent(inout)::co_virial_tensor(3,3)
    !    integer::i,j
    !    real(8)::ab,ba
    !    
    !    do i=1,3
    !        do j=i+1,3
    !            ab=co_virial_tensor(i,j)
    !            ba=co_virial_tensor(j,i)
    !            co_virial_tensor(i,j)=(ab+ba)/2
    !            co_virial_tensor(j,i)=(ab+ba)/2
    !        end do
    !    end do
    !end subroutine
    !
    !
    !subroutine calc_coms(x,particles,l,coms)
    !    real(8),intent(in)::x(num_images(),mp*N/num_images(),3),particles(mp*N,2),l
    !    real(8),intent(out):: coms(num_images(),mp*N/num_images(),3)
    !    real(8)::mol(mp,3),cntr_mass(3)
    !    integer::i,j,npi,imag
    !    
    !    imag=this_image()
    !    npi=mp*N/num_images()
    !    do i=1,npi,mp
    !        mol(:,:)=x(imag,i:i+mp-1,:)
    !        do j=2,mp
    !            mol(j,:)=mol(j,:)+l*int((mol(1,:)-mol(j,:))*2/l)
    !        end do
    !        cntr_mass=particles(i+(mp-1)*npi,1)*mol(1,:)+particles(i+1+(mp-1)*npi,1)*mol(2,:)+particles(i+2+(mp-1)*npi,1)*mol(3,:)
    !        cntr_mass=cntr_mass/(particles(i+(mp-1)*npi,1)+particles(i+(mp-1)*npi+1,1)+particles(i+(mp-1)*npi+2,1))
    !        do j=1,mp
    !            coms(imag,i+j-1,:)=cntr_mass
    !        end do
    !        
    !    end do
    !    call synh_arrays(coms)
    !end subroutine
    !
    !subroutine synh_arrays(array)
    !    real(8),intent(inout)::array(:,:,:)
    !    integer::i
    !    
    !    SYNC ALL
    !    do i=1,num_images()
    !        call CO_BROADCAST(array(i,:,:),i)
    !    end do
    !    SYNC ALL
    !end subroutine
    !
    !
    !subroutine get_k_bound(start,ending)
    !    implicit none
    !    integer, intent(out)::start,ending
    !    integer::i=0,base,remainder,image
    !    image=this_image()
    !    
    !    base=(2*k+1)**3/num_images()
    !    remainder=mod((2*k+1)**3,num_images())
    !    if (this_image() .le. remainder) then
    !        start=1+(image-1)*(base+1)
    !        ending=image*(base+1)
    !    else
    !        start=1+(image-1)*(base)+remainder
    !        ending=base*image+remainder
    !    end if
    !end subroutine
    !
    !subroutine get_ks(k1,k2,k3,k_index)
    !    implicit none
    !    integer, intent(in)::k_index
    !    integer, intent(out)::k1,k2,k3
    !    integer::i,allk
    !    
    !    i=k_index-1
    !    allk=2*k+1
    !    k1=i/allk**2-k
    !    k2=mod(i,allk**2)/allk-k
    !    k3=mod(mod(i,allk**2),allk)-k
    !end subroutine
    !    
    !subroutine invert_33_matrix(a)
    !    implicit none
    !    real(8),intent(inout)::a(3,3)
    !    
    !    real(8)::det,a_in(3,3)
    !    
    !    det=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) - a(1,2)*(a(2,1)*a(3,3)-a(3,1)*a(2,3)) + a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
    !    
    !    
    !    a_in(1,1)=(a(2,2)*a(3,3)-a(2,3)*a(3,2))
    !    a_in(1,2)=-(a(1,2)*a(3,3)-a(1,3)*a(3,2))
    !    a_in(1,3)=(a(1,2)*a(2,3)-a(1,3)*a(2,2))
    !    
    !    a_in(2,1)=-(a(2,1)*a(3,3)-a(3,1)*a(2,3))
    !    a_in(2,2)=(a(1,1)*a(3,3)-a(1,3)*a(3,1))
    !    a_in(2,3)=-(a(1,1)*a(2,3)-a(1,3)*a(2,1))
    !    
    !    a_in(3,1)=(a(2,1)*a(3,2)-a(2,2)*a(3,1))
    !    a_in(3,2)=-(a(1,1)*a(3,2)-a(1,2)*a(3,1))
    !    a_in(3,3)=(a(1,1)*a(2,2)-a(1,2)*a(2,1))
    !    
    !    a=a_in/det
    !    
    !    
    !end subroutine        
end module