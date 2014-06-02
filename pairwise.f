  subroutine mp_dir_dot(lmax, r_inv)
    integer, intent(in) :: lmax
    real*8, intent(in) :: r_inv !1/R
    !
    integer :: i
    real*8 :: eax, eay, eaz, ebx, eby, fz
    !
    e_pair = qta(1)*vab(1)
    fz  = qta(1)*vabr(1)
    eax = 0d0
    eay = 0d0
    eaz = 0d0
    ebx = 0d0
    eby = 0d0
    do i=2,(lmax+1)*(lmax+1)
       e_pair = e_pair + qta(i)*vab(i)
       fz  = fz + qta(i)*vabr(i)
       eax = eax + qtax(i)*vab(i)
       eay = eay + qtay(i)*vab(i)
       eaz = eaz + qtaz(i)*vab(i)
       ebx = ebx + qtbx(i)*vba(i)
       eby = eby + qtby(i)*vba(i)
       !Ebz = -Eax
    enddo
    f_pair = (/ (eay+eby)*r_inv, -(eax+ebx)*r_inv, -fz/)
    qa_pair = (/ -eax, -eay, -eaz/)
    qb_pair = (/ -ebx, -eby, eaz/)
    !DEBUGPRINT
    write(*,*)
    write(*,*) "====================================================================================================================================="
    write(*,*)
    write(*,*) "Computing energy, forces and torques in the QI frame for this pair"
    write(*,'(A15F16.10)') "Energy = ", e_pair
    write(*,'(A15F16.10)') "E_AB^Ax = ", Eax
    write(*,'(A15F16.10)') "E_AB^Ay = ", Eay
    write(*,'(A15F16.10)') "E_AB^Az = ", Eaz
    write(*,'(A15F16.10)') "E_AB^Bx = ", Ebx
    write(*,'(A15F16.10)') "E_AB^By = ", Eby
    write(*,'(A15F16.10)') "Fz = ", fz
    write(*,*)
    write(*,'(A15(3F16.10))') "QI force = ", f_pair
    write(*,'(A15(3F16.10))') "QI torque A = ", qa_pair
    write(*,'(A15(3F16.10))') "QI torque B = ", qb_pair
    write(*,*)
    write(*,*) "====================================================================================================================================="
    write(*,*)
  end subroutine mp_dir_dot


  subroutine mp_dir_accumulate(atom_a, atom_b, e, f, q)
    !
    integer, intent(in) :: atom_a, atom_b
    real*8, intent(in out) :: e, f(3,*), q(3,*)
    !
    integer :: xyz
    real*8 :: f_pair_lab(3), qa_pair_lab(3), qb_pair_lab(3), tempf, tempa, tempb
    ! rotate from QI to lab
    tempf = f_pair(1)
    tempa = qa_pair(1)
    tempb = qb_pair(1)
    f_pair_lab(1) = tempf*Uz(1)
    qa_pair_lab(1) = tempa*Uz(1)
    qb_pair_lab(1) = tempb*Uz(1)
    f_pair_lab(2) = tempf*Uz(2)
    qa_pair_lab(2) = tempa*Uz(2)
    qb_pair_lab(2) = tempb*Uz(2)
    f_pair_lab(3) = tempf*Uz(3)
    qa_pair_lab(3) = tempa*Uz(3)
    qb_pair_lab(3) = tempb*Uz(3)
    tempf = f_pair(2)
    tempa = qa_pair(2)
    tempb = qb_pair(2)
    f_pair_lab(1) = f_pair_lab(1) + tempf*Uz(4)
    qa_pair_lab(1) = qa_pair_lab(1) + tempa*Uz(4)
    qb_pair_lab(1) = qb_pair_lab(1) + tempb*Uz(4)
    f_pair_lab(2) = f_pair_lab(2) + tempf*Uz(5)
    qa_pair_lab(2) = qa_pair_lab(2) + tempa*Uz(5)
    qb_pair_lab(2) = qb_pair_lab(2) + tempb*Uz(5)
    f_pair_lab(3) = f_pair_lab(3) + tempf*Uz(6)
    qa_pair_lab(3) = qa_pair_lab(3) + tempa*Uz(6)
    qb_pair_lab(3) = qb_pair_lab(3) + tempb*Uz(6)
    tempf = f_pair(3)
    tempa = qa_pair(3)
    tempb = qb_pair(3)
    f_pair_lab(1) = f_pair_lab(1) + tempf*Uz(7)
    qa_pair_lab(1) = qa_pair_lab(1) + tempa*Uz(7)
    qb_pair_lab(1) = qb_pair_lab(1) + tempb*Uz(7)
    f_pair_lab(2) = f_pair_lab(2) + tempf*Uz(8)
    qa_pair_lab(2) = qa_pair_lab(2) + tempa*Uz(8)
    qb_pair_lab(2) = qb_pair_lab(2) + tempb*Uz(8)
    f_pair_lab(3) = f_pair_lab(3) + tempf*Uz(9)
    qa_pair_lab(3) = qa_pair_lab(3) + tempa*Uz(9)
    qb_pair_lab(3) = qb_pair_lab(3) + tempb*Uz(9)
    ! accumulate
    e = e + e_pair
    do xyz=1,3
      f(xyz,atom_a) = f(xyz,atom_a) + f_pair_lab(xyz)
      f(xyz,atom_b) = f(xyz,atom_b) - f_pair_lab(xyz)
      q(xyz,atom_a) = q(xyz,atom_a) + qa_pair_lab(xyz)
      q(xyz,atom_b) = q(xyz,atom_b) + qb_pair_lab(xyz)
    enddo
    !DEBUGPRINT
    write(*,*)
    write(*,*) "====================================================================================================================================="
    write(*,*)
    write(*,*) "Rotating forces and torques for this pair back to the lab frame, using the Uab matrix printed above"
    write(*,'(A15(3F16.10))') "lab force = ", f_pair_lab
    write(*,'(A15(3F16.10))') "lab torque A = ", qa_pair_lab
    write(*,'(A15(3F16.10))') "lab torque B = ", qb_pair_lab
    write(*,*) "====================================================================================================================================="
    return
  end subroutine mp_dir_accumulate


  !> Forms the QI frame multipoles and torque intermediates.
  subroutine form_qi(atom_a, atom_b, lmax)
    use mp_rot, only: mp_rot_form_d
    use psf, only: mp_lab_q
    !
    integer, intent(in) :: atom_a, atom_b, lmax
    !
    real*8 :: d1(9), d2(25), d3(49), d4(81), temp
    integer :: i, j, ij, kk
    ! L = 0
    qta(1) = mp_lab_q(1, atom_a)
    qtb(1) = mp_lab_q(1, atom_b)
    if (lmax .lt. 1) return
    ! L = 1
    ! Uz is the Uab rotation matrix that orients the internuclear vector
    d1(5) = Uz(1)
    d1(8) = Uz(2)
    d1(2) = Uz(3)
    d1(6) = Uz(4)
    d1(9) = Uz(5)
    d1(3) = Uz(6)
    d1(4) = Uz(7)
    d1(7) = Uz(8)
    d1(1) = Uz(9)
    temp = mp_lab_q(2, atom_a)
    qta(2) = temp*d1(1)
    qta(3) = temp*d1(2)
    qta(4) = temp*d1(3)
    temp = mp_lab_q(3, atom_a)
    qta(2) = qta(2) + temp*d1(4)
    qta(3) = qta(3) + temp*d1(5)
    qta(4) = qta(4) + temp*d1(6)
    temp = mp_lab_q(4, atom_a)
    qta(2) = qta(2) + temp*d1(7)
    qta(3) = qta(3) + temp*d1(8)
    qta(4) = qta(4) + temp*d1(9)
    qtax(2) = qta(4)
    qtax(3) = 0d0
    qtax(4) = -qta(2)
    qtay(2) = -qta(3)
    qtay(3) = qta(2)
    qtay(4) = 0d0
    qtaz(2) = 0d0
    qtaz(3) = -qta(4)
    qtaz(4) = qta(3)
    temp = mp_lab_q(2, atom_b)
    qtb(2) = temp*d1(1)
    qtb(3) = temp*d1(2)
    qtb(4) = temp*d1(3)
    temp = mp_lab_q(3, atom_b)
    qtb(2) = qtb(2) + temp*d1(4)
    qtb(3) = qtb(3) + temp*d1(5)
    qtb(4) = qtb(4) + temp*d1(6)
    temp = mp_lab_q(4, atom_b)
    qtb(2) = qtb(2) + temp*d1(7)
    qtb(3) = qtb(3) + temp*d1(8)
    qtb(4) = qtb(4) + temp*d1(9)
    qtbx(2) = qtb(4)
    qtbx(3) = 0d0
    qtbx(4) = -qtb(2)
    qtby(2) = -qtb(3)
    qtby(3) = qtb(2)
    qtby(4) = 0d0
    if (lmax .lt. 2) return
    ! L = 2
    call mp_rot_form_d(lmax, d1, d2, d3, d4)
    temp = mp_lab_q(5, atom_a)
    qta(5) = temp*d2(1)
    qta(6) = temp*d2(2)
    qta(7) = temp*d2(3)
    qta(8) = temp*d2(4)
    qta(9) = temp*d2(5)
    temp = mp_lab_q(6, atom_a)
    qta(5) = qta(5) + temp*d2(6)
    qta(6) = qta(6) + temp*d2(7)
    qta(7) = qta(7) + temp*d2(8)
    qta(8) = qta(8) + temp*d2(9)
    qta(9) = qta(9) + temp*d2(10)
    temp = mp_lab_q(7, atom_a)
    qta(5) = qta(5) + temp*d2(11)
    qta(6) = qta(6) + temp*d2(12)
    qta(7) = qta(7) + temp*d2(13)
    qta(8) = qta(8) + temp*d2(14)
    qta(9) = qta(9) + temp*d2(15)
    temp = mp_lab_q(8, atom_a)
    qta(5) = qta(5) + temp*d2(16)
    qta(6) = qta(6) + temp*d2(17)
    qta(7) = qta(7) + temp*d2(18)
    qta(8) = qta(8) + temp*d2(19)
    qta(9) = qta(9) + temp*d2(20)
    temp = mp_lab_q(9, atom_a)
    qta(5) = qta(5) + temp*d2(21)
    qta(6) = qta(6) + temp*d2(22)
    qta(7) = qta(7) + temp*d2(23)
    qta(8) = qta(8) + temp*d2(24)
    qta(9) = qta(9) + temp*d2(25)
    qtax(5) = 1.7320508075688772d0*qta(7)
    qtax(6) = qta(9)
    qtax(7) = -1.7320508075688772d0*qta(5) - qta(8)
    qtax(8) = qta(7)
    qtax(9) = -qta(6)
    qtay(5) = -1.7320508075688772d0*qta(6)
    qtay(6) = 1.7320508075688772d0*qta(5) - qta(8)
    qtay(7) = -qta(9)
    qtay(8) = qta(6)
    qtay(9) = qta(7)
    qtaz(5) = 0.d0
    qtaz(6) = -qta(7)
    qtaz(7) = qta(6)
    qtaz(8) = -2.d0*qta(9)
    qtaz(9) = 2.d0*qta(8)
    temp = mp_lab_q(5, atom_b)
    qtb(5) = temp*d2(1)
    qtb(6) = temp*d2(2)
    qtb(7) = temp*d2(3)
    qtb(8) = temp*d2(4)
    qtb(9) = temp*d2(5)
    temp = mp_lab_q(6, atom_b)
    qtb(5) = qtb(5) + temp*d2(6)
    qtb(6) = qtb(6) + temp*d2(7)
    qtb(7) = qtb(7) + temp*d2(8)
    qtb(8) = qtb(8) + temp*d2(9)
    qtb(9) = qtb(9) + temp*d2(10)
    temp = mp_lab_q(7, atom_b)
    qtb(5) = qtb(5) + temp*d2(11)
    qtb(6) = qtb(6) + temp*d2(12)
    qtb(7) = qtb(7) + temp*d2(13)
    qtb(8) = qtb(8) + temp*d2(14)
    qtb(9) = qtb(9) + temp*d2(15)
    temp = mp_lab_q(8, atom_b)
    qtb(5) = qtb(5) + temp*d2(16)
    qtb(6) = qtb(6) + temp*d2(17)
    qtb(7) = qtb(7) + temp*d2(18)
    qtb(8) = qtb(8) + temp*d2(19)
    qtb(9) = qtb(9) + temp*d2(20)
    temp = mp_lab_q(9, atom_b)
    qtb(5) = qtb(5) + temp*d2(21)
    qtb(6) = qtb(6) + temp*d2(22)
    qtb(7) = qtb(7) + temp*d2(23)
    qtb(8) = qtb(8) + temp*d2(24)
    qtb(9) = qtb(9) + temp*d2(25)
    qtbx(5) = 1.7320508075688772d0*qtb(7)
    qtbx(6) = qtb(9)
    qtbx(7) = -1.7320508075688772d0*qtb(5) - qtb(8)
    qtbx(8) = qtb(7)
    qtbx(9) = -qtb(6)
    qtby(5) = -1.7320508075688772d0*qtb(6)
    qtby(6) = 1.7320508075688772d0*qtb(5) - qtb(8)
    qtby(7) = -qtb(9)
    qtby(8) = qtb(6)
    qtby(9) = qtb(7)
    if (lmax .lt. 3) return
    ! L = 3
    !call DGEMV('n', 7, 7, ONE, d3, 7, mp_lab_Q(10,atom_a), 1, ZERO, qta(10), 1)
    temp = mp_lab_q(10, atom_a)
    qta(10) = temp*d3(1)
    qta(11) = temp*d3(2)
    qta(12) = temp*d3(3)
    qta(13) = temp*d3(4)
    qta(14) = temp*d3(5)
    qta(15) = temp*d3(6)
    qta(16) = temp*d3(7)
    temp = mp_lab_q(11, atom_a)
    qta(10) = qta(10) + temp*d3(8)
    qta(11) = qta(11) + temp*d3(9)
    qta(12) = qta(12) + temp*d3(10)
    qta(13) = qta(13) + temp*d3(11)
    qta(14) = qta(14) + temp*d3(12)
    qta(15) = qta(15) + temp*d3(13)
    qta(16) = qta(16) + temp*d3(14)
    temp = mp_lab_q(12, atom_a)
    qta(10) = qta(10) + temp*d3(15)
    qta(11) = qta(11) + temp*d3(16)
    qta(12) = qta(12) + temp*d3(17)
    qta(13) = qta(13) + temp*d3(18)
    qta(14) = qta(14) + temp*d3(19)
    qta(15) = qta(15) + temp*d3(20)
    qta(16) = qta(16) + temp*d3(21)
    temp = mp_lab_q(13, atom_a)
    qta(10) = qta(10) + temp*d3(22)
    qta(11) = qta(11) + temp*d3(23)
    qta(12) = qta(12) + temp*d3(24)
    qta(13) = qta(13) + temp*d3(25)
    qta(14) = qta(14) + temp*d3(26)
    qta(15) = qta(15) + temp*d3(27)
    qta(16) = qta(16) + temp*d3(28)
    temp = mp_lab_q(14, atom_a)
    qta(10) = qta(10) + temp*d3(29)
    qta(11) = qta(11) + temp*d3(30)
    qta(12) = qta(12) + temp*d3(31)
    qta(13) = qta(13) + temp*d3(32)
    qta(14) = qta(14) + temp*d3(33)
    qta(15) = qta(15) + temp*d3(34)
    qta(16) = qta(16) + temp*d3(35)
    temp = mp_lab_q(15, atom_a)
    qta(10) = qta(10) + temp*d3(36)
    qta(11) = qta(11) + temp*d3(37)
    qta(12) = qta(12) + temp*d3(38)
    qta(13) = qta(13) + temp*d3(39)
    qta(14) = qta(14) + temp*d3(40)
    qta(15) = qta(15) + temp*d3(41)
    qta(16) = qta(16) + temp*d3(42)
    temp = mp_lab_q(16, atom_a)
    qta(10) = qta(10) + temp*d3(43)
    qta(11) = qta(11) + temp*d3(44)
    qta(12) = qta(12) + temp*d3(45)
    qta(13) = qta(13) + temp*d3(46)
    qta(14) = qta(14) + temp*d3(47)
    qta(15) = qta(15) + temp*d3(48)
    qta(16) = qta(16) + temp*d3(49)
    qtax(10) = 2.449489742783178d0*qta(12)
    qtax(11) = 1.5811388300841898d0*qta(14)
    qtax(12) = -2.449489742783178d0*qta(10) - 1.5811388300841898d0*qta(13)
    qtax(13) = 1.5811388300841898d0*qta(12) + 1.224744871391589d0*qta(16)
    qtax(14) = -1.5811388300841898d0*qta(11) - 1.224744871391589d0*qta(15)
    qtax(15) = 1.224744871391589d0*qta(14)
    qtax(16) = -1.224744871391589d0*qta(13)
    qtay(10) = -2.449489742783178d0*qta(11)
    qtay(11) = 2.449489742783178d0*qta(10) - 1.5811388300841898d0*qta(13)
    qtay(12) = -1.5811388300841898d0*qta(14)
    qtay(13) = 1.5811388300841898d0*qta(11) - 1.224744871391589d0*qta(15)
    qtay(14) = 1.5811388300841898d0*qta(12) - 1.224744871391589d0*qta(16)
    qtay(15) = 1.224744871391589d0*qta(13)
    qtay(16) = 1.224744871391589d0*qta(14)
    qtaz(10) = 0.d0
    qtaz(11) = -qta(12)
    qtaz(12) = qta(11)
    qtaz(13) = -2.d0*qta(14)
    qtaz(14) = 2.d0*qta(13)
    qtaz(15) = -3.d0*qta(16)
    qtaz(16) = 3.d0*qta(15)
    temp = mp_lab_q(10, atom_b)
    qtb(10) = temp*d3(1)
    qtb(11) = temp*d3(2)
    qtb(12) = temp*d3(3)
    qtb(13) = temp*d3(4)
    qtb(14) = temp*d3(5)
    qtb(15) = temp*d3(6)
    qtb(16) = temp*d3(7)
    temp = mp_lab_q(11, atom_b)
    qtb(10) = qtb(10) + temp*d3(8)
    qtb(11) = qtb(11) + temp*d3(9)
    qtb(12) = qtb(12) + temp*d3(10)
    qtb(13) = qtb(13) + temp*d3(11)
    qtb(14) = qtb(14) + temp*d3(12)
    qtb(15) = qtb(15) + temp*d3(13)
    qtb(16) = qtb(16) + temp*d3(14)
    temp = mp_lab_q(12, atom_b)
    qtb(10) = qtb(10) + temp*d3(15)
    qtb(11) = qtb(11) + temp*d3(16)
    qtb(12) = qtb(12) + temp*d3(17)
    qtb(13) = qtb(13) + temp*d3(18)
    qtb(14) = qtb(14) + temp*d3(19)
    qtb(15) = qtb(15) + temp*d3(20)
    qtb(16) = qtb(16) + temp*d3(21)
    temp = mp_lab_q(13, atom_b)
    qtb(10) = qtb(10) + temp*d3(22)
    qtb(11) = qtb(11) + temp*d3(23)
    qtb(12) = qtb(12) + temp*d3(24)
    qtb(13) = qtb(13) + temp*d3(25)
    qtb(14) = qtb(14) + temp*d3(26)
    qtb(15) = qtb(15) + temp*d3(27)
    qtb(16) = qtb(16) + temp*d3(28)
    temp = mp_lab_q(14, atom_b)
    qtb(10) = qtb(10) + temp*d3(29)
    qtb(11) = qtb(11) + temp*d3(30)
    qtb(12) = qtb(12) + temp*d3(31)
    qtb(13) = qtb(13) + temp*d3(32)
    qtb(14) = qtb(14) + temp*d3(33)
    qtb(15) = qtb(15) + temp*d3(34)
    qtb(16) = qtb(16) + temp*d3(35)
    temp = mp_lab_q(15, atom_b)
    qtb(10) = qtb(10) + temp*d3(36)
    qtb(11) = qtb(11) + temp*d3(37)
    qtb(12) = qtb(12) + temp*d3(38)
    qtb(13) = qtb(13) + temp*d3(39)
    qtb(14) = qtb(14) + temp*d3(40)
    qtb(15) = qtb(15) + temp*d3(41)
    qtb(16) = qtb(16) + temp*d3(42)
    temp = mp_lab_q(16, atom_b)
    qtb(10) = qtb(10) + temp*d3(43)
    qtb(11) = qtb(11) + temp*d3(44)
    qtb(12) = qtb(12) + temp*d3(45)
    qtb(13) = qtb(13) + temp*d3(46)
    qtb(14) = qtb(14) + temp*d3(47)
    qtb(15) = qtb(15) + temp*d3(48)
    qtb(16) = qtb(16) + temp*d3(49)
    qtbx(10) = 2.449489742783178d0*qtb(12)
    qtbx(11) = 1.5811388300841898d0*qtb(14)
    qtbx(12) = -2.449489742783178d0*qtb(10) - 1.5811388300841898d0*qtb(13)
    qtbx(13) = 1.5811388300841898d0*qtb(12) + 1.224744871391589d0*qtb(16)
    qtbx(14) = -1.5811388300841898d0*qtb(11) - 1.224744871391589d0*qtb(15)
    qtbx(15) = 1.224744871391589d0*qtb(14)
    qtbx(16) = -1.224744871391589d0*qtb(13)
    qtby(10) = -2.449489742783178d0*qtb(11)
    qtby(11) = 2.449489742783178d0*qtb(10) - 1.5811388300841898d0*qtb(13)
    qtby(12) = -1.5811388300841898d0*qtb(14)
    qtby(13) = 1.5811388300841898d0*qtb(11) - 1.224744871391589d0*qtb(15)
    qtby(14) = 1.5811388300841898d0*qtb(12) - 1.224744871391589d0*qtb(16)
    qtby(15) = 1.224744871391589d0*qtb(13)
    qtby(16) = 1.224744871391589d0*qtb(14)
    if (lmax .lt. 4) return
    ! L = 4
    temp = mp_lab_q(17, atom_a)
    qta(17) = temp*d4(1)
    qta(18) = temp*d4(2)
    qta(19) = temp*d4(3)
    qta(20) = temp*d4(4)
    qta(21) = temp*d4(5)
    qta(22) = temp*d4(6)
    qta(23) = temp*d4(7)
    qta(24) = temp*d4(8)
    qta(25) = temp*d4(9)
    temp = mp_lab_q(18, atom_a)
    qta(17) = qta(17) + temp*d4(10)
    qta(18) = qta(18) + temp*d4(11)
    qta(19) = qta(19) + temp*d4(12)
    qta(20) = qta(20) + temp*d4(13)
    qta(21) = qta(21) + temp*d4(14)
    qta(22) = qta(22) + temp*d4(15)
    qta(23) = qta(23) + temp*d4(16)
    qta(24) = qta(24) + temp*d4(17)
    qta(25) = qta(25) + temp*d4(18)
    temp = mp_lab_q(19, atom_a)
    qta(17) = qta(17) + temp*d4(19)
    qta(18) = qta(18) + temp*d4(20)
    qta(19) = qta(19) + temp*d4(21)
    qta(20) = qta(20) + temp*d4(22)
    qta(21) = qta(21) + temp*d4(23)
    qta(22) = qta(22) + temp*d4(24)
    qta(23) = qta(23) + temp*d4(25)
    qta(24) = qta(24) + temp*d4(26)
    qta(25) = qta(25) + temp*d4(27)
    temp = mp_lab_q(20, atom_a)
    qta(17) = qta(17) + temp*d4(28)
    qta(18) = qta(18) + temp*d4(29)
    qta(19) = qta(19) + temp*d4(30)
    qta(20) = qta(20) + temp*d4(31)
    qta(21) = qta(21) + temp*d4(32)
    qta(22) = qta(22) + temp*d4(33)
    qta(23) = qta(23) + temp*d4(34)
    qta(24) = qta(24) + temp*d4(35)
    qta(25) = qta(25) + temp*d4(36)
    temp = mp_lab_q(21, atom_a)
    qta(17) = qta(17) + temp*d4(37)
    qta(18) = qta(18) + temp*d4(38)
    qta(19) = qta(19) + temp*d4(39)
    qta(20) = qta(20) + temp*d4(40)
    qta(21) = qta(21) + temp*d4(41)
    qta(22) = qta(22) + temp*d4(42)
    qta(23) = qta(23) + temp*d4(43)
    qta(24) = qta(24) + temp*d4(44)
    qta(25) = qta(25) + temp*d4(45)
    temp = mp_lab_q(22, atom_a)
    qta(17) = qta(17) + temp*d4(46)
    qta(18) = qta(18) + temp*d4(47)
    qta(19) = qta(19) + temp*d4(48)
    qta(20) = qta(20) + temp*d4(49)
    qta(21) = qta(21) + temp*d4(50)
    qta(22) = qta(22) + temp*d4(51)
    qta(23) = qta(23) + temp*d4(52)
    qta(24) = qta(24) + temp*d4(53)
    qta(25) = qta(25) + temp*d4(54)
    temp = mp_lab_q(23, atom_a)
    qta(17) = qta(17) + temp*d4(55)
    qta(18) = qta(18) + temp*d4(56)
    qta(19) = qta(19) + temp*d4(57)
    qta(20) = qta(20) + temp*d4(58)
    qta(21) = qta(21) + temp*d4(59)
    qta(22) = qta(22) + temp*d4(60)
    qta(23) = qta(23) + temp*d4(61)
    qta(24) = qta(24) + temp*d4(62)
    qta(25) = qta(25) + temp*d4(63)
    temp = mp_lab_q(24, atom_a)
    qta(17) = qta(17) + temp*d4(64)
    qta(18) = qta(18) + temp*d4(65)
    qta(19) = qta(19) + temp*d4(66)
    qta(20) = qta(20) + temp*d4(67)
    qta(21) = qta(21) + temp*d4(68)
    qta(22) = qta(22) + temp*d4(69)
    qta(23) = qta(23) + temp*d4(70)
    qta(24) = qta(24) + temp*d4(71)
    qta(25) = qta(25) + temp*d4(72)
    temp = mp_lab_q(25, atom_a)
    qta(17) = qta(17) + temp*d4(73)
    qta(18) = qta(18) + temp*d4(74)
    qta(19) = qta(19) + temp*d4(75)
    qta(20) = qta(20) + temp*d4(76)
    qta(21) = qta(21) + temp*d4(77)
    qta(22) = qta(22) + temp*d4(78)
    qta(23) = qta(23) + temp*d4(79)
    qta(24) = qta(24) + temp*d4(80)
    qta(25) = qta(25) + temp*d4(81)
    qtax(17) = 3.1622776601683795d0*qta(19)
    qtax(18) = 2.1213203435596424d0*qta(21)
    qtax(19) = -3.1622776601683795d0*qta(17) - 2.1213203435596424d0*qta(20)
    qtax(20) = 2.1213203435596424d0*qta(19) + 1.8708286933869707d0*qta(23)
    qtax(21) = -2.1213203435596424d0*qta(18) - 1.8708286933869707d0*qta(22)
    qtax(22) = 1.8708286933869707d0*qta(21) + 1.4142135623730951d0*qta(25)
    qtax(23) = -1.8708286933869707d0*qta(20) - 1.4142135623730951d0*qta(24)
    qtax(24) = 1.4142135623730951d0*qta(23)
    qtax(25) = -1.4142135623730951d0*qta(22)
    qtay(17) = -3.1622776601683795d0*qta(18)
    qtay(18) = 3.1622776601683795d0*qta(17) - 2.1213203435596424d0*qta(20)
    qtay(19) = -2.1213203435596424d0*qta(21)
    qtay(20) = 2.1213203435596424d0*qta(18) - 1.8708286933869707d0*qta(22)
    qtay(21) = 2.1213203435596424d0*qta(19) - 1.8708286933869707d0*qta(23)
    qtay(22) = 1.8708286933869707d0*qta(20) - 1.4142135623730951d0*qta(24)
    qtay(23) = 1.8708286933869707d0*qta(21) - 1.4142135623730951d0*qta(25)
    qtay(24) = 1.4142135623730951d0*qta(22)
    qtay(25) = 1.4142135623730951d0*qta(23)
    qtaz(17) = 0.d0
    qtaz(18) = -qta(19)
    qtaz(19) = qta(18)
    qtaz(20) = -2.d0*qta(21)
    qtaz(21) = 2.d0*qta(20)
    qtaz(22) = -3.d0*qta(23)
    qtaz(23) = 3.d0*qta(22)
    qtaz(24) = -4.d0*qta(25)
    qtaz(25) = 4.d0*qta(24)
    !call DGEMV('n', 9, 9, ONE, d4, 9, mp_lab_Q(17,atom_b), 1, ZERO, qtb(17), 1)
    temp = mp_lab_q(17, atom_b)
    qtb(17) = temp*d4(1)
    qtb(18) = temp*d4(2)
    qtb(19) = temp*d4(3)
    qtb(20) = temp*d4(4)
    qtb(21) = temp*d4(5)
    qtb(22) = temp*d4(6)
    qtb(23) = temp*d4(7)
    qtb(24) = temp*d4(8)
    qtb(25) = temp*d4(9)
    temp = mp_lab_q(18, atom_b)
    qtb(17) = qtb(17) + temp*d4(10)
    qtb(18) = qtb(18) + temp*d4(11)
    qtb(19) = qtb(19) + temp*d4(12)
    qtb(20) = qtb(20) + temp*d4(13)
    qtb(21) = qtb(21) + temp*d4(14)
    qtb(22) = qtb(22) + temp*d4(15)
    qtb(23) = qtb(23) + temp*d4(16)
    qtb(24) = qtb(24) + temp*d4(17)
    qtb(25) = qtb(25) + temp*d4(18)
    temp = mp_lab_q(19, atom_b)
    qtb(17) = qtb(17) + temp*d4(19)
    qtb(18) = qtb(18) + temp*d4(20)
    qtb(19) = qtb(19) + temp*d4(21)
    qtb(20) = qtb(20) + temp*d4(22)
    qtb(21) = qtb(21) + temp*d4(23)
    qtb(22) = qtb(22) + temp*d4(24)
    qtb(23) = qtb(23) + temp*d4(25)
    qtb(24) = qtb(24) + temp*d4(26)
    qtb(25) = qtb(25) + temp*d4(27)
    temp = mp_lab_q(20, atom_b)
    qtb(17) = qtb(17) + temp*d4(28)
    qtb(18) = qtb(18) + temp*d4(29)
    qtb(19) = qtb(19) + temp*d4(30)
    qtb(20) = qtb(20) + temp*d4(31)
    qtb(21) = qtb(21) + temp*d4(32)
    qtb(22) = qtb(22) + temp*d4(33)
    qtb(23) = qtb(23) + temp*d4(34)
    qtb(24) = qtb(24) + temp*d4(35)
    qtb(25) = qtb(25) + temp*d4(36)
    temp = mp_lab_q(21, atom_b)
    qtb(17) = qtb(17) + temp*d4(37)
    qtb(18) = qtb(18) + temp*d4(38)
    qtb(19) = qtb(19) + temp*d4(39)
    qtb(20) = qtb(20) + temp*d4(40)
    qtb(21) = qtb(21) + temp*d4(41)
    qtb(22) = qtb(22) + temp*d4(42)
    qtb(23) = qtb(23) + temp*d4(43)
    qtb(24) = qtb(24) + temp*d4(44)
    qtb(25) = qtb(25) + temp*d4(45)
    temp = mp_lab_q(22, atom_b)
    qtb(17) = qtb(17) + temp*d4(46)
    qtb(18) = qtb(18) + temp*d4(47)
    qtb(19) = qtb(19) + temp*d4(48)
    qtb(20) = qtb(20) + temp*d4(49)
    qtb(21) = qtb(21) + temp*d4(50)
    qtb(22) = qtb(22) + temp*d4(51)
    qtb(23) = qtb(23) + temp*d4(52)
    qtb(24) = qtb(24) + temp*d4(53)
    qtb(25) = qtb(25) + temp*d4(54)
    temp = mp_lab_q(23, atom_b)
    qtb(17) = qtb(17) + temp*d4(55)
    qtb(18) = qtb(18) + temp*d4(56)
    qtb(19) = qtb(19) + temp*d4(57)
    qtb(20) = qtb(20) + temp*d4(58)
    qtb(21) = qtb(21) + temp*d4(59)
    qtb(22) = qtb(22) + temp*d4(60)
    qtb(23) = qtb(23) + temp*d4(61)
    qtb(24) = qtb(24) + temp*d4(62)
    qtb(25) = qtb(25) + temp*d4(63)
    temp = mp_lab_q(24, atom_b)
    qtb(17) = qtb(17) + temp*d4(64)
    qtb(18) = qtb(18) + temp*d4(65)
    qtb(19) = qtb(19) + temp*d4(66)
    qtb(20) = qtb(20) + temp*d4(67)
    qtb(21) = qtb(21) + temp*d4(68)
    qtb(22) = qtb(22) + temp*d4(69)
    qtb(23) = qtb(23) + temp*d4(70)
    qtb(24) = qtb(24) + temp*d4(71)
    qtb(25) = qtb(25) + temp*d4(72)
    temp = mp_lab_q(25, atom_b)
    qtb(17) = qtb(17) + temp*d4(73)
    qtb(18) = qtb(18) + temp*d4(74)
    qtb(19) = qtb(19) + temp*d4(75)
    qtb(20) = qtb(20) + temp*d4(76)
    qtb(21) = qtb(21) + temp*d4(77)
    qtb(22) = qtb(22) + temp*d4(78)
    qtb(23) = qtb(23) + temp*d4(79)
    qtb(24) = qtb(24) + temp*d4(80)
    qtb(25) = qtb(25) + temp*d4(81)
    qtbx(17) = 3.1622776601683795d0*qtb(19)
    qtbx(18) = 2.1213203435596424d0*qtb(21)
    qtbx(19) = -3.1622776601683795d0*qtb(17) - 2.1213203435596424d0*qtb(20)
    qtbx(20) = 2.1213203435596424d0*qtb(19) + 1.8708286933869707d0*qtb(23)
    qtbx(21) = -2.1213203435596424d0*qtb(18) - 1.8708286933869707d0*qtb(22)
    qtbx(22) = 1.8708286933869707d0*qtb(21) + 1.4142135623730951d0*qtb(25)
    qtbx(23) = -1.8708286933869707d0*qtb(20) - 1.4142135623730951d0*qtb(24)
    qtbx(24) = 1.4142135623730951d0*qtb(23)
    qtbx(25) = -1.4142135623730951d0*qtb(22)
    qtby(17) = -3.1622776601683795d0*qtb(18)
    qtby(18) = 3.1622776601683795d0*qtb(17) - 2.1213203435596424d0*qtb(20)
    qtby(19) = -2.1213203435596424d0*qtb(21)
    qtby(20) = 2.1213203435596424d0*qtb(18) - 1.8708286933869707d0*qtb(22)
    qtby(21) = 2.1213203435596424d0*qtb(19) - 1.8708286933869707d0*qtb(23)
    qtby(22) = 1.8708286933869707d0*qtb(20) - 1.4142135623730951d0*qtb(24)
    qtby(23) = 1.8708286933869707d0*qtb(21) - 1.4142135623730951d0*qtb(25)
    qtby(24) = 1.4142135623730951d0*qtb(22)
    qtby(25) = 1.4142135623730951d0*qtb(23)
    !!DEBUGPRINT
    write(*,*)
    write(*,*) "====================================================================================================================================="
    write(*,*)
    write(*,*) "Building QI intermediates for pair", atom_a, ",", atom_b
    write(*,*)
    write(*,*) "The Uab orientation matrix"
    write(*,'(3F16.10)') Uz(1), Uz(4), Uz(7)
    write(*,'(3F16.10)') Uz(2), Uz(5), Uz(7)
    write(*,'(3F16.10)') Uz(3), Uz(6), Uz(7)
    write(*,*)
    write(*,*)
    write(*,*) "The D1 rotation matrix for dipoles (Uab transpose, reordered for spherical harmonics)"
    write(*,'(3F16.10)') d1(1), d1(4), d1(7)
    write(*,'(3F16.10)') d1(2), d1(5), d1(7)
    write(*,'(3F16.10)') d1(3), d1(6), d1(7)
    write(*,*)
    write(*,*)
    write(*,*) "The D2 rotation matrix for quadrupoles"
    write(*,'(5F16.10)') d2(1), d2(6), d2(11), d2(16), d2(21)
    write(*,'(5F16.10)') d2(2), d2(7), d2(12), d2(17), d2(22)
    write(*,'(5F16.10)') d2(3), d2(8), d2(13), d2(18), d2(23)
    write(*,'(5F16.10)') d2(4), d2(9), d2(14), d2(19), d2(24)
    write(*,'(5F16.10)') d2(5), d2(10), d2(15), d2(20), d2(25)
    write(*,*)
    write(*,*)
    write(*,*) "The D3 rotation matrix for octopoles"
    write(*,'(7F16.10)') d3(1), d3(8), d3(15), d3(22), d3(29), d3(36), d3(43)
    write(*,'(7F16.10)') d3(2), d3(9), d3(16), d3(23), d3(30), d3(37), d3(44)
    write(*,'(7F16.10)') d3(3), d3(10), d3(17), d3(24), d3(31), d3(38), d3(45)
    write(*,'(7F16.10)') d3(4), d3(11), d3(18), d3(25), d3(32), d3(39), d3(46)
    write(*,'(7F16.10)') d3(5), d3(12), d3(19), d3(26), d3(33), d3(40), d3(47)
    write(*,'(7F16.10)') d3(6), d3(13), d3(20), d3(27), d3(34), d3(41), d3(48)
    write(*,'(7F16.10)') d3(7), d3(14), d3(21), d3(28), d3(35), d3(42), d3(49)
    write(*,*)
    write(*,*)
    write(*,*) "The D4 rotation matrix for hexadecapoles"
    write(*,'(9F16.10)') d4(1), d4(10), d4(19), d4(28), d4(37), d4(46), d4(55), d4(64), d4(73)
    write(*,'(9F16.10)') d4(2), d4(11), d4(20), d4(29), d4(38), d4(47), d4(56), d4(65), d4(74)
    write(*,'(9F16.10)') d4(3), d4(12), d4(21), d4(30), d4(39), d4(48), d4(57), d4(66), d4(75)
    write(*,'(9F16.10)') d4(4), d4(13), d4(22), d4(31), d4(40), d4(49), d4(58), d4(67), d4(76)
    write(*,'(9F16.10)') d4(5), d4(14), d4(23), d4(32), d4(41), d4(50), d4(59), d4(68), d4(77)
    write(*,'(9F16.10)') d4(6), d4(15), d4(24), d4(33), d4(42), d4(51), d4(60), d4(69), d4(78)
    write(*,'(9F16.10)') d4(7), d4(16), d4(25), d4(34), d4(43), d4(52), d4(61), d4(70), d4(79)
    write(*,'(9F16.10)') d4(8), d4(17), d4(26), d4(35), d4(44), d4(53), d4(62), d4(71), d4(80)
    write(*,'(9F16.10)') d4(9), d4(18), d4(27), d4(36), d4(45), d4(54), d4(63), d4(72), d4(81)
    write(*,*)
    write(*,*) "The multipoles on A in the QI frame."
    do kk = 1, 25
        write(*,'(A10I2A5F16.10)') "Qta(",kk,") = ", qta(kk)
    enddo
    write(*,*)
    write(*,*) "The multipoles on B in the QI frame."
    do kk = 1, 25
        write(*,'(A10I2A5F16.10)') "Qtb(",kk,") = ", qtb(kk)
    enddo
    write(*,*)
    write(*,*) "The x torque intermediates for A in the QI frame."
    do kk = 1, 25
        write(*,'(A10I2A5F16.10)') "Qtax(",kk,") = ", qtax(kk)
    enddo
    write(*,*)
    write(*,*) "The y torque intermediates for A in the QI frame."
    do kk = 1, 25
        write(*,'(A10I2A5F16.10)') "Qtay(",kk,") = ", qtay(kk)
    enddo
    write(*,*)
    write(*,*) "The z torque intermediates for A in the QI frame."
    do kk = 1, 25
        write(*,'(A10I2A5F16.10)') "Qtaz(",kk,") = ", qtaz(kk)
    enddo
    write(*,*)
    write(*,*) "The x torque intermediates for B in the QI frame."
    do kk = 1, 25
        write(*,'(A10I2A5F16.10)') "Qtbx(",kk,") = ", qtbx(kk)
    enddo
    write(*,*)
    write(*,*) "The y torque intermediates for B in the QI frame."
    do kk = 1, 25
        write(*,'(A10I2A5F16.10)') "Qtby(",kk,") = ", qtby(kk)
    enddo
    write(*,*)
    write(*,*) "====================================================================================================================================="
    write(*,*)
    return
  end subroutine form_qi
