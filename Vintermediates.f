    !
    ! N.B. The intermediate refered to as B herein is really the D_l(k R) intermediate from eq 29 of our paper.
    !
    integer, intent(in) :: lmax ! this allows us to shortcut cases where only low angular momentum is needed
    real*8, intent(in) :: X, rvec(10), krvec(11), b(6), qta(25), qtb(25)
    real*8, intent(out) :: vab(25), vba(25), vabR(25)
    integer kk
    !
    real*8 :: e, d ! just some local variables
    ! L = 0
    ! C-C
    ! m=0
    e = rvec(1)*(B(2) - krvec(1)*X)
    d = -B(2)*rvec(2)
    vab(1) = e * qtb(1)
    vba(1) = e * qta(1)
    vabR(1) = d * qtb(1)
    if(lmax .lt. 1) return
    ! L = 1
    ! C-D
    ! m=0
    e = -B(2)*rvec(2)
    d = 2.d0*rvec(3)*(B(2) + krvec(3)*X)
    vab(1) = vab(1) + e * qtb(2)
    vba(2) = e * qta(1)
    vabR(1) = vabR(1) + d * qtb(2)
    ! D-C
    ! m=0
    e = B(2)*rvec(2)
    d = -2.d0*rvec(3)*(B(2) + krvec(3)*X)
    vab(2) = e * qtb(1)
    vba(1) = vba(1) + e * qta(2)
    vabR(2) = d * qtb(1)
    ! D-D
    ! m=0
    e = -0.6666666666666667d0*rvec(3)*(3.d0*B(3) + krvec(3)*X)
    d = rvec(4)*(6.d0*B(3) + 4.d0*krvec(5)*X)
    vab(2) = vab(2) + e * qtb(2)
    vba(2) = vba(2) + e * qta(2)
    vabR(2) = vabR(2) + d * qtb(2)
    ! m=1
    e = rvec(3)*(B(3) - 0.6666666666666667d0*krvec(3)*X)
    d = -3.d0*B(3)*rvec(4)
    vab(3) =  e * qtb(3)
    vba(3) =  e * qta(3)
    vabR(3) =  d * qtb(3)
    vab(4) =  e * qtb(4)
    vba(4) =  e * qta(4)
    vabR(4) =  d * qtb(4)
    if(lmax .lt. 2) return

    ! L = 2
    ! C-Q
    ! m=0
    e = B(3)*rvec(3)
    d = -0.3333333333333333d0*rvec(4)*(9.d0*B(3) + 4.d0*krvec(5)*X)
    vab(1) = vab(1) + e * qtb(5)
    vba(5) = e * qta(1)
    vabR(1) = vabR(1) + d * qtb(5)
    ! D-Q
    ! m=0
    e = rvec(4)*(3.d0*B(3) + 1.333333333333333d0*krvec(5)*X)
    d = -1.333333333333333d0*rvec(5)*(9.d0*B(3) + 2.d0*(1.d0 + krvec(2))*krvec(5)*&
               X)
    vab(2) = vab(2) + e * qtb(5)
    vba(5) = vba(5) + e * qta(2)
    vabR(2) = vabR(2) + d * qtb(5)
    ! m=1
    e = -1.732050807568877d0*B(3)*rvec(4)
    d = 2.309401076758503d0*rvec(5)*(3.d0*B(3) + krvec(5)*X)
    vab(3) = vab(3) +  e * qtb(6)
    vba(6) =  e * qta(3)
    vabR(3) = vabR(3) +  d * qtb(6)
    vab(4) = vab(4) +  e * qtb(7)
    vba(7) =  e * qta(4)
    vabR(4) = vabR(4) +  d * qtb(7)
    ! Q-C
    ! m=0
    e = B(3)*rvec(3)
    d = -0.3333333333333333d0*rvec(4)*(9.d0*B(3) + 4.d0*krvec(5)*X)
    vab(5) = e * qtb(1)
    vba(1) = vba(1) + e * qta(5)
    vabR(5) = d * qtb(1)
    ! Q-D
    ! m=0
    e = -3.d0*rvec(4)*(B(3) + 0.4444444444444444d0*krvec(5)*X)
    d = 1.333333333333333d0*rvec(5)*(9.d0*B(3) + 2.d0*(1.d0 + krvec(2))*krvec(5)*&
               X)
    vab(5) = vab(5) + e * qtb(2)
    vba(2) = vba(2) + e * qta(5)
    vabR(5) = vabR(5) + d * qtb(2)
    ! m=1
    e = 1.732050807568877d0*B(3)*rvec(4)
    d = -2.309401076758503d0*rvec(5)*(3.d0*B(3) + krvec(5)*X)
    vab(6) =  e * qtb(3)
    vba(3) = vba(3) +  e * qta(6)
    vabR(6) =  d * qtb(3)
    vab(7) =  e * qtb(4)
    vba(4) = vba(4) +  e * qta(7)
    vabR(7) =  d * qtb(4)
    ! Q-Q
    ! m=0
    e = rvec(5)*(6.d0*B(4) + 0.08888888888888889d0*(-3.d0 + 10.d0*krvec(2))*&
               krvec(5)*X)
    d = -0.2222222222222222d0*rvec(6)*(135.d0*B(4) + 4.d0*(1.d0 + 2.d0*krvec(2))*&
               krvec(7)*X)
    vab(5) = vab(5) + e * qtb(5)
    vba(5) = vba(5) + e * qta(5)
    vabR(5) = vabR(5) + d * qtb(5)
    ! m=1
    e = -0.2666666666666667d0*rvec(5)*(15.d0*B(4) + krvec(5)*X)
    d = rvec(6)*(20.d0*B(4) + 2.666666666666667d0*krvec(7)*X)
    vab(6) = vab(6) +  e * qtb(6)
    vba(6) = vba(6) +  e * qta(6)
    vabR(6) = vabR(6) +  d * qtb(6)
    vab(7) = vab(7) +  e * qtb(7)
    vba(7) = vba(7) +  e * qta(7)
    vabR(7) = vabR(7) +  d * qtb(7)
    ! m=2
    e = rvec(5)*(B(4) - 0.2666666666666667d0*krvec(5)*X)
    d = -5.d0*B(4)*rvec(6)
    vab(8) =  e * qtb(8)
    vba(8) =  e * qta(8)
    vabR(8) =  d * qtb(8)
    vab(9) =  e * qtb(9)
    vba(9) =  e * qta(9)
    vabR(9) =  d * qtb(9)
    if(lmax .lt. 3) return

    ! L = 3
    ! C-O
    ! m=0
    e = rvec(4)*(-B(3) - 0.2666666666666667d0*krvec(5)*X)
    d = 0.2666666666666667d0*rvec(5)*(15.d0*B(3) + 2.d0*(2.d0*krvec(5) + krvec(7))*&
               X)
    vab(1) = vab(1) + e * qtb(10)
    vba(10) = e * qta(1)
    vabR(1) = vabR(1) + d * qtb(10)
    ! D-O
    ! m=0
    e = -4.d0*rvec(5)*(B(4) + 0.1333333333333333d0*krvec(7)*X)
    d = 0.2666666666666667d0*rvec(6)*(75.d0*B(4) + 4.d0*(1.d0 + krvec(2))*krvec(7)*&
               X)
    vab(2) = vab(2) + e * qtb(10)
    vba(10) = vba(10) + e * qta(2)
    vabR(2) = vabR(2) + d * qtb(10)
    ! m=1
    e = 2.449489742783178d0*B(4)*rvec(5)
    d = -0.1632993161855452d0*rvec(6)*(75.d0*B(4) + 8.d0*krvec(7)*X)
    vab(3) = vab(3) +  e * qtb(11)
    vba(11) =  e * qta(3)
    vabR(3) = vabR(3) +  d * qtb(11)
    vab(4) = vab(4) +  e * qtb(12)
    vba(12) =  e * qta(4)
    vabR(4) = vabR(4) +  d * qtb(12)
    ! Q-O
    ! m=0
    e = rvec(6)*(-10.d0*B(4) - 0.1777777777777778d0*(3.d0 + 2.d0*krvec(2))*&
               krvec(7)*X)
    d = 0.08888888888888889d0*rvec(7)*(675.d0*B(4) + 2.d0*(27.d0 + 4.d0*krvec(2)**&
               2)*krvec(7)*X)
    vab(5) = vab(5) + e * qtb(10)
    vba(10) = vba(10) + e * qta(5)
    vabR(5) = vabR(5) + d * qtb(10)
    ! m=1
    e = 7.071067811865475d0*rvec(6)*(B(4) + 0.1066666666666667d0*krvec(7)*X)
    d = -0.1885618083164127d0*rvec(7)*(225.d0*B(4) + 8.d0*(2.d0 + krvec(2))*&
               krvec(7)*X)
    vab(6) = vab(6) +  e * qtb(11)
    vba(11) = vba(11) +  e * qta(6)
    vabR(6) = vabR(6) +  d * qtb(11)
    vab(7) = vab(7) +  e * qtb(12)
    vba(12) = vba(12) +  e * qta(7)
    vabR(7) = vabR(7) +  d * qtb(12)
    ! m=2
    e = -2.23606797749979d0*B(4)*rvec(6)
    d = 0.298142396999972d0*rvec(7)*(45.d0*B(4) + 4.d0*krvec(7)*X)
    vab(8) = vab(8) +  e * qtb(13)
    vba(13) =  e * qta(8)
    vabR(8) = vabR(8) +  d * qtb(13)
    vab(9) = vab(9) +  e * qtb(14)
    vba(14) =  e * qta(9)
    vabR(9) = vabR(9) +  d * qtb(14)
    ! O-C
    ! m=0
    e = rvec(4)*(B(3) + 0.2666666666666667d0*krvec(5)*X)
    d = -0.2666666666666667d0*rvec(5)*(15.d0*B(3) + 2.d0*(2.d0*krvec(5) +&
                krvec(7))*X)
    vab(10) = e * qtb(1)
    vba(1) = vba(1) + e * qta(10)
    vabR(10) = d * qtb(1)
    ! O-D
    ! m=0
    e = -4.d0*rvec(5)*(B(4) + 0.1333333333333333d0*krvec(7)*X)
    d = 0.2666666666666667d0*rvec(6)*(75.d0*B(4) + 4.d0*(1.d0 + krvec(2))*krvec(7)*&
               X)
    vab(10) = vab(10) + e * qtb(2)
    vba(2) = vba(2) + e * qta(10)
    vabR(10) = vabR(10) + d * qtb(2)
    ! m=1
    e = 2.449489742783178d0*B(4)*rvec(5)
    d = -0.1632993161855452d0*rvec(6)*(75.d0*B(4) + 8.d0*krvec(7)*X)
    vab(11) =  e * qtb(3)
    vba(3) = vba(3) +  e * qta(11)
    vabR(11) =  d * qtb(3)
    vab(12) =  e * qtb(4)
    vba(4) = vba(4) +  e * qta(12)
    vabR(12) =  d * qtb(4)
    ! O-Q
    ! m=0
    e = rvec(6)*(10.d0*B(4) + 0.1777777777777778d0*(3.d0 + 2.d0*krvec(2))*krvec(7)*&
               X)
    d = -0.08888888888888889d0*rvec(7)*(675.d0*B(4) + 2.d0*(27.d0 + 4.d0*krvec(2)**&
               2)*krvec(7)*X)
    vab(10) = vab(10) + e * qtb(5)
    vba(5) = vba(5) + e * qta(10)
    vabR(10) = vabR(10) + d * qtb(5)
    ! m=1
    e = -7.071067811865475d0*rvec(6)*(B(4) + 0.1066666666666667d0*krvec(7)*X)
    d = 0.1885618083164127d0*rvec(7)*(225.d0*B(4) + 8.d0*(2.d0 + krvec(2))*&
               krvec(7)*X)
    vab(11) = vab(11) +  e * qtb(6)
    vba(6) = vba(6) +  e * qta(11)
    vabR(11) = vabR(11) +  d * qtb(6)
    vab(12) = vab(12) +  e * qtb(7)
    vba(7) = vba(7) +  e * qta(12)
    vabR(12) = vabR(12) +  d * qtb(7)
    ! m=2
    e = 2.23606797749979d0*B(4)*rvec(6)
    d = -0.298142396999972d0*rvec(7)*(45.d0*B(4) + 4.d0*krvec(7)*X)
    vab(13) =  e * qtb(8)
    vba(8) = vba(8) +  e * qta(13)
    vabR(13) =  d * qtb(8)
    vab(14) =  e * qtb(9)
    vba(9) = vba(9) +  e * qta(14)
    vabR(14) =  d * qtb(9)
    ! O-O
    ! m=0
    e = rvec(7)*(-20.d0*B(5) - 0.005079365079365079d0*(15.d0 + 28.d0*krvec(2) +&
                28.d0*krvec(4))*krvec(7)*X)
    d = 0.01777777777777778d0*rvec(8)*(7875.d0*B(5) + 4.d0*(41.d0 - 4.d0*krvec(2) +&
                4.d0*krvec(4))*krvec(9)*X)
    vab(10) = vab(10) + e * qtb(10)
    vba(10) = vba(10) + e * qta(10)
    vabR(10) = vabR(10) + d * qtb(10)
    ! m=1
    e = rvec(7)*(15.d0*B(5) + 0.01523809523809524d0*(-5.d0 + 28.d0*krvec(2))*&
               krvec(7)*X)
    d = -0.01333333333333333d0*rvec(8)*(7875.d0*B(5) + 32.d0*(3.d0 + 2.d0*&
               krvec(2))*krvec(9)*X)
    vab(11) = vab(11) +  e * qtb(11)
    vba(11) = vba(11) +  e * qta(11)
    vabR(11) = vabR(11) +  d * qtb(11)
    vab(12) = vab(12) +  e * qtb(12)
    vba(12) = vba(12) +  e * qta(12)
    vabR(12) = vabR(12) +  d * qtb(12)
    ! m=2
    e = rvec(7)*(-6.d0*B(5) - 0.07619047619047619d0*krvec(7)*X)
    d = rvec(8)*(42.d0*B(5) + 1.066666666666667d0*krvec(9)*X)
    vab(13) = vab(13) +  e * qtb(13)
    vba(13) = vba(13) +  e * qta(13)
    vabR(13) = vabR(13) +  d * qtb(13)
    vab(14) = vab(14) +  e * qtb(14)
    vba(14) = vba(14) +  e * qta(14)
    vabR(14) = vabR(14) +  d * qtb(14)
    ! m=3
    e = rvec(7)*(B(5) - 0.07619047619047619d0*krvec(7)*X)
    d = -7.d0*B(5)*rvec(8)
    vab(15) =  e * qtb(15)
    vba(15) =  e * qta(15)
    vabR(15) =  d * qtb(15)
    vab(16) =  e * qtb(16)
    vba(16) =  e * qta(16)
    vabR(16) =  d * qtb(16)
    if(lmax .lt. 4) return

    ! L = 4
    ! C-H
    ! m=0
    e = rvec(5)*(B(4) + 0.07619047619047619d0*krvec(7)*X)
    d = -0.009523809523809524d0*rvec(6)*(525.d0*B(4) + 8.d0*(5.d0*krvec(7) + 2.d0*&
               krvec(9))*X)
    vab(1) = vab(1) + e * qtb(17)
    vba(17) = e * qta(1)
    vabR(1) = vabR(1) + d * qtb(17)
    ! D-H
    ! m=0
    e = rvec(6)*(5.d0*B(4) + 0.07619047619047619d0*(5.d0 + 2.d0*krvec(2))*krvec(7)*&
               X)
    d = -0.01904761904761905d0*rvec(7)*(1575.d0*B(4) + 8.d0*(3.d0*(5.d0 + 2.d0*&
               krvec(2))*krvec(7) + 2.d0*(-2.d0 + krvec(2))*krvec(9))*X)
    vab(2) = vab(2) + e * qtb(17)
    vba(17) = vba(17) + e * qta(2)
    vabR(2) = vabR(2) + d * qtb(17)
    ! m=1
    e = -3.162277660168379d0*rvec(6)*(B(4) + 0.07619047619047619d0*krvec(7)*X)
    d = 0.06023386019368342d0*rvec(7)*(315.d0*B(4) + 8.d0*(3.d0*krvec(7) +&
                krvec(9))*X)
    vab(3) = vab(3) +  e * qtb(18)
    vba(18) =  e * qta(3)
    vabR(3) = vabR(3) +  d * qtb(18)
    vab(4) = vab(4) +  e * qtb(19)
    vba(19) =  e * qta(4)
    vabR(4) = vabR(4) +  d * qtb(19)
    ! Q-H
    ! m=0
    e = rvec(7)*(15.d0*B(5) + 0.05079365079365079d0*(3.d0 + 2.d0*krvec(2))*&
               krvec(9)*X)
    d = -0.003174603174603175d0*rvec(8)*(33075.d0*B(5) + 16.d0*(39.d0 - 2.d0*&
               krvec(2) + 4.d0*krvec(4))*krvec(9)*X)
    vab(5) = vab(5) + e * qtb(17)
    vba(17) = vba(17) + e * qta(5)
    vabR(5) = vabR(5) + d * qtb(17)
    ! m=1
    e = -10.95445115010332d0*rvec(7)*(B(5) + 0.0253968253968254d0*krvec(9)*X)
    d = 0.0347760353971534d0*rvec(8)*(2205.d0*B(5) + 16.d0*(2.d0 + krvec(2))*&
               krvec(9)*X)
    vab(6) = vab(6) +  e * qtb(18)
    vba(18) = vba(18) +  e * qta(6)
    vabR(6) = vabR(6) +  d * qtb(18)
    vab(7) = vab(7) +  e * qtb(19)
    vba(19) = vba(19) +  e * qta(7)
    vabR(7) = vabR(7) +  d * qtb(19)
    ! m=2
    e = 3.872983346207417d0*B(5)*rvec(7)
    d = -0.03688555567816588d0*rvec(8)*(735.d0*B(5) + 16.d0*krvec(9)*X)
    vab(8) = vab(8) +  e * qtb(20)
    vba(20) =  e * qta(8)
    vabR(8) = vabR(8) +  d * qtb(20)
    vab(9) = vab(9) +  e * qtb(21)
    vba(21) =  e * qta(9)
    vabR(9) = vabR(9) +  d * qtb(21)
    ! O-H
    ! m=0
    e = rvec(8)*(35.d0*B(5) + 0.01015873015873016d0*(65.d0 + 2.d0*krvec(2) + 4.d0*&
               krvec(4))*krvec(9)*X)
    d = -0.005079365079365079d0*rvec(9)*(55125.d0*B(5) + 8.d0*(115.d0 + 31.d0*&
               krvec(2) - 4.d0*krvec(4) + 2.d0*krvec(6))*krvec(9)*X)
    vab(10) = vab(10) + e * qtb(17)
    vba(17) = vba(17) + e * qta(10)
    vabR(10) = vabR(10) + d * qtb(17)
    ! m=1
    e = -0.002459037045211058d0*rvec(8)*(11025.d0*B(5) + 32.d0*(5.d0 + 2.d0*&
               krvec(2))*krvec(9)*X)
    d = 0.01967229636168847d0*rvec(9)*(11025.d0*B(5) + 2.d0*(95.d0 + 8.d0*&
               krvec(2) + 8.d0*krvec(4))*krvec(9)*X)
    vab(11) = vab(11) +  e * qtb(18)
    vba(18) = vba(18) +  e * qta(11)
    vabR(11) = vabR(11) +  d * qtb(18)
    vab(12) = vab(12) +  e * qtb(19)
    vba(19) = vba(19) +  e * qta(12)
    vabR(12) = vabR(12) +  d * qtb(19)
    ! m=2
    e = 12.12435565298214d0*rvec(8)*(B(5) + 0.0217687074829932d0*krvec(9)*X)
    d = -0.1319657758147716d0*rvec(9)*(735.d0*B(5) + 4.d0*(3.d0 + krvec(2))*&
               krvec(9)*X)
    vab(13) = vab(13) +  e * qtb(20)
    vba(20) = vba(20) +  e * qta(13)
    vabR(13) = vabR(13) +  d * qtb(20)
    vab(14) = vab(14) +  e * qtb(21)
    vba(21) = vba(21) +  e * qta(14)
    vabR(14) = vabR(14) +  d * qtb(21)
    ! m=3
    e = -2.645751311064591d0*B(5)*rvec(8)
    d = 0.2015810522715879d0*rvec(9)*(105.d0*B(5) + 2.d0*krvec(9)*X)
    vab(15) = vab(15) +  e * qtb(22)
    vba(22) =  e * qta(15)
    vabR(15) = vabR(15) +  d * qtb(22)
    vab(16) = vab(16) +  e * qtb(23)
    vba(23) =  e * qta(16)
    vabR(16) = vabR(16) +  d * qtb(23)
    ! H-C
    ! m=0
    e = rvec(5)*(B(4) + 0.07619047619047619d0*krvec(7)*X)
    d = -0.009523809523809524d0*rvec(6)*(525.d0*B(4) + 8.d0*(5.d0*krvec(7) + 2.d0*&
               krvec(9))*X)
    vab(17) = e * qtb(1)
    vba(1) = vba(1) + e * qta(17)
    vabR(17) = d * qtb(1)
    ! H-D
    ! m=0
    e = -5.d0*rvec(6)*(B(4) + 0.01523809523809524d0*(5.d0 + 2.d0*krvec(2))*&
               krvec(7)*X)
    d = rvec(7)*(30.d0*B(4) + 0.1523809523809524d0*(3.d0*(5.d0 + 2.d0*krvec(2))*&
               krvec(7) + 2.d0*(-2.d0 + krvec(2))*krvec(9))*X)
    vab(17) = vab(17) + e * qtb(2)
    vba(2) = vba(2) + e * qta(17)
    vabR(17) = vabR(17) + d * qtb(2)
    ! m=1
    e = 3.162277660168379d0*rvec(6)*(B(4) + 0.07619047619047619d0*krvec(7)*X)
    d = -0.06023386019368342d0*rvec(7)*(315.d0*B(4) + 8.d0*(3.d0*krvec(7) +&
                krvec(9))*X)
    vab(18) =  e * qtb(3)
    vba(3) = vba(3) +  e * qta(18)
    vabR(18) =  d * qtb(3)
    vab(19) =  e * qtb(4)
    vba(4) = vba(4) +  e * qta(19)
    vabR(19) =  d * qtb(4)
    ! H-Q
    ! m=0
    e = rvec(7)*(15.d0*B(5) + 0.05079365079365079d0*(3.d0 + 2.d0*krvec(2))*&
               krvec(9)*X)
    d = -0.003174603174603175d0*rvec(8)*(33075.d0*B(5) + 16.d0*(39.d0 - 2.d0*&
               krvec(2) + 4.d0*krvec(4))*krvec(9)*X)
    vab(17) = vab(17) + e * qtb(5)
    vba(5) = vba(5) + e * qta(17)
    vabR(17) = vabR(17) + d * qtb(5)
    ! m=1
    e = -10.95445115010332d0*rvec(7)*(B(5) + 0.0253968253968254d0*krvec(9)*X)
    d = 0.0347760353971534d0*rvec(8)*(2205.d0*B(5) + 16.d0*(2.d0 + krvec(2))*&
               krvec(9)*X)
    vab(18) = vab(18) +  e * qtb(6)
    vba(6) = vba(6) +  e * qta(18)
    vabR(18) = vabR(18) +  d * qtb(6)
    vab(19) = vab(19) +  e * qtb(7)
    vba(7) = vba(7) +  e * qta(19)
    vabR(19) = vabR(19) +  d * qtb(7)
    ! m=2
    e = 3.872983346207417d0*B(5)*rvec(7)
    d = -0.03688555567816588d0*rvec(8)*(735.d0*B(5) + 16.d0*krvec(9)*X)
    vab(20) =  e * qtb(8)
    vba(8) = vba(8) +  e * qta(20)
    vabR(20) =  d * qtb(8)
    vab(21) =  e * qtb(9)
    vba(9) = vba(9) +  e * qta(21)
    vabR(21) =  d * qtb(9)
    ! H-O
    ! m=0
    e = rvec(8)*(-35.d0*B(5) - 0.01015873015873016d0*(65.d0 + 2.d0*krvec(2) + 4.d0*&
               krvec(4))*krvec(9)*X)
    d = 0.005079365079365079d0*rvec(9)*(55125.d0*B(5) + 8.d0*(115.d0 + 31.d0*&
               krvec(2) - 4.d0*krvec(4) + 2.d0*krvec(6))*krvec(9)*X)
    vab(17) = vab(17) + e * qtb(10)
    vba(10) = vba(10) + e * qta(17)
    vabR(17) = vabR(17) + d * qtb(10)
    ! m=1
    e = 0.002459037045211058d0*rvec(8)*(11025.d0*B(5) + 32.d0*(5.d0 + 2.d0*&
               krvec(2))*krvec(9)*X)
    d = -0.01967229636168847d0*rvec(9)*(11025.d0*B(5) + 2.d0*(95.d0 + 8.d0*&
               krvec(2) + 8.d0*krvec(4))*krvec(9)*X)
    vab(18) = vab(18) +  e * qtb(11)
    vba(11) = vba(11) +  e * qta(18)
    vabR(18) = vabR(18) +  d * qtb(11)
    vab(19) = vab(19) +  e * qtb(12)
    vba(12) = vba(12) +  e * qta(19)
    vabR(19) = vabR(19) +  d * qtb(12)
    ! m=2
    e = -12.12435565298214d0*rvec(8)*(B(5) + 0.0217687074829932d0*krvec(9)*X)
    d = 0.1319657758147716d0*rvec(9)*(735.d0*B(5) + 4.d0*(3.d0 + krvec(2))*&
               krvec(9)*X)
    vab(20) = vab(20) +  e * qtb(13)
    vba(13) = vba(13) +  e * qta(20)
    vabR(20) = vabR(20) +  d * qtb(13)
    vab(21) = vab(21) +  e * qtb(14)
    vba(14) = vba(14) +  e * qta(21)
    vabR(21) = vabR(21) +  d * qtb(14)
    ! m=3
    e = 2.645751311064591d0*B(5)*rvec(8)
    d = -0.2015810522715879d0*rvec(9)*(105.d0*B(5) + 2.d0*krvec(9)*X)
    vab(22) =  e * qtb(15)
    vba(15) = vba(15) +  e * qta(22)
    vabR(22) =  d * qtb(15)
    vab(23) =  e * qtb(16)
    vba(16) = vba(16) +  e * qta(23)
    vabR(23) =  d * qtb(16)
    ! H-H
    ! m=0
    e = rvec(9)*(70.d0*B(6) + 0.0004837490551776266d0*(-35.d0 + 6.d0*krvec(2)*&
               (95.d0 - 2.d0*krvec(2) + 4.d0*krvec(4)))*krvec(9)*X)
    d = -0.00018140589569161d0*rvec(10)*(3.472875d6*B(6) + 16.d0*krvec(11)*&
               (615.d0 + 198.d0*krvec(2) - 28.d0*krvec(4) + 8.d0*krvec(6))*X)
    vab(17) = vab(17) + e * qtb(17)
    vba(17) = vba(17) + e * qta(17)
    vabR(17) = vabR(17) + d * qtb(17)
    ! m=1
    e = rvec(9)*(-56.d0*B(6) - 0.002418745275888133d0*(7.d0 + 48.d0*krvec(2) +&
                24.d0*krvec(4))*krvec(9)*X)
    d = 0.0036281179138322d0*rvec(10)*(138915.d0*B(6) + 4.d0*krvec(11)*(117.d0 +&
                8.d0*krvec(4))*X)
    vab(18) = vab(18) +  e * qtb(18)
    vba(18) = vba(18) +  e * qta(18)
    vabR(18) = vabR(18) +  d * qtb(18)
    vab(19) = vab(19) +  e * qtb(19)
    vba(19) = vba(19) +  e * qta(19)
    vabR(19) = vabR(19) +  d * qtb(19)
    ! m=2
    e = rvec(9)*(28.d0*B(6) + 0.002418745275888133d0*(-7.d0 + 54.d0*krvec(2))*&
               krvec(9)*X)
    d = -0.0163265306122449d0*rvec(10)*(15435.d0*B(6) + 8.d0*krvec(11)*(5.d0 +&
                2.d0*krvec(2))*X)
    vab(20) = vab(20) +  e * qtb(20)
    vba(20) = vba(20) +  e * qta(20)
    vabR(20) = vabR(20) +  d * qtb(20)
    vab(21) = vab(21) +  e * qtb(21)
    vba(21) = vba(21) +  e * qta(21)
    vabR(21) = vabR(21) +  d * qtb(21)
    ! m=3
    e = rvec(9)*(-8.d0*B(6) - 0.01693121693121693d0*krvec(9)*X)
    d = rvec(10)*(72.d0*B(6) + 0.3047619047619048d0*krvec(11)*X)
    vab(22) = vab(22) +  e * qtb(22)
    vba(22) = vba(22) +  e * qta(22)
    vabR(22) = vabR(22) +  d * qtb(22)
    vab(23) = vab(23) +  e * qtb(23)
    vba(23) = vba(23) +  e * qta(23)
    vabR(23) = vabR(23) +  d * qtb(23)
    ! m=4
    e = rvec(9)*(B(6) - 0.01693121693121693d0*krvec(9)*X)
    d = -9.d0*B(6)*rvec(10)
    vab(24) =  e * qtb(24)
    vba(24) =  e * qta(24)
    vabR(24) =  d * qtb(24)
    vab(25) =  e * qtb(25)
    vba(25) =  e * qta(25)
    vabR(25) =  d * qtb(25)
    !DEBUGPRINT
    write(*,*)
    write(*,*) "====================================================================================================================================="
    write(*,*)
    write(*,*) "T tensor construction"
    write(*,*)
    write(*,'(A15F16.10)') "X intermediate", X
    write(*,*)
    write(*,*) "The D_l intermediate (eq. 29 in the JCP paper)"
    do kk = 1, 6
        write(*,'(A10I2A5F16.10)') "D_(",kk,") = ", B(kk)
    enddo
    write(*,*)
    write(*,*) "The rvec intermediate 1/(4 pi espsilon_0 R^n)"
    do kk = 1, 10
        write(*,'(A10I2A5F16.10)') "rvec(",kk,") = ", rvec(kk)
    enddo
    write(*,*)
    write(*,*) "The krvec intermediate (kappa*R)^n"
    do kk = 1, 11
        write(*,'(A10I2A5F16.10)') "krvec(",kk,") = ", krvec(kk)
    enddo
    write(*,*)
    write(*,*) "The Vab intermediate"
    do kk = 1, 25
        write(*,'(A10I2A5F16.10)') "Vab(",kk,") = ", Vab(kk)
    enddo
    write(*,*)
    write(*,*) "The Vba intermediate"
    do kk = 1, 25
        write(*,'(A10I2A5F16.10)') "Vba(",kk,") = ", Vba(kk)
    enddo
    write(*,*)
    write(*,*) "The VabR intermediate"
    do kk = 1, 25
        write(*,'(A10I2A5F16.10)') "VabR(",kk,") = ", VabR(kk)
    enddo
    write(*,*) "====================================================================================================================================="
    write(*,*)

    return
