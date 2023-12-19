module m_zrand
implicit none
contains
  ! Wrapper to generate random uniform value, it is here to quickly substitute
  ! F77 drand function.
  function zrand() result (rand_numb)
    real(selected_real_kind(8)) :: rand_numb
    call random_number( rand_numb )
  end function zrand
end module m_zrand