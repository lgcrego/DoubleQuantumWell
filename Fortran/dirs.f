module dirs 
contains
  function dirname(number)
    integer,intent(in) :: number
    character(len=8)   :: dirname

    ! Cast the (rounded) number to string using 6 digits and
    ! leading zeros
    ! This is the same w/o leading zeros  
    write (dirname, '(I8)')  number

    ! This is for one digit (no rounding)
    !write (dirname, '(F4.1)')  number
  end function
end module

