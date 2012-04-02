!------------------------
! FORTRAN unit test utility
!
! Author: Andrew H. Chen meihome @at@ gmail.com
!------------------------
!
! Unit test framework for FORTRAN.  (FoRtran UnIT)
!
! This package is to perform unit test for FORTRAN subroutines
!
! The method used most are: assert_true, assert_equals
! 
! Coding convention:  
!   1) All private subroutines end with underscore.  i.e. init_fruit_
!   2) All methods must be exposed by interface.  i.e. interface init_fruit
!   3) Variable and methods are lower case connected with underscores.  i.e. init_fruit_, and failed_assert_count
!
! Modified by James Wookey to allow compilation in gfortran (see lines marked
! JW2012)
!

module fruit_util
  interface equals
     module procedure equalEpsilon
     module procedure floatEqual
     module procedure integerEqual
     module procedure doublePrecisionEqual
     module procedure stringEqual
     module procedure logicalEqual
  end interface

  interface to_s
     module procedure to_s_int_
     module procedure to_s_real_
     module procedure to_s_logical_
     module procedure to_s_double_
     module procedure to_s_complex_
     module procedure to_s_double_complex_
     module procedure to_s_string_
  end interface

  interface strip
     module procedure strip_
  end interface

  private:: &
       to_s_int_, to_s_real_, to_s_logical_, to_s_double_, to_s_complex_, to_s_double_complex_, to_s_string_, &
       strip_

contains

  function to_s_int_ (value)
    implicit none
    character(len=500):: to_s_int_
    integer, intent(in) :: value
    character(len=500) :: result
    write (result, *) value
    to_s_int_ = adjustl(trim(result))
  end function to_s_int_

  function to_s_real_ (value)
    implicit none
    character(len=500):: to_s_real_
    real, intent(in) :: value
    character(len=500) :: result
    write (result, *) value
    to_s_real_ = adjustl(trim(result))
  end function to_s_real_

  function to_s_double_ (value)
    implicit none
    character(len=500):: to_s_double_
    double precision, intent(in) :: value
    character(len=500) :: result
    write (result, *) value
    to_s_double_ = adjustl(trim(result))
  end function to_s_double_

  function to_s_complex_ (value)
    implicit none
    character(len=500):: to_s_complex_
    complex, intent(in) :: value
    character(len=500) :: result
    write (result, *) value
    to_s_complex_ = adjustl(trim(result))
  end function to_s_complex_

  function to_s_double_complex_ (value)
    implicit none
    character(len=500):: to_s_double_complex_
    double complex, intent(in) :: value
    character(len=500) :: result
    write (result, *) value
    to_s_double_complex_ = adjustl(trim(result))
  end function to_s_double_complex_

  function to_s_logical_ (value)
    implicit none
    character(len=500):: to_s_logical_
    logical, intent(in) :: value
    character(len=500) :: result
    write (result, *) value
    to_s_logical_ = adjustl(trim(result))
  end function to_s_logical_

  function to_s_string_ (value)
    implicit none
    character(len=500):: to_s_string_
    character(len=*), intent(in) :: value
    to_s_string_ = value
  end function to_s_string_

  function strip_(value)
    implicit none
    character(len=500):: strip_
    character(len=*), intent(in) :: value
    strip_ = trim(adjustl(value))
  end function strip_

  !------------------------
  ! test if 2 values are close
  !------------------------
  !logical function equals (number1, number2) 
  !  real,  intent (in) :: number1, number2
  !  
  !  return equalEpsilon (number1, number2, epsilon(number1))
  !
  !end function equals


  function equalEpsilon (number1, number2, epsilon ) result (resultValue)
    real , intent (in) :: number1, number2, epsilon 
    logical :: resultValue 

    resultValue = .false.

    ! test very small number1
    if ( abs(number1) < epsilon .and.  abs(number1 - number2) < epsilon ) then
       resultValue = .true.
    else 
       if ((abs(( number1 - number2)) / number1) < epsilon ) then
          resultValue = .true.
       else
          resultValue = .false.
       end if
    end if

  end function equalEpsilon

  function floatEqual (number1, number2 ) result (resultValue)
    real , intent (in) :: number1, number2
    real :: epsilon 
    logical :: resultValue 

    resultValue = .false.
    epsilon = 1E-6

    ! test very small number1
    if ( abs(number1) < epsilon .and.  abs(number1 - number2) < epsilon ) then
       resultValue = .true.
    else 
       if ((abs(( number1 - number2)) / number1) < epsilon ) then
          resultValue = .true.
       else
          resultValue = .false.
       end if
    end if
  end function floatEqual

  function doublePrecisionEqual (number1, number2 ) result (resultValue)
    double precision , intent (in) :: number1, number2
    real :: epsilon 
    logical :: resultValue 

    resultValue = .false.
    epsilon = 1E-6
    !epsilon = epsilon (number1)

    ! test very small number1
    if ( abs(number1) < epsilon .and.  abs(number1 - number2) < epsilon ) then
       resultValue = .true.
    else 
       if ((abs(( number1 - number2)) / number1) < epsilon ) then
          resultValue = .true.
       else
          resultValue = .false.
       end if
    end if
  end function doublePrecisionEqual

  function integerEqual (number1, number2 ) result (resultValue)
    integer , intent (in) :: number1, number2
    logical :: resultValue 

    resultValue = .false.

    if ( number1 .eq. number2 ) then
       resultValue = .true.
    else 
       resultValue = .false.
    end if
  end function integerEqual

  function stringEqual (str1, str2 ) result (resultValue)
    character(*) , intent (in) :: str1, str2
    logical :: resultValue 

    resultValue = .false.

    if ( str1 .eq. str2 ) then
       resultValue = .true.
    end if
  end function stringEqual

  function logicalEqual (l1, l2 ) result (resultValue)
    logical, intent (in) :: l1, l2
    logical              :: resultValue 

    resultValue = .false.

    if ( l1 .eqv. l2 ) then
       resultValue = .true.
    end if
  end function logicalEqual

end module fruit_util


module fruit
  use fruit_util
  implicit none

  integer, parameter :: MSG_LENGTH = 1500
  integer, parameter :: MSG_STACK_SIZE = 300000

  integer, private, save :: successful_assert_count = 0
  integer, private, save :: failed_assert_count = 0
  character (len = MSG_LENGTH), private, dimension (MSG_STACK_SIZE), save :: messageArray
  character (len = MSG_LENGTH), private, save :: msg = '[unit name not set from set_name]: '
  character (len = MSG_LENGTH), private, save :: unit_name  = '_not_set_'
  integer, private, save :: messageIndex = 1

  integer, private, save :: successful_case_count = 0
  integer, private, save :: failed_case_count = 0
  integer, private, save :: testCaseIndex = 1
  logical, private, save :: last_passed = .false.

  interface init_fruit
     module procedure init_fruit_
  end interface

  interface initializeFruit
     module procedure obsolete_subroutine_delete_later_initializeFruit_
  end interface

  interface fruit_summary
     module procedure fruit_summary_
  end interface

  interface getTestSummary
     module procedure obsolete_subroutine_delete_later_getTestSummary_  
  end interface

  interface get_last_message
     module procedure get_last_message_
  end interface

  interface is_last_passed
     module procedure is_last_passed_
  end interface

  interface assert_true
     module procedure assert_true_logical_
  end interface

  interface assertTrue
     module procedure obsolete_delete_this_later_assert_true_logical_
  end interface

  interface assert_equals
     module procedure assert_equals_int_
     module procedure assert_equals_double_
     module procedure assert_equals_real_
     module procedure assert_equals_logical_
     module procedure assert_equals_string_
     module procedure assert_equals_complex_
     module procedure assert_equals_real_within_range_
     module procedure assert_equals_double_within_range_

     module procedure assert_equals_1d_int_
     module procedure assert_equals_1d_double_
     module procedure assert_equals_1d_real_
     module procedure assert_equals_1d_string_
     module procedure assert_equals_1d_complex_
     module procedure assert_equals_1d_real_within_range_
     module procedure assert_equals_1d_double_within_range_

     module procedure assert_equals_2d_int_
     module procedure assert_equals_2d_double_
     module procedure assert_equals_2d_real_
     module procedure assert_equals_2d_complex_

  end interface

  interface assertEquals
     module procedure assert_equals_int_
     module procedure assert_equals_double_
     module procedure assert_equals_real_
     module procedure assert_equals_logical_
     module procedure assert_equals_string_
     module procedure assert_equals_complex_
     module procedure assert_equals_real_within_range_
     module procedure assert_equals_double_within_range_

     module procedure assert_equals_1d_int_
     module procedure assert_equals_1d_double_
     module procedure assert_equals_1d_real_
     module procedure assert_equals_1d_string_
     module procedure assert_equals_1d_complex_
     module procedure assert_equals_1d_real_within_range_
     module procedure assert_equals_1d_double_within_range_

     module procedure assert_equals_2d_int_
     module procedure assert_equals_2d_double_
     module procedure assert_equals_2d_real_
     module procedure assert_equals_2d_complex_

  end interface

  interface assert_not_equals
     module procedure assert_not_equals_real_
     module procedure assert_not_equals_1d_real_
     module procedure assert_not_equals_double_
  end interface

  interface assertNotEquals
     module procedure assert_not_equals_real_
     module procedure assert_not_equals_1d_real_
     module procedure assert_not_equals_double_
  end interface

  interface add_success
     module procedure add_success_
  end interface

  interface addSuccess
     module procedure obsolete_subroutine_delete_later_addSuccess_
  end interface

  interface add_fail
     module procedure add_fail_
     module procedure add_fail_unit_
  end interface

  interface addFail
     module procedure add_fail_
     module procedure add_fail_unit_
  end interface

  interface set_unit_name
     module procedure set_unit_name_
  end interface

  interface get_unit_name
     module procedure get_unit_name_
  end interface

  interface getTotalCount
     module procedure obsolete_subroutine_delete_later_getTotalCount_
  end interface

  interface get_total_count
     module procedure get_total_count_
  end interface

  interface getFailedCount
     module procedure obsolete_subroutine_delete_later_getFailedCount_
  end interface

  interface get_failed_count
     module procedure get_failed_count_
  end interface

  interface success_assert_action
     module procedure success_assert_action_
  end interface

  interface failed_assert_action
     module procedure failed_assert_action_
  end interface

  interface isAllSuccessful
     module procedure obsolete_subroutine_delete_later_isAllSuccessful_
  end interface

  interface is_all_successful
     module procedure is_all_successful_
  end interface

  private :: &
       init_fruit_, obsolete_subroutine_delete_later_initializeFruit_, &
       fruit_summary_, obsolete_subroutine_delete_later_getTestSummary_, &
       get_last_message_, obsolete_, make_error_msg_, &
       get_total_count_, obsolete_subroutine_delete_later_getTotalCount_,&
       get_failed_count_, obsolete_subroutine_delete_later_getFailedCount_, &
       add_fail_, add_fail_unit_, increase_message_stack_, &
       success_assert_action_, failed_assert_action_, success_mark_, failed_mark_, &
       assert_true_logical_, obsolete_delete_this_later_assert_true_logical_, &
       assert_equals_int_, &
       assert_equals_double_, &
       assert_equals_real_, &
       assert_equals_logical_, &
       assert_equals_string_, &
       assert_equals_complex_, &
       assert_equals_real_within_range_, &
       assert_equals_double_within_range_, &
       assert_equals_1d_int_, &
       assert_equals_1d_double_, &
       assert_equals_1d_real_, &
       assert_equals_1d_string_, &
       assert_equals_1d_complex_, &
       assert_equals_2d_int_, &
       assert_equals_2d_double_, &
       assert_equals_2d_real_, &
       assert_equals_2d_complex_, &
       assert_equals_1d_real_within_range_, &
       assert_equals_1d_double_within_range_, &
       assert_not_equals_real_, &
       assert_not_equals_double_, &
       assert_not_equals_1d_real_, &
       add_success_, obsolete_subroutine_delete_later_addSuccess_, &
       is_all_successful_, obsolete_subroutine_delete_later_isAllSuccessful_, &
       set_unit_name_, get_unit_name_, is_last_passed_

contains

  subroutine init_fruit_
    successful_assert_count = 0
    failed_assert_count = 0
    messageIndex = 1
    write (*,*)
    write (*,*) "Test module initialized"
    write (*,*)
    write (*,*) "   . : successful assert,   F : failed assert "
    write (*,*)
  end subroutine init_fruit_

  subroutine obsolete_subroutine_delete_later_initializeFruit_
    call obsolete_ ("initializeFruit is OBSOLETE.  replaced by init_fruit")
    call init_fruit
  end subroutine obsolete_subroutine_delete_later_initializeFruit_

  subroutine obsolete_subroutine_delete_later_getTestSummary_
    call obsolete_ ( "getTestSummary is OBSOLETE.  replaced by fruit_summary")
    call fruit_summary
  end subroutine obsolete_subroutine_delete_later_getTestSummary_

  subroutine fruit_summary_
    integer :: i

    write (*,*)
    write (*,*)
    write (*,*) '    Start of FRUIT summary: '
    write (*,*)

    if (failed_assert_count > 0) then
       write (*,*) 'Some tests failed!'
    else
       write (*,*) 'SUCCESSFUL!'
    end if

    write (*,*)
    if ( messageIndex > 1) then
       write (*,*) '  -- Failed assertion messages:'

       do i = 1, messageIndex - 1
          write (*,"(A)") trim(strip(messageArray(i)))
       end do

       write (*,*) '  -- end of failed assertion messages.'
       write (*,*)
    else
       write (*,*) '  No messages '
    end if

    if (successful_assert_count + failed_assert_count /= 0) then

       write (*,*) 'Total asserts :   ', successful_assert_count + failed_assert_count
       write (*,*) 'Successful    :   ', successful_assert_count
       write (*,*) 'Failed        :   ', failed_assert_count
       write (*,'("Successful rate:   ",f6.2,"%")')  real(successful_assert_count) * 100.0 / &
            real (successful_assert_count + failed_assert_count)
       write (*, *)
       write (*,*) 'Successful asserts / total asserts : [ ',&
            successful_assert_count, '/', successful_assert_count + failed_assert_count, ' ]'
       write (*,*) 'Successful cases   / total cases   : [ ', successful_case_count, '/', &
            successful_case_count + failed_case_count, ' ]'
       write (*, *) '  -- end of FRUIT summary'

    end if
  end subroutine fruit_summary_

  subroutine obsolete_subroutine_delete_later_addSuccess_
    call obsolete_ ("addSuccess is OBSOLETE.  replaced by add_success")
    call add_success_
  end subroutine obsolete_subroutine_delete_later_addSuccess_

  subroutine add_success_
    call success_assert_action_
  end subroutine add_success_

  subroutine add_fail_ (message)
    character (*), intent (in), optional :: message
    call failed_assert_action_('none', 'none', message)
  end subroutine add_fail_

  subroutine add_fail_unit_ (unitName, message)
    character (*), intent (in) :: unitName
    character (*), intent (in) :: message

    call add_fail_ ("[in " //  unitName // "(fail)]: " //  message)
  end subroutine add_fail_unit_

  subroutine obsolete_subroutine_delete_later_isAllSuccessful_(result)
    logical, intent(out) :: result
    call obsolete_ ('subroutine isAllSuccessful is changed to function is_all_successful.')
    result = (failed_assert_count .eq. 0 )
  end subroutine obsolete_subroutine_delete_later_isAllSuccessful_

  subroutine is_all_successful_(result)
    logical, intent(out) :: result
    result= (failed_assert_count .eq. 0 )
  end subroutine is_all_successful_

  subroutine success_mark_
    write(*,"(A1)",ADVANCE='NO') '.'
  end subroutine success_mark_

  subroutine failed_mark_
    write(*,"(A1)",ADVANCE='NO') 'F'
  end subroutine failed_mark_

  subroutine increase_message_stack_
    if (messageIndex > MSG_STACK_SIZE) then
       write (*, *) "Too many errors to put into stack"
       call getTestSummary ()
       stop 1
    end if

    messageArray (messageIndex) = msg
    messageIndex = messageIndex + 1
  end subroutine increase_message_stack_

!  function get_last_message_
  function get_last_message_() ! JW2012
    character(len=MSG_LENGTH) :: get_last_message_
    if (messageIndex > 1) then
       get_last_message_ = strip(messageArray(messageIndex-1))
    else
       get_last_message_ = ''
    end if
  end function get_last_message_

  subroutine obsolete_subroutine_delete_later_getTotalCount_ (count)
    integer, intent (out) :: count
    call obsolete_ (' getTotalCount subroutine is replaced by function get_total_count')
    call get_total_count_(count)
  end subroutine obsolete_subroutine_delete_later_getTotalCount_

  subroutine get_total_count_(count) 
    integer, intent(out) :: count

    count = successful_assert_count + failed_assert_count
  end subroutine get_total_count_

  subroutine obsolete_subroutine_delete_later_getFailedCount_ (count)
    integer, intent (out) :: count

    call obsolete_ (' getFailedCount subroutine is replaced by function get_failed_count')
    call get_failed_count_ (count)

  end subroutine obsolete_subroutine_delete_later_getFailedCount_

  subroutine get_failed_count_ (count)
    integer, intent(out) :: count
    count = failed_assert_count
  end subroutine get_failed_count_

  subroutine obsolete_ (message)
    character (*), intent (in), optional :: message
    write (*,*) 
    write (*,*) "<<<<<<<<<<<<<<<<<<<<<<<<<< WARNING from FRUIT  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    write (*,*) message
    write (*,*) 
    write (*,*) " old calls will be replaced in the next release in Jan 2009"
    write (*,*) " Naming convention for all the method calls are changed to: first_name from firstName"
    write (*,*) " Subroutines will be deleted: assertEquals, assertNotEquals, assertTrue, addSuccessful, addFail, etc."
    write (*,*) "<<<<<<<<<<<<<<<<<<<<<<<<<< WARNING from FRUIT  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    write (*,*) 
  end subroutine obsolete_

  subroutine success_assert_action_
    successful_assert_count = successful_assert_count + 1
    last_passed = .true.
    call success_mark_  
  end subroutine success_assert_action_

  subroutine failed_assert_action_ (expected, got, message)
    character(*), intent(in) :: expected, got
    character(*), intent(in), optional :: message

    call make_error_msg_ (expected, got, message)
    call increase_message_stack_
    failed_assert_count = failed_assert_count + 1
    last_passed = .false.
    call failed_mark_
  end subroutine failed_assert_action_

  subroutine set_unit_name_(value)
    character(*), intent(in) :: value
    unit_name = strip(value)
  end subroutine set_unit_name_

  subroutine get_unit_name_(value)
    character(*), intent(out) :: value
    value = strip(unit_name)
  end subroutine get_unit_name_

  subroutine make_error_msg_ (var1, var2, message)
    character(*), intent(in) :: var1, var2
    character(*), intent(in), optional :: message
    if (present(message)) then
       msg = '[' // trim(strip(unit_name)) // ']: Expected [' // trim(strip(var1)) // '], Got [' // trim(strip(var2)) // ']; '// &
            ' User message: [' // message // ']'
    else
       msg = '[' // trim(strip(unit_name)) // ']: Expected [' // trim(strip(var1)) // '], Got [' // trim(strip(var2)) // ']' 
    endif
  end subroutine make_error_msg_

!  function is_last_passed_
  function is_last_passed_() ! JW2012
    logical:: is_last_passed_
    is_last_passed_ = last_passed 
  end function is_last_passed_

  !--------------------------------------------------------------------------------
  ! all assertions
  !--------------------------------------------------------------------------------
  subroutine obsolete_delete_this_later_assert_true_logical_(var1, message)
    logical, intent (in) :: var1
    character (*), intent (in), optional :: message

    call obsolete_ ('assertTrue subroutine is replaced by function assert_true')
    call assert_true_logical_(var1, message)
  end subroutine obsolete_delete_this_later_assert_true_logical_

  subroutine assert_true_logical_ (var1, message)
    logical, intent (in) :: var1
    character (*), intent (in), optional :: message

    if ( var1 .eqv. .true.) then
       call success_assert_action_
    else
       call failed_assert_action_(to_s(.true.), to_s(var1), message)
    end if
  end subroutine assert_true_logical_

  subroutine assert_equals_int_ (var1, var2, message)
    integer, intent(in) :: var1, var2
    character (*), intent(in), optional :: message

    if ( var1 .eq. var2) then
       call success_assert_action_
    else
       call failed_assert_action_ (to_s(var1), to_s(var2), message)
    end if
  end subroutine assert_equals_int_

  subroutine assert_equals_logical_ (var1, var2, message)
    logical, intent (in)  :: var1, var2
    character (*), intent (in), optional :: message

    if ( var1 .eqv. var2 ) then
       call success_assert_action_
    else
       call failed_assert_action_(to_s(var1), to_s(var2), message)
    end if
  end subroutine assert_equals_logical_

  subroutine assert_equals_string_ (var1, var2, message)
    character(*), intent (in)  :: var1, var2
    character (*), intent (in), optional :: message

    if ( trim(strip(var1)) == trim(strip(var2))) then
       call success_assert_action_
    else
       call failed_assert_action_(var1, var2, message)
    end if
  end subroutine assert_equals_string_

  subroutine assert_equals_real_ (var1, var2, message)
    real, intent (in) :: var1, var2
    character (*), intent (in), optional :: message

    if ( var1 .eq. var2) then
       call success_assert_action_
    else
7      call failed_assert_action_(to_s(var1), to_s(var2), message)
    end if
  end subroutine assert_equals_real_

  subroutine assert_equals_double_ (var1, var2, message)
    double precision, intent (in) :: var1, var2
    character(*), intent(in), optional :: message

    if ( var1 .eq. var2) then
       call success_assert_action_
    else
       call failed_assert_action_(to_s(var1), to_s(var2), message)
    end if
  end subroutine assert_equals_double_

  subroutine assert_equals_complex_ (var1, var2, message)
    double complex,   intent(IN) :: var1, var2
    character (*),    intent(IN), optional :: message
    integer count

    if ( var1 .ne. var2) then
       call failed_assert_action_(to_s(var1), to_s(var2), message)
    else
       call success_assert_action_
    end if

  end subroutine assert_equals_complex_

  subroutine assert_equals_real_within_range_(var1, var2, var3, message)
    real, intent (in) :: var1, var2, var3
    character(*), intent(in), optional :: message

    if ( abs( var1 - var2) .le. var3) then
       call success_assert_action_
    else
       call failed_assert_action(to_s(var1), to_s(var2), message)
    end if

  end subroutine assert_equals_real_within_range_

  subroutine assert_equals_double_within_range_(var1, var2, var3, message)
    double precision, intent (in) :: var1, var2, var3
    character(*), intent(in), optional :: message

    if ( abs( var1 - var2) .le. var3) then
       call success_assert_action_
    else
       call failed_assert_action_(to_s(var1), to_s(var2), message)
    end if
  end subroutine assert_equals_double_within_range_

  subroutine assert_equals_1d_int_ (var1, var2, n, message)
    integer, intent (in) :: n
    integer, intent (in) :: var1(n), var2(n)
    character (*), intent (in), optional :: message

    integer count

    loop_dim1: do count = 1, n
       if ( var1(count) .ne. var2(count)) then
          call failed_assert_action_(to_s(var1(count)), to_s(var2(count)), message)
          return
       end if
    end do loop_dim1

    call success_assert_action_
  end subroutine assert_equals_1d_int_

  subroutine assert_equals_1d_string_ (var1, var2, n, message)
    integer, intent (in) :: n
    character(*), intent (in) :: var1(n), var2(n)
    character (*), intent (in), optional :: message
    integer count

    loop_dim1: do count = 1, n
       if ( strip(var1(count)) .ne. strip(var2(count))) then
          call failed_assert_action_(var1(count), var2(count), message)
          return
       end if
    end do loop_dim1

    call success_assert_action_
  end subroutine assert_equals_1d_string_

  subroutine assert_equals_1d_real_within_range_(var1, var2, n, var3, message)
    integer, intent(in) :: n
    real, intent (in) :: var1(n), var2(n), var3
    character(*), intent(in), optional :: message

    if ( maxval( abs( var1 - var2)) .le. var3) then
       call success_assert_action_
    else
       call failed_assert_action_(to_s(var1(1)), to_s(var2(1)), '1D array real has difference' // ' ' // message)
    end if
  end subroutine assert_equals_1d_real_within_range_

  subroutine assert_equals_1d_double_within_range_(var1, var2, n, var3, message)
    integer, intent(in) :: n
    double precision, intent (in) :: var1(n), var2(n), var3
    character(*), intent(in), optional :: message

    if ( maxval( abs( var1 - var2)) .le. var3) then
       call success_assert_action_
    else
       call failed_assert_action_(to_s(var1(1)), to_s(var2(1)), message)
    end if
  end subroutine assert_equals_1d_double_within_range_

  subroutine assert_equals_1d_double (var1, var2, n, message)
    integer, intent (in) :: n
    double precision, intent (in) :: var1(n), var2(n)
    character(*), intent(in), optional :: message

    integer count

    loop_dim1: do count = 1, n
       if ( var1(count) .ne. var2(count)) then
          call failed_assert_action_(to_s(var1(count)), to_s(var2(count)), &
               'Array different at count: ' // to_s(count) // ' ' // message)
          return
       end if
    end do loop_dim1

    call success_assert_action_
  end subroutine assert_equals_1d_double

  subroutine assert_equals_2d_real (var1, var2, n, m)
    integer, intent (in) :: n, m
    real, intent (in) :: var1(n,m), var2(n,m)

    integer count1, count2

    loop_dim2: do count2 = 1, m
       loop_dim1: do count1 = 1, n
          if ( var1(count1,count2) .ne. var2(count1,count2)) then
             call failed_assert_action(to_s(var1(count1, count2)), to_s(var2(count1, count2)),&
                  'Array (' // to_s(count1) // ',' // to_s( count2) //')')
             return
          end if
       end do loop_dim1
    end do loop_dim2

    call success_assert_action_
  end subroutine assert_equals_2d_real

  subroutine assert_equals_2d_double (var1, var2, n, m)
    integer, intent (in) :: n, m
    double precision, intent (in) :: var1(n,m), var2(n,m)

    integer count1, count2

    loop_dim2: do count2 = 1, m
       loop_dim1: do count1 = 1, n
          if ( var1(count1,count2) .ne. var2(count1,count2)) then
             call failed_assert_action_(to_s(var1(count1, count2)), to_s(var2(count1, count2)), &
                  'Array difference at (' // to_s(count1) // ',' // to_s(count2) // ')')
             return
          end if
       end do loop_dim1
    end do loop_dim2

    call success_assert_action_
  end subroutine assert_equals_2d_double

  subroutine assert_equals_2d_int_ (var1, var2, n, m, message)
    integer, intent (in) :: n, m
    integer, intent (in) :: var1(n,m), var2(n,m)
    character (*), intent (in), optional :: message

    integer count1, count2

    loop_dim2: do count2 = 1, m
       loop_dim1: do count1 = 1, n
          if ( var1(count1,count2) .ne. var2(count1,count2)) then
             call failed_assert_action_(to_s(var1(count1, count2)), to_s(var2(count1, count2)), message)
             return
          end if
       end do loop_dim1
    end do loop_dim2

    call success_assert_action_
  end subroutine assert_equals_2d_int_

  subroutine assert_equals_1d_real_ (var1, var2, n, message)
    integer, intent (in) :: n
    real, intent (in) :: var1(n), var2(n)
    character (*), intent (in), optional :: message

    integer count

    loop_dim1: do count = 1, n
       if ( var1(count) .ne. var2(count)) then
          call failed_assert_action_(to_s(var1(count)), to_s(var2(count)), message)
          return
       end if
    end do loop_dim1
    call success_assert_action_
  end subroutine assert_equals_1d_real_

  subroutine assert_equals_2d_real_ (var1, var2, n, m, message)
    integer, intent (in) :: n, m
    real, intent (in) :: var1(n,m), var2(n,m)
    character (*), intent(in), optional :: message

    integer count1, count2

    loop_dim2: do count2 = 1, m
       loop_dim1: do count1 = 1, n
          if ( var1(count1,count2) .ne. var2(count1,count2)) then
             call failed_assert_action_(to_s(var1(count1, count2)), to_s(var2(count1, count2)), message)
             return
          end if
       end do loop_dim1
    end do loop_dim2

    call success_assert_action_
  end subroutine assert_equals_2d_real_

  subroutine assert_equals_1d_double_ (var1, var2, n, message)
    integer, intent (in) :: n
    double precision, intent (in) :: var1(n), var2(n)
    character (*), intent (in), optional :: message
    integer count

    loop_dim1: do count = 1, n
       if ( var1(count) .ne. var2(count)) then
          call failed_assert_action_(to_s(var1(count)), to_s(var2(count)), message)
          return
       end if
    end do loop_dim1

    call success_assert_action_
  end subroutine assert_equals_1d_double_

  subroutine assert_equals_2d_double_ (var1, var2, n, m, message)
    integer, intent (in) :: n, m
    double precision, intent (in) :: var1(n,m), var2(n,m)
    character (*), intent (in), optional :: message
    integer count1, count2

    loop_dim2: do count2 = 1, m
       loop_dim1: do count1 = 1, n
          if ( var1(count1,count2) .ne. var2(count1,count2)) then
             call failed_assert_action_(to_s(var1(count1, count2)), to_s(var2(count1, count2)), message)
             return
          end if
       end do loop_dim1
    end do loop_dim2

    call success_assert_action_
  end subroutine assert_equals_2d_double_

  subroutine assert_equals_1d_complex_ (var1, var2, n, message)
    integer,          intent(IN) :: n
    double complex,   intent(IN) :: var1(n), var2(n)
    character (*),    intent(IN), optional :: message
    integer count

    loop_dim1: do count = 1, n
       if ( var1(count) .ne. var2(count)) then
          call failed_assert_action_(to_s(var1(count)), to_s(var2(count)), message)
          return
       end if
    enddo loop_dim1

    call success_assert_action_
  end subroutine assert_equals_1d_complex_

  subroutine assert_equals_2d_complex_ (var1, var2, n, m, message)
    integer,        intent(IN) :: n, m
    double complex, intent(IN) :: var1(n,m), var2(n,m)
    character (*),    intent(IN), optional :: message
    integer count1, count2

    loop_dim2: do count2 = 1, m
       loop_dim1: do count1 = 1, n
          if ( var1(count1,count2) .ne. var2(count1,count2)) then
             call failed_assert_action_(to_s(var1(count1, count2)), to_s(var2(count1, count2)), message)
             return
          endif
       enddo loop_dim1
    enddo loop_dim2

    call success_assert_action_
  end subroutine assert_equals_2d_complex_

  subroutine assert_not_equals_real_ (var1, var2, message)
    real, intent (in) :: var1, var2
    character (*), intent (in), optional :: message

    if ( var1 .ne. var2) then
       call success_assert_action_
    else
       call failed_assert_action_(to_s(var1), to_s(var2), message)
    end if
  end subroutine assert_not_equals_real_

  subroutine assert_not_equals_double_ (var1, var2, message)
    double precision, intent (in) :: var1, var2
    character(*), intent(in), optional :: message

    if ( var1 .ne. var2) then
       call success_assert_action_
    else
       call failed_assert_action_(to_s(var1), to_s(var2), message)
    end if
  end subroutine assert_not_equals_double_

  subroutine assert_not_equals_1d_real_ (var1, var2, n)
    integer, intent (in) :: n
    real, intent (in) :: var1(n), var2(n)

    integer count

    loop_dim1: do count = 1, n
       if ( var1(count) .ne. var2(count)) then
          call failed_assert_action(to_s(var1(count)), to_s(var2(count)),&
               'Array (' // to_s(count)//')')
          return
       end if
    end do loop_dim1

    call success_assert_action_

  end subroutine assert_not_equals_1d_real_

end module fruit
